# pdm5b_simultest5.py - Testing JAGS fits of Hierarchical Diffusion Models with
#         trial-to-trial variability in visual encoding time (VET) and drift
#
# Copyright (C) 2019 Michael D. Nunez, <mdnunez1@uci.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 12/10/18      Michael Nunez              Converted from pdm5b_simultest4.py
# 12/19/18      Michael Nunez      Fix lapse-trial generation in simuldiffn200()
# 01/03/19      Michael Nunez          Increase range of residual motor response time based on 
#                                      empirical NDT data (bottom left subfigure 6 in Nunez et al. 2018)

# Modules
from __future__ import division
import numpy as np
import pyjags
import scipy.io as sio
from scipy.stats import truncnorm
import warnings
import random

# Simulate diffusion models
def simuldiffn200(N=100,Alpha=1,Vet=.2,Rmr=.2,Nu=1,Zeta=None,rangeVet=0,rangeRmr=0,rangeZeta=0,Eta=.3,Varsigma=1):
    """
    SIMULDIFF  Generates data according to a diffusion model


    Reference:
    Tuerlinckx, F., Maris, E.,
    Ratcliff, R., & De Boeck, P. (2001). A comparison of four methods for
    simulating the diffusion process. Behavior Research Methods,
    Instruments, & Computers, 33, 443-456.

    Parameters
    ----------
    N: a integer denoting the size of the output vector
    (defaults to 100 experimental trials)

    Alpha: the mean boundary separation across trials  in evidence units
    (defaults to 1 evidence unit)

    Vet: the mean visual encoding time across trials in seconds
    (defaults to .2 seconds)

    Rmr: the mean residual motor response time across trials in seconds
    (defaults to .2 seconds)

    Nu: the mean drift rate across trials in evidence units per second
    (defaults to 1 evidence units per second, restricted to -5 to 5 units)

    Zeta: the initial bias in the evidence process for choice A
    (defaults to 50% of total evidence units given by Alpha)

    rangeVet: Visual encoding time across trials is generated from a uniform
    distribution of Vet - rangeVet/2 to  Vet + rangeVet/2 across trials
    (defaults to 0 seconds)

    rangeRmr: Residual motor response time across trials is generated from a uniform
    distribution of Rmr - rangeRmr/2 to  Rmr + rangeRmr/2 across trials
    (defaults to 0 seconds)

    rangeZeta: Bias across trials is generated from a uniform distribution
    of Zeta - rangeZeta/2 to Zeta + rangeZeta/2 across trials
    (defaults to 0 evidence units)

    Eta: Standard deviation of the drift rate across trials
    (defaults to 3 evidence units per second, restricted to less than 3 evidence units)

    Varsigma: The diffusion coefficient, the standard deviation of the
    evidence accumulation process within one trial. It is recommended that
    this parameter be kept fixed unless you have reason to explore this parameter
    (defaults to 1 evidence unit per second)

    Returns
    -------
    Numpy complex vector with 1) Real component ( np.real(x) ): reaction times (in seconds) multiplied by the response vector
    such that negative reaction times encode response B and positive reaction times
    encode response A  and 2) Imaginary component ( np.imag(x) ): N200 peak-latencies in seconds
    
    
    Converted from simuldiff.m MATLAB script by Joachim Vandekerckhove
    See also http://ppw.kuleuven.be/okp/dmatoolbox.
    """

    if Zeta is None:
        Zeta = .5*Alpha

    if (Nu < -5) or (Nu > 5):
        Nu = np.sign(Nu)*5
        warnings.warn('Nu is not in the range [-5 5], bounding drift rate to %.1f...' % (Nu))

    if (Eta > 3):
        warning.warn('Standard deviation of drift rate is out of bounds, bounding drift rate to 3')
        eta = 3

    if (Eta == 0):
        Eta = 1e-16

    #Initialize output vectors
    result = np.zeros(N)
    T = np.zeros(N)
    XX = np.zeros(N)
    N200 = np.zeros(N)

    #Called sigma in 2001 paper
    D = np.power(Varsigma,2)/2

    #Program specifications
    eps = 2.220446049250313e-16 #precision from 1.0 to next double-precision number
    delta=eps

    for n in range(0,N):
        r1 = np.random.normal()
        mu = Nu + r1*Eta
        zz = Zeta - rangeZeta/2 + rangeZeta*np.random.uniform()
        finish = 0
        totaltime = 0
        startpos = 0
        Aupper = Alpha - zz
        Alower = -zz
        radius = np.min(np.array([np.abs(Aupper), np.abs(Alower)]))
        while (finish==0):
            lambda_ = 0.25*np.power(mu,2)/D + 0.25*D*np.power(np.pi,2)/np.power(radius,2)
            # eq. formula (13) in 2001 paper with D = sigma^2/2 and radius = Alpha/2
            F = D*np.pi/(radius*mu)
            F = np.power(F,2)/(1 + np.power(F,2) )
            # formula p447 in 2001 paper
            prob = np.exp(radius*mu/D)
            prob = prob/(1 + prob)
            dir_ = 2*(np.random.uniform() < prob) - 1
            l = -1
            s2 = 0
            while (s2>l):
                s2=np.random.uniform()
                s1=np.random.uniform()
                tnew=0
                told=0
                uu=0
                while (np.abs(tnew-told)>eps) or (uu==0):
                    told=tnew
                    uu=uu+1
                    tnew = told + (2*uu+1) * np.power(-1,uu) * np.power(s1,(F*np.power(2*uu+1,2)));
                    # infinite sum in formula (16) in BRMIC,2001
                l = 1 + np.power(s1,(-F)) * tnew;
            # rest of formula (16)
            t = np.abs(np.log(s1))/lambda_;
            # is the negative of t* in (14) in BRMIC,2001
            totaltime=totaltime+t
            dir_=startpos+dir_*radius
            vetime = Vet - rangeVet/2 + rangeVet*np.random.uniform()
            rmrt = Rmr - rangeRmr/2 + rangeRmr*np.random.uniform()
            if ( (dir_ + delta) > Aupper):
                T[n]=vetime+totaltime+rmrt
                XX[n]=1
                N200[n]=vetime
                finish=1
            elif ( (dir_-delta) < Alower ):
                T[n]=vetime+totaltime+rmrt
                XX[n]=-1
                N200[n]=vetime
                finish=1
            else:
                startpos=dir_
                radius=np.min(np.abs([Aupper, Alower]-startpos))

    result = T*XX + N200*1.j
    return result



# Number of simulations
nsims = 30

# Number of simulated subjects
nsubs = 100

# Number of trials per subject
ntrials = 100

# Generate samples from the joint-model of reactiont time, choice, and N200 latencies

# Set random seed
random.seed(11)

vetsub = np.random.uniform(.15, .275, size=(nsims, nsubs)) # Uniform from .15 to .275 seconds
rmrsub = np.random.uniform(.15, .475, size=(nsims, nsubs)) # Uniform from .15 to .475 seconds
alphasub = np.random.uniform(.8, 1.4, size=(nsims, nsubs)) # Uniform from .8 to 1.4 evidence units
deltasub = np.random.uniform(-4, 4, size=(nsims, nsubs)) # Uniform from -4 to 4 evidence units per second
vettrialrange = np.random.uniform(0,.1, size=(nsims,nsubs)) # Uniform from 0 to .1 seconds
rmrtrialrange = np.random.uniform(0,.1, size=(nsims,nsubs)) # Uniform from 0 to .1 seconds
deltatrialsd = np.random.uniform(0, 2, size=(nsims,nsubs)) # Uniform from 0 to 2 evidence units per second
prob_mindwander = np.matlib.repmat(np.linspace(0, .1, num=nsims),nsubs,1).T # From 0 to 10 percent of trials
y = np.zeros((nsims, nsubs, ntrials))
rt = np.zeros((nsims, nsubs, ntrials))
acc = np.zeros((nsims, nsubs, ntrials))
n200 = np.zeros((nsims, nsubs, ntrials))
for n in range(0,nsims):
    for s in range(0,nsubs):
        tempout = simuldiffn200(N=ntrials, Alpha= alphasub[n,s], Vet= vetsub[n,s], Rmr= rmrsub[n,s], 
            Nu= deltasub[n,s], Eta= deltatrialsd[n,s], rangeVet=vettrialrange[n,s], rangeRmr=rmrtrialrange[n,s])
        tempx = np.sign(np.real(tempout))
        tempt = np.abs(np.real(tempout))
        mindwanderx = np.random.randint(low=0,high=2,size=ntrials)*2 -1
        mindwandert = np.random.uniform(low=0,high=2,size=ntrials) # Randomly distributed from 0 to 2 seconds

        mindwander_trials = np.random.choice(ntrials, size=np.round(ntrials*prob_mindwander[n,s]), replace=False)
        tempx[mindwander_trials] = mindwanderx[mindwander_trials]
        tempt[mindwander_trials] = mindwandert[mindwander_trials]
        y[n,s,:] = tempx*tempt
        rt[n,s,:] = tempt
        acc[n,s,:] = (tempx + 1)/2
        #Assume N200s are unaffect by the mind wandering process
        n200[n,s,:] = np.imag(tempout)

genparam = dict()
genparam['vetsub'] = vetsub
genparam['rmrsub'] = rmrsub
genparam['alphasub'] = alphasub
genparam['deltasub'] = deltasub
genparam['vettrialrange'] = vettrialrange
genparam['rmrtrialrange'] = rmrtrialrange
genparam['deltatrialsd'] = deltatrialsd
genparam['prob_mindwander'] = prob_mindwander
genparam['rt'] = rt
genparam['acc'] = acc
genparam['y'] = y
genparam['n200'] = n200
sio.savemat('genparam_test4.mat', genparam)

genparam = sio.loadmat('genparam_test4.mat')

# Send to JAGS
nchains = 6
burnin = 2000  # Note that scientific notation breaks pyjags
nsamps = 10000

pyjags.modules.load_module('wiener')
pyjags.modules.load_module('dic')
pyjags.modules.list_modules()

tojags = '''
model {
    ##########
    # Fixed Parameters
    beta <- .5
    ##########
    # Between-subject variability in non-decision time
    tersubsd ~ dgamma(.2,1)

    # Between-subject variability in drift
    deltasubsd ~ dgamma(1,1)

    # Between-subject variability in boundary separation
    alphasubsd ~ dgamma(1,1)

    ##########
    # Block-level parameters
    ##########

    # Condition-level non-decision time
    tercond ~ dnorm(.3, pow(.25,-2))

    # Condition-level drift rate
    deltacond ~ dnorm(1, pow(2, -2))

    # Condition-level boundary separation
    alphacond ~ dnorm(1, pow(.5,-2))

    # Subject-level parameters
    for (sub in 1:nsubs) {
        # Subject-level non-decision time
        tersub[sub] ~ dnorm(tercond, pow(tersubsd, -2))T(0,1)

        # Subject-level drift rate
        deltasub[sub] ~ dnorm(deltacond, pow(deltasubsd, -2))T(-9, 9)

        # Subject-level boundary separation
        alphasub[sub] ~ dnorm(alphacond, pow(alphasubsd, -2))T(.1,3)

        # Subject-level mind wandering
        probsub[sub,1:2] ~ ddirch(c(1,1))
    }

    ##########
    # Wiener likelihoods
    for (i in 1:N) {
        # Log density for DDM process
        ld_comp[i, 1] <- dlogwiener(y[i],alphasub[subject[i]],
        tersub[subject[i]],
        beta,
        deltasub[subject[i]])

        # Log density for lapse trials (negative max RT to positive max RT)
        ld_comp[i, 2] <- logdensity.unif(y[i], -maxrt[subject[i]], maxrt[subject[i]])

        # Select one of these two densities (Mixture of nonlapse and lapse trials)
        density[i] <- exp(ld_comp[i, component_chosen[i]] - Constant)
        
        # Generate a likelihood for the MCMC sampler using a trick to maximize density value
        Ones[i] ~ dbern(density[i])

        # Probability of mind wandering trials (lapse trials)
        component_chosen[i] ~ dcat(probsub[subject[i],1:2])
    }
}
'''

modelfile = 'simultrialparams2.jags'
f = open(modelfile, 'w')
f.write(tojags)
f.close()

# Track these variables
trackvars = ['tersubsd', 'deltasubsd', 'alphasubsd',
             'tercond', 'deltacond', 'alphacond',
             'tersub', 'deltasub', 'alphasub',
             'probsub']


for n in range(0, nsims):
    #Initialize vectors
    N = np.sum(genparam['rt'][n, :, :] > .01) #Use 1 ms cutoff
    subject = np.zeros(N)
    rt = np.zeros(N)
    acc = np.zeros(N)
    y = np.zeros(N)
    minrt = np.zeros(nsubs)
    maxrt = np.zeros(nsubs)
    Ones = np.ones(N)
    indextrack = 0
    for s in range(1, nsubs + 1):
        # Use only very small cutoff when fitting a mixture model
        # "Invalid parent" error occurs without this cutoff
        whereindx = np.where(genparam['rt'][n, s - 1, :] > .01) #Use 1 ms cutoff
        subn = whereindx[0].shape[0]
        subject[indextrack:(indextrack+subn)] = np.ones(subn) * s
        rt[indextrack:(indextrack+subn)] = genparam['rt'][n, s - 1, whereindx]
        acc[indextrack:(indextrack+subn)] = genparam['acc'][n, s - 1, whereindx]
        y[indextrack:(indextrack+subn)] = genparam['y'][n, s - 1, whereindx]
        minrt[s - 1] = np.min(rt[indextrack:(indextrack+subn)])
        maxrt[s - 1] = np.max(rt[indextrack:(indextrack+subn)])
        indextrack += subn
    # Create dictionary of initial values
    initials = []
    for c in range(0, nchains):
        chaininit = {
            'alphasub': np.random.uniform(.5, 2., size=(nsubs)),
            'deltasub': np.random.uniform(-4., 4., size=(nsubs)),
            'tersub': np.zeros((nsubs)),
            'alphacond': np.random.uniform(.5, 2.),
            'deltacond': np.random.uniform(-4., 4.),
            'tercond': np.random.uniform(0., .3),
            'alphasubsd': np.random.uniform(.01, 1.),
            'deltasubsd': np.random.uniform(.01, 3.),
            'tersubsd': np.random.uniform(.01, .1)
        }
        for j in range(0, nsubs):
            chaininit['tersub'][j] = np.random.uniform(0., minrt[j]/2)
        initials.append(chaininit)
    print 'Fitting model %i ...' % n
    threaded = pyjags.Model(file=modelfile, init=initials,
                            data=dict(y=y, N=N, nsubs=nsubs, maxrt=maxrt,
                                      subject=subject, Ones=Ones, Constant=10),
                            chains=nchains, adapt=burnin, threads=6,
                            progress_bar=True)
    samples = threaded.sample(nsamps, vars=trackvars, thin=10)
    savestring = ('modelfits/trialparam5_test_model%i.mat') % (n + 1)
    print 'Saving results to: \n %s' % savestring
    sio.savemat(savestring, samples)