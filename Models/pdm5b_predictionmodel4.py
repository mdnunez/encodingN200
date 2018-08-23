# pdm5b_predictionmodel4.py - Finds and evaluates posterior predictive distributions
#                           for RTs and accuracies
#
# Copyright (C) 2018 Michael D. Nunez, <mdnunez1@uci.edu>
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
# 08/08/18      Michael Nunez                               Original code
#
#

#To do: Add Eta parameter to simulations

# Modules
from __future__ import division
import numpy as np
import pyjags
import scipy.io as sio
from scipy.stats import truncnorm
import warnings
import random

# Simulate diffusion models
def simuldiff(N=100,Alpha=1,Ter=.35,Nu=1,Zeta=None,rangeTer=0,rangeZeta=0,Eta=.3,Varsigma=1):
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

    Ter: the mean non-decision time across trials in seconds
    (defaults to .35 seconds)

    Nu: the mean drift rate across trials in evidence units per second
    (defaults to 1 evidence units per second, restricted to -5 to 5 units)

    Zeta: the initial bias in the evidence process for choice A
    (defaults to 50% of total evidence units given by Alpha)

    rangeTer: Non-decision time across trials is generated from a uniform
    distribution of Ter - rangeTer/2 to  Ter + rangeTer/2 across trials
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
    Numpy vector with reaction times (in seconds) multiplied by the response vector
    such that negative reaction times encode response B and positive reaction times
    encode response A  
    
    
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

    #Called sigma in 2001 paper
    D = np.power(Varsigma,2)/2

    #Program specifications
    eps = 2.220446049250313e-16 #precision from 1.0 to next double-precision number
    delta=eps

    for j in range(0,N):
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
            ndrt = Ter - rangeTer/2 + rangeTer*np.random.uniform()
            if ( (dir_ + delta) > Aupper):
                T[j]=totaltime+ndrt
                XX[j]=1
                finish=1
            elif ( (dir_-delta) < Alower ):
                T[j]=totaltime+ndrt
                XX[j]=-1
                finish=1
            else:
                startpos=dir_
                radius=np.min(np.abs([Aupper, Alower]-startpos))

    result = T*XX
    return result


#Combine chains of samples from AJGS
def combinechains(possamps):

    # Number of chains
    nchains = possamps.shape[-1]

    # Number of samples per chain
    nsamps = possamps.shape[-2]

    # Number of dimensions
    ndims = possamps.ndim - 2

    reshapeit = ()
    for i in possamps.shape[0:ndims]:
        reshapeit += (i,)
    reshapeit += (nchains*nsamps,)

    # Reshape data
    reshaped = np.reshape(possamps, reshapeit)

    return reshaped


#Percentage of variance explained by prediction
def rsquaredpred(datastats,simstats):
    
    datastats = datastats.flatten()
    simstats = simstats.flatten()
    sspred = np.nansum(np.power(datastats-simstats,2))/np.sum(np.isfinite(datastats))
    ssinsamp = np.nanvar(datastats,ddof=1)
    varexplained = 1 - sspred/ssinsamp

    return varexplained


# Set random seed
random.seed(12)

# Set up training data
traindata = np.genfromtxt('../Data/N200_rt_window_150_275.csv', delimiter=',')
y = (traindata[:, 2]/1000. ) * (traindata[:, 3] * 2 - 1)
condition = np.array(traindata[: , 4], 'int')
condition += 1
experiment = np.array(traindata[:, 6], 'int')
nconds = np.unique(condition).shape[0]
_, sessioncount = np.unique(traindata[:, 5], return_inverse=True)  # Unique EEG session index
sessioncount = sessioncount + 1
_, expindex = np.unique(traindata[:, 5], return_index=True)
experiment = experiment[expindex] - 1 # Index experiment by session and convert to 0 and 1
nses = np.unique(sessioncount).shape[0]
N = y.shape[0]
n200lat = traindata[:, 0]/1000. #convert from ms to seconds
maxrt = np.empty((nconds, nses))
for k in range(0, nconds):
    for j in range(0, nses):
        where = (sessioncount == (j + 1)) & (condition == (k + 1))
        maxrt[k, j] = np.max(np.abs(y[where]))


# Set up test data
testdata = np.genfromtxt('../Data/TEST_N200_rt_window_150_275.csv', delimiter=',')
y_test = (testdata[:, 2]/1000. ) * (testdata[:, 3] * 2 - 1)
condition_test = np.array(testdata[: , 4], 'int')
condition_test += 1
experiment_test = np.array(testdata[:, 6], 'int')
nconds_test = np.unique(condition_test).shape[0]
_, sessioncount_test = np.unique(testdata[:, 5], return_inverse=True)  # Unique EEG session index
sessioncount_test = sessioncount_test + 1
_, expindex_test = np.unique(testdata[:, 5], return_index=True)
experiment_test = experiment_test[expindex_test] - 1 # Index experiment by session and convert to 0 and 1
nses_test = np.unique(sessioncount_test).shape[0]
N_test = y_test.shape[0]
n200lat_test = testdata[:, 0]/1000. #convert from ms to seconds

#Load model samples
jagsmodel = 'jagsmodel_all_n1lat_random_lapseJul_11_18_10_26.mat'

samples = sio.loadmat(jagsmodel)

# Get posterior distributions
tersub_samps = combinechains(samples['tersub'])
alphasub_samps = combinechains(samples['alphasub'])
deltasub_samps = combinechains(samples['deltasub'])
probsub_samps = combinechains(samples['probsub'])


# Generate posterior predictive distributions
npossamps = np.shape(tersub_samps)[-1]
ndraws = 1000
nses = 49
nconds = 3
simulated_trainn200_y = np.zeros((nconds,nses,ndraws,120))
simulated_testn200_y = np.zeros((nconds,nses,ndraws,120))
draw = np.random.choice(npossamps,ndraws)
drawcount = 0
for d in draw:
    print("Simulation draw %d" % drawcount)
    for s in range(0,nses):
        for c in range(0,nconds):
            traintrials = np.where( (condition==(c+1)) & (sessioncount==(s+1)) )[0]
            testtrials = np.where( (condition_test==(c+1)) & (sessioncount_test==(s+1)) )[0]
            #####
            tempdiffsamps = simuldiff(N=np.size(traintrials),Ter=tersub_samps[c,s,d],Alpha=alphasub_samps[c,s,d],Nu=deltasub_samps[c,s,d],Eta=0)
            estimate_lapse_proportion = probsub_samps[c,s,1,d] #Estimated proportion of lapse trials
            lapsetrials = np.random.choice(np.size(traintrials),np.round(np.size(traintrials)*estimate_lapse_proportion),replace=False)
            tempdiffsamps[lapsetrials] = np.random.uniform(-maxrt[c,s],maxrt[c,s],size=np.size(lapsetrials)) # Lapse process
            simulated_trainn200_y[c,s,drawcount,0:np.size(traintrials)] = tempdiffsamps
            simulated_trainn200_y[c,s,drawcount,np.size(traintrials):] = np.nan
            #####
            tempdiffsamps_test = simuldiff(N=np.size(testtrials),Ter=tersub_samps[c,s,d],Alpha=alphasub_samps[c,s,d],Nu=deltasub_samps[c,s,d],Eta=0)
            lapsetrials_test = np.random.choice(np.size(testtrials),np.round(np.size(testtrials)*estimate_lapse_proportion),replace=False)
            tempdiffsamps_test[lapsetrials_test] = np.random.uniform(-maxrt[c,s],maxrt[c,s],size=np.size(lapsetrials_test)) # Lapse process
            simulated_testn200_y[c,s,drawcount,0:np.size(testtrials)] = tempdiffsamps_test
            simulated_testn200_y[c,s,drawcount,np.size(testtrials):] = np.nan
                    
    drawcount += 1


# Calculate in-sample and out-of-sample prediction metrics
train_cRT_25percentiles = np.zeros((nconds,nses))
train_cRT_medians = np.zeros((nconds,nses))
train_cRT_75percentiles = np.zeros((nconds,nses))
train_iRT_25percentiles = np.zeros((nconds,nses))
train_iRT_medians = np.zeros((nconds,nses))
train_iRT_75percentiles = np.zeros((nconds,nses))
train_accuracies = np.zeros((nconds,nses))
simtrain_cRT_25percentiles = np.zeros((nconds,nses))
simtrain_cRT_medians = np.zeros((nconds,nses))
simtrain_cRT_75percentiles = np.zeros((nconds,nses))
simtrain_iRT_25percentiles = np.zeros((nconds,nses))
simtrain_iRT_medians = np.zeros((nconds,nses))
simtrain_iRT_75percentiles = np.zeros((nconds,nses))
simtrain_accuracies = np.zeros((nconds,nses))
test_cRT_25percentiles = np.zeros((nconds,nses))
test_cRT_medians = np.zeros((nconds,nses))
test_cRT_75percentiles = np.zeros((nconds,nses))
test_iRT_25percentiles = np.zeros((nconds,nses))
test_iRT_medians = np.zeros((nconds,nses))
test_iRT_75percentiles = np.zeros((nconds,nses))
test_accuracies = np.zeros((nconds,nses))
simtest_cRT_25percentiles = np.zeros((nconds,nses))
simtest_cRT_medians = np.zeros((nconds,nses))
simtest_cRT_75percentiles = np.zeros((nconds,nses))
simtest_iRT_25percentiles = np.zeros((nconds,nses))
simtest_iRT_medians = np.zeros((nconds,nses))
simtest_iRT_75percentiles = np.zeros((nconds,nses))
simtest_accuracies = np.zeros((nconds,nses))
for s in range(0,nses):
    for c in range(0,nconds):
        train_cRT_25percentiles[c,s] = np.nanpercentile(np.abs(y[(condition==(c+1)) & (sessioncount==(s+1)) & (y > 0)]),25)
        train_cRT_medians[c,s] = np.nanpercentile(np.abs(y[(condition==(c+1)) & (sessioncount==(s+1)) & (y > 0)]),50)
        train_cRT_75percentiles[c,s] = np.nanpercentile(np.abs(y[(condition==(c+1)) & (sessioncount==(s+1)) & (y > 0)]),75)
        train_iRT_25percentiles[c,s] = np.nanpercentile(np.abs(y[(condition==(c+1)) & (sessioncount==(s+1)) & (y < 0)]),25)
        train_iRT_medians[c,s] = np.nanpercentile(np.abs(y[(condition==(c+1)) & (sessioncount==(s+1)) & (y < 0)]),50)
        train_iRT_75percentiles[c,s] = np.nanpercentile(np.abs(y[(condition==(c+1)) & (sessioncount==(s+1)) & (y < 0)]),75)
        train_accuracies[c,s] = np.nanmean((np.sign(y[(condition==(c+1)) & (sessioncount==(s+1))]) + 1)/2.)
        simtrain_cRT_25percentiles[c,s] = np.nanpercentile(np.abs(simulated_trainn200_y[c,s,(simulated_trainn200_y[c,s,:])>0]),25)
        simtrain_cRT_medians[c,s] = np.nanpercentile(np.abs(simulated_trainn200_y[c,s,(simulated_trainn200_y[c,s,:])>0]),50)
        simtrain_cRT_75percentiles[c,s] = np.nanpercentile(np.abs(simulated_trainn200_y[c,s,(simulated_trainn200_y[c,s,:])>0]),75)
        simtrain_iRT_25percentiles[c,s] = np.nanpercentile(np.abs(simulated_trainn200_y[c,s,(simulated_trainn200_y[c,s,:])<0]),25)
        simtrain_iRT_medians[c,s] = np.nanpercentile(np.abs(simulated_trainn200_y[c,s,(simulated_trainn200_y[c,s,:])<0]),50)
        simtrain_iRT_75percentiles[c,s] = np.nanpercentile(np.abs(simulated_trainn200_y[c,s,(simulated_trainn200_y[c,s,:])<0]),75)
        simtrain_accuracies[c,s] = np.nanmean((np.sign(simulated_testn200_y[c,s,:]) + 1)/2.)
        test_cRT_25percentiles[c,s] = np.nanpercentile(np.abs(y_test[(condition_test==(c+1)) & (sessioncount_test==(s+1)) & (y_test > 0)]),25)
        test_cRT_medians[c,s] = np.nanpercentile(np.abs(y_test[(condition_test==(c+1)) & (sessioncount_test==(s+1)) & (y_test > 0)]),50)
        test_cRT_75percentiles[c,s] = np.nanpercentile(np.abs(y_test[(condition_test==(c+1)) & (sessioncount_test==(s+1)) & (y_test > 0)]),75)
        test_iRT_25percentiles[c,s] = np.nanpercentile(np.abs(y_test[(condition_test==(c+1)) & (sessioncount_test==(s+1)) & (y_test < 0)]),25)
        test_iRT_medians[c,s] = np.nanpercentile(np.abs(y_test[(condition_test==(c+1)) & (sessioncount_test==(s+1)) & (y_test < 0)]),50)
        test_iRT_75percentiles[c,s] = np.nanpercentile(np.abs(y_test[(condition_test==(c+1)) & (sessioncount_test==(s+1)) & (y_test < 0)]),75)
        test_accuracies[c,s] = np.nanmean((np.sign(y_test[(condition_test==(c+1)) & (sessioncount_test==(s+1))]) + 1)/2.)
        simtest_cRT_25percentiles[c,s] = np.nanpercentile(np.abs(simulated_testn200_y[c,s,(simulated_testn200_y[c,s,:])>0]),25)
        simtest_cRT_medians[c,s] = np.nanpercentile(np.abs(simulated_testn200_y[c,s,(simulated_testn200_y[c,s,:])>0]),50)
        simtest_cRT_75percentiles[c,s] = np.nanpercentile(np.abs(simulated_testn200_y[c,s,(simulated_testn200_y[c,s,:])>0]),75)
        simtest_iRT_25percentiles[c,s] = np.nanpercentile(np.abs(simulated_testn200_y[c,s,(simulated_testn200_y[c,s,:])<0]),25)
        simtest_iRT_medians[c,s] = np.nanpercentile(np.abs(simulated_testn200_y[c,s,(simulated_testn200_y[c,s,:])<0]),50)
        simtest_iRT_75percentiles[c,s] = np.nanpercentile(np.abs(simulated_testn200_y[c,s,(simulated_testn200_y[c,s,:])<0]),75)
        simtest_accuracies[c,s] = np.nanmean((np.sign(simulated_testn200_y[c,s,:]) + 1)/2.)



insample_cRT_25per = rsquaredpred(train_cRT_25percentiles,simtrain_cRT_25percentiles)
print "The in-sample prediction of 25th correct RT percentiles explains %.3f %% of the variance" % (insample_cRT_25per*100)

insample_cRT_50per = rsquaredpred(train_cRT_medians,simtrain_cRT_medians)
print "The in-sample prediction of correct RT medians explains %.3f %% of the variance" % (insample_cRT_50per*100)

insample_cRT_75per = rsquaredpred(train_cRT_75percentiles,simtrain_cRT_75percentiles)
print "The in-sample prediction of 75th correct RT percentiles explains %.3f %% of the variance" % (insample_cRT_75per*100)

insample_iRT_25per = rsquaredpred(train_iRT_25percentiles,simtrain_iRT_25percentiles)
print "The in-sample prediction of 25th incorrect RT percentiles explains %.3f %% of the variance" % (insample_iRT_25per*100)

insample_iRT_50per = rsquaredpred(train_iRT_medians,simtrain_iRT_medians)
print "The in-sample prediction of incorrect RT medians explains %.3f %% of the variance" % (insample_iRT_50per*100)

insample_iRT_75per = rsquaredpred(train_iRT_75percentiles,simtrain_iRT_75percentiles)
print "The in-sample prediction of 75th incorrect RT percentiles explains %.3f %% of the variance" % (insample_iRT_75per*100)

insample_accuracies = rsquaredpred(train_accuracies,simtrain_accuracies)
print "The in-sample prediction of EEG session accuracies explains %.3f %% of the variance" % (insample_accuracies*100)

outsample_cRT_25per = rsquaredpred(test_cRT_25percentiles,simtest_cRT_25percentiles)
print "The out-of-sample prediction of 25th correct RT percentiles explains %.3f %% of the variance" % (outsample_cRT_25per*100)

outsample_cRT_50per = rsquaredpred(test_cRT_medians,simtest_cRT_medians)
print "The out-of-sample prediction of correct RT medians explains %.3f %% of the variance" % (outsample_cRT_50per*100)

outsample_cRT_75per = rsquaredpred(test_cRT_75percentiles,simtest_cRT_75percentiles)
print "The out-of-sample prediction of 75th correct RT percentiles explains %.3f %% of the variance" % (outsample_cRT_75per*100)

outsample_iRT_25per = rsquaredpred(test_iRT_25percentiles,simtest_iRT_25percentiles)
print "The out-of-sample prediction of 25th incorrect RT percentiles explains %.3f %% of the variance" % (outsample_iRT_25per*100)

outsample_iRT_50per = rsquaredpred(test_iRT_medians,simtest_iRT_medians)
print "The out-of-sample prediction of incorrect RT medians explains %.3f %% of the variance" % (outsample_iRT_50per*100)

outsample_iRT_75per = rsquaredpred(test_iRT_75percentiles,simtest_iRT_75percentiles)
print "The out-of-sample prediction of 75th incorrect RT percentiles explains %.3f %% of the variance" % (outsample_iRT_75per*100)

outsample_accuracies = rsquaredpred(test_accuracies,simtest_accuracies)
print "The out-of-sample prediction of EEG session accuracies explains %.3f %% of the variance" % (outsample_accuracies*100)