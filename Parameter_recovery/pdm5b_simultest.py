# pdm5b_simultest.py - Testing JAGS fits of Hierarchical Diffusion Models with
#         trial-to-trial variability in non-decision time and drift
#
# Copyright (C) 2017 Michael D. Nunez, <mdnunez1@uci.edu>
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
# 12/10/17      Michael Nunez              Converted from hgs_ratcliffhdm.py
# 12/14/17      Michael Nunez         Add mind wandering process, use 350 ms cutoff
# 12/21/17      Michael Nunez         Systematically very the mind wandering process
# 12/28/17      Michael Nunez        Increase the number of simulations
# Notes
# To install Pyjags from github on Anaconda:
# which pip [make sure 'pip' is connected to Anaconda]
# git clone --recursive https://github.com/tmiasko/pyjags.git
# pip install ./pyjags
#
# Modules
import numpy as np
import pyjags
from pymatbridge import Matlab
import scipy.io as sio

# Number of simulations
nsims = 30

# Number of simulated subjects
nsubs = 100

# Number of trials per subject
ntrials = 100

# MATLAB code
matlabcode = '''
%%Generate samples from the drift-diffusion model with trial-to-trial drift variability
rng(11);
tersub = rand({0},{1})*.15 + .35; % Uniform .35 to .5
alphasub = rand({0},{1})*.06 + .08; % Uniform from .08 to .14
deltasub = rand({0},{1})*.8 - .4; % Uniform from  -.4 to .4
%tertrialsd = repmat(linspace(0,.1,{0})',1,{1}); % From 0 to .1
tertrialsd = rand({0},{1})*.1; % Uniform from 0 to .1
deltatrialsd = rand({0},{1})*.2; % Uniform from 0 to .2
%prob_mindwander = rand({0},{1})*.10; % Uniform from 0 to 10 percent of trials
prob_mindwander = repmat(linspace(0,.10,{0})',1,{1}); % From 0 to 10 percent of trials
y=zeros({0},{1},{2});
rt=zeros({0},{1},{2});
acc=zeros({0},{1},{2});
for n=1:{0},
	for s=1:{1},
		[tempt, tempx] = simuldiff([alphasub(n, s), tersub(n, s),...
			deltatrialsd(n, s), alphasub(n, s)*.5, 0, tertrialsd(n, s), deltasub(n, s)],{2});
		mindwanderx = randi(2,1,{2}) - 1; 
		mindwandert = rand(1,{2})*2; % Randomly distributed from 0 to 2 seconds

		mindwander_trials = randperm({2},round({2}*prob_mindwander(n, s)));
		tempx(mindwander_trials) = mindwanderx(mindwander_trials);
		tempt(mindwander_trials) = mindwandert(mindwander_trials);

		y(n, s, :) = tempt.*(tempx*2 - 1);
		rt(n, s, :) = tempt;
		acc(n, s, :) = tempx;
	end
end
'''.format(nsims, nsubs, ntrials)

# Write out code
modelfile = 'simultrialparams.m'
f = open(modelfile, 'w')
f.write(matlabcode)
f.close()

# Simulate some data from simuldiff in MATLAB (find a way to do this in Python)
mlab = Matlab(maxtime=240)  # Increase max start time
mlab.start()
results2 = mlab.run_code('simultrialparams;')
genparam = dict()
genparam['tersub'] = mlab.get_variable('tersub')
genparam['alphasub'] = mlab.get_variable('alphasub')
genparam['deltasub'] = mlab.get_variable('deltasub')
genparam['tertrialsd'] = mlab.get_variable('tertrialsd')
genparam['deltatrialsd'] = mlab.get_variable('deltatrialsd')
genparam['prob_mindwander'] = mlab.get_variable('prob_mindwander')
genparam['rt'] = mlab.get_variable('rt')
genparam['acc'] = mlab.get_variable('acc')
genparam['y'] = mlab.get_variable('y')
mlab.stop()
sio.savemat('genparam_test.mat', genparam)

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
	}

	##########
	# Wiener likelihoods
	for (i in 1:N) {
		y[i] ~ dwiener(alphasub[subject[i]],
		tersub[subject[i]],
		beta,
		deltasub[subject[i]])
	}
}
'''

modelfile = 'simultrialparams.jags'
f = open(modelfile, 'w')
f.write(tojags)
f.close()

# Track these variables
trackvars = ['tersubsd', 'deltasubsd', 'alphasubsd',
             'tercond', 'deltacond', 'alphacond',
             'tersub', 'deltasub', 'alphasub']


for n in range(0, nsims):
    success = False
    #Initialize vectors
    N = np.sum(genparam['rt'][n, :, :] > .35)
    subject = np.zeros(N)
    rt = np.zeros(N)
    acc = np.zeros(N)
    y = np.zeros(N)
    minrt = np.zeros(nsubs)
    indextrack = 0
    for s in range(1, nsubs + 1):
        whereindx = np.where(genparam['rt'][n, s - 1, :] > .35) #Use 350 ms cutoff
        subn = whereindx[0].shape[0]
        subject[indextrack:(indextrack+subn)] = np.ones(subn) * s
        rt[indextrack:(indextrack+subn)] = genparam['rt'][n, s - 1, whereindx]
        acc[indextrack:(indextrack+subn)] = genparam['acc'][n, s - 1, whereindx]
        y[indextrack:(indextrack+subn)] = genparam['y'][n, s - 1, whereindx]
        minrt[s - 1] = np.min(rt[indextrack:(indextrack+subn)])
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
            chaininit['tersub'][j] = np.random.uniform(0., minrt[j] - .1)
        initials.append(chaininit)
    print 'Fitting model %i ...' % n
    threaded = pyjags.Model(file=modelfile, init=initials,
                            data=dict(y=y, N=N, nsubs=nsubs,
                                      subject=subject),
                            chains=nchains, adapt=burnin, threads=6,
                            progress_bar=True)
    samples = threaded.sample(nsamps, vars=trackvars, thin=10)
    savestring = ('modelfits/trialparam_test_model%i.mat') % (n + 1)
    print 'Saving results to: \n %s' % savestring
    sio.savemat(savestring, samples)
