# pdm5b_eegbehavmodel4.py - Creates models for training data of all subjects
#                         with subject-level N1 latencies, lapse trials explicitly modeled
#
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
# 06/13/18      Michael Nunez            Converted from pdm5b_eegbehavmodel1.py

# Imports
from __future__ import division
import numpy as np
import numpy.ma as ma
import scipy.io as sio
import pyjags
from scipy import stats
from time import strftime
import os
from pdm5b_papermodels import *

# Set up reaction time data
trialdata = np.genfromtxt('../Data/N200_rt_window_150_275.csv', delimiter=',')

y = (trialdata[:, 2]/1000. ) * (trialdata[:, 3] * 2 - 1)
condition = np.array(trialdata[: , 4], 'int')
condition += 1
experiment = np.array(trialdata[:, 6], 'int')
nconds = np.unique(condition).shape[0]
_, sessioncount = np.unique(trialdata[:, 5], return_inverse=True)  # Unique EEG session index
sessioncount = sessioncount + 1
_, expindex = np.unique(trialdata[:, 5], return_index=True)
experiment = experiment[expindex] - 1 # Index experiment by session and convert to 0 and 1
nses = np.unique(sessioncount).shape[0]
N = y.shape[0]

# Set up N1 latency matrix
sesdata = np.genfromtxt('../Data/N1deflec2_allSNR_window_150_275.csv', delimiter=',')

# Set up N1 latency matrix
n1lat = sesdata[:, 0]
n1onset= sesdata[:, 1]
n1mag = sesdata[:, 2]
conds = np.array(sesdata[:, 3], dtype='int')
conds += 1
_, sescount = np.unique(sesdata[:, 4], return_inverse=True) # Unique EEG session index
sescount = sescount + 1
exp = np.array(sesdata[:, 5], dtype='int')
trueses = np.array(sesdata[:, 6], dtype='int')

n1deflec = np.empty((nconds, nses))
n1sub = np.empty((nconds, nses))
condmat = np.empty((nconds, nses))
submat = np.empty((nconds, nses))
n1sub[:] = np.NAN
n1deflec[:] = np.NAN
condmat[:] = np.NAN
submat[:] = np.NAN
for n in range(0, n1lat.shape[0]):
    n1sub[conds[n] - 1, sescount[n] - 1] = n1lat[n]
    n1deflec[conds[n] - 1, sescount[n] - 1] = n1onset[n]
    condmat[conds[n] - 1, sescount[n] - 1] = conds[n]
    submat[conds[n] - 1, sescount[n] - 1] = sescount[n]
n1deflec /= 1000 #Convert from ms to seconds
n1deflec = ma.masked_invalid(n1deflec)  # Remove NANs for pyjags
n1sub /= 1000 #Convert from ms to seconds
n1sub = ma.masked_invalid(n1sub)  # Remove NANs for pyjags


# Initialize non-decision time with mininum RT for each subject and condition
# Use maximum RT for the bounds on the lapse process, modeled by a uniform distribution
minrt = np.empty((nconds, nses))
maxrt = np.empty((nconds, nses))
for k in range(0, nconds):
    for j in range(0, nses):
        where = (sessioncount == (j + 1)) & (condition == (k + 1))
        minrt[k, j] = np.min(np.abs(y[where]))
        maxrt[k, j] = np.max(np.abs(y[where]))


# Input for mixture modeling
Ones = np.ones(N)
Constant = 10


# pyjags code

# Make sure $LD_LIBRARY_PATH sees /usr/local/lib
pyjags.modules.load_module('wiener')
pyjags.modules.load_module('dic')
pyjags.modules.list_modules()

nchains = 6
burnin = 2000  # Note that scientific notation breaks pyjags
nsamps = 50000

trackvars = ['alphasub', 'deltasub', 'tersub', 'n1sub',
             'n1gammacond', 'n1cond', 'alphacond', 'deltacond', 'tercond',
             'n1subsd', 'alphasubsd', 'deltasubsd', 'tersubsd', 'n1gammault',
             'n1gammasd', 'probsub']

initials = []
for c in range(0, nchains):
    chaininit = {
        'alphasub': np.random.uniform(.5, 2., size=(nconds, nses)),
        'deltasub': np.random.uniform(-4., 4., size=(nconds, nses)),
        'tersub': np.random.uniform(0., .3, size=(nconds, nses)),
        'n1gammacond': np.random.uniform(-2., 2., size=(3, 2, nconds)),
        'n1gammasd': np.random.uniform(.01, 3., size=(1,3)),
        'n1gammault': np.random.uniform(-2., 2., size=(1,3)),
        'n1cond': np.random.uniform(.16, .273, size=(2, nconds)),
        'alphacond': np.random.uniform(.5, 2., size=(2, nconds)),
        'deltacond': np.random.uniform(-4., 4., size=(2, nconds)),
        'tercond': np.random.uniform(0., .3, size=(2, nconds)),
        'n1subsd': np.random.uniform(.01, .1),
        'alphasubsd': np.random.uniform(.01, 1.),
        'deltasubsd': np.random.uniform(.01, 3.),
        'tersubsd': np.random.uniform(.01, .1)
    }
    for k in range(0, nconds):
        for j in range(0, nses):
            chaininit['tersub'][k, j] = np.random.uniform(0., minrt[k, j]/2)
    initials.append(chaininit)

# Run JAGS model

# Choose JAGS model type
modelname = 'all_n1lat_random_lapse'

thismodel = jagsmodels[modelname]

# Save model
timestart = strftime('%b') + '_' + strftime('%d') + '_' + \
    strftime('%y') + '_' + strftime('%H') + '_' + strftime('%M')
modelfile = 'jagsmodel_' + modelname + timestart + '.jags'
f = open(modelfile, 'w')
f.write(thismodel)
f.close()
print 'Fitting model %s ...' % (modelname + timestart)


indata = dict(N=N, y=y, nses=nses, nconds=nconds, maxrt=maxrt, Ones=Ones, Constant=Constant,
                                  EEGsession=sessioncount, condition=condition,
                                  n1sub=n1sub, experiment=experiment, nexps=2)

threaded = pyjags.Model(file=modelfile, init=initials,
                        data=indata,
                        chains=nchains, adapt=burnin, threads=6,
                        progress_bar=True)


samples = threaded.sample(nsamps, vars=trackvars, thin=10)

savestring = 'jagsmodel_' + \
    modelname + timestart + ".mat"

print 'Saving results to: \n %s' % savestring

sio.savemat(savestring, samples)