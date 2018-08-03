# pdm5b_eegbehavmodel7.py - Fits a four parameter generative model with 
#                   trial-by-trial N200 latencies, lapse trials explicitly modeled
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
# 07/19/18      Michael Nunez            Converted from pdm5b_eegbehavmodel6.py
# 08/01/18      Michael Nunez                       Use of better predictive model

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

n200lat = trialdata[:, 0]/1000. #convert from ms to seconds

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

# Track these variables
trackvars = ['rmrsub', 'vetsub', 'alphasub', 'deltasub', 'probsub',
            'rmrcond', 'vetcond', 'deltacond', 'alphacond',
             'rmrsubsd', 'vetsubsd', 'alphasubsd', 'deltasubsd', 'n200trialsd']

#Note: initialize so that all probability is on the lapse trial process, this avoids the "invalid parent" errors from JAGS
initials = []
for c in range(0, nchains):
    chaininit = {
        'rmrsub': np.random.uniform(.1, .3, size=(nconds, nses)),
        'vetsub': np.random.uniform(.16, .273, size=(nconds, nses)),
        'alphasub': np.random.uniform(.5, 2., size=(nconds, nses)),
        'deltasub': np.random.uniform(-4., 4., size=(nconds, nses)),
        'rmrcond': np.random.uniform(.1, .3, size=(2, nconds)),
        'vetcond': np.random.uniform(.16, .273, size=(2, nconds)),
        'alphacond': np.random.uniform(.5, 2., size=(2, nconds)),
        'deltacond': np.random.uniform(-4., 4., size=(2, nconds)),
        'alphasubsd': np.random.uniform(.01, 1.),
        'deltasubsd': np.random.uniform(.01, 3.),
        'rmrsubsd': np.random.uniform(.01, .1),
        'vetsubsd': np.random.uniform(.01, .1),
        'n200trialsd': np.random.uniform(.01, .1),
        'probsub': np.stack((np.zeros((nconds,nses)),np.ones((nconds,nses))),axis=2)
    }
    # for k in range(0, nconds):
    #     for j in range(0, nses):
    #         chaininit['vetsub'][k, j] = np.random.uniform(0., minrt[k, j]/3)
    #         chaininit['rmrsub'][k, j] = np.random.uniform(0., minrt[k, j]/3)
    initials.append(chaininit)

# Run JAGS model

# Choose JAGS model type
modelname = '4parameter_lapse_alt'

thismodel = jagsmodels[modelname]

# Save model
timestart = strftime('%b') + '_' + strftime('%d') + '_' + \
    strftime('%y') + '_' + strftime('%H') + '_' + strftime('%M')
modelfile = 'jagsmodel_' + modelname + timestart + '.jags'
f = open(modelfile, 'w')
f.write(thismodel)
f.close()
print 'Fitting model %s ...' % (modelname + timestart)


indata = data=dict(N=N, y=y, nses=nses, nconds=nconds, maxrt=maxrt, Ones=Ones, Constant=Constant,
                                  EEGsession=sessioncount, condition=condition,
                                  n200lat=n200lat,  experiment=experiment,
                                  nexps=2)

threaded = pyjags.Model(file=modelfile, init=initials,
                        data=indata,
                        chains=nchains, adapt=burnin, threads=6,
                        progress_bar=True)


samples = threaded.sample(nsamps, vars=trackvars, thin=10)

savestring = 'jagsmodel_' + \
    modelname + timestart + ".mat"

print 'Saving results to: \n %s' % savestring

sio.savemat(savestring, samples)