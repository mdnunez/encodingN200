# pdm5b_prepinputs.py - Creates data matrix for statistical models
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
# 11/20/17     Michael Nunez               Coverted from pdm5b_prepinputs2.py
#                                            Modified deflection calculation
# 12/01/17                             Change deflection detection window
# 12/04/17                  Fix finddeflection output, use last found sample
# 06/01/18                                Create index for N200 waveforms
# 01/04/19      Michael Nunez          Export data after lowpass filtering with different parameters
#                                          Revert to original script


# Imports
from __future__ import division
import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt
from IPython import get_ipython  # Run magic functions from script
# Initialize ipython matplotlib plotting graphics
get_ipython().magic('pylab')

# Initial
dataloc = '/data10/michael/pdm/{1}/{0}/{0}_allcleaned.npz'
indxloc = '/data10/michael/pdm/{1}/{0}/{0}_traintestindx.npz'
svdloc = '/data10/michael/pdm/{2}/{0}/{1}/erp_svd_{0}_{1}_v5.mat'
saveloc = '/data10/michael/pdm/exp5data/jagsin/behav_strint7'
# svdloc = '/data10/michael/pdm/{2}/{0}/{1}/erp_svd_{0}_{1}_v6.mat'
# saveloc = '/data10/michael/pdm/exp5data/jagsin/behav_strint8'

# Definitions


def finddeflection(erp):
    '''
    Find where the derivative begins to be negative before the minimum value
    '''
    endsamp = np.shape(erp)[0]
    derivative = erp[1:endsamp] - erp[0:(endsamp - 1)]
    wheremin = np.argmin(erp)
    # Find the derivative up until the minimum
    posdeflection = derivative[:(wheremin)]
    # Find where the derivative is positive
    testtime = np.where(posdeflection > 0)[0]
    if np.size(testtime) > 0:
        onset = testtime[-1]  # Take the last positive derivative value
    else:
        onset = np.nan
    return onset


# Code

dataout = dict()
dataout['experiments'] = ['exp4data/subjects', 'exp5data/subjects/training']
dataout['sessions'] = {'exp4b': ['s1', 's2'], 'exp5b': ['ses1']}
for i in range(2, 8):
    dataout['sessions']['exp5b'].extend(['ses%i' % i])
dataout['subjects'] = ['s59', 's64', 's68', 's80', 's82', 's93',
                       's94', 's95', 's96', 's97', 's100', 's101', 's109',
                       's110']

# Initialize data vectors
fields = ['rt', 'correct', 'randrots', 'block', 'artifact',
          'subject', 'experiment', 'n200lat', 'n1lat', 'n200', 'n1',
          'n1deflec', 'n1deflec_slope']
for f in fields:
    dataout[f] = np.empty((480 * (12 * 2 + 4 * 7)))
dataout['session'] = np.empty((480 * (12 * 2 + 4 * 7)))
dataout['missing'] = np.zeros((480 * (12 * 2 + 4 * 7)))
dataout['train'] = np.ones((480 * (12 * 2 + 4 * 7)))
dataout['condition'] = np.zeros((480 * (12 * 2 + 4 * 7)))
dataout['n1data'] = np.zeros((375, 3 * (12 * 2 + 4 * 7)))
dataout['n1datacond'] = np.zeros((3 * (12 * 2 + 4 * 7)))
dataout['n1dataexp'] = np.zeros((3 * (12 * 2 + 4 * 7)))

# Conditions
snrlabels = ['low', 'med', 'high']

# Extract data
exptrack = 1
startindex = 0
n1dataindx = 0
for exp in dataout['experiments']:
    if exptrack == 1:
        experiment = 'exp4b'
        nses = 2
        posconds = np.array([0.2, 2.0, 20.0])
    else:
        experiment = 'exp5b'
        nses = 7
        posconds = np.array([0.5, 1.0, 2.0])
    subtrack = 1
    for subdes in dataout['subjects']:
        # These missing subjects are already accounted for
        if os.path.isfile(dataloc.format(subdes, exp)):
            print 'Loading data for %s in %s...' % (subdes, exp)
            alldata = np.load(dataloc.format(subdes, exp))
            index = np.load(indxloc.format(subdes, exp))
            alltrials = np.arange(0, 480 * nses) + startindex
            dataout['train'][alltrials[index['test']]] = 0
            # _, uindex = np.unique(alldata['snrvec'], return_inverse=True)
            dataout['condition'][alltrials] += 1 * \
                (alldata['snrvec'] == posconds[1])
            dataout['condition'][alltrials] += 2 * \
                (alldata['snrvec'] == posconds[2])
            dataout['experiment'][alltrials] = exptrack
            dataout['subject'][alltrials] = subtrack
            for f in fields:
                if f in alldata:
                    dataout[f][alltrials] = alldata[f]
            sestrack = 1
            # These missing sessions are not accounted for in the
            # initialization
            for ses in dataout['sessions'][experiment]:
                sestrials = np.arange(0, 480) + startindex
                if os.path.isfile(svdloc.format(subdes, ses, exp)):
                    condtrack = 0
                    svds = sio.loadmat(svdloc.format(subdes, ses, exp))
                    for snr in posconds:
                        thesetrials = np.intersect1d(sestrials,
                                                     np.where(dataout['condition'] == condtrack))
                        # Note that this window exists from -100 to 1000 ms time-locked
                        # to the onset of the response interval
                        # Minimum of PCA component of traditional ERP during
                        # the response interval
                        dataout['n1data'][:, n1dataindx] = svds['svd_rint_%s' % (snrlabels[condtrack])][
                            'u'][0][0][0:375, 0]
                        dataout['n1datacond'][n1dataindx] = condtrack
                        dataout['n1dataexp'][n1dataindx] = exptrack
                        n1dataindx += 1
                        dataout['n1lat'][thesetrials] = np.argmin(
                            svds['svd_rint_%s' % (snrlabels[condtrack])]['u'][0][0][250:375, 0]) + 150
                        deflectdata = svds['svd_rint_%s' % (snrlabels[condtrack])][
                            'u'][0][0][100:375, 0]
                        dataout['n1deflec'][thesetrials] = finddeflection(
                            deflectdata)
                        dataout['n1'][thesetrials] = np.min(
                            svds['svd_rint_%s' % (snrlabels[condtrack])]['u'][0][0][250:375, 0] *
                            svds['svd_rint_%s' % (snrlabels[condtrack])]['s'][0][0][0, 0])
                        # Minimum of single-trial ERP estimates
                        dataout['n200lat'][thesetrials] = np.argmin(
                            svds['st_erp_rint_%s' % (snrlabels[condtrack])][250:375], axis=0) + 150
                        dataout['n200'][thesetrials] = np.min(
                            svds['st_erp_rint_%s' % (snrlabels[condtrack])][250:375], axis=0) / (
                            dataout['n1'][thesetrials])
                        condtrack += 1

                else:
                    n1dataindx += 3
                    dataout['missing'][sestrials] = np.ones(480)
                    dataout['artifact'][sestrials] = np.ones(480)
                    # Generate random conditions for missing sessions
                    dataout['condition'][sestrials] = np.concatenate((
                        np.ones(160) * 0, np.ones(160) * 1, np.ones(160) * 2))
                dataout['session'][sestrials] = sestrack
                sestrack += 1
                startindex += 480
        subtrack += 1
    exptrack += 1

# Save all data
print 'Saving all data as .npz'
np.savez(saveloc + '.npz', **dataout)
print 'Saving all data as .mat'
sio.savemat(saveloc + '.mat', dataout)

# Plot all N1 waveforms
plt.figure()
plt.plot(np.arange(-100, 275), dataout['n1data'])
zeroline = plt.plot(np.array([0., 0.]),np.array([-.15, .10]))
plt.setp(zeroline, linewidth=3, linestyle='--')