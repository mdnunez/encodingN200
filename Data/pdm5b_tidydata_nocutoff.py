# pdm5b_tidydata_nocutoff.py - Creates tidy data csv file
#                              without 350 ms cutoff
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
# 06/13/18      Michael Nunez               Converted from pdm5b_tidydata.py


# Imports
import numpy as np
import pandas as pd
import scipy.io as sio

# Initial
dataloc = '/home/michael/data10/michael/pdm/exp5data/jagsin/behav_strint7.npz'

# Code
data = np.load(dataloc)
n200 = data['n200']  # Single trial amplitude
n200lat = data['n200lat']  # Single trial latency
n1amp = data['n1']  # Trial-averaged amplitude
n1lat = data['n1lat']  # Trial-averaged latency
n1deflec = data['n1deflec']  # Trial-averaged deflection
conds = data['condition']  # Condition index
truesubs = data['subject']  # True subject index
exps = data['experiment']  # Experiment (visual noise) type index
ses = data['session']  # Run session index
tempses = (exps - 1) * 2 + ses
tempsubs = truesubs + .1 * tempses
_, EEGses = np.unique(tempsubs, return_inverse=True)  # Unique EEG session index


# Remove repeating data
_, alln1indx = np.unique(data['n1'], return_index=True)
n1indx = alln1indx[data['n1lat'][alln1indx] >= 100]


# Calculate reaction time percentiles
rt1per = np.zeros(n1indx.shape[0])
rt10per = np.zeros(n1indx.shape[0])
rt30per = np.zeros(n1indx.shape[0])
rt50per = np.zeros(n1indx.shape[0])
rt70per = np.zeros(n1indx.shape[0])
rt90per = np.zeros(n1indx.shape[0])
rtmean = np.zeros(n1indx.shape[0])
accmean = np.zeros(n1indx.shape[0])


cutoffs = np.empty(n1indx.shape[0])
perrmed = np.empty(n1indx.shape[0])
for n in range(0, n1indx.shape[0]):
    thosetrials = np.intersect1d(
        np.where(data['train']), np.where(data['condition'] == conds[n1indx][n]))
    thosetrials = np.intersect1d(
        thosetrials, np.where(data['subject'] == truesubs[n1indx][n]))
    thosetrials = np.intersect1d(
        thosetrials, np.where(data['session'] == ses[n1indx][n]))
    thosetrials = np.intersect1d(
        thosetrials, np.where(data['experiment'] == exps[n1indx][n]))
    thosetrials = np.intersect1d(
        thosetrials, np.where(np.isfinite(data['rt'])))
    thosetrials = np.intersect1d(thosetrials, np.where(data['rt'] > 0))
    # Use fixed cutoff
    cutoffs[n] = 0
    thesetrials = np.intersect1d(
        thosetrials, np.where(data['rt'] > cutoffs[n]))
    perrmed[n] = (1 - thesetrials.shape[0] /
                  float(thosetrials.shape[0])) * 100  # Percent removed
    rt10per[n] = np.percentile(data['rt'][thesetrials], 10)
    rt1per[n] = np.percentile(data['rt'][thesetrials], 1)
    rt30per[n] = np.percentile(data['rt'][thesetrials], 30)
    rt50per[n] = np.percentile(data['rt'][thesetrials], 50)
    rt70per[n] = np.percentile(data['rt'][thesetrials], 70)
    rt90per[n] = np.percentile(data['rt'][thesetrials], 90)
    rtmean[n] = np.nanmean(data['rt'][thesetrials])
    accmean[n] = np.nanmean(data['correct'][thesetrials])


N1rtdata = np.vstack((n1lat[n1indx], n1deflec[n1indx], n1amp[n1indx], conds[n1indx], EEGses[n1indx], exps[n1indx], ses[n1indx],
                      rt1per, rt10per, rt30per, rt50per, rt70per, rt90per,
                      rtmean, cutoffs, perrmed, truesubs[n1indx], accmean)).T


# Save out EEG session data
# See JASP for GUI analysis program: https://jasp-stats.org/

np.savetxt(
    'N1deflec2_allSNR_window_150_275.csv',    # file name
    N1rtdata,               # array to save
    fmt='%.4f',             # formatting, 4 digits in this case
    delimiter=',',          # column delimiter
    newline='\n',           # new line character
    header='N1 latencies, N1 deflection, N1 Magnitude, SNR condition, EEG Session Counter, Experiment, Session, 1st RT percentiles, 10th RT percentiles, 30th RT percentiles, 50th RT percentiles, 70th RT percentiles, 90th RT percentiles, RT means, Cutoffs, Percent Removed, True Subject Index, Accuracy'
)      # file header

# Save out single-trial data
alltrials = np.intersect1d(np.where(np.logical_not(
    data['missing'])), np.where(data['train']))

n200lat = data['n200lat'][alltrials]
rt = data['rt'][alltrials]
acc = data['correct'][alltrials]
n200 = data['n200'][alltrials]
conds = data['condition'][alltrials]
truesubs = data['subject'][alltrials]
exps = data['experiment'][alltrials]
ses = data['session'][alltrials]
EEGses = EEGses[alltrials]

# Remove bad data
keeptrials = (n200lat != np.min(n200lat)) & (n200lat != np.max(n200lat)) & np.isfinite(n200lat) & np.isfinite(rt) & (rt > 0)
n200lat = n200lat[keeptrials]
rt = rt[keeptrials]
acc = acc[keeptrials]
n200 = n200[keeptrials]
conds = conds[keeptrials]
truesubs = truesubs[keeptrials]
exps = exps[keeptrials]
ses = ses[keeptrials]
EEGses = EEGses[keeptrials]


# Export data to JASP
N200rtdata = np.vstack((n200lat, n200, rt, acc, conds, EEGses, exps, ses, truesubs)).T

np.savetxt(
    'N200_rt_window_150_275.csv',    # file name
    N200rtdata,               # array to save
    fmt='%.4f',             # formatting, 4 digits in this case
    delimiter=',',          # column delimiter
    newline='\n',           # new line character
    header='Single-trial N200 latencies, N200 amplitudes, RT, Accuracy, SNR condition, EEG Session Counter, Experiment, Session, True Subject Index'
)      # file header

