# pdm5b_plotN200.py - Creates data matrix for statistical models
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
# 01/03/18     Michael Nunez               Coverted from pdm5b_prepinputs.py
# 01/05/18     Michael Nunez                  Print summary statistics
# 02/27/18     Michael Nunez               Change fontsize, print 10th percentiles
# 04/13/18     Michael Nunez                Change figure size

# Imports
from __future__ import division
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from IPython import get_ipython  # Run magic functions from script
# Initialize ipython matplotlib plotting graphics
get_ipython().magic('pylab')

# Parameters

fontsize = 16
teal = "#009E73"
blue = "#0072B2"
gold = "#E69F00"

# Initial
loadloc = '/home/michael/data10/michael/pdm/exp5data/jagsin/behav_strint7'


# Load all data
print 'Loading all N200 data...'
dataout = np.load(loadloc + '.npz')

# Set up EEG session data
print 'Loading all tidy data...'
trialdata = np.genfromtxt(
    '../Data/N200_rt_window_150_275_fixed350cutoff.csv', delimiter=',')
sesdata = np.genfromtxt(
    '../Data/N1deflec2_cutoffs_allSNR_window_150_275_fixed350cutoff.csv', delimiter=',')

n1lat = sesdata[:, 0]
n1deflec = sesdata[:, 1]
rt1per = sesdata[:, 7]
rt10per = sesdata[:, 8]
rtmeanper = sesdata[:, 13]
rtmedian = sesdata[:, 11]
accmean = sesdata[:, 17]

n200lat = trialdata[:, 0]
rt = trialdata[:, 2]

# Plot all N200 waveforms
print 'Plotting the N200 data'
plt.figure(figsize=(12, 6))
plt.plot(np.arange(-100, 275), dataout['n1data'])
zeroline = plt.plot(np.array([0., 0.]), np.array([-.15, .15]))
plt.setp(zeroline, linewidth=3, linestyle='--')
plt.xlim(-100, 275)
plt.ylim(-.15, .15)

# Plot N200 deflection and peak latency distributions
kde_deflec = stats.gaussian_kde(n1deflec)
x_deflec = np.linspace(0, 275, 100)
p_deflec = kde_deflec(x_deflec)
p_deflec = .15 * p_deflec / np.max(p_deflec)
# deflec_dist = plt.plot(x_deflec, p_deflec)
# plt.setp(deflec_dist, color=teal, linewidth=3)
plt.fill_between(x_deflec, p_deflec, np.zeros((100)), color=teal)

kde_lat = stats.gaussian_kde(n1lat)
x_lat = np.linspace(150, 275, 100)
p_lat = kde_lat(x_lat)
p_lat = .15 * p_lat / np.max(p_lat)
# lat_dist = plt.plot(x_lat, p_lat)
# plt.setp(lat_dist, color=blue, linewidth=3)
plt.fill_between(x_lat, p_lat, np.zeros((100)), color=blue)

# Format plot
frame = plt.gca()
frame.axes.get_yaxis().set_visible(False)
plt.tick_params(axis='both', which='major', labelsize=fontsize)
plt.xlabel('Time (ms) after Gabor onset', fontsize=fontsize)

# Save the figure
print 'Saving the figure...'
plt.savefig('N200_waveforms.png', dpi=300, format='png')

# Generate values for the maunscript text
print n1lat.shape[0]
print np.mean(n1lat)
print np.median(n1lat)
print np.percentile(n1lat, [5, 50, 95])
print np.mean(n1deflec)
print np.median(n1deflec)
print np.percentile(n1deflec, [5, 50, 95])
print stats.pearsonr(n1lat, n1deflec)
print np.mean(rt10per)
print np.median(rt10per)
print np.percentile(rt10per, [5, 50, 95])
print stats.pearsonr(n1lat, rt10per)
print stats.pearsonr(n1deflec, rt10per)
print np.mean(rtmedian)
print np.median(rtmedian)
print np.percentile(rtmedian, [5, 50, 95])
print np.mean(accmean)
print np.median(accmean)
print np.percentile(accmean, [5, 50, 95])

sixconds = sesdata[:, 3] + (sesdata[:, 5] - 1) * 3
for cond in np.unique(sixconds):
    print 'Condition %d:' % cond
    print np.mean(n1lat[sixconds == cond])
    print np.mean(n1deflec[sixconds == cond])
    print np.mean(rt10per[sixconds == cond])
    print np.mean(rtmedian[sixconds == cond])
    print np.mean(accmean[sixconds == cond])

truesubject = sesdata[:, 16]
subn1lat = np.zeros((np.unique(truesubject).shape[0]))
subn1deflec = np.zeros((np.unique(truesubject).shape[0]))
subrt10 = np.zeros((np.unique(truesubject).shape[0]))
subrtmed = np.zeros((np.unique(truesubject).shape[0]))
subacc = np.zeros((np.unique(truesubject).shape[0]))
for truesub in np.int8(np.unique(truesubject)):
    subn1lat[truesub - 1] = np.mean(n1lat[truesubject == truesub])
    subn1deflec[truesub - 1] = np.mean(n1deflec[truesubject == truesub])
    subrt10[truesub - 1] = np.mean(rt10per[truesubject == truesub])
    subrtmed[truesub - 1] = np.mean(rtmedian[truesubject == truesub])
    subacc[truesub - 1] = np.mean(accmean[truesubject == truesub])

print np.min(subn1lat)
print np.mean(subn1lat)
print np.max(subn1lat)
print np.min(subn1deflec)
print np.mean(subn1deflec)
print np.max(subn1deflec)
print np.min(subrt10)
print np.mean(subrt10)
print np.max(subrt10)
print np.min(subrtmed)
print np.mean(subrtmed)
print np.max(subrtmed)
print np.min(subacc)
print np.mean(subacc)
print np.max(subacc)
print stats.pearsonr(subn1lat, subn1deflec)
print stats.pearsonr(subn1lat, subrt10)
print stats.pearsonr(subn1deflec, subrt10)
print stats.pearsonr(subrtmed, subacc)
