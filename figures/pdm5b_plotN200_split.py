# pdm5b_plotN200_split.py - Plots N200 waveforms split by experiment and condition
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
# 06/01/18     Michael Nunez               Coverted from pdm5b_plotN200.py
# 06/15/18     Michael Nunez                      Data from both experiments
# 06/18/18     Michael Nunez                    Plot deflection distributions
# 08/02/18     Michael Nunez                        Use of data without cutoffs
# 01/08/19     Michael Nunez          Add stars to differentiate experiment 1 vs 2


## References:
# http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=6
# https://matplotlib.org/gallery/text_labels_and_annotations/custom_legends.html

# Imports
from __future__ import division
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from IPython import get_ipython  # Run magic functions from script
# Initialize ipython matplotlib plotting graphics
get_ipython().magic('pylab')

# Parameters

fontsize = 16
color1 = "#D95F02"
color2 = "#1B9E77"
color3 = "#7570B3"
color4 = "#E6AB02"
color5 = "#66A61E"
color6 = "#E7298A"
color7 = "#A6761D"
colororder = [(color1, color2, color3),(color4, color5, color6)]

# Initial
loadloc = '/home/michael/data10/michael/pdm/exp5data/jagsin/behav_strint7'


# Load all data
print 'Loading all N200 data...'
dataout = np.load(loadloc + '.npz')

# Set up EEG session data
sesdata = np.genfromtxt(
    '/home/michael/Dropbox/Research/PDM/encodingN200/Data/N1deflec2_allSNR_window_150_275.csv', delimiter=',')

n1lat = sesdata[:, 0]
n1deflec = sesdata[:, 1]
condition = sesdata[:, 3]
experiment = sesdata[:, 5]
rt1per = sesdata[:, 7]
rt10per = sesdata[:, 8]
rtmeanper = sesdata[:, 13]
rtmedian = sesdata[:, 11]
accmean = sesdata[:, 17]


# Plot all N200 waveforms in Experiment 2
print 'Plotting the N200 data'
plt.figure(figsize=(18, 9))
for n in range(0,dataout['n1data'].shape[1]):
    if (np.min(dataout['n1data'][:,n]) != 0):
        plt.plot(np.arange(-100, 275), dataout['n1data'][:,n],color=colororder[int(dataout['n1dataexp'][n]) -1][int(dataout['n1datacond'][n])])
zeroline = plt.plot(np.array([0., 0.]), np.array([-.15, .15]))
plt.setp(zeroline, linewidth=3, linestyle='--')
plt.xlim(-50, 275)
plt.ylim(-.15, .15)

# Plot N200 deflection distributions
kde_deflec = stats.gaussian_kde(n1deflec)
x_deflec = np.linspace(0, 275, 100)
p_deflec = kde_deflec(x_deflec)
p_deflec = .15 * p_deflec / np.max(p_deflec)
plt.fill_between(x_deflec, p_deflec, np.zeros((100)), color=color7, alpha=0.25)
p_deflec[(p_deflec < .001)] = np.nan
plt.plot(x_deflec, p_deflec, color='k', alpha=0.25)

# Plot peak latency distributions
for e in range(1,3):
    for n in range(0,3):
        thisindex = ((experiment == e) & (condition == n)) 
        kde_lat = stats.gaussian_kde(n1lat[thisindex])
        x_lat = np.linspace(150, 275, 100)
        p_lat = kde_lat(x_lat)
        p_lat = .15 * p_lat / np.max(p_lat)
        if e==1:
            thishatch='0'
        else:
            thishatch='*'
        plt.fill_between(x_lat, p_lat, np.zeros((100)), color=colororder[e-1][n], alpha=0.5,hatch=thishatch)
        p_lat[(p_lat < .001)] = np.nan
        plt.plot(x_lat, p_lat, color='k', alpha=0.25)

# Format plot
frame = plt.gca()
frame.axes.get_yaxis().set_visible(False)
plt.tick_params(axis='both', which='major', labelsize=fontsize)
plt.xlabel('Time (ms) after Gabor onset', fontsize=fontsize)

# Custom legend
custom_lines = [Line2D([0], [0], color=color3, lw=4),
                Line2D([0], [0], color=color2, lw=4),
                Line2D([0], [0], color=color1, lw=4),
                Line2D([0], [0], color=color6, lw=4, marker='*',markersize=12),
                Line2D([0], [0], color=color5, lw=4, marker='*',markersize=12),
                Line2D([0], [0], color=color4, lw=4, marker='*',markersize=12)]
plt.legend(custom_lines, ['Exp. 1 Low Noise', 'Exp. 1 Med Noise', 'Exp. 1 High Noise',
    'Exp. 2 Low Noise', 'Exp. 2 Med Noise', 'Exp. 2 High Noise'],loc=2)

# Save the figure
print 'Saving the figure...'
plt.savefig('N200_waveforms_split.png', dpi=300, format='png')