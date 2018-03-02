# pdm5b_scatterplot.py - Creates scatter plots of N200 latencies versus RT percentiles
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
# 06/24/16      Michael Nunez                             Original code
# 12/06/17      Michael Nunez                                Cleanup

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Parameters

fontsize = 20

# Initial

# Set up EEG session data
trialdata = np.genfromtxt(
    '../Data/N200_rt_window_150_275_fixed350cutoff.csv', delimiter=',')
sesdata = np.genfromtxt(
    '../Data/N1deflec2_cutoffs_allSNR_window_150_275_fixed350cutoff.csv', delimiter=',')

n1lat = sesdata[:, 0]
rt1per = sesdata[:, 7]
rt10per = sesdata[:, 8]
rtmeanper = sesdata[:, 13]

n200lat = trialdata[:, 0]
rt = trialdata[:, 2]

# Plot data

f, (ax1, ax2) = plt.subplots(1, 2)

index = np.isfinite(n1lat)
ultb = np.polyfit(n1lat[index], rt10per[index], deg=1)
ax1.scatter(n1lat[index], rt10per[index])
ax1.plot(n1lat[index], ultb[0] * n1lat[index] +
         ultb[1], '--', color='red', linewidth=3)
ax1.set_xlabel('Trial-averaged N200 Latency (ms)', fontsize=fontsize)
ax1.set_xticks(np.arange(150, 300, 25))
ax1.set_xlim((150, 275))
ax1.set_ylim((400, 1400))
ax1.set_ylabel('10th RT Percentile (ms)', fontsize=fontsize)

ax2.scatter(n200lat, rt)
tb = np.polyfit(n200lat, rt, deg=1)
ax2.plot(n200lat, tb[0] *
         n200lat + tb[1], '--', color='red', linewidth=3)
ax2.set_xlabel('Single-trial N200 Latency (ms)', fontsize=fontsize)
ax2.set_xticks(np.arange(150, 300, 25))
ax2.set_xlim((150, 275))
ax2.set_ylim((400, 1400))
ax2.set_ylabel('Reaction Time (ms)', fontsize=fontsize)


plt.tight_layout()
fig = plt.gcf()
fig.set_size_inches(20, 12)
figManager = plt.get_current_fig_manager()  # Maximize screen
figManager.window.showMaximized()
plt.savefig('n1lat_rt10per_n200_rt.png', dpi=300)
