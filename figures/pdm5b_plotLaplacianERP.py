# pdm5b_plotLaplacianERP.py - Plots Laplacian transform of the average ERP waveform
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
# 08/28/18     Michael Nunez                              Original code

## References:
# https://matplotlib.org/api/_as_gen/matplotlib.pyplot.subplots.html
# https://stackoverflow.com/questions/22607444/pyplot-shared-axes-and-no-space-between-subplots
# https://matplotlib.org/gallery/pyplots/align_ylabels.html#sphx-glr-gallery-pyplots-align-ylabels-py
# https://matplotlib.org/users/usetex.html

# Imports
from __future__ import division
import numpy as np
from scipy import stats
import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from IPython import get_ipython  # Run magic functions from script
# Initialize ipython matplotlib plotting graphics
get_ipython().magic('pylab')


# Parameters

fontsize = 18
color6 = "#E7298A"
color7 = "#A6761D"
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Load all data
ERPdata = sio.loadmat('ERPchanpos.mat')
Lapdata = sio.loadmat('ERPLap.mat')
timelocked = np.arange(-99,601)

plotindex = 101

# Create two subplots that are stacked on top of each other

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(24, 12))
fig.subplots_adjust(hspace=0)
plt.xlim(-100,600)
ax1.plot(timelocked, np.squeeze(ERPdata['scalperp'][plotindex,0:700,:]),color='k')
ax2.plot(timelocked, np.squeeze(Lapdata['laperp'][plotindex,0:700,:]), color='k')



# Create zero lines and topography markers
zeroline1 = ax1.plot(np.array([0., 0.]), np.array([-70, 30]))
plt.setp(zeroline1, linewidth=3, linestyle='--', color='b')
plotline11 = ax1.plot(np.array([92., 92.]), np.array([-70, 30]))
plt.setp(plotline11, linewidth=3, linestyle='--', color=color7)
plotline12 = ax1.plot(np.array([181., 181.]), np.array([-70, 30]))
plt.setp(plotline12, linewidth=3, linestyle='--', color=color6)
ax1.set_ylim(-70,30)
zeroline2 = ax2.plot(np.array([0., 0.]), np.array([-.125, .075]))
plt.setp(zeroline2, linewidth=3, linestyle='--', color='b')
plotline21 = ax2.plot(np.array([92., 92.]), np.array([-.125, .075]))
plt.setp(plotline21, linewidth=3, linestyle='--', color=color7)
plotline22 = ax2.plot(np.array([181., 181.]), np.array([-.125, .075]))
plt.setp(plotline22, linewidth=3, linestyle='--', color=color6)
ax2.set_ylim(-.125,.075)


# Axis labels and font size changes

ax1.tick_params(axis='both', which='major', labelsize=fontsize)
ax2.tick_params(axis='both', which='major', labelsize=fontsize)
plt.xlabel(r'Time (ms) after Gabor onset', fontsize=fontsize)


ylabelx= -.05

ax1.set_ylabel(r'ERP Amplitude ($\mu V$)', fontsize=fontsize)
ax1.yaxis.set_label_coords(ylabelx, 0.5)
ax2.set_ylabel(r'ERP spline Laplacian ($\nabla \cdot \nabla$)', fontsize=fontsize)
ax2.yaxis.set_label_coords(ylabelx, 0.5)

# Save the figure
print 'Saving the figure...'
plt.savefig('Example_ERP_Laplacian.png', dpi=300, format='png')
