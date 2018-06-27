# pdm5b_parameterrecovery.py - Evaluates simulation results of fitting
#                         models with trial-to-trial variability
#                         in non-decision time and drift-rate
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
# 01/30/18      Michael Nunez              Converted from pdm5b_simulresults.py
# 06/18/18      Michael Nunez             Addition of results with modeled lapse trials
# 06/20/18      Michael Nunez           Addition of results without 350 ms cutoffs

# Imports
import numpy as np
import scipy.io as sio
from scipy import stats
import matplotlib.pyplot as plt
from IPython import get_ipython  # Run magic functions from script
get_ipython().magic('pylab')  # Initialize ipython matplotlib plotting graphics


def recovery(possamps, truevals):  # Parameter recovery plots
    """Plots true parameters versus 99% and 95% credible intervals of recovered
    parameters. Also plotted are the median and mean of the posterior
    distributions

    Parameters
    ----------
    possamps : ndarray of posterior chains where the last dimension is the
    number of chains, the second to last dimension is the number of samples in
    each chain, all other dimensions must match the dimensions of truevals

    truevals : ndarray of true parameter values
    """

    # Number of chains
    nchains = possamps.shape[-1]

    # Number of samples per chain
    nsamps = possamps.shape[-2]

    # Number of variables to plot
    nvars = np.prod(possamps.shape[0:-2])

    # Reshape data
    alldata = np.reshape(possamps, (nvars, nchains, nsamps))
    alldata = np.reshape(alldata, (nvars, nchains * nsamps))
    truevals = np.reshape(truevals, (nvars))

    # Plot properties
    LineWidths = np.array([2, 5])
    teal = np.array([0, .7, .7])
    blue = np.array([0, 0, 1])
    orange = np.array([1, .3, 0])
    Colors = [teal, blue]

    for v in range(0, nvars):
        # Compute percentiles
        bounds = stats.scoreatpercentile(alldata[v, :], (.5, 2.5, 97.5, 99.5))
        for b in range(0, 2):
            # Plot credible intervals
            credint = np.ones(100) * truevals[v]
            y = np.linspace(bounds[b], bounds[-1 - b], 100)
            lines = plt.plot(credint, y)
            plt.setp(lines, color=Colors[b], linewidth=LineWidths[b])
            if b == 1:
                # Mark median
                mmedian = plt.plot(truevals[v], np.median(alldata[v, :]), 'o')
                plt.setp(mmedian, markersize=10, color=[0., 0., 0.])
                # Mark mean
                mmean = plt.plot(truevals[v], np.mean(alldata[v, :]), '*')
                plt.setp(mmean, markersize=10, color=teal)
    # Plot line y = x
    tempx = np.linspace(np.min(truevals), np.max(
        truevals), num=100)
    recoverline = plt.plot(tempx, tempx)
    plt.setp(recoverline, linewidth=3, color=orange)


def diagnostic(insamples):  # Diagnostics
    """
    Returns Rhat (measure of convergence, less is better with an approximate
    1.10 cutoff) and Neff, number of effective samples).

    Reference: Gelman, A., Carlin, J., Stern, H., & Rubin D., (2004).
              Bayesian Data Analysis (Second Edition). Chapman & Hall/CRC:
              Boca Raton, FL.


    Parameters
    ----------
    insamples: dic
        Sampled values of monitored variables as a dictionary where keys
        are variable names and values are numpy arrays with shape:
        (dim_1, dim_n, iterations, chains). dim_1, ..., dim_n describe the
        shape of variable in JAGS model.

    Returns
    -------
    dict:
        Rhat for each variable. Prints Maximum Rhat
    """

    result = {}  # Initialize dictionary
    maxrhats = np.zeros((len(insamples.keys())), dtype=float)
    keyindx = 0
    for key in insamples.keys():
        if key[0] != '_':
            result[key] = {}

            possamps = insamples[key]

            nchains = possamps.shape[-1]
            nsamps = possamps.shape[-2]
            # Mean of each chain
            chainmeans = np.mean(possamps, axis=-2)
            # Global mean of each parameter
            globalmean = np.mean(chainmeans, axis=-1)
            result[key]['mean'] = globalmean
            globalmeanext = np.expand_dims(
                globalmean, axis=-1)  # Expand the last dimension
            globalmeanext = np.repeat(
                globalmeanext, nchains, axis=-1)  # For differencing
            # Between-chain variance
            between = np.sum(np.square(chainmeans - globalmeanext),
                             axis=-1) * nsamps / (nchains - 1.)
            # Mean of the variances of each chain
            within = np.mean(np.var(possamps, axis=-2), axis=-1)
            # Total estimated variance
            totalestvar = (1. - (1. / nsamps)) * \
                within + (1. / nsamps) * between
            # Rhat (Gelman-Rubin statistic)
            temprhat = np.sqrt(totalestvar / within)
            maxrhats[keyindx] = np.nanmax(temprhat)  # Ignore NANs
            result[key]['rhat'] = temprhat
            keyindx += 1
            # Possible number of effective samples?
            # Geweke statistic?
    print "Maximum Rhat: %3.2f" % (np.max(maxrhats))
    return result

nsims = 30
nsubs = 100

genparam = sio.loadmat('../Parameter_recovery/genparam_test.mat')
rtpercentiles_cutoff = np.zeros((nsims, nsubs))
rtpercentiles = np.zeros((nsims, nsubs))
for n in range(0, nsims):
    for s in range(0, nsubs):
        # Use no cutoff
        rtpercentiles[n, s] = np.percentile(
            genparam['rt'][n, s, :], 10, axis=-1)
        # Use 350 ms cutoff
        whereindx_cutoff = np.where(genparam['rt'][n, s, :] > .35)
        rtpercentiles_cutoff[n, s] = np.percentile(
            genparam['rt'][n, s, whereindx_cutoff], 10, axis=-1)

linearmodel1 = '../Parameter_recovery/modelfits/trialparam_test_model%i.mat'
linearmodel2 = '../Parameter_recovery/modelfits/trialparam2_test_model%i.mat'
linearmodel3 = '../Parameter_recovery/modelfits/trialparam3_test_model%i.mat'


diags = dict()
samples = dict()
diags2 = dict()
samples2 = dict()
diags3 = dict()
samples3 = dict()
betas = np.empty((nsims, 5))

for m in range(1, nsims + 1):
    samples[m] = sio.loadmat(linearmodel1 % (m))
    samples2[m] = sio.loadmat(linearmodel2 % (m))
    samples3[m] = sio.loadmat(linearmodel3 % (m))
    diags[m] = diagnostic(samples[m])
    diags2[m] = diagnostic(samples2[m])
    diags3[m] = diagnostic(samples3[m])
    tempbetas0 = np.polyfit(
        genparam['tersub'][(m - 1), :].T, rtpercentiles[(m - 1), :], deg=1)
    betas[m - 1, 0] = tempbetas0[0]
    tempbetas1 = np.polyfit(genparam['tersub'][(m - 1), :] * 1000, np.median(
        samples[m]['tersub'][:].reshape((nsubs, 6000)), axis=1) * 1000, deg=1)
    betas[m - 1, 1] = tempbetas1[0]
    tempbetas2 = np.polyfit(genparam['tersub'][(m - 1), :] * 1000, np.median(
        samples2[m]['tersub'][:].reshape((nsubs, 6000)), axis=1) * 1000, deg=1)
    betas[m - 1, 2] = tempbetas2[0]
    tempbetas3 = np.polyfit(genparam['tersub'][(m - 1), :] * 1000, np.median(
        samples3[m]['tersub'][:].reshape((nsubs, 6000)), axis=1) * 1000, deg=1)
    betas[m - 1, 3] = tempbetas3[0]
    tempbetas4 = np.polyfit(
        genparam['tersub'][(m - 1), :].T, rtpercentiles_cutoff[(m - 1), :], deg=1)
    betas[m - 1, 4] = tempbetas4[0]


plt.figure()
plt.scatter(genparam['tersub'][0, :].T *
            1000, rtpercentiles[0, :] * 1000)
tempbetas5 = np.polyfit(genparam['tersub'][0, :].T *
            1000, rtpercentiles[0, :] * 1000, deg=1)
regress5 = plt.plot(genparam['tersub'][0, :].T *
            1000, tempbetas5[0] * genparam['tersub'][0, :].T *
            1000 + tempbetas5[1], linewidth=3)
plt.plot(genparam['tersub'][0, :].T * 1000, genparam['tersub']
         [0, :].T * 1000, linewidth=3, color=np.array([1, .3, 0]))
plt.xlabel('Real non-decision time (ms)', fontsize=16)
plt.ylabel('10th reaction time percentile (ms)', fontsize=16)
plt.figure()
recovery(samples[1]['tersub'][:] * 1000,
         genparam['tersub'][0, :].T * 1000)
plt.xlabel('Real non-decision time (ms)', fontsize=16)
plt.ylabel('Non-decision time posteriors (ms)', fontsize=16)
# Note that simuldiff() is on a different evidence scale since the
# diffusion parameter is .1 and dwiener assumes the diffusion parameter =1
plt.figure()
recovery(samples[1]['deltasub'][:],
         genparam['deltasub'][0, :].T * 10)
plt.xlabel('Real evidence accumulation rate (evidence units)', fontsize=16)
plt.ylabel('Evidence accumulation posteriors (evidence units)', fontsize=16)
plt.figure()
recovery(samples[1]['alphasub'][:],
         genparam['alphasub'][0, :].T * 10)
plt.xlabel('Real evidence boundaries (evidence units)', fontsize=16)
plt.ylabel('Evidence boundary posteriors (evidence units)', fontsize=16)




percentcontam = np.array([0, 2, 4, 6, 8, 10])
converttoindex = percentcontam*2.9 + 1 #Convert to range [1, 30]

plt.figure()
plt.scatter(np.arange(1, nsims + 1), betas[:, 0])
tempbetas5 = np.polyfit(np.arange(1, nsims + 1), betas[:, 0], deg=1)
regress5 = plt.plot(np.arange(1, nsims + 1), tempbetas5[0] * np.arange(
    1, nsims + 1) + tempbetas5[1], linewidth=3)
ax1 = plt.gca()
ax1.set_xticks(converttoindex)
ax1.set_xticklabels(percentcontam)
plt.xlim([1, 30])
plt.ylim([0.5, 1.5])
plt.title('Changing slopes of real NDT versus 10th percentile', fontsize=16)
plt.xlabel('Percentage of contaminant trials (%)', fontsize=16)
plt.ylabel('Beta value', fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), betas[:, 1])
tempbetas6 = np.polyfit(np.arange(1, nsims + 1), betas[:, 1], deg=1)
regress6 = plt.plot(np.arange(1, nsims + 1), tempbetas6[0] * np.arange(
    1, nsims + 1) + tempbetas6[1], linewidth=3)
ax2 = plt.gca()
ax2.set_xticks(converttoindex)
ax2.set_xticklabels(percentcontam)
plt.xlim([1, 30])
plt.ylim([0.5, 1.5])
plt.title('Changing slopes of real NDT versus posterior medians', fontsize=16)
plt.xlabel('Percentage of contaminant trials (%)', fontsize=16)
plt.ylabel('Beta value', fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), betas[:, 2])
tempbetas7 = np.polyfit(np.arange(1, nsims + 1), betas[:, 2], deg=1)
regress7 = plt.plot(np.arange(1, nsims + 1), tempbetas7[0] * np.arange(
    1, nsims + 1) + tempbetas7[1], linewidth=3)
ax2 = plt.gca()
ax2.set_xticks(converttoindex)
ax2.set_xticklabels(percentcontam)
plt.xlim([1, 30])
plt.ylim([0.5, 1.5])
plt.title('Changing slopes of real NDT versus posterior medians with lapse trials and 350 ms cutoffs', fontsize=16)
plt.xlabel('Percentage of contaminant trials (%)', fontsize=16)
plt.ylabel('Beta value', fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), betas[:, 3])
tempbetas8 = np.polyfit(np.arange(1, nsims + 1), betas[:, 3], deg=1)
regress8 = plt.plot(np.arange(1, nsims + 1), tempbetas8[0] * np.arange(
    1, nsims + 1) + tempbetas8[1], linewidth=3)
ax2 = plt.gca()
ax2.set_xticks(converttoindex)
ax2.set_xticklabels(percentcontam)
plt.xlim([1, 30])
plt.ylim([0.5, 1.5])
plt.title('Changing slopes of real NDT versus posterior medians with lapse trials', fontsize=16)
plt.xlabel('Percentage of contaminant trials (%)', fontsize=16)
plt.ylabel('Beta value', fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), betas[:, 4])
tempbetas9 = np.polyfit(np.arange(1, nsims + 1), betas[:, 4], deg=1)
regress9 = plt.plot(np.arange(1, nsims + 1), tempbetas9[0] * np.arange(
    1, nsims + 1) + tempbetas9[1], linewidth=3)
ax1 = plt.gca()
ax1.set_xticks(converttoindex)
ax1.set_xticklabels(percentcontam)
plt.xlim([1, 30])
plt.ylim([0.5, 1.5])
plt.title('Changing slopes of real NDT versus 10th percentile after 350 ms cutoffs', fontsize=16)
plt.xlabel('Percentage of contaminant trials (%)', fontsize=16)
plt.ylabel('Beta value', fontsize=16)

changingbetas = dict()
changingbetas['percentilebetas'] = betas[:, 0]
changingbetas['posteriorbetas'] = betas[:, 1]
changingbetas['posteriorbetaslapse350'] = betas[:, 2]
changingbetas['posteriorbetaslapse'] = betas[:, 3]
changingbetas['percentilebetas_cutoffs'] = betas[:, 4]
changingbetas['converttoindex'] = converttoindex
changingbetas['percentcontam'] = percentcontam
sio.savemat('recovery_changingbetas.mat', changingbetas)
