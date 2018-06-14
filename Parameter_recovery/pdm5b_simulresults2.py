# pdm5b_simulresults.py - Evaluates simulation results of fitting
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
# 06/12/18      Michael Nunez           Converted from pdm5b_simulresults.py

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

genparam = sio.loadmat('genparam_test.mat')
rtpercentiles = np.zeros((nsims, nsubs))
for n in range(0, nsims):
    for s in range(0, nsubs):
        # Use 350 ms cutoff
        whereindx = np.where(genparam['rt'][n, s, :] > .35)
        rtpercentiles[n, s] = np.percentile(
            genparam['rt'][n, s, whereindx], 10, axis=-1)

linearmodel = 'modelfits/trialparam2_test_model%i.mat'

diags = dict()
samples = dict()
betas = np.empty((nsims, 2))

for m in range(1, nsims + 1):
    samples[m] = sio.loadmat(linearmodel % (m))
    diags[m] = diagnostic(samples[m])
    plt.figure()
    plt.scatter(genparam['tersub'][(m - 1), :].T *
                1000, rtpercentiles[(m - 1), :] * 1000)
    tempbetas0 = np.polyfit(
        genparam['tersub'][(m - 1), :].T, rtpercentiles[(m - 1), :], deg=1)
    betas[m - 1, 0] = tempbetas0[0]
    plt.plot(genparam['tersub'][(m - 1), :].T * 1000, genparam['tersub']
             [(m - 1), :].T * 1000, linewidth=3, color=np.array([1, .3, 0]))
    plt.xlabel('Real non-decision time (ms)', fontsize=16)
    plt.ylabel('10th RT percentile (ms)', fontsize=16)
    plt.title('Recovery of Non-decision times from 10th RT percentiles, Model %i' %
              (m), fontsize=16)
    plt.figure()
    recovery(samples[m]['tersub'][:] * 1000,
             genparam['tersub'][(m - 1), :].T * 1000)
    plt.title('Recovery of Non-decision times from model fit, Model %i' %
              (m), fontsize=16)
    plt.xlabel('Real non-decision time (ms)', fontsize=16)
    plt.ylabel('Non-decision time Posteriors (ms)', fontsize=16)
    tempbetas1 = np.polyfit(genparam['tersub'][(m - 1), :] * 1000, np.median(
        samples[m]['tersub'][:].reshape((nsubs, 6000)), axis=1) * 1000, deg=1)
    betas[m - 1, 1] = tempbetas1[0]
    # Note that simuldiff() is on a different evidence scale since the
    # diffusion parameter is .1 and dwiener assumes the diffusion parameter =1
    plt.figure()
    recovery(samples[m]['deltasub'][:],
             genparam['deltasub'][(m - 1), :].T * 10)
    plt.title('Recovery of Mean Drift Rates from model fit, Model %i' %
              (m), fontsize=16)
    plt.figure()
    recovery(samples[m]['alphasub'][:],
             genparam['alphasub'][(m - 1), :].T * 10)
    plt.title('Recovery of Evidence bounds from model fit, Model %i' %
              (m), fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), betas[:, 0])
tempbetas2 = np.polyfit(np.arange(1, nsims + 1), betas[:, 0], deg=1)
regress2 = plt.plot(np.arange(1, nsims + 1), tempbetas2[0] * np.arange(
    1, nsims + 1) + tempbetas2[1], linewidth=3, color=np.array([1, .3, 0]))
plt.title('Changing slopes of real ndt versus 10th percentile', fontsize=16)
plt.xlabel('Simulation # (Increasing wind wandering)', fontsize=16)
plt.ylabel('Beta value', fontsize=16)
plt.figure()
plt.scatter(np.arange(1, nsims + 1), betas[:, 1])
tempbetas3 = np.polyfit(np.arange(1, nsims + 1), betas[:, 1], deg=1)
regress3 = plt.plot(np.arange(1, nsims + 1), tempbetas3[0] * np.arange(
    1, nsims + 1) + tempbetas3[1], linewidth=3, color=np.array([1, .3, 0]))
plt.title('Changing slopes of real ndt versus posterior medians', fontsize=16)
plt.xlabel('Simulation # (Increasing wind wandering)', fontsize=16)
plt.ylabel('Beta value', fontsize=16)
