# pdm5b_simulresults5.py - Evaluates simulation results of 
#                       Hierarchical Diffusion Models with
#         trial-to-trial variability in visual encoding time (VET) and drift
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
# 12/10/18      Michael Nunez           Converted from pdm5b_simulresults4.py
# 01/02/19      Michael Nunez            Calculate variance explained (R^2 adj) for each regression


#References:
# https://stackoverflow.com/questions/893657/how-do-i-calculate-r-squared-using-python-and-numpy

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

nsims = 25
nsubs = 100

genparam = sio.loadmat('genparam_test4.mat')
rtpercentiles = np.zeros((nsims, nsubs))
for n in range(0, nsims):
    for s in range(0, nsubs):
        # Use no cutoff
        whereindx = np.where(genparam['rt'][n, s, :] > .01)
        rtpercentiles[n, s] = np.percentile(
            genparam['rt'][n, s, whereindx], 10, axis=-1)

linearmodel = 'modelfits/trialparam5_test_model%i.mat'

diags = dict()
samples = dict()
betas = np.empty((nsims, 4))
r2adj = np.empty((nsims, 4))

for m in range(1, nsims + 1):
    print 'Loading samples for simulation %d...' % m
    samples[m] = sio.loadmat(linearmodel % (m))
    diags[m] = diagnostic(samples[m])
    # plt.figure()
    # plt.scatter(np.mean(genparam['n200'][(m - 1), :,:],axis=1).T *
    #             1000, rtpercentiles[(m - 1), :] * 1000)
    VET = np.mean(genparam['n200'][(m - 1), :,:],axis=1).T * 1000
    measured_N200 = VET + np.random.normal(loc=0,scale=10,size=np.shape(VET))
    y0 = rtpercentiles[(m - 1), :] * 1000
    tempbetas0 = np.polyfit(VET,y0,deg=1)
    regression_equation0 = np.poly1d(tempbetas0)
    yhat0 = regression_equation0(VET)
    VARResidual0 = np.sum((y0 - yhat0)**2)/float(len(y0)-2)
    VARTotal0 = np.var(y0,ddof=1)
    r2adj[m - 1, 0] = 1 - VARResidual0/VARTotal0
    betas[m - 1, 0] = tempbetas0[0]
    # plt.plot(np.mean(genparam['n200'][(m - 1), :,:],axis=1).T* 1000, 
    #          np.mean(genparam['n200'][(m - 1), :,:],axis=1).T* 1000, linewidth=3, color=np.array([1, .3, 0]))
    # plt.xlabel('Real visual encoding time (ms)', fontsize=16)
    # plt.ylabel('10th RT percentile (ms)', fontsize=16)
    # plt.title('Recovery of VET from 10th RT percentiles, Model %i' %
    #           (m), fontsize=16)
    # plt.figure()
    # recovery(samples[m]['tersub'][:] * 1000,
    #          np.mean(genparam['n200'][(m - 1), :,:],axis=1).T * 1000)
    # plt.title('Recovery of VET from model fit NDT, Model %i' %
    #           (m), fontsize=16)
    # plt.xlabel('Real visual encoding time (ms)', fontsize=16)
    # plt.ylabel('Non-decision time Posteriors (ms)', fontsize=16)
    y1 = np.median(samples[m]['tersub'][:].reshape((nsubs, 6000)), axis=1) * 1000
    tempbetas1 = np.polyfit(VET,y1,deg=1)
    regression_equation1 = np.poly1d(tempbetas1)
    yhat1 = regression_equation1(VET)
    VARResidual1 = np.sum((y1 - yhat1)**2)/float(len(y1)-2)
    VARTotal1 = np.var(y1,ddof=1)
    r2adj[m - 1, 1] = 1 - VARResidual1/VARTotal1
    betas[m - 1, 1] = tempbetas1[0]
    # plt.figure()
    # recovery(samples[m]['deltasub'][:],
    #          genparam['deltasub'][(m - 1), :].T)
    # plt.title('Recovery of Mean Drift Rates from model fit, Model %i' %
    #           (m), fontsize=16)
    # plt.figure()
    # recovery(samples[m]['alphasub'][:],
    #          genparam['alphasub'][(m - 1), :].T)
    # plt.title('Recovery of Evidence bounds from model fit, Model %i' %
    #           (m), fontsize=16)
    tempbetas2 = np.polyfit(measured_N200,y0,deg=1)
    regression_equation2 = np.poly1d(tempbetas2)
    yhat2 = regression_equation2(measured_N200)
    VARResidual2 = np.sum((y0 - yhat2)**2)/float(len(y0)-2)
    VARTotal2 = np.var(y0,ddof=1)
    r2adj[m - 1, 2] = 1 - VARResidual2/VARTotal2
    betas[m - 1, 2] = tempbetas2[0]
    tempbetas3 = np.polyfit(measured_N200,y1,deg=1)
    regression_equation3 = np.poly1d(tempbetas3)
    yhat3 = regression_equation3(measured_N200)
    VARResidual3 = np.sum((y1 - yhat3)**2)/float(len(y1)-2)
    VARTotal3 = np.var(y1,ddof=1)
    r2adj[m - 1, 3] = 1 - VARResidual3/VARTotal3
    betas[m - 1, 3] = tempbetas3[0]




plt.figure()
plt.scatter(np.arange(1, nsims + 1), betas[:, 0])
tempbetas4 = np.polyfit(np.arange(1, nsims + 1), betas[:, 0], deg=1)
regress4 = plt.plot(np.arange(1, nsims + 1), tempbetas4[0] * np.arange(
    1, nsims + 1) + tempbetas4[1], linewidth=3, color=np.array([1, .3, 0]))
plt.title('Changing slopes of real VET versus 10th percentile', fontsize=16)
plt.xlabel('Simulation # (Increasing wind wandering)', fontsize=16)
plt.ylabel('Beta value', fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), betas[:, 1])
tempbetas3 = np.polyfit(np.arange(1, nsims + 1), betas[:, 1], deg=1)
regress3 = plt.plot(np.arange(1, nsims + 1), tempbetas3[0] * np.arange(
    1, nsims + 1) + tempbetas3[1], linewidth=3, color=np.array([1, .3, 0]))
plt.title('Changing slopes of real VET versus NDT posterior medians', fontsize=16)
plt.xlabel('Simulation # (Increasing wind wandering)', fontsize=16)
plt.ylabel('Beta value', fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), betas[:, 2])
tempbetas5 = np.polyfit(np.arange(1, nsims + 1), betas[:, 2], deg=1)
regress5 = plt.plot(np.arange(1, nsims + 1), tempbetas5[0] * np.arange(
    1, nsims + 1) + tempbetas5[1], linewidth=3, color=np.array([1, .3, 0]))
plt.title('Changing slopes of measured N200 versus 10th percentile', fontsize=16)
plt.xlabel('Simulation # (Increasing wind wandering)', fontsize=16)
plt.ylabel('Beta value', fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), betas[:, 3])
tempbetas6 = np.polyfit(np.arange(1, nsims + 1), betas[:, 3], deg=1)
regress6 = plt.plot(np.arange(1, nsims + 1), tempbetas6[0] * np.arange(
    1, nsims + 1) + tempbetas6[1], linewidth=3, color=np.array([1, .3, 0]))
plt.title('Changing slopes of measured N200 versus NDT posterior medians', fontsize=16)
plt.xlabel('Simulation # (Increasing wind wandering)', fontsize=16)
plt.ylabel('Beta value', fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), r2adj[:, 0])
temprsqu2 = np.polyfit(np.arange(1, nsims + 1), r2adj[:, 0], deg=1)
regress7 = plt.plot(np.arange(1, nsims + 1), temprsqu2[0] * np.arange(
    1, nsims + 1) + temprsqu2[1], linewidth=3, color=np.array([1, .3, 0]))
plt.title('Changing R^2_adj of real VET versus 10th percentile', fontsize=16)
plt.xlabel('Simulation # (Increasing wind wandering)', fontsize=16)
plt.ylabel('R^2_adj', fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), r2adj[:, 1])
temprsqu3 = np.polyfit(np.arange(1, nsims + 1), r2adj[:, 1], deg=1)
regress8 = plt.plot(np.arange(1, nsims + 1), temprsqu3[0] * np.arange(
    1, nsims + 1) + temprsqu3[1], linewidth=3, color=np.array([1, .3, 0]))
plt.title('Changing R^2_adj of real VET versus NDT posterior medians', fontsize=16)
plt.xlabel('Simulation # (Increasing wind wandering)', fontsize=16)
plt.ylabel('R^2_adj', fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), r2adj[:, 2])
temprsqu4 = np.polyfit(np.arange(1, nsims + 1), r2adj[:, 2], deg=1)
regress9 = plt.plot(np.arange(1, nsims + 1), temprsqu4[0] * np.arange(
    1, nsims + 1) + temprsqu4[1], linewidth=3, color=np.array([1, .3, 0]))
plt.title('Changing R^2_adj of measured N200 versus 10th percentile', fontsize=16)
plt.xlabel('Simulation # (Increasing wind wandering)', fontsize=16)
plt.ylabel('R^2_adj', fontsize=16)

plt.figure()
plt.scatter(np.arange(1, nsims + 1), r2adj[:, 3])
temprsqu5 = np.polyfit(np.arange(1, nsims + 1), r2adj[:, 3], deg=1)
regress10 = plt.plot(np.arange(1, nsims + 1), temprsqu5[0] * np.arange(
    1, nsims + 1) + temprsqu5[1], linewidth=3, color=np.array([1, .3, 0]))
plt.title('Changing R^2_adj of measured N200 versus NDT posterior medians', fontsize=16)
plt.xlabel('Simulation # (Increasing wind wandering)', fontsize=16)
plt.ylabel('R^2_adj', fontsize=16)


percentcontam = np.array([0, 2, 4, 6, 8, 10])
converttoindex = percentcontam*2.9 + 1 #Convert to range [1, 30]

# VETrecovery = dict()
# VETrecovery['percentilebetas'] = betas[:, 0]
# VETrecovery['posteriorbetas'] = betas[:, 1]
# VETrecovery['posteriorbetaslapse350'] = betas[:, 2]
# VETrecovery['posteriorbetaslapse'] = betas[:, 3]
# VETrecovery['percentilebetas_cutoffs'] = betas[:, 4]
# VETrecovery['converttoindex'] = converttoindex
# VETrecovery['percentcontam'] = percentcontam
# sio.savemat('VETrecovery.mat', VETrecovery)