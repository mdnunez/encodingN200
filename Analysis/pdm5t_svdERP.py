# pdm5t_svdERP.py - Finds SVD decompositions of
#                   stimulus- and response-locked event related potentials
#                   SVD decomposition per condition, less aggressive filtering
#                   aggressive windowing, specific decomposition per condition
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
# 11/20/17      Michael Nunez                        Converted from pdm5b_svdERP2.py
# 01/03/19      Michael Nunez          Run lowpass filtering with different parameters
# 01/04/19      Michael Nunez            Use attenuation parameter 5.0 dB for more accurate filter
#                                          Revert to original script



# Imports
import numpy as np
import scipy.io as sio
import os
from scipy.signal import buttord, butter, filtfilt, freqz
import matplotlib.pyplot as plt
from IPython import get_ipython  # Run magic functions from script
ipython = get_ipython()

# Butterworth and filtfilt wrapper


def butterfilt(data, sr, passband=(1.0, 10.0), stopband=(0.25, 20.0),
               attenuation=(1.0, 10.0), plotfreqz=False, plotlim=100):
    """Wrapper for Butterworth IIR filter of sample*channel*trial EEG data

    See: https://dsp.stackexchange.com/questions/31420/butterworth-band-pass-filter

    Inputs:
    data - sample*channel*trial EEG data
    sr - sample rate (samps/sec)

    Optional Inputs:
    passband - low and high cutoffs for Butterworth passband
    stopband - low and high cuttoffs for Butterworth stopband
    attenuation - Maximum loss in the passband and stopband respectively (dB)
    plotfreqz - Flag for plotting frequency responses of both filters
    plotlim - Upper limit of frequencies to plot

    Outputs:
    filtdata - Filtered sample*channel*trial EEG data
    """

    nyquist = .5 * float(sr)

    # Butterworth filter order selection
    N, Wn = buttord(wp=[float(passband[0]) / nyquist, float(passband[1]) / nyquist], ws=[float(stopband[0]) / nyquist, float(stopband[1]) / nyquist],
                    gpass=attenuation[0], gstop=attenuation[1], analog=False)

    # Filter the data
    b, a = butter(N, Wn, 'band', False)
    filtdata = filtfilt(b, a, data, axis=0)

    if plotfreqz:
        w, h = freqz(b, a)
        plt.plot((nyquist / np.pi) * w, abs(h))
        plt.setp(plt.gca(), XLim=[0, plotlim], YLim=[0, 1.1])
        plt.plot([0, nyquist], [np.sqrt(0.5), np.sqrt(0.5)], '--')
        plt.title('Butterworth Passband Frequency Response')

    return filtdata


def baseline(data, wind=range(0, 50)):
    """Simple function to baseline EEG data for ERP calculations

    Inputs:
    data - sample*channel*trial or sample*trial*channel EEG data
    wind - subtracting the mean of this window to re-center EEG data

    Outputs:
    recentered - Re-centered EEG data
    """

    baselines = np.squeeze(np.mean(data[wind, :, :], axis=0))
    recentered = data - np.tile(baselines, (data.shape[0], 1, 1))

    return recentered


def epochsubset(data, newindex, lockindex=None):
    """Reepochs each epoch of EEG data by timelocking each epoch
    "i" to "newindex[i]" with maximum window size available

    Inputs:  
    data - sample*channel*trial EEG data
    newindex - Vector of length "trial"

    Optional inputs:
    lockindex - Sample in which newdata is timelocked
                Default: nanmin(newindex)

    Outputs:  
    newdata - Re-timelocked EEG data          
    lockindex - Sample in which newdata is timelocked
    badtrials - Index of trials where newindex contained nans
    """

    if lockindex is None:
        lockindex = int(np.nanmin(newindex))

    windsize = (np.shape(data)[0] - int(np.nanmax(newindex))) + int(lockindex)

    newdata = np.zeros(
        (int(windsize), int(np.shape(data)[1]), int(np.shape(data)[2])))

    for t in range(0, np.size(newindex)):
        if np.isfinite(newindex[t]):
            begin_index = int(newindex[t]) - lockindex
            end_index = windsize + begin_index
            newdata[:, :, t] = data[begin_index: end_index, :, t]

    badtrials = np.where(np.isnan(newindex))

    return newdata, lockindex, badtrials


# Data save and load locations
dataloc = '/data10/michael/pdm/exp5data/subjects/training/{0}/{0}_allcleaned.npz'
svd_saveloc = '/data10/michael/pdm/exp5data/subjects/training/{0}/{1}/erp_svd_{0}_{1}_v5.mat'
# svd_saveloc = '/data10/michael/pdm/exp5data/subjects/training/{0}/{1}/erp_svd_{0}_{1}_v6.mat'
indxloc = '/data10/michael/pdm/exp5data/subjects/training/{0}/{0}_traintestindx.npz'

# Subjects
subjects = ['s100', 's59', 's109', 's110']

# Conditions
snrlabels = ['low', 'med', 'high']

# Sessions
sessions = []
for i in range(1, 8):
    sessions.extend(['ses%i' % i])

# Parietal channels to keep similar signs (+ or -) of weights across subjects
poschans = np.array([30, 36, 41, 46])
poschans = np.concatenate((poschans, np.arange(50, 55)))
poschans = np.concatenate((poschans, np.arange(57, 101)))


# Percentage variance explained
def perexp(s):
    pexp = np.square(s) / float(np.sum(np.square(s)))
    return pexp

for subdes in subjects:
    print 'Loading data for subject %s...' % (subdes)
    # Load EEG Data
    data = np.load(dataloc.format(subdes))

    sr = 1000.  # Make sure it is a floating point number

    print 'Filtering data for subject %s...' % (subdes)
    filtered = butterfilt(data['eeg'], sr, passband=(
        1.0, 10.0), stopband=(0.25, 20.0))
    # filtered = butterfilt(data['eeg'], sr, passband=(
    #     1.0, 30.0), stopband=(0.25, 40.0), attenuation=(1.0,5.0))
    filtered2 = butterfilt(data['eeg'], sr, passband=(
        0.1, 4.0), stopband=(0.01, 8.0))

    noiseonsets = data['noiseonsets']
    noiseonsets[noiseonsets < 100] = np.nan
    noiseonsets[noiseonsets > 1000] = np.nan
    datart = data['rt']
    datart[datart < 0] = np.nan

    if any(data['session'] == 'missing'):
        sestrack = np.where(data['session'] == 'missing')[0][0]
        nantrials = np.zeros((480 * np.shape(data['session'])[0]), dtype=bool)
        nantrials[np.arange(0, 480) + 480 * (sestrack)] = True
        noiseonsets[nantrials] = np.nan
        datart[nantrials] = np.nan

    print 'Finding cue-locked and response-locked EEG for subject %s...' % (subdes)
    cue_eeg, cue_lockindex, cue_badtrials = epochsubset(
        filtered, noiseonsets)
    resp_eeg, resp_lockindex, resp_badtrials = epochsubset(
        filtered2, datart + 1250)

    print 'Removing baselines for subject %s...' % (subdes)
    # Response interval, stimulus-locked
    rint_eeg = baseline(filtered[1150:2250, :, :], wind=range(0, 100))

    # Cue interval, stimulus-locked
    cue_eeg = baseline(
        cue_eeg[(cue_lockindex - 100):(cue_lockindex + 500), :, :], wind=range(0, 100))

    # Response interval, response-locked
    resp_eeg = baseline(
        resp_eeg[(resp_lockindex - 1250):(resp_lockindex), :, :], wind=range(0, 100))

    # Load training/test index
    indx = np.load(indxloc.format(subdes))

    # Trials to keep
    keeptrials = np.union1d(np.where(np.isfinite(data['noiseonsets']))[0],
                            np.where(np.sum(data['artifact'], axis=0) != 129)[0])
    keeptrials = np.setdiff1d(keeptrials, cue_badtrials)
    keeptrials = np.setdiff1d(keeptrials, resp_badtrials)

    # Non-artifact training trials
    traintrials = np.union1d(indx['train'], keeptrials)

    sestrack = 0
    for ses in sessions:
        temptrials = np.intersect1d(
            traintrials, np.arange(0, 480) + 480 * (sestrack))
        temptrials2 = np.arange(0, 480) + 480 * (sestrack)
        snrtrack = 0
        svds = dict()  # Run SVD to extract single-trial estimates
        for snr in np.array([0.5, 1.0, 2.0]):
            print 'Calculating weights for subject %s, session %s, %0.1f SNR...' % (subdes, ses, snr)
            # Extract trials
            snrtrials = np.where(data['snrvec'] == snr)

            thesetrials = np.intersect1d(temptrials, snrtrials)
            alltrials = np.intersect1d(temptrials2, snrtrials)

            # Channels to keep
            keepchans = np.squeeze(np.array(
                np.where(np.logical_not(np.all(np.squeeze(data['artifact'][:, alltrials]), 1)))))

            # Maintain same extraction across conditions, do this for 3 time-locked
            # intervals
            windows = ['rint', 'cue', 'resp']
            for w in np.arange(0, 3):
                exec('tempeeg = %s_eeg' % windows[w])

                ixgrid1 = np.ix_(np.arange(0, tempeeg.shape[0]),
                                 keepchans, thesetrials)  # Use only training data for ERP
                ixgrid2 = np.ix_(np.arange(0, tempeeg.shape[0]),
                                 keepchans, alltrials)  # Calculate single-trial ERPs in all data

                erp = np.squeeze(np.mean(tempeeg[ixgrid1], axis=-1))

                # Non-statistical SVD
                print 'Computing singular value decomposition of %s %s for window %s' % (
                    subdes, ses, windows[w])
                u, s, v = np.linalg.svd(erp, full_matrices=0)
                # Flip weights in the first component
                if np.sum(v[0, poschans]) < 0.:
                    v[0, :] = -v[0, :]
                    u[:, 0] = -u[:, 0]
                svds['svd_%s_%s' %
                     (windows[w], snrlabels[snrtrack])] = dict()
                print 'svd_%s_%s' % (windows[w], snrlabels[snrtrack])
                svds['svd_%s_%s' % (windows[w], snrlabels[snrtrack])][
                    'erpinput'] = erp
                svds['svd_%s_%s' %
                     (windows[w], snrlabels[snrtrack])]['u'] = u
                svds['svd_%s_%s' %
                     (windows[w], snrlabels[snrtrack])]['s'] = s
                svds['svd_%s_%s' % (windows[w], snrlabels[snrtrack])][
                    'perexp'] = perexp(s)
                svds['svd_%s_%s' %
                     (windows[w], snrlabels[snrtrack])]['v'] = v
                svds['goodchans'] = keepchans + 1
                svds['st_erp_%s_%s' % (windows[w], snrlabels[snrtrack])] = np.tensordot(
                    tempeeg[ixgrid2], v[0, :], axes=(1, 0))  # Note the superiority of np.tensordot() to matrix algebra
                print 'Component 1 of svd_%s_%s explained %0.3f%% variance' % (
                    windows[w], snrlabels[snrtrack], perexp(s)[0])
            snrtrack += 1
        sio.savemat(svd_saveloc.format(subdes, ses), svds)
        sestrack += 1
