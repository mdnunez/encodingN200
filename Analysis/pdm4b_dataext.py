# pdm4b_dataext.py - Extracts cleaned EEG from each subject, saves data
#
# Copyright (C) 2023 Michael D. Nunez, <m.d.nunez@uva.nl>
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
# 06/22/16     Michael Nunez                   Converted from pdm5b_dataext.py
# 11-Aug-2023  Michael Nunez         Extract spatial frequency information

# Imports
import numpy as np
import scipy.io as sio
import os.path

dataloc = '/media/michael/My Book/data4/pdm/exp4data/subjects/{0}/{1}/{0}_cleaned.mat'
saveloc = '/media/michael/My Book/data4/pdm/exp4data/subjects/{0}/{0}_allcleaned.npz'

subjects = ['s100', 's101', 's64', 's68', 's80',
            's82', 's93', 's94', 's95', 's96', 's97', 's59']

subtrack = 0
for subdes in subjects:
    sestrack = 0
    data = dict()
    modelins = dict()

    data['eeg'] = np.empty((3500, 129, 480 * 2))
    data['rt'] = np.empty((480 * 2))
    data['correct'] = np.empty((480 * 2))
    data['noiseonsets'] = np.empty((480 * 2))
    data['noisetimes'] = np.empty((480 * 2))
    data['snrvec'] = np.empty((480 * 2))
    data['randrots'] = np.empty((480 * 2))
    data['spfs'] = np.empty((480 * 2)) # Added 11-Aug-2023
    data['block'] = np.empty((480 * 2))
    data['artifact'] = np.empty((129, 480 * 2), dtype='bool')
    data['session'] = []
    for i in range(1, 3):
        data['session'].extend(['s%i' % i])
    for session in data['session']:
        if os.path.isfile(dataloc.format(subdes, session)):
            print 'Loading data from subject %s, session %s...' % (subdes, session)
            dataout = sio.loadmat(dataloc.format(subdes, session), variable_names=[
                'data', 'sr', 'photo', 'artifact', 'expinfo'])

            thesetrials = np.arange(0, 480) + 480 * sestrack
            data['eeg'][:, :, thesetrials] = dataout['data']
            data['artifact'][:, thesetrials] = np.array(dataout[
                'artifact'], dtype=bool)
            data['rt'][thesetrials] = np.squeeze(dataout[
                'expinfo'][0][0]['rt'])
            data['correct'][thesetrials] = np.squeeze(dataout[
                'expinfo'][0][0]['correct'])
            data['snrvec'][thesetrials] = np.squeeze(dataout[
                'expinfo'][0][0]['snrvec'])
            data['noisetimes'][thesetrials] = np.squeeze(dataout[
                'expinfo'][0][0]['noisetimes'])
            data['randrots'][thesetrials] = np.squeeze(dataout[
                'expinfo'][0][0]['randrots'])
            data['spfs'][thesetrials] = np.squeeze(dataout[
                'expinfo'][0][0]['spfs']) # Added 11-Aug-2023
            data['block'][thesetrials] = np.concatenate((np.ones(60)*1,
                np.ones(60)*2, np.ones(60)*3, np.ones(60)*4, np.ones(60)*5,
                np.ones(60)*6, np.ones(60)*7, np.ones(60)*8))
            for t in range(0, 480):
                tempfound = np.where(
                    (dataout['photo'][:, 1, t] > 1) & (np.cumsum(
                        dataout['photo'][:, 1, t]) > 0) & (np.cumsum(
                            dataout['photo'][:, 1, t]) < 90))
                if (tempfound[0].shape[0] == 0):  # If empty'
                    data['artifact'][:, thesetrials[t]] = True
                else:
                    data['noiseonsets'][thesetrials[t]] = np.array(
                        tempfound[0][0], dtype='int32')

            dataout = None
        else:
            print 'Could not find data from subject %s, session %s ...' % (subdes, session)
            data['session'][sestrack] = 'missing'

        sestrack += 1
    print 'Saving all %s''s EEG data...' % subdes
    np.savez(saveloc.format(subdes), **data)
    subtrack += 1
