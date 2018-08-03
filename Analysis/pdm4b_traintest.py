# pdm4b_traintest.py - Saves indicies for training and test
#                      trial indices
#
# Copyright (C) 2016 Michael D. Nunez, <mdnunez1@uci.edu>
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
# 06/22/16      Michael Nunez                   Converted from pdm5t_traintest.py


# Imports
import numpy as np

# Save location
indxloc = '/data10/michael/pdm/exp4data/subjects/{0}/{0}_traintestindx.npz'

# Subjects
subjects = ['s100', 's101', 's59', 's64', 's68', 's80',
            's82', 's93', 's94', 's95', 's96', 's97']


for subdes in subjects:
    np.random.seed(int(subdes[1:]))  # Use this seed for training / test split
    shuffled = np.empty([480 * 2], dtype=int)
    index = dict()
    index['train'] = np.empty([320 * 2], dtype=int)  # 2/3 training
    index['test'] = np.empty([160 * 2], dtype=int)  # 1/3 test
    for n in range(0, 2):
        shuffled[(480 * n):(480 + 480 * n)
                 ] = np.random.permutation(480) + 480 * n
        index['train'][(320 * n):(320 + 320 * n)
                       ] = shuffled[(480 * n):(320 + 480 * n)]
        index['test'][(160 * n):(160 + 160 * n)
                      ] = shuffled[(320 + 480 * n):(480 + 480 * n)]
    np.savez(indxloc.format(subdes), **index)
