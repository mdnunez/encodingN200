%PDM5_BATCH - Script performs data extraction in batch
%
% Copyright (C) 2018 Michael D. Nunez, <mdnunez1@uci.edu>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  06/14/18       Michael Nunez                 Original code

%% Initial
subjects = {'s100', 's59', 's109', 's110'};

sessions = {'ses1', 'ses2', 'ses3', 'ses4', 'ses5', 'ses6', 'ses7'};

fail = zeros(length(subjects),length(sessions));

% %% Code

for m = 1:numel(subjects),
	for n=1:numel(sessions),
		sub = subjects{m};
		ses = sessions{n};
		cd(sprintf('/data10/michael/pdm/exp5data/subjects/training/%s/%s',sub,ses));

		if exist(sprintf('%s_%s_rawfilt.mat',sub,ses)) ~= 2, %If file doesn't already exist
			fprintf('Extracting data for %s %s...\n',sub,ses);
			%Set up for artscreenEEG scripts
			try
				eeg = pdm5_dataext;
			catch me
				fprintf('Error in extracting data from %s %s!\n',sub,ses);
				fprintf('Reason: \n');
				fprintf('%s\n',me.message);
				fail(m,n) = 1;
				break;
			end

			%Save raw data
			fprintf('%%Saving raw data...\n');
			save(sprintf('%s_%s_rawfilt.mat',sub,ses),'-struct','eeg');
		end
	end
end