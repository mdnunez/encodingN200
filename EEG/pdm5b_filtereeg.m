function data=pdm5b_filtereeg(data,sr)
% function fdata=pdm5b_filtereeg(data,sr)
%
%Filter specifications for experiment type exp5
%
%   data - time by chan (by trial)
%   sr - sampling rate


% Copyright (C) 2017 Michael D. Nunez, <mdnunez1@uci.edu>
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

%% Code

% Detrend the data
data=ndetrend(data);

% Highpass the data
Fstop = 0.25;        % Stopband Frequency
Fpass = 1;           % Passband Frequency
Astop = 10;          % Stopband Attenuation (dB)
Apass = 1;           % Passband Ripple (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, sr);
Hd = design(h, 'butter', 'MatchExactly', match);

% Carry out the filtering
data=filter(Hd,data);
data=flipdim(data,1);
data=filter(Hd,data);
data=flipdim(data,1);

% Notch at 60 Hz
Fpass1 = 59;          % First Passband Frequency
Fstop1 = 59.9;        % First Stopband Frequency
Fstop2 = 60.1;        % Second Stopband Frequency
Fpass2 = 61;          % Second Passband Frequency
Apass1 = 1;           % First Passband Ripple (dB)
Astop  = 10;          % Stopband Attenuation (dB)
Apass2 = 1;           % Second Passband Ripple (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, sr);
Hd = design(h, 'butter', 'MatchExactly', match);

% Carry out the filtering
data=filter(Hd,data);
data=flipdim(data,1);
data=filter(Hd,data);
data=flipdim(data,1);

% Lowpass the data
Fpass = 90;          % Passband Frequency
Fstop = 100;          % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 10;          % Stopband Attenuation (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, sr);
Hd = design(h, 'butter', 'MatchExactly', match);

% Carry out the filtering
data=filter(Hd,data);
data=flipdim(data,1);
data=filter(Hd,data);
data=flipdim(data,1);
