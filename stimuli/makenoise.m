function noisemat = makenoise(radiuscm,numsine,stdpasscm,fixationcm,freqbndscm,graphit)
%MAKENOISE - Builds a circular visual noise patch
%
%
%Useage: 
%  >> noisemat = makenoise(radius,freqbnds,rotation,fixation,,graphit)
%
%Inputs:
%   radius - length of circular field raidus in centimeters 
%                                        (default= 500 pixels)
%
%   numsine - number of frequency bands to generate noise (default = 2)
%
%   stdpass - standard deviation of bandpass Guassian filters 
%                                             (default = .01 cycle/pixel)
%
%   fixation - length of fixation spot radius in centimeters 
%                                        (default = 10 pixels)
%
%   freqbnds - spatial frequencies for random noise generation in 
%                cycle/cm == cycle/degree visual angle at 57 cm
%               (should be same length as numsine; 
%                 default = uniform draw from [1 3] cycle/cm)
%
%   graphit (optional): Graphs the image in grayscale if = 1 (default = 0)
%
%See also MAKEGABOR, MAKEFIXATION.

%% Record of revisions:
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  11/11/14       Michael Nunez                 Original Code
%  11/17/14       Michael Nunez                 Fixation spot
%                                               Pixel density
%  01/25/15       Michael Nunez       Noise now sum of sinusoids + Guassian
%                                               Different inputs
%  01/27/15       Michael Nunez             Output standardized to [-1,1]
%                                        Input in cycle per degree
%  02/06/15       "                 Frequency inputs same length as numsine 
%                                    Changed type of noise generation
%  02/10/15       "                 Changed defaults
%  02/12/15       "
%  04/09/15       Michael Nunez             Fixation spot now zero
%  04/29/15       Michael Nunez             Fixation spot now -1
%  01/22/16       Michael Nunez                 Remove fixation

% Copyright (C) 2016 Michael D. Nunez, <mdnunez1@uci.edu>
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

%% Arguments

%Size inputs in degrees assume a 24" == 60.96 cm at 1920 * 1280 resolution
ppcm = 37.8536;

if nargin < 1 || isempty(radiuscm)
    radius = 500;
else
    radius = round(ppcm*radiuscm);
end
if nargin < 2 || isempty(numsine)
    numsine = 2;
end
if nargin < 3 || isempty(stdpasscm)
    stdpass = .01;
else
    stdpass = stdpasscm/ppcm;
end
if nargin < 4 || isempty(fixationcm)
    fixation = 10;
else
    fixation = round(ppcm*fixationcm);
end
if nargin < 5 || isempty(freqbndscm)
    freqs = (1 + 2*rand(numsine,1))/ppcm;
else
    if length(freqbndscm) ~= numsine
       error('The number of frequency bands should equal to the number of sinusoids.');
    end
    freqs = freqbndscm/ppcm; %Convert to cycles per pixel
end
if nargin < 6 || isempty(graphit)
    graphit = 0;
end

%% Code
bigimage = 2*radius;

whereput = (-radius+1):radius;
%Make grids
[xg,yg]=meshgrid(whereput,whereput);

% %Guassian random noise
% unfiltnoise = randn(bigimage);
% %Fourier transform of noise
% fouriermat = fftshift(fft2(unfiltnoise - mean2(unfiltnoise))); %cycles per radius?

%Guassian random noise
unfiltnoise = randn(bigimage);
%Fourier transform of noise
fouriermat = fftshift(fft2(unfiltnoise - mean2(unfiltnoise))); %cycles per radius?

xcmgrid = (xg./radius); %Cycles per pixel
ycmgrid = (yg./radius); %Cycles per pixel
radii = sqrt(xcmgrid.^2 + ycmgrid.^2); %radii in cycles per pixel
bandpass = zeros(bigimage);
for m=1:numsine
    bandpass = bandpass + exp(-((radii-2*freqs(m)).^2)/2/stdpass^2);
end

%Bandpass spectrum
specmat = fouriermat.*bandpass;

%Bandpass filtered noise matrix
noisemat = real(ifft2(fftshift(specmat) + mean2(unfiltnoise)));


%% Cut noisemat to circular display

[xc,yc] = meshgrid(1:bigimage);
z = sqrt((xc-radius).^2 + (yc-radius).^2);

noisemat(z > radius) = NaN;

%Standardize to [-1, 1]
noisemat = 2*((noisemat - min(noisemat(:)) )/(max(noisemat(:))-min(noisemat(:))) - .5);

% %Make fixation spot
% noisemat(z < fixation) = -1;
% 
% %Make border ring
% noisemat(z > (radius - fixation) & z <= radius) = -1;

%% Graph the image

if graphit == 1
figure;
imagesc(noisemat);
colormap(gray);
screensize = get(0,'ScreenSize');
set(gcf,'Position', [1 1 screensize(3) screensize(4)]);
axis('square');
end