function gabormat = makegabor(radiuscm,gaborsizecm,rotation,fixationcm,freqcm,graphit)
%MAKEGABOR - Builds a gabor patch with specific orientation in a circluar
%field
%
%   This function creates an image matrix of a gabor with a specific
%   rotation
%
%Useage: 
%  >> gabormat = makegabor(radius,gaborsize,rotation,fixation,freq,graphit)
%
%Inputs:
%   radius - length of circular field raidus in centimeters 
%                                        (default= 500 pixels)
%
%   gaborsize*gaborsize is the Gabor box size 
%                                        (default = 3*radius/5)
%
%   rotation - in degrees, rotation clockwise from vertical (default = 45)
%
%   fixation - length of fixation spot radius in centimeters
%                                        (default = 10 pixels)
%
%   freq - spatial frequency of gabor in centimeters 
%                                        (default = .5/gaborsize)
%
%   graphit (optional): Graphs the image in grayscale if = 1 (default = 0)
%
%See also MAKEGRATING, MAKENOISE.

%% Record of revisions:
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  10/31/14       Michael Nunez                 Original Code
%  11/11/14       Michael Nunez                Circular display
%  11/17/14       Michael Nunez                 Fixation spot
%  01/21/15       Michael Nunez              Spatial frequency input
%  01/23/15       Michael Nunez    Changed frequency input behav, randphase
%  01/25/15       Michael Nunez              Fixation spot is NaN
%  01/27/14       Michael Nunez             Output standardized to [-1,1]
%                                        Input in cycle per degree
%  05/08/15       Michael Nunez             Fixation spot now -1
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
if nargin < 2 || isempty(gaborsizecm)
    gaborsize = 3*radius/5;
else
    gaborsize = round(ppcm*gaborsizecm);
end
if nargin < 3 || isempty(rotation)
    rotation = 45;
end
if nargin < 4 || isempty(fixationcm)
    fixation = 10;
else
    fixation = round(ppcm*fixationcm);
end
if nargin < 5 || isempty(freqcm)
    freq = .5/gaborsize;
else
    freq = freqcm/ppcm;
end
if nargin < 6 || isempty(graphit)
    graphit = 0;
end

%% Code
halfgabor = floor(gaborsize/2);
bigimage = 2*radius;
whereput = (-radius+1):radius;

%Rotation to radians
%Randomly flip the gabor patch
% randflip = pi*(2-randi(2));
% figangle = rotation*(pi/180) + randflip;
randphase = rand*2*pi;
figangle = rotation*(pi/180);

%Make grids
[xg,yg]=meshgrid(whereput,whereput);

% %Put some random placement into the gabors
% randplace = randi(floor((2*radius-gaborsize)/3)-1,[1,2]).*[cos(90*rand) sin(90*rand)];
% %randplace = (floor((2*radius-gaborsize)/3)-1).*[cos(90*rand) sin(90*rand)];
% %randplace = (floor((2*radius-gaborsize)/3)-1)*ones(1,2);
% randdir = 3 - 2*randi(2,[1,2]);
% %shift = randplace.*randdir;
shift = [0 0];

%%%Draw the figure-gabor
%Rotation matrix on meshgrid & shifts
xp = xg.*cos(figangle) + yg.*sin(figangle) + shift(1);
yp = -xg.*sin(figangle) + yg.*cos(figangle) + shift(2);

justgabor = exp(-((xp/(halfgabor/2)).^2)-(.4*(yp/(halfgabor/2)).^2)).*sin(2*pi*freq*(xp) + randphase);

%% Cut gabormat to circular display

[xc,yc] = meshgrid(1:bigimage);
z = sqrt((xc-radius).^2 + (yc-radius).^2);
justgabor(z > radius) = NaN;

%Standardize to [-1, 1]
gabormat = 2*((justgabor - min(justgabor(:)) )/(max(justgabor(:))-min(justgabor(:))) - .5);

% %Make fixation spot
% gabormat(z < fixation) = -1;
% 
% %Make border ring
% gabormat(z > (radius - fixation) & z <= radius) = -1;

%% Graph the image

if graphit == 1
%figure;
imagesc(gabormat);
colormap(gray);
screensize = get(0,'ScreenSize');
set(gcf,'Position', [1 1 screensize(3) screensize(4)]);
axis('square');
end

