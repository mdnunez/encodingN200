function fixmat = makefixation(radiuscm,fixationcm,graphit)
%MAKEFIXATION - Builds fixation dot
%
%Useage: 
%  >> fixmat = makefixation(radius,fixation,graphit)
%
%Inputs:
%   radius - large radius to match makenoise and makegabor functions, in
%   centimeters (default = 500 pixels)
%
%   fixation - length of fixation spot radius in centimeters 
%                                        (default = 10 pixels)
%
%   graphit (optional): Graphs the image in grayscale if = 1 (default = 0)
%
%SEE ALSO: MAKENOISE, MAKEGABOR

%% Record of revisions:
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  11/14/14       Michael Nunez                 Original Code
%  02/10/15       Michael Nunez            Input in cycle per degree
%  04/09/15       Michael Nunez             Fixation spot now zero
%  05/28/15       Michael Nunez             Fixed input arguments

% Copyright (C) 2015 Michael D. Nunez, <mdnunez1@uci.edu>
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
if nargin < 2 || isempty(fixationcm)
    fixation = 10;
else
    fixation = round(ppcm*fixationcm);
end
if nargin < 3 || isempty(graphit)
    graphit = 0;
end

%% Code
bigimage = 2*radius;

fixmat = ones(bigimage);

%% Cut fixmat to circular display

[xc,yc] = meshgrid(1:bigimage);
z = sqrt((xc-radius).^2 + (yc-radius).^2);
fixmat(z > radius) = NaN;

%Make fixation spot
fixmat(z < fixation) = 0;

%% Graph the image

if graphit == 1
%figure;
imagesc(fixmat);
colormap(gray);
screensize = get(0,'ScreenSize');
set(gcf,'Position', [1 1 screensize(3) screensize(4)]);
axis('square');
end