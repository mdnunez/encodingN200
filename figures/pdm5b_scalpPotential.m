function pdm5b_scalpPotential(sp)

%pdm5b_scalpPotential - Creates scalp potential images
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
%  04/13/18        Michael Nunez             Adapted from pdm5b_corticallapfig.m


%% Code

fontsize = 40;

setpainters;
load ScalpN200Potential.mat;
load ~/data10/michael/headmodels/s05.mat;

muVtoV = 1e-6; %Conversion from muVs tp Vs
%Note that the Gain matrix is amps/mm^2

viewpos = {{-60,8},{60,8},{0,0}};


f2=figure('units','normalized','outerposition',[0 0 1 1]);
h2 = drawmesh(Scalp);
setmesh(h2,'interp',scalpPotential);
set(gca,'CLim', [-30 30]);
newmap = (bone + jet)/2;
colormap(newmap);
view(gca,viewpos{sp}{:});
