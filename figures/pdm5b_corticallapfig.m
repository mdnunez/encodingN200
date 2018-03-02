function pdm5b_corticallapfig(sp)

%pdm5b_corticallapfig - Creates cortical Laplacian images on brains 
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
%  12/29/17        Michael Nunez             Adapted from pdm3b_corticallapfig.m
%  02/27/18 	   Michael Nunez           Plot scalp Current Source Density


%% Code

fontsize = 40;

setpainters;
load ScalpN200CSD.mat;
load InverseN200Laplacian.mat;
load ~/data10/michael/headmodels/s05.mat;

muVtoV = 1e-6; %Conversion from muVs tp Vs
%Note that the Gain matrix is amps/mm^2

viewpos = {{-60,8},{60,8},{0,0}};

% f1 = figure('units','normalized','outerposition',[0 0 1 1]);
% h = drawmesh(Brain);
% setmesh(h,'interp',inversesolution);
% set(gca,'CLim', [-20 20]); %Note that this scale is muAmps because the muVtoV was not applied
% newmap = (bone + jet)/2;
% colormap(newmap);

% %Color bar
% cb = colorbar('EastOutside','Fontsize',fontsize);
% yl = ylabel(cb,'\muA/mm^2','Fontsize',fontsize);
% set(yl,'Rotation',270);

% %h2 = drawmesh(Scalp);
% %setmesh(h2,'glassy');
% %alpha(h2,.1);
% %Change viewpoint (azimuth and elevation) of the subplots
% view(gca,viewpos{sp}{:});


f2=figure('units','normalized','outerposition',[0 0 1 1]);
h2 = drawmesh(Scalp);
setmesh(h2,'interp',scalpCSD);
set(gca,'CLim', [-.03 .03]); %Note that this scale is muAmps because the muVtoV was not applied
newmap = (bone + jet)/2;
colormap(newmap);
view(gca,viewpos{sp}{:});


%% Save figure

% export_fig(f1,sprintf('corticalLap_%i',sp),'-opengl','-jpg','-r200');