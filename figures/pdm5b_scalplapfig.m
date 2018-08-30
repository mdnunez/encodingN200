function pdm5b_scalplapfig(sp)

%pdm5b_scalplapfig - Creates ERP Laplacian images on MNI scalp 
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
%  08/28/18        Michael Nunez             Adapted from pdm5b_corticallapfig.m


%% Code

setpainters;
load ERPLap.mat;
load ERPchanpos.mat;
load ~/data10/michael/headmodels/s05.mat;
load ~/data10/finalheadmodel/HNL128_ELECTRODE.mat; %This contains bad channels for Laplacian

plotindex = 101;
timeindex = 181 +100;
viewpos = {{-60,8},{60,8},{0,0}};

temperp = squeeze(scalperp(plotindex,timeindex,1:128))';
goodchan = setdiff(1:128,find(isnan(temperp)));
Electrode.directspline3D = splineinterp(.1,Electrode.Coordinate(goodchan,:),Scalp.Vertex);
erpCSD = temperp(goodchan)*Electrode.directspline3D';


Electrode.spline3D = splineinterp(.1,Electrode.Coordinate,Scalp.Vertex);
lapCSD = squeeze(laperp(plotindex,timeindex,1:128))'*Electrode.spline3D';

%Zero out periphery values for the Laplacian
u1=find(Scalp.Vertex(:,3)<-40);
u3=find(Scalp.Vertex(:,2)>39 & Scalp.Vertex(:,3)<14);
lapCSD(u1)=0;
lapCSD(u3)=0;

% f1=figure('units','normalized','outerposition',[0 0 1 1]);
% h1 = drawmesh(Scalp);
% setmesh(h1,'interp',erpCSD);
% set(gca,'CLim', [-60 60]);
% newmap = (bone + jet)/2;
% colormap(newmap);
% view(gca,viewpos{sp}{:});

f2=figure('units','normalized','outerposition',[0 0 1 1]);
h2 = drawmesh(Scalp);
setmesh(h2,'interp',lapCSD);
set(gca,'CLim', [-.1 .1]);
newmap = (bone + jet)/2;
colormap(newmap);
view(gca,viewpos{sp}{:});
