%PDM5B_N200LOC - Script that calculates spline-surface Laplacians of average N200 response
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
%  12/12/17        Michael Nunez             Adapted from pdm4_n200loc
%  12/15/17        Michael Nunez           Save out N200chanpos.mat
%  12/29/17        Michael Nunez          Updating new load locations
%  02/27/18 	   Michael Nunez        Save out scalp Current Source Density
%  04/13/18        Michael Nunez            Change percentage valid

% To do:
% 1) Create simulated dura surface
% 2) Find labeled locations
% 3) Create movie with network leadup to peak

% %% Preliminary 

% experiments = {'exp4data/subjects', 'exp5data/subjects/training'};
% sessions = {{'s1', 's2'}, {'ses1', 'ses2', 'ses3', 'ses4', 'ses5', 'ses6', 'ses7'}};
% subjects = {'s59', 's64', 's68', 's80', 's82', 's93', ...
%                        's94', 's95', 's96', 's97', 's100', 's101', 's109', 's110'};
% conds = {'high', 'med', 'low'};

% svdloc = '/home/michael/data10/michael/pdm/%s/%s/%s/erp_svd_%s_%s_v5.mat';
% eegloc = '/home/michael/data10/michael/pdm/%s/%s/%s/%s_cleaned.mat';

% scalpdata = zeros(147,129);

% trackscalp = 1;
% for expe = 1:2,
%     nses = length(sessions{expe});
%     for sub = 1:length(subjects),
%         for ses = 1:nses,
%             tempsvdloc = sprintf(svdloc,experiments{expe},subjects{sub},sessions{expe}{ses},subjects{sub},sessions{expe}{ses});
%             if expe == 1
%                 tempeegloc = sprintf(eegloc,experiments{expe},subjects{sub},sessions{expe}{ses},subjects{sub});
%             else
%                 tempeegloc = sprintf(eegloc,experiments{expe},subjects{sub},sessions{expe}{ses},[subjects{sub},'_',sessions{expe}{ses}]);
%             end
%             if exist(tempsvdloc)==2,
%                 fprintf('Loading EEG data %s \n',tempeegloc);
%                 eeg = load(tempeegloc);
%                 fprintf('Loading SVD data %s \n',tempsvdloc);
%                 svddata = load(tempsvdloc);
%                 [goodchans,goodtrials,badchans,badtrials] = goodbad(eeg); %artscreenEEG function
%                 for cond = 1:3,
%                     this_svd = svddata.(sprintf('svd_rint_%s',conds{cond}));
%                     [~,wheremin] = min(this_svd.u(251:375));
%                     wheremin = wheremin + 250;
%                     tempN200 = this_svd.u(wheremin,1)*this_svd.v(1,:)*this_svd.s(1);
%                     % Calculate spline interpolated data
%                     chanpos = eeg.hm.Electrode.CoordOnSphere;
%                     splinemat = splineinterp(.1, chanpos, chanpos(goodchans,:));
%                     N200projected = tempN200*splinemat;
%                     N200projected(eeg.hm.Electrode.NoInterp) = NaN;
%                     scalpdata(trackscalp, :) = N200projected;
%                     trackscalp = trackscalp + 1;
%                     % lapmat = splinenlap(.1, chanpos, chanpos);
%                     % N200projectedsph = N200projected*lapmat;
%                     % N200projectedsph(eeg.hm.Electrode.NoInterp) = NaN;
%                 end
%             else
%                 fprintf('Did not find file %s \n',tempsvdloc);
%             end
%         end
%     end
% end

% save('N200chanpos.mat','scalpdata');

% Experiment 1: 21 Unique EEG sessions
% Experiment 2: 28 Unique EEG sessions

fprintf('Loading previously generated N200 scalp topographies...\n');
load('N200chanpos.mat');

fprintf('Loading head models...\n');
load ~/data10/finalheadmodel/HNL128_ELECTRODE.mat %This contains bad channels for Laplacian
load ~/data10/michael/headmodels/s05.mat
load eginn128hm.mat

Electrode.badchan.eeg = union([19],Electrode.badchan.eeg); %Add dead channel 19 to bad channel list
goodchan = setdiff(1:size(Electrode.Coordinate,1),Electrode.badchan.eeg);

Electrode = makegeolap(Electrode,Scalp); %Get Laplacian matrix
Electrode.directspline3D = splineinterp(.1,Electrode.Coordinate(goodchan,:),Scalp.Vertex); %Constant is based on Bill's topo.m
Electrode.spline3D = splineinterp(.1,Electrode.Coordinate,Scalp.Vertex);

lapdata = [scalpdata(:,goodchan)*Electrode.lap' zeros(size(scalpdata,1),1)];

%Plot the mean response with Cort Horton's function
fprintf('Plotting the mean response on a 2D scalp...\n');
figure; hortontopo(mean(scalpdata,1),EGINN128,'channumbers',0,'drawelectrodes',0,...
    'chanfontsize',8,'cmap','jet','badchans',[Electrode.badchan.eeg 129]);

%Plot the Laplacian with Cort Horton's function
fprintf('Plotting the mean Laplacian on a 2D scalp...\n');
figure; hortontopo(mean(lapdata,1),EGINN128,'channumbers',0,'drawelectrodes',0,...
    'chanfontsize',8,'cmap','jet','badchans',[129]);

%Plot a random sample of single subjects
randsamp = randperm(size(scalpdata,1),9);

fprintf('Plotting the a random sample of EEG topographies and Laplacians on a 2D scalp...\n');
eegplot = figure;
lapplot = figure;
for n=1:9
figure(eegplot);
subplot(3,3,n);
hortontopo(scalpdata(randsamp(n),:),EGINN128,'channumbers',0,'drawelectrodes',0,...
'chanfontsize',8,'cmap','jet','badchans',[Electrode.badchan.eeg 129]);
figure(lapplot);
subplot(3,3,n);
hortontopo(lapdata(randsamp(n),:),EGINN128,'channumbers',0,'drawelectrodes',0,...
'chanfontsize',8,'cmap','jet','badchans',[129]);
end

%% Average plots

%Plot the mean response
fprintf('Plotting the mean response on a 3D scalp...\n');
scalpPotential = mean(scalpdata(:,goodchan),1)*Electrode.directspline3D';

%Zero out periphery values
u1=find(Scalp.Vertex(:,3)<-40);
u3=find(Scalp.Vertex(:,2)>39 & Scalp.Vertex(:,3)<14);
scalpPotential(u1)=0;
scalpPotential(u3)=0;

f1 = figure('units','normalized','outerposition',[0 0 1 1]);
h = drawmesh(Scalp);
setmesh(h,'interp',scalpPotential);
% set(gca,'CLim', [-.0025 .0025]);
colormap jet
alpha(.8);
h2 = drawmesh(Brain);
set(h2,'facecolor',[1 1 1],'edgecolor',[0 0 0]);
%Change viewpoint (azimuth and elevation) of the figure
view(gca,0,0);

fprintf('Saving the scalp potential...\n');
save('ScalpN200Potential.mat','scalpPotential');



%Plot the mean Laplacian
fprintf('Plotting the Laplacian response on a 3D scalp...\n');
scalpCSD = mean(lapdata(:,1:128),1)*Electrode.spline3D';

%Zero out periphery values
u1=find(Scalp.Vertex(:,3)<-40);
u3=find(Scalp.Vertex(:,2)>39 & Scalp.Vertex(:,3)<14);
scalpCSD(u1)=0;
scalpCSD(u3)=0;

f2 = figure('units','normalized','outerposition',[0 0 1 1]);
h = drawmesh(Scalp);
setmesh(h,'interp',scalpCSD);
colormap jet
alpha(.8);
h2 = drawmesh(Brain);
set(h2,'facecolor',[1 1 1],'edgecolor',[0 0 0]);
%Change viewpoint (azimuth and elevation) of the subplots
view(gca,0,0);

fprintf('Saving the scalp CSD...\n');
save('ScalpN200CSD.mat','scalpCSD');

%Plot the mean Laplacian
fprintf('Plotting the Laplacian response on a simulated dura mater surface...\n');
Dura = Scalp;
Dura.Vertex = Scalp.Vertex*(9/10);
mask = mean(lapdata(:,1:128),1)*Electrode.spline3D';

%Zero out periphery values
u1=find(Scalp.Vertex(:,3)<-40);
u3=find(Scalp.Vertex(:,2)>39 & Scalp.Vertex(:,3)<14);
mask(u1)=NaN;
mask(u3)=NaN;

f2 = figure('units','normalized','outerposition',[0 0 1 1]);
h = drawmesh(Dura);
setmesh(h,'interp',mask);
colormap jet
alpha(.8);
h2 = drawmesh(Brain);
set(h2,'facecolor',[1 1 1],'edgecolor',[0 0 0]);
%Change viewpoint (azimuth and elevation) of the subplots
view(gca,0,0);

%% Inverse model

%Find vertices of gain matrix
theseverts = [];
for e=1:128
    [~,minind] = min(sqrt(sum((ones(10241,1)*Electrode.Coordinate(e,:) - Scalp.Vertex).^2,2))); %Min of norm of difference
    theseverts = [theseverts minind];
end
Electrode.newGain = Gain(theseverts,:);

Electrode.LapGain = Electrode.lap'*Electrode.newGain;
%Use a small regularization parameter. Larger eigenvalues correspond to components that cover more space because they explain a larger percentage of the data?
[inverse, stat] = inversemodel(Electrode.LapGain,'prctile',40);

Electrode.lapinverse = inverse;

inversesolution = mean(scalpdata(:,goodchan),1)*Electrode.lapinverse';

fprintf('Saving the inverse solution...\n');
save('InverseN200Laplacian.mat','inversesolution');

fprintf('Generating one possible cortical surface projected Laplacian...\n');
f3 = figure('units','normalized','outerposition',[0 0 1 1]);
h = drawmesh(Brain);
setmesh(h,'interp',inversesolution);
set(gca,'CLim', [-20 20]);
colormap jet
%h2 = drawmesh(Scalp);
%setmesh(h2,'glassy');
%alpha(h2,.1);
%Change viewpoint (azimuth and elevation) of the subplots
view(gca,0,0);

% Find coritcal locations, these locations were found with point and click methods
labelindex = (inversesolution < -10);
theselocs = Brain.Vertex(inversesolution < -10,:);


labmask = zeros(1,size(Brain.Vertex,1));
labmask(labelindex) = 1;

fprintf('Generating figure with labeled cortical areas...\n');
f4 = figure('units','normalized','outerposition',[0 0 1 1]);
h = drawmesh(Brain);
setmesh(h,'interp',labmask);
colormap jet
%h2 = drawmesh(Scalp);
%setmesh(h2,'glassy');
%alpha(h2,.1);
%Change viewpoint (azimuth and elevation) of the subplots
view(gca,0,0);

locinds = [];
for l=1:size(theselocs,1)
    [~,locind] = min(sqrt(sum((ones(size(Brain.Vertex,1),1)*theselocs(l,:) - Brain.Vertex).^2,2)));
    locinds(l) = locind;
end

foundlocs = unique(Dest.lab(Dest.ind(locinds)));
fprintf('Locations found from averaged localized Laplacian:\n');
foundlocs(:)

% Find coritcal locations by performing localization for each EEG session and noise condition
fprintf('Finding coritcal locations by performing localization for each EEG session and noise condition...\n')
allsubfoundlocs = cell(0);
allsubfoundmins = cell(0);
for s=1:size(scalpdata,1)
	localizedvals = scalpdata(s,goodchan)*Electrode.lapinverse';
	theselocs = Brain.Vertex(localizedvals < -10,:);
	minlocs = Brain.Vertex(argmin(localizedvals), :);

	locinds = [];
	for l=1:size(theselocs,1)
	    [~,locind] = min(sqrt(sum((ones(size(Brain.Vertex,1),1)*theselocs(l,:) - Brain.Vertex).^2,2)));
	    locinds(l) = locind;
	end
	tempfoundlocs = unique(Dest.lab(Dest.ind(locinds)));
	allsubfoundlocs = {allsubfoundlocs{:} tempfoundlocs{:}};

	minlocinds = [];
	for l=1:size(minlocs,1)
	    [~,locind] = min(sqrt(sum((ones(size(Brain.Vertex,1),1)*minlocs(l,:) - Brain.Vertex).^2,2)));
	    minlocinds(l) = locind;
	end

	tempfoundmins = unique(Dest.lab(Dest.ind(minlocinds)));
	allsubfoundmins = {allsubfoundmins{:} tempfoundmins{:}};
end

uniquefoundlocs = unique(allsubfoundlocs);
totaloffoundlocs = zeros(1,length(uniquefoundlocs));
for l=1:length(uniquefoundlocs),
	totaloffoundlocs(l) = sum(strcmp(allsubfoundlocs,uniquefoundlocs{l}));
end

uniquefoundmins = unique(allsubfoundmins);
totaloffoundmins = zeros(1,length(uniquefoundmins));
for l=1:length(uniquefoundmins),
	totaloffoundmins(l) = sum(strcmp(allsubfoundmins,uniquefoundmins{l}));
end


percentagevalid = .7; % 1 corresponds to 100 percent
foundlocs2 = uniquefoundlocs(find(totaloffoundlocs > percentagevalid*size(scalpdata,1)));
fprintf('All locations found from localized Laplacian of each EEG session and condition with cutoff of %d %%:\n',percentagevalid*100);
foundlocs2(:)
[sortedcounts,sortbytotal] = sort(totaloffoundmins,2,'descend');
for c=1:length(uniquefoundmins),
	fprintf('%s was the minimum of %.1f%% observations \n',uniquefoundmins{sortbytotal(c)},(sortedcounts(c)/size(scalpdata,1))*100);
end