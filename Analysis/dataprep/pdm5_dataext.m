function pdmdata = pdm5_dataext(varargin)
%PDM5_DATAEXT - Extracts cnt data and cleans photocells from PDM experiment4
%
%Usage: pdmdata = pdm5_dataext(nblocks, MaxTrials, MinTrials, cntdirectory);
%
%
%Inputs:
%  nblock: Number of blocks
%  MaxTrials: Maximum number of trials in each block
%  MinTrials: Minimum number of trials in each block
%  cntdirectory: Directory with only 'nblocks' .cnt files corresponding to 
%'nblocks' blocks
%
%Outputs:
% pdmdata: A structure to be used in Cort's functions
%
%To save output (example):
%save s00.mat -struct pdmdata
%
%To load output (example):
%pdmdata = load('s00.mat')

% Copyright (C) 2016 Michael D. Nunez, <mdnunez1@uci.edu>, Cort Horton, Bill Winter
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
%  1/12/16        Michael Nunez             Adapted from pdm4_dataext
%  6/22/16        Michael Nunez                  Clarify filtering

%% Initial
if nargin < 1
    nblocks = 8;
    MaxTrials = 60;
    MinTrials = 60;
    cntdir = pwd;
elseif nargin == 1
    nblocks = varargin{1};
    MaxTrials = 60;
    MinTrials = 60;
    cntdir = pwd;
elseif nargin == 2
    nblocks = varargin{1};
    MaxTrials = varargin{2};
    MinTrials = 60;
    cntdir = pwd;
elseif nargin == 3
    nblocks = varargin{1};
    MaxTrials = varargin{2};
    MinTrials = varargin{3};
    cntdir = pwd;
elseif nargin == 3
    nblocks = varargin{1};
    MaxTrials = varargin{2};
    MinTrials = varargin{3};
    cntdir = varargin{4};
else
    error('Too many inputs.');
end

cnt = dir([cntdir,'/*.RAW']);

if numel(cnt) >= 9
    eyes = 1;
else
    eyes = 0;
end

fprintf('File load order:\n');
for bnum = 1:nblocks
    if eyes
        fprintf('%s ',cnt(bnum+2).name);
    else
        fprintf('%s ',cnt(bnum).name);
    end
end
fprintf('\n');

behavdir = dir([cntdir,'/behav/*.mat']);
behavfile = behavdir.name;
fprintf('Behavior file name: %s\n',behavfile);

pdmdata.expinfo = load([cntdir,'/behav/',behavfile]);
flickerhz = pdmdata.expinfo.flickerhz;
noisehz = pdmdata.expinfo.noisehz;

%% Code
badchans = [44 49 56 63 100 108 114 120 129];

try
    for bnum = 1:nblocks
        disp(['Cleaning block ',num2str(bnum)]);

        if eyes
            [raw, srate, info] = rdata([cntdir,'/',cnt(bnum+2).name],'egi');
        else
            [raw, srate, info] = rdata([cntdir,'/',cnt(bnum).name],'egi');
        end

        %Filtering step ([1 100] bandpass and detrending-baselining) using Cort's filtereeg function
        [tempdata, ~] = rmchannel(raw',badchans); %Remove photocell channels and average reference the data
        block = pdm5_filtereeg(tempdata,srate);
        [block, ~] = rmchannel(block,[]); %Average reference the remaining data
        
        nbounds = round(srate/(noisehz/2) +[-10 10]);
        sbounds = round(srate/flickerhz +[-10 10]);

        %The monitor is drawing from the top, so these three photocells
        %should be slightly delayed from one another
%         pnoise = getspike(block([56 100 114],:),round(srate/noisehz -10)); %Photocell of 30 Hz Noise
%         psignal = getspike(block([63 108 120],:),round(srate/flickerhz -10)); %Photocell of 20 Hz Signal
%         

        spikenoise = getspike(raw([56 114],:),nbounds(1)); %Photocell of 30 Hz Noise
        pnoise = zeros(size(spikenoise));
        pnoise(spikenoise > nbounds(1) & spikenoise < nbounds(2)) = ...
            spikenoise(spikenoise > nbounds(1) & spikenoise < nbounds(2));
        spikesignal = getspike(raw([63 120],:),sbounds(1)); %Photocell of 20 Hz Signal
        psignal = zeros(size(spikesignal));
        psignal(spikesignal > sbounds(1) & spikesignal < sbounds(2)) = ...
            spikesignal(spikesignal > sbounds(1) & spikesignal < sbounds(2));
        
%         pdmdata.pherror.ndiffs = diff(pnoise')';
%         pdmdata.pherror.sdiffs = diff(psignal')';
        
        %Use photocell 1 to find beginning of the trial, note that
        %photocell 1 and photocell 5 are delayed by a few samples
        nonzero = (pnoise(2,:) > 0);
        trialstart = strfind(nonzero,[zeros(1,1500) 1]) + 1500;
        
        %Use photocell 6 to find the beginning of the stimulus, note that
        %photocell 2 and photocell 6 are also delayed
        %nonzero = (psignal(1,:) > 0);
        nonzero = (psignal(2,:) > 0);
        respstart = strfind(nonzero,[zeros(1,1500) 1]) + 1500;
        
        %ntrials = length(trialstart);
        ntrials = length(respstart);
        
        fprintf('%i trials were found in block %i \n',ntrials,bnum);
        
        missing(bnum) = 0;
        if ntrials < MinTrials
            missing(bnum) = MinTrials - ntrials;
            respstart = [nan(1,missing(bnum)) respstart];
            trialstart = [nan(1,missing(bnum)) trialstart];
        end
        if ntrials > MaxTrials
            error(sprintf('Number of trials found is larger the number of trials desired in block %i \n',bnum));
        end
        
        %Put data into equal trial windows (3500 sample points to cover all possible window lengths of
        %.25+1+2+.25 seconds) locked on first SIGNAL photocell and approximately on the last
        %photocell (exact photocell placement differences by one or two samples
        %between trials)
        %Note that there may be some EEG overlap
        resplen = 2250;
        data{bnum} = zeros(resplen+1250,129,ntrials);
        photo{bnum} = zeros(resplen+1250,4,ntrials);
        button{bnum} = zeros(resplen+1250,2,ntrials); %Extract behavioral data
        correct{bnum} = nan(1,ntrials);
        rt{bnum} = nan(1,ntrials);
        flickon{bnum} = nan(1,ntrials);
        for i=(missing(bnum)+1):(ntrials+missing(bnum))
            data{bnum}(:,:,i) = block((respstart(i)-1250):(respstart(i)+resplen-1),:); %Segment 250 ms before trial start
            photo{bnum}(:,:,i) = [spikenoise(:,(respstart(i)-1250):(respstart(i)+resplen-1))' , ...
                spikesignal(:,(respstart(i)-1250):(respstart(i)+resplen-1))'];
            button{bnum}(:,:,i) = [info.event.left((respstart(i)-1250):(respstart(i)+resplen-1))' , ...
                info.event.righ((respstart(i)-1250):(respstart(i)+resplen-1))'];
            %%Get RT and accuracy
            %Correct, incorrect or answered inappropriately?
            thistrial = i + MinTrials*(bnum-1) + missing(bnum);
            if pdmdata.expinfo.highlow(thistrial) == 2 && any(button{bnum}(:,2,i)) && ~any(button{bnum}(:,1,i))
                correct{bnum}(i) = 1;
                tempsamps = find(button{bnum}(:,2,i));
            elseif pdmdata.expinfo.highlow(thistrial) == 1 && any(button{bnum}(:,1,i)) && ~any(button{bnum}(:,2,i))
                correct{bnum}(i) = 1;
                tempsamps = find(button{bnum}(:,1,i));
            elseif pdmdata.expinfo.highlow(thistrial) == 2 && any(button{bnum}(:,1,i)) && ~any(button{bnum}(:,2,i))
                correct{bnum}(i) = 0;
                tempsamps = find(button{bnum}(:,1,i));
            elseif pdmdata.expinfo.highlow(thistrial) == 1 && any(button{bnum}(:,2,i)) && ~any(button{bnum}(:,1,i))
                correct{bnum}(i) = 0;
                tempsamps = find(button{bnum}(:,2,i));
            else
                tempsamps = NaN;
            end
            %%Check this
            %rt{bnum}(i) = tempsamps(1) - respstart(i) + 1;
            rt{bnum}(i) = tempsamps(1) - 1250;
        end
        flickon{bnum} = 1251 - (respstart-trialstart);
        flickon{bnum} = flickon{bnum}((missing(bnum)+1):(ntrials+missing(bnum)));

    end
catch me
    rethrow(me)
    if bnum>1
        disp('Saving all other data...');
        alldata = [];
        allphoto = [];
        allbutton = [];
        allcorrect = [];
        allrt = [];
        allflickon = [];
        for i=1:(bnum-1)
            alldata = cat(3,alldata,data{i});
            allphoto = cat(3,allphoto,photo{i});
            allbutton = cat(3,allbutton,button{i});
            allcorrect = [allcorrect,correct{i}];
            allrt = [allrt,rt{i}];
            allflickon = [allflickon,flickon{i}];
        end
        pdmdata.data = alldata(:,1:129,:);
        pdmdata.photo = allphoto;
        pdmdata.button = allbutton;
        pdmdata.flickon = allflickon;
        pdmdata.expinfo.correct = allcorrect;
        pdmdata.expinfo.rt = allrt;
        pdmdata.sr = srate;
        pdmdata = addhm(pdmdata,'eginn128');
        pdmdata.error = me;
        pdmdata.missing = missing;
    end
    return
end

alldata = [];
allphoto = [];
allbutton = [];
allcorrect = [];
allrt = [];
allflickon = [];
for i=1:nblocks
    alldata = cat(3,alldata,data{i});
    allphoto = cat(3,allphoto,photo{i});
    allbutton = cat(3,allbutton,button{i});
    allcorrect = [allcorrect,correct{i}];
    allrt = [allrt,rt{i}];
    allflickon = [allflickon,flickon{i}];
end

pdmdata.data = alldata(:,1:129,:);
pdmdata.photo = allphoto;
pdmdata.button = allbutton;
pdmdata.flickon = allflickon;
pdmdata.expinfo.correct = allcorrect;
pdmdata.expinfo.rt = allrt;
pdmdata.sr = srate;
pdmdata = addhm(pdmdata,'eginn128');
pdmdata.missing = missing;