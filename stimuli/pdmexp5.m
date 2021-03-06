function pdmexp5(demotest)
%% Record of revisions:
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  10/28/15       Michael Nunez               Converted from pdmexp4
%                    Training session blocked, truncated normal spf draws
%  12/18/15       Michael Nunez      Noise is now has two spatial
%                                  frequencies that mask the signal equally
%  01/06/16       Michael Nunez            Cap size text box
%  01/07/16       Michael Nunez  Block appropriately depending upon session
%                                        Fixed signal spatial frequencies
%  01/11/16       Michael Nunez    Changed timing during messages slightly
%                                           Changed save location
%  01/12/16       Michael Nunez     Changed slum and lower bound of SNR
%  01/21/16       Michael Nunez    Use "makenoise" instead of "makepixeled"
%                        Signal luminance changes based on lower SNR bound
%  01/22/16       Michael Nunez    Default SNR change, fixation radius=.1

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

%% Initial PTB3 Code

PsychJavaTrouble;   %make GetChar work (hack fix; needed?)

AssertOpenGL; %Issue warning if PTB3 with non-openGL used

%InitializePsychSound; % Initialize the sound driver

% if strcmp(whichcmp,'k')
%     Screen('Preference', 'SkipSyncTests', 1);
% end

if ~IsLinux
    error('This program was written to run on Ubuntu Linux.');
end

if nargin < 1
    rundemo = 0;
    runtest = 0;
elseif demotest == 1
    rundemo = 1;
    runtest = 0;   
else
    rundemo = 1;
    runtest = 1;   
end

if runtest
    devs(1).name = 'This is a test';
else
    %Find port for reponse box
    devs = dir('/dev/ttyUSB*');
    if isempty(devs)
        error('The Cedrus Response Box could not be found!');
    end
end
    
if length(devs) > 1
    !dmesg | grep tty
    warning('There are multiple devices that could be the Cedrus Response Box!\n Find "ttyUSB*" in the above output on the same line as "FTDI USB Serial Device"');
end

if nargin > 1
    error('Too many function inputs.');
end
%% Experimenter Prompt

%Inputs Prompt and Output Setup
%Experimenter Prompt
Screenres = get(0,'Screensize');

prompt1={'Subject Number (must begin with letter):','Session Number:','Window Pointer:',...
    'Screen Length (x-axis):','Screen Width (y-axis):','Refresh Rate (fps):',...
    'Noise Frequency:','Gabor Flicker Frequency:','Number of Blocks:','Trials per Block:','SNR:','Cedrus Port [ttyUSB0 ...]:'};
if rundemo && runtest
    def1={'SZZ_test','100','0',num2str(Screenres(3)),num2str(Screenres(4)),'120','40','30','1','6','[2 1 .5]',sprintf('/dev/%s',devs(1).name)};
    studytitle='PDM Experiment 4b Test';
elseif rundemo
    def1={'SZZ_test','1','0',num2str(Screenres(3)),num2str(Screenres(4)),'120','40','30','8','6','[2 1 .5]',sprintf('/dev/%s',devs(1).name)};
    studytitle='PDM Experiment 4b DEMO';
else
    def1={'SZZ_test','','0',num2str(Screenres(3)),num2str(Screenres(4)),'120','40','30','8','60','[2 1 .5]',sprintf('/dev/%s',devs(1).name)};
    studytitle='PDM Experiment 4b';
end

lineNo=1;
answer=inputdlg(prompt1,studytitle,lineNo,def1);
%Subject Number
subnum = answer{1};
%ExpSession Number
sesnum = str2num(answer{2});
if isempty(sesnum)
    error('Please enter an appropriate session number (1 or greater)!');
end
output.sesnum = sesnum;
%Window Pointer / 'Home Screen'.  0 - the primary monitor; 1 - the secondary monitor.
whichScreen = str2num(answer{3});
%Screen resolution on the x-axis
xres = str2num(answer{4});
output.xres = xres;
%Screen resolution on the y-axis
yres = str2num(answer{5});
output.yres = yres;
%This should be the same as the Refresh Rate shown in the Display
%Properties on the computer.  Always check before running the experiment to
%match flicker frequency.
%This code is currently set up to only handle multiples of 60 fps.
refrate = str2num(answer{6});
realrefrate = Screen(0,'FrameRate');
if refrate ~= Screen(0,'FrameRate')
    error(['The real screen refresh rate is set to ',num2str(realrefrate),...
       'Hz while the proposed screen refresh rate is ',num2str(refrate),'Hz.']);
end
output.refrate = refrate;
%Noise frequency (Hz)
noisehz = str2num(answer{7});
output.noisehz = noisehz;
if round(refrate/noisehz) ~= refrate/noisehz
    error('The noise frequency should be divisor of the refresh rate.');
end

%Gabor flicker frequency (Hz)
flickerhz = str2num(answer{8});
output.flickerhz = flickerhz;
if round(refrate/2/flickerhz) ~= refrate/2/flickerhz
    error('The gabor flicker frequency should be divisor of half the refresh rate.');
end

%Number of blocks
blocknum = str2num(answer{9});

%Trials per block
output.tperb = str2num(answer{10});
if round(output.tperb/2) ~= output.tperb/2
    error('There should be an even number of trials per block.');
end

%Number of Trials
trialnum = blocknum*output.tperb;

%Initialize block
block = 1;

%Noise levels
snrs = str2num(answer{11});
output.snrs = snrs;

%signal luminance
lowboundsnr = min(snrs);
slum = lowboundsnr/(1 + lowboundsnr);

%Gabor spatial frequencies (cycles per degree at 57 cm)
gaborspf = [2.4 2.6];

%Noise spatial frequency (cycles per degree at 57 cm)
noisespf = 10;

%Radius of fixation spot (degrees visual angle at 57 cm)
fixrad = .10;

%Number of trials should be a multiple of the number of cells
ncells = length(snrs);
if round(trialnum/ncells) ~= trialnum/ncells
    error('The number of trials should be a multiple of %i.\n',ncells);
end

%Cedrus Handle
if ~runtest
    cport = answer{12};
    try
        chandle = CedrusResponseBox('Open',cport);
    catch me1
        rethrow(me1);
        fprintf('Cedrus port may need a ''chmod 777'' if you''re getting permission issues.\n');
    end
end

%% Subject Prompt
prompt2={'What is your gender?',...
    'Age:','Do you consider yourself right handed, left handed, or both?  (''r'',''l'', or''b'')',...
    'Visual acuity result? (''20/30'' or ''20/35'' or ''20/20'')',...
    'What is your EEG cap size? (''Large'', ''Medium'', or ''Small'')',...
    'Do you have any personal or family history of epilepsy? (''y'' or ''n'')'
    };
demographtitle='Subject Demographics';
lineNo=1;
subdemo=inputdlg(prompt2,demographtitle,lineNo);
switch subdemo{6}
    case 'n'
    otherwise
    error('You have indicated that you have a personal or family history of epilepsy. This experiment involves a fast flickering image. It is recommended that you NOT participate in this study due to a possible risk of seizure.  Please discuss your options with the experimenters.');
end
output.gender = subdemo{1};
output.age = str2num(subdemo{2});
output.hand = subdemo{3};
output.vision = subdemo{4};
output.capsize = subdemo{5};

%% Code

%Get date and time that the session begins
output.date = date;
output.start_time = clock;
    
%number of rows and columns of image
nCols = 1000;
nRows = 1000;

%Initialize estimated accuracy vector, for speed
estcorrect = zeros(1,trialnum);

%Keyboard keypress variables
advancechar = ' ';
escapechar = 27;

%Colors
%txtcolor = [125 125 255]; %Light blue
txtcolor = round([0 .6 .6]*255); %Teal
black = [0 0 0];
white = [255 255 255];
gray = 255*sqrt([.5 .5 .5]);
%darkgray = 255*sqrt([.25 .25 .25]);
blackwhite{1} = black;
blackwhite{2} = white;

% Load fonts
myfont = '-bitstream-courier 14 pitch-bold-i-normal--0-0-0-0-m-0-ascii-0';
fontsize = 26;

%Define photocell placement
% photorect = [0 0 100 90];
% Photo{2} = [0 0 100 90];
% Photo{1} = [0 180 100 270];
k = 0;
photorect = [0 0 100 90];
for m = 1:6
    pRect(m,:) = CenterRectOnPoint(photorect,50,50+k);
    k = k + 160;
end


%Load sounds

%Sound is public domain: http://creativecommons.org/publicdomain/zero/1.0/
%Sound found on Freesound.org, made by "pan14" 
[ygood, Fsgood] = audioread('263133_tone-beep.wav');

%Sound under the Attribution creative commons license: http://creativecommons.org/licenses/by/3.0/
%Sound found on Freesound.org, made by "Autistic Lucario" 
[ybad,  Fsbad ] = audioread('142608_error.wav');

goodsound = audioplayer(ygood,Fsgood);
badsound = audioplayer(ybad,Fsbad);

%Flush Cedrus Events
if ~runtest
    CedrusResponseBox('FlushEvents',chandle);
end


%The following TRY, CATCH, END statement ends psychtoolbox if an error
%occurs
try
    %Open a fullscreen window on the first screen with black background
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'UseVirtualFramebuffer');
    %PsychImaging('AddTask','FinalFormatting','DisplayColorCorrection','SimpleGamma');
    [wptr,windowRect] = PsychImaging('OpenWindow', whichScreen,gray);
    %PsychColorCorrection('SetEncodingGamma',win,1/displayGamma);
    PsychGPUControl('FullScreenWindowDisablesCompositor', 1);

    %sets size of gabor field that will be pasted onto Screen
    imageRect=SetRect(0,0,nCols,nRows);
    destRect=CenterRect(imageRect,windowRect);

    %Creates a window of a black screen with gray circle and fixation spot
    fiximage = makefixation([],fixrad);
    
    
    fiximage(fiximage == 1) = gray(1); %Gamma correction, gray
    fiximage(isnan(fiximage)) = gray(1);
    fiximage(fiximage == 0) = black(1);
    fixglind = Screen('MakeTexture',wptr,fiximage);
    
    %This vector defines the noise frequency for our image
    noiseflic = [];
    for i=1:ceil(4*noisehz)
        noiseflic = [noiseflic 1 zeros(1,(round(refrate/noisehz)- 1))];
    end
    noiseonfind = find(noiseflic);
    
    %This vector defines the Gabor flicker frequency for our image
    gaborflic = [];
    for i=1:ceil(4*flickerhz)
        gaborflic = [gaborflic 2*ones(1,round(refrate/2/flickerhz)) ones(1,round(refrate/2/flickerhz))];
    end
    
    %Set seed based on the time. Backwards compatible with older MATLAB
    %versions
    output.seed = round(sum(100*clock));
    rng('default');
    if runtest
        rng(15);
    else
        rng(output.seed);
    end
    
    %Gabor will not be shown for the first 500ms to 1000ms of the trial
    numframes = round(refrate/noisehz);
    minframe = round(.5*refrate/numframes)*numframes;
    maxframe = round(refrate/numframes)*numframes;
    posframes = minframe:numframes:maxframe;
    trialnframes = posframes(randi(length(posframes),1,trialnum));
    output.noisetimes = trialnframes/refrate;
    
    %Inter-trial interval, 1500ms to 2000ms
    output.intertrial = 1.5 + rand(1,trialnum)*.5;
    
    
    %Define SNR vector, ensure even cell counts
    snrvec = [];
    for b=1:blocknum
        %Show Mix, High, Med, Low, High, Med, Low, Mix if it is the training session
        if sesnum == 1
            blkvec = [];
            if b==1 || b==8
                for n=1:length(snrs)
                    blkvec = [blkvec snrs(n)*ones(1,output.tperb/length(snrs))];
                end
            else
                snrindx = mod(b-1,3)+3*(mod(b-1,3)==0);
                blkvec = [snrs(snrindx)*ones(1,output.tperb)];
            end
        else
            blkvec = [];
            for n=1:length(snrs)
                blkvec = [blkvec snrs(n)*ones(1,output.tperb/length(snrs))];
            end
        end
        snrvec = [snrvec blkvec(randperm(numel(blkvec)))];
    end
    output.snrvec = snrvec;
    
    cut = 0; %Counter for ESC
    
    %Calculate the number of frames in a cycle of an image flicker
    numCycleFrames = trialnframes + ceil(refrate*1.2) + ceil(refrate*rand(1,trialnum)*.8);

    %Output stimulus display time in seconds
    output.stimtime = (numCycleFrames/refrate);
    
    %Initialize recording of trialflic
    output.trialflic = cell(1,trialnum);
    
    output.nlum = slum./snrvec; %noise luminance
    
    Screen('TextFont',wptr,'Arial');
    Screen('TextSize',wptr,18);
    ShowCursor(0);	% arrow cursor
    sessiontext = 'Loading images...';
    sessiontext2 = 'The experiment will begin shortly';
    sessiontext3 = sprintf('Session %i of the experiment has started! Good luck!',sesnum);
    trialtext2 = 'Please wait for the experimenter';
    
    HideCursor;
    Screen('DrawTexture',wptr,fixglind,[],destRect);
    Screen('TextSize', wptr, fontsize);
    Screen('TextFont', wptr, myfont);
    Screen('DrawText',wptr, sessiontext,(xres - length(sessiontext)*9)/2,yres*(5/12),txtcolor);
    Screen(wptr,'FillRect',black,pRect'); 
    Screen('Flip',wptr);
    
    %Track generation time
    tic;
    
    %Generate Gabor and noise images for all blocks
    noisecount = noisehz*4;
    %pregenexp5(noisecount,noisespf,trialnum,gaborspf,fixrad);
    [noises, allgabors, allspfs, allrandrot] = genexp5(noisecount,noisespf,trialnum,gaborspf,fixrad);
    
    %Concatenate Gabor patches from high and low distributions
    gabors = cat(3,allgabors{2},allgabors{1});
    spfs = [allspfs{2} allspfs{1}];
    randrots = [allrandrot{2} allrandrot{1}];
    
    %Index of spatial frequencies for Gabor images
    temp = [ones(1,(trialnum/2))*2 ones(1,(trialnum/2))*1];
%     if sesnum <= 1
%         repthisperm = randperm(trialnum,output.tperb); %Create first random permulation of length: trials per block
%         spfperm = [];
%         for b=1:blocknum
%             spfperm = [spfperm repthisperm(randperm(output.tperb))]; %If it is the training session, repeat exact distribution draws
%         end
%     else
        spfperm = randperm(trialnum);
%     end
    output.highlow = temp(spfperm); %Index if it was draw from high or low distribution
    clear temp;
    gabors = gabors(:,:,spfperm);
    output.spfs = spfs(spfperm);
    
    altperm = randperm(trialnum);
    output.randrots = randrots(altperm);
    
    %Locate fixation spot
    blackfix = (makefixation([],fixrad) == 0);
    
    %Locate border spot
    %darkfix = (makefixborder([],fixrad) == 0);
    
    %Change NaN to gray
    noises(isnan(noises)) = 0;
    gabors(isnan(gabors)) = 0;
    
    %Random order of noise images
    whichnoises = nan(trialnum,noisehz*4);
    for t=1:trialnum
        whichnoises(t,:) = randperm(noisehz*4);
    end
    
    %Index for unique SNRs
    [thesesnrs,~,usnrs] = unique(snrvec); 
    
    %Multiply each noise matrix by its luminance
    trialnoises = nan(size(noises,1),size(noises,2),numel(snrs),noisehz*4);
    for b=1:numel(snrs)
        splum = slum./thesesnrs(b); %noise luminance
        trialnoises(:,:,b,:) = splum.*noises;
    end
    thesenoises = 255*sqrt(trialnoises/2 + .5);
    %thesenoises(repmat(repmat(darkfix,1,1,b),1,1,1,noisehz*4)) = darkgray(1); %Recreate border
    thesenoises(repmat(repmat(blackfix,1,1,b),1,1,1,noisehz*4)) = black(1); %Recreate fixation spot
    noisescreen = nan(1,size(noises,3));
    
    %Create noise images as OpenGL textures
    for n=1:(noisehz*4)
        for b=1:numel(snrs)
            noisescreen(b,n) = Screen('MakeTexture',wptr,thesenoises(:,:,b,n)); %The square root is in order to account for monitor gamma. That is, the monitor approximately squares the input stimulus color value
        end
    end
    clear thesenoises noises
    
    %Save generation time
    output.gentime = toc;
    
    %Display second text screen
    Screen('DrawTexture',wptr,fixglind,[],destRect);
    Screen('TextSize', wptr, fontsize);
    Screen('TextFont', wptr, myfont);
    Screen('DrawText',wptr, sessiontext2,(xres - length(sessiontext2)*9)/2,yres*(5/12),txtcolor);
    Screen(wptr,'FillRect',black,pRect'); 
    Screen('Flip',wptr);
    
    %Wait for spacebar
    FlushEvents('keyDown');
    [char,when] = GetChar; %Wait for keypress to continue
    notspace=1;
    while notspace
        switch char
            case ' '
                notspace =0;
            otherwise
                [char,when] = GetChar; %Wait for keypress to continue
                notspace =1;
        end
    end
    
    %Display third text screen
    Screen('DrawTexture',wptr,fixglind,[],destRect);
    Screen('TextSize', wptr, fontsize);
    Screen('TextFont', wptr, myfont);
    Screen('DrawText',wptr, sessiontext3,(xres - length(sessiontext3)*9)/2,yres*(5/12),txtcolor);
    Screen(wptr,'FillRect',black,pRect'); 
    Screen('Flip',wptr);
    pause(2);
    
    %Pause before beginning of the block, pause for a second after text
    Screen('DrawTexture',wptr,fixglind,[],destRect);
    Screen(wptr,'FillRect',black,pRect'); 
    Screen('Flip',wptr);
    pause(1);
    
    %Initialize timer
    tic;
    for trials = 1:trialnum
      %printiter(trials);
      if ~cut %ESC key track
        %trialinbl = trials - output.tperb*(block-1);
        
        %Create image with specified snr
        trialgabor = slum*gabors(:,:,trials);
        %Shift domain to [0 1], convert to appropriate luminance
        %transformation for gaborimage color values, takes into account
        %monitor gamma
        if trials == 1
            bothscreen = nan(1,noisehz*4);
        end
        for n=1:noisehz*4
            thisimage = 255*sqrt( (trialgabor(:,:)+trialnoises(:,:,usnrs(trials),whichnoises(trials,n)) ...
               )/2 + .5);
            %thisimage(darkfix) = darkgray(1);
            thisimage(blackfix) = black(1);
            if trials > 1 %Clear former textures to save memory
                Screen('Close',bothscreen(n));
            end
            bothscreen(n) = Screen('MakeTexture',wptr,thisimage);
        end
        clear trialgabor fixx fixy trialflic
        
        trialflic = [ones(1,trialnframes(trials)) gaborflic];
        output.trialflic{trials} = trialflic;
          
        %Wait at least lboundwait seconds between trials
        lboundwait = output.intertrial(trials);
        output.elapsedtime(trials) = toc;
        if output.elapsedtime(trials) < lboundwait
            pause(lboundwait-output.elapsedtime(trials));
        end
        output.fixedtime(trials) = toc;
        
        if ~runtest
            CedrusResponseBox('FlushEvents',chandle);
        end
        
        %Display rush loops (Rush is apparently obsolete in PTB3, test this)
        Priority(MaxPriority(wptr)); %New to PTB3
        
        %Loop 1: Noise interval for 500ms - 1000ms
        %Response interval for 1000ms - 2000ms, accept responses
        noisenum = 0;
        bwswitch = 1; 
        for i = 1:numCycleFrames(trials)
            if noiseflic(i)
                noisenum = noisenum + 1;
                bwswitch = mod(bwswitch,2) + 1; %Changes 1 to 2 and vica versa
            end
            if trialflic(i) == 2
                Screen('DrawTexture',wptr,bothscreen(noisenum),[],destRect);
            else
                Screen('DrawTexture',wptr,noisescreen(usnrs(trials),whichnoises(trials,noisenum)),[],destRect);
            end
            if i == 1
                Screen(wptr,'FillRect',white,pRect([3 4],:)');
            else
                Screen(wptr,'FillRect',black,pRect([3 4],:)');
            end
            Screen(wptr,'FillRect',blackwhite{bwswitch},pRect([1 5],:)');
            Screen(wptr,'FillRect',blackwhite{trialflic(i)},pRect([2 6],:)');
            Screen('Flip',wptr);
        end
    
        %Loop 2: Keep displaying black fixation spot (only) for 250ms to collect responses
        for frames = 1:round(refrate/4)
            Screen('DrawTexture',wptr,fixglind,[],destRect);
            if frames == 1
                Screen(wptr,'FillRect',white,pRect([3 4],:)');
                Screen(wptr,'FillRect',black,pRect([1 2 5 6],:)'); 
            else
                Screen(wptr,'FillRect',black,pRect'); 
            end
            Screen('Flip',wptr);
        end

        %Timer to calculate time between the last trial and the next
        tic;
        
        %Play feedback sound
        if ~runtest
            evt = CedrusResponseBox('GetButtons',chandle);
            if isempty(evt)
                correct = 0;
            elseif (evt.button == 6 && output.highlow(trials) == 2) || ... %High frequency
                    (evt.button == 4 && output.highlow(trials) == 1) %Low frequency
                correct = 1;
            else
                correct = 0;
            end
        else
            correct = randi(2) - 1;
        end
        estcorrect(trials) = correct;
        stop(goodsound);
        stop(badsound);
        if correct
            %sound(ygood, Fsgood);
            %PsychPortAudio('Start', goodhand , 1, 0, 1);
            play(goodsound,[1 length(ygood)]);
        else
            %sound(ybad, Fsbad);
            %PsychPortAudio('Start', badhand , 1, 0, 1);
            play(badsound,[1 length(ybad)]);
        end
        
        % Trial Number Display
        Screen('DrawTexture',wptr,fixglind,[],destRect);
        Screen('TextSize', wptr, 36);
        Screen('TextFont', wptr, myfont);
        Screen('DrawText', wptr, sprintf('B%i/%i',block,blocknum), 10, 930, black);
        Screen('DrawText', wptr, sprintf('T%i/%i',trials+1,trialnum), 10, 1030, black);
        Screen(wptr,'FillRect',black,pRect');
        Screen('Flip',wptr);

        if ~cut
            if trials == trialnum
                %Show ending screen for 1 second
                percorrect = sum(estcorrect((trials-output.tperb+1):trials))/output.tperb;
                endtext = ['Done!  ',...
                    num2str(round(percorrect*100)),'% correct this block. Thank you for participating!'];
                Screen('DrawTexture',wptr,fixglind,[],destRect);
                Screen('TextSize', wptr, fontsize);
                Screen('TextFont', wptr, myfont);
                Screen('DrawText',wptr, endtext,(xres - length(endtext)*9)/2,yres*(5/12),txtcolor);
                Screen(wptr,'FillRect',black,pRect'); 
                %Pause for 1 second
                pause(1);
                Screen('Flip',wptr);
                %Pause for 1 second
                pause(1);
                %Wait for spacebar to end program
                FlushEvents('keyDown');
                [char,~] = GetChar; %Wait for keypress to continue
                notspace=1;
                while notspace
                    switch char
                        case advancechar
                            notspace =0;
                        otherwise
                            [char,~] = GetChar; %Wait for keypress to continue
                            notspace =1;
                    end
                end
            elseif trials/output.tperb == round(trials/output.tperb)
                 %Take a break every 'output.tperb' trials and show ending Screens
                percorrect = sum(estcorrect((trials-output.tperb+1):trials))/output.tperb;
                trialtext = ['Block ',num2str(block),' complete!  ',...
                    num2str(round(percorrect*100)),'% correct this block. You may now take a break!'];
                
                block = block + 1;
                
                Screen('DrawTexture',wptr,fixglind,[],destRect);
                Screen('TextSize', wptr, fontsize);
                Screen('TextFont', wptr, myfont);
                Screen('DrawText',wptr,trialtext,(xres - length(trialtext)*9)/2,yres*(5/12),txtcolor);
                Screen('DrawText',wptr,trialtext2,(xres - length(trialtext2)*9)/2,yres*(5/12) + 32,txtcolor);
                Screen(wptr,'FillRect',black,pRect'); 
                %Pause for 1 second
                pause(1);
                Screen('Flip',wptr);
                %Pause for 10 seconds
                pause(10);
                
                %Wait for spacebar
                FlushEvents('keyDown');
                [char,~] = GetChar; %Wait for keypress to continue
                notspace=1;
                while notspace
                    switch char
                        case advancechar
                            notspace =0;
                            Screen('DrawTexture',wptr,fixglind,[],destRect);
                            Screen('TextSize', wptr, fontsize);
                            Screen('TextFont', wptr, myfont);
                            Screen('DrawText',wptr,sprintf('Block %i of the experiment has started! Good luck!',block),(xres - length(sessiontext3)*9)/2,yres*(5/12),txtcolor);
                            Screen(wptr,'FillRect',black,pRect'); 
                            Screen('Flip',wptr);
                            %Timer to calculate time between the last trial and the next
                            %Pause for 1 second
                            pause(2);
                            tic;
                        case escapechar %Escape from experiment and save current data (for experimenter)
                            notspace =0;
                            %RestoreScreen(whichScreen);
                            ShowCursor;
                            Screen('CloseAll');
                            output.ESC_time = clock;
                            output.estcorrect = estcorrect;
                            eval([subnum,'_ExpSession',num2str(sesnum),'=output;']);
                            if ~exist([pwd,'/exp5behav'],'dir')
                                mkdir('exp5behav');
                            end
                            if ~runtest
                                eval(['save(''exp5behav/',subnum,'_Exp_',num2str(sesnum),'_',date,'.mat'',''-struct'', ''',subnum,'_ExpSession',num2str(sesnum),''');']);
                            end
                            warning on all;
                            if ~runtest
                                CedrusResponseBox('Close',chandle);
                            end
                            %PsychPortAudio('Close', goodhand);
                            %PsychPortAudio('Close', badhand);
                            return
                        otherwise
                            [char,when] = GetChar; %Wait for keypress to continue
                            notspace =1;
                    end
                end
                %Pause before beginning of the block, pause for a second after text
                Screen('DrawTexture',wptr,fixglind,[],destRect);
                Screen(wptr,'FillRect',black,pRect'); 
                Screen('Flip',wptr);
                pause(1);
            end
            
        end
      end  
    end
catch me
    fprintf('\n');
    %RestoreScreen(whichScreen);
    ShowCursor;
    Screen('CloseAll');
    output.error_time = clock;
    output.estcorrect = estcorrect;
    %Organize data and time in a string
    rightnow = clock;
    rightnow = num2cell(rightnow)';
    timestr = sprintf('_%i',rightnow{1:5});
    if ~exist([pwd,'/exp5behav'],'dir')
        mkdir('exp5behav');
    end
    if ~runtest
        eval(['save(''exp5behav/',subnum,'_Exp_',num2str(sesnum),timestr,'.mat'',''-struct'', ''output'');']);
    end
    if ~runtest
        CedrusResponseBox('Close',chandle);
    end
    %PsychPortAudio('Close', goodhand);
    %PsychPortAudio('Close', badhand);
    rethrow(me); %rethrow reproduces the original error, stored in the object 'me'
end

fprintf('\n');

%RestoreScreen(whichScreen);
ShowCursor;
Screen('CloseAll');

%Output time finished
output.finish_time = clock;

%Estimated accuracy
output.estcorrect = estcorrect;

%Organize data and time in a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

if ~exist([pwd,'/exp5behav'],'dir')
    mkdir('exp5behav');
end
if ~runtest
    eval(['save(''exp5behav/',subnum,'_Exp_',num2str(output.sesnum),timestr,'.mat'',''-struct'', ''output'');']);
end
warning on all;
if ~runtest
    CedrusResponseBox('Close',chandle);
end
%PsychPortAudio('Close', goodhand);
%PsychPortAudio('Close', badhand);

%% ----------------------------------------------------------------------%%
function [noises, allgabors, allspfs, allrandrot] = genexp5(numnoise,noisespfrq,numgabor,gaborspfrq,radius)
%GENEXP% - generates images for pdmexp5
%
%Useage: 
%  >> pregenexp5(numnoise,noisespfrq,ngabor,gaborspfrq)
%
%Inputs:
%   numnoise - Number of noise images to generate
%
%   noisespfrq - Spatial frequency of pixelated visual noise (cycles per cm)
%
%   numgabor - Number of gabor images to pregenerate
%
%   gaborspfrq - Spatial frequencies of gabors (cycles per cm)
%
%   radius - Radius size of fixation spot (cycles per cm)

%% Code

%Gabor size
gaborsize = 10;

if round(numgabor/length(gaborspfrq)) ~= numgabor/length(gaborspfrq)
    error('Number of Gabor spatial frequencies must be a divisor of the number of Gabor images');
end

if numnoise > 0
    fprintf('Building noise images...\n');
end
noises = nan(1000,1000,numnoise);
% for m=1:numnoise
%     noises(:,:,m) = makepixeled([],radius,noisespfrq,[]);
% end

%%Combine two spatial frequencies to equally mask both high and low signal
%%stimuli
for m=1:numnoise
%     temphn = makepixeled([],radius,3,[]);
%     templn = makepixeled([],radius,2,[]);
%     noises(:,:,m) = (temphn + templn)/2;
      noises(:,:,m) = makenoise([],2,.1,radius,[2 3]);
end


%%Make Gabor images
if ~exist([pwd,'/pregen/gabor'],'dir')
    mkdir('pregen/gabor');
end

if numgabor > 0
    fprintf('Building Gabor images...\n');
    for f=1:length(gaborspfrq)
        numimages = numgabor/length(gaborspfrq);
        gabors = nan(1000,1000,numimages);
        spf = nan(numimages,1);
        randrot = nan(numimages,1);
        for m=1:numimages
            %Draw from truncated normal distributions
            if f==1
%                 draw = 100;
%                 while draw > 0 %There is a 2.28% chance of the draw being less than 0
%                     draw = normrnd(-.1,.05);
%                 end
            draw = -.1;
            elseif f==2
%                 draw = -100;
%                 while draw < 0 %There is a 2.28% chance of the draw being more than 0
%                     draw = normrnd(.1,.05);
%                 end
            draw = .1;
            end
            spf(m) = draw + 2.5; %modes = 2.4 | 2.6
            randrot(m) = rand(1)*360;
            gabors(:,:,m) = makegabor([],gaborsize,randrot(m)*360,[],spf(m));
        end
        allgabors{f} = gabors;
        allspfs{f} = spf;
        allrandrot{f} = randrot;
    end
end

timeels = toc;
fprintf('The images took %3.2f minutes to generate and save.\n',timeels/60);