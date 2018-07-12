%pdm5b_simplemodels.m - Calculates Bayes Factors for simple regression models
%
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


%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  01/09/18        Michael Nunez                Original code
%  07/12/18        Michael Nunez           Run models on data without RT cutoffs

%% Initial
singletrial = readtable('../Data/N200_rt_window_150_275_fixed350cutoff.csv');
sesdata = readtable('../Data/N1deflec2_cutoffs_allSNR_window_150_275_fixed350cutoff.csv');


%% JAGS code for simple linear regression
model = {
'model {'
    '##########'
    '#Prior of Regression intercept (remaining reaction time)'
    'tercond ~ dnorm(.3, pow(.25,-2))'
    '#Prior of effect of single-trial N200 latency on reaction time'
    'n1gammault ~ dnorm(1,pow(3,-2)) #''Informative'' prior for Bayes Factor calculation'
    '#Prior of Variability in reaction time'
    'tersubsd ~ dgamma(.2,1)'

    '##########'
    '#Linear regression (normal likelihood)'
    'for (i in 1:N) {'
        'rtpercentile[i] ~ dnorm(tercond + n200lat[i]*n1gammault, pow(tersubsd,-2) )'
    '}'
'}'
};

%% Code for Trinity

rawrt = singletrial.RT/1000; %convert to secs from ms
stN200 = singletrial.x_Single_trialN200Latencies/1000; %convert to secs from ms

rt10per = sesdata.x10thRTPercentiles/1000; %convert to secs from ms
taN200 = sesdata.x_N1Latencies/1000; %convert to secs from ms

tadeflec = sesdata.N1Deflection/1000; %convert to secs from ms

datast.rtpercentile = rawrt;
datast.n200lat = stN200;
datast.N = length(stN200);

datata.rtpercentile = rt10per;
datata.n200lat = taN200;
datata.N = length(taN200);

datadf.rtpercentile = rt10per;
datadf.n200lat = tadeflec;
datadf.N = length(tadeflec);

% Track these variables
params = {'tercond', 'n1gammault', 'tersubsd'};

initstruct = @()struct(...
'n1gammault', rand*4 - 2, ...
'tercond', rand, ...
'tersubsd', rand*.09 + .01);

%% Run single-trial model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF_singletrial';

modelname = [modeltype,timestr];

fprintf('Building JAGS model %s and saving output...\n',modelname);

nsamples = 5e3;
nburnin = 2e3;
nchains =3;
thin =10;
verbosity =1;
parallelit = 0; %Set this to 1 if GNU Parallel is installed
maxcores = 3;
modules = {'wiener' 'dic'};

tic
[stats, chains, diagnostics, info] = callbayes('jags', ...
    'model', model, ...
    'data', datast, ...
    'nsamples', nsamples, ...
    'nburnin', nburnin, ...
    'nchains', nchains, ...
    'thin',thin,...
    'verbosity', verbosity, ...
    'monitorparams', params, ...
    'parallel',parallelit, ...
    'maxcores',maxcores, ...
    'modules',modules, ...
    'init', initstruct); 

info.comptime = toc/60;
fprintf('JAGS took %f minutes!\n', info.comptime)

%Calculate Bayes Factor using Savage-Dickey ratio
xi = -4:.01:4;
[density_st] = ksdensity(chains.n1gammault(:),xi);
figure;
plot(xi,density_st,'r');
line(xi,normpdf(xi,1,3));
title('Single-trial N200 versus RT Slope Posterior')
numerator = density_st(xi==1);
denominator = normpdf(1,1,3);
bf_st = numerator/denominator;
fprintf('The Bayes Factor of the single trial data is %.3f ! \n',bf_st);

save(sprintf('jagsmodel%s.mat',modelname),'stats', 'chains', 'diagnostics','info','params','datast','model','bf_st');
system('rm -r wdir');

%% Run trial-averaged model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF_trialaveraged';

modelname = [modeltype,timestr];

fprintf('Building JAGS model %s and saving output...\n',modelname);

nsamples = 5e3;
nburnin = 2e3;
nchains =3;
thin =10;
verbosity =1;
parallelit = 0; %Set this to 1 if GNU Parallel is installed
maxcores = 3;
modules = {'wiener' 'dic'};

tic
[stats, chains, diagnostics, info] = callbayes('jags', ...
    'model', model, ...
    'data', datata, ...
    'nsamples', nsamples, ...
    'nburnin', nburnin, ...
    'nchains', nchains, ...
    'thin',thin,...
    'verbosity', verbosity, ...
    'monitorparams', params, ...
    'parallel',parallelit, ...
    'maxcores',maxcores, ...
    'modules',modules, ...
    'init', initstruct); 

info.comptime = toc/60;
fprintf('JAGS took %f minutes!\n', info.comptime)

%Calculate Bayes Factor using Savage-Dickey ratio
xi = -4:.01:4;
[density_ta] = ksdensity(chains.n1gammault(:),xi);
figure;
plot(xi,density_ta,'r');
line(xi,normpdf(xi,1,3));
title('Single-trial N200 versus RT Slope Posterior')
numerator = density_ta(xi==1);
denominator = normpdf(1,1,3);
bf_ta = numerator/denominator;
fprintf('The Bayes Factor of the trial-averaged data is %.3f ! \n',bf_ta);

save(sprintf('jagsmodel%s.mat',modelname),'stats', 'chains', 'diagnostics','info','params','datata','model','bf_ta');
system('rm -r wdir');

%% Run trial-averaged model with deflection data

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF_deflection';

modelname = [modeltype,timestr];

fprintf('Building JAGS model %s and saving output...\n',modelname);

nsamples = 5e3;
nburnin = 2e3;
nchains =3;
thin =10;
verbosity =1;
parallelit = 0; %Set this to 1 if GNU Parallel is installed
maxcores = 3;
modules = {'wiener' 'dic'};

tic
[stats, chains, diagnostics, info] = callbayes('jags', ...
    'model', model, ...
    'data', datadf, ...
    'nsamples', nsamples, ...
    'nburnin', nburnin, ...
    'nchains', nchains, ...
    'thin',thin,...
    'verbosity', verbosity, ...
    'monitorparams', params, ...
    'parallel',parallelit, ...
    'maxcores',maxcores, ...
    'modules',modules, ...
    'init', initstruct); 

info.comptime = toc/60;
fprintf('JAGS took %f minutes!\n', info.comptime)

%Calculate Bayes Factor using Savage-Dickey ratio
xi = -4:.01:4;
[density_df] = ksdensity(chains.n1gammault(:),xi);
figure;
plot(xi,density_df,'r');
line(xi,normpdf(xi,1,3));
title('Single-trial N200 versus RT Slope Posterior')
numerator = density_df(xi==1);
denominator = normpdf(1,1,3);
bf_df = numerator/denominator;
fprintf('The Bayes Factor of the trial-averaged data is %.3f ! \n',bf_df);
fprintf('The Bayes Factor for the null hypothesis is %.3f ! \n',1/bf_df);

save(sprintf('jagsmodel%s.mat',modelname),'stats', 'chains', 'diagnostics','info','params','datadf','model','bf_df');
system('rm -r wdir');