
%%Generate samples from the drift-diffusion model with trial-to-trial drift variability
rng(11);
tersub = rand(30,100)*.15 + .35; % Uniform .35 to .5
alphasub = rand(30,100)*.06 + .08; % Uniform from .08 to .14
deltasub = rand(30,100)*.8 - .4; % Uniform from  -.4 to .4
%tertrialsd = repmat(linspace(0,.1,30)',1,100); % From 0 to .1
tertrialsd = rand(30,100)*.1; % Uniform from 0 to .1
deltatrialsd = rand(30,100)*.2; % Uniform from 0 to .2
%prob_mindwander = rand(30,100)*.10; % Uniform from 0 to 10 percent of trials
prob_mindwander = repmat(linspace(0,.10,30)',1,100); % From 0 to 10 percent of trials
y=zeros(30,100,100);
rt=zeros(30,100,100);
acc=zeros(30,100,100);
for n=1:30,
	for s=1:100,
		[tempt, tempx] = simuldiff([alphasub(n, s), tersub(n, s),...
			deltatrialsd(n, s), alphasub(n, s)*.5, 0, tertrialsd(n, s), deltasub(n, s)],100);
		mindwanderx = randi(2,1,100) - 1; 
		mindwandert = rand(1,100)*2; % Randomly distributed from 0 to 2 seconds

		mindwander_trials = randperm(100,round(100*prob_mindwander(n, s)));
		tempx(mindwander_trials) = mindwanderx(mindwander_trials);
		tempt(mindwander_trials) = mindwandert(mindwander_trials);

		y(n, s, :) = tempt.*(tempx*2 - 1);
		rt(n, s, :) = tempt;
		acc(n, s, :) = tempx;
	end
end
