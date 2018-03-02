# pdm5b_posteriordistributions1.R - Plots overlapping posterior distributions from Model 1
#
# Copyright (C) 2018 Michael D. Nunez, <mdnunez1@uci.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Record of Revisions
#   Date           Programmers               Description of change
#   ====        =================            =====================
#  01/26/18        Michael Nunez                 Original code
#  01/29/18        Michael Nunez                 Remove title

## Necessary packages
library(ggplot2)
library(ggjoy)
library(viridis)
library(R.matlab)

loadloc = '../Models'

## Code

# Read in the reaction times
samples = readMat(paste(loadloc,
  '/jagsmodel_all_n1lat_randomJan_23_18_12_37.mat',sep=""))

mainN200effect = as.vector(samples[1]$n1gammault[1,1,,])
print(sprintf('Length of N200 effect posterior samples is %d',length(mainN200effect)))
mainN200effectCI = quantile(mainN200effect, prob=c(.025,.5,.975))
print(sprintf('Effect (median posterior and 95%% credible interval) of trial-averaged N200 latency on non-decision time: %.2f, CI: [%.2f, %.2f]',mainN200effectCI[2],mainN200effectCI[1],mainN200effectCI[3]));
mainN200effect_density = density(mainN200effect)
x11()
plot(N200effect_density)


nsamps = length(mainN200effect)
N200effects <- data.frame(matrix(ncol = 2, nrow = nsamps*7))
colnames(N200effects) <- c("Effect", "Condition")

N200effects$Effect[1:nsamps] = as.vector(samples[1]$n1gammault[1,1,,])
N200effects$Condition[1:nsamps] = sprintf('G')
N200effects$Effect[(nsamps+1):(nsamps*2)] = as.vector(samples[4]$n1gammacond[1,1,1,,])
N200effects$Condition[(nsamps+1):(nsamps*2)] = sprintf('A')
N200effects$Effect[(nsamps*2+1):(nsamps*3)] = as.vector(samples[4]$n1gammacond[1,1,2,,])
N200effects$Condition[(nsamps*2+1):(nsamps*3)] = sprintf('B')
N200effects$Effect[(nsamps*3+1):(nsamps*4)] = as.vector(samples[4]$n1gammacond[1,1,3,,])
N200effects$Condition[(nsamps*3+1):(nsamps*4)] = sprintf('C')
N200effects$Effect[(nsamps*4+1):(nsamps*5)] = as.vector(samples[4]$n1gammacond[1,2,1,,])
N200effects$Condition[(nsamps*4+1):(nsamps*5)] = sprintf('D')
N200effects$Effect[(nsamps*5+1):(nsamps*6)] = as.vector(samples[4]$n1gammacond[1,2,2,,])
N200effects$Condition[(nsamps*5+1):(nsamps*6)] = sprintf('E')
N200effects$Effect[(nsamps*6+1):(nsamps*7)] = as.vector(samples[4]$n1gammacond[1,2,3,,])
N200effects$Condition[(nsamps*6+1):(nsamps*7)] = sprintf('F')

condlabs = c("Exp. 1 High Noise", "Exp. 1 Med Noise", "Exp. 1 Low Noise", "Exp. 2 High Noise", "Exp. 2 Med Noise", "Exp. 2 Low Noise", "Overall (hierarchical) effect") 

#Color-blind friendly palette (http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## Plot
png(paste('Model1_N200effects.png',sep=""),units="in",width=10,height=10,res=300)
plot1 = ggplot(N200effects,aes(x=Effect,y=Condition,fill=Condition)) +
  geom_joy(scale=3, alpha=.67) + xlim(-3,3) + 
  xlab('Effect of trial-averaged N200 latency on non-decision time \n (ms increase per ms increase in N200 latency)') +
  ylab('') +
  # ggtitle('ERP latency effect posterior distributions') +
  annotate("segment", x=0, xend=0, y=1, yend=10, colour = 'red',size=3,alpha=.5) +
  annotate("segment", x=1, xend=1, y=1, yend=10, colour = 'blue',size=3,alpha=.5) +
  theme(text = element_text(size=20), legend.position="none") +
  scale_fill_manual(values = tail(cbbPalette,n=7)) +
  scale_y_discrete(labels=condlabs)
#Print Bayes Factors
BF = vector()
for (n in seq(1,7)) {
	temp_density = density(N200effects$Effect[(nsamps*(n-1)+1):(nsamps*n)])
	numerator = approx(temp_density$x,temp_density$y,xout=1)
	denominator = dnorm(1,mean=1,sd=3)
	BF[n] = numerator$y / denominator
	if (n==1) {
		yplace = 7.5
	} else {
		yplace = n - 0.5
	}
	plot1 = plot1 + annotate("text",x=2,y=yplace,label=sprintf('BF1: %3.2f',BF[n]),size=5,colour = 'blue')
}
plot(plot1)
dev.off()
