# pdm5b_scatterplots.R - Creates scatter plots of N200 latencies versus RT percentiles
#
# Copyright (C) 2017 Michael D. Nunez, <mdnunez1@uci.edu>
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
#  12/06/17        Michael Nunez                 Original code
#  02/23/18        Michael Nunez          Addition of deflection time
#  07/10/18        Michael Nunez          Adding  N200 peak-latency vs NDT
#                                     Fix manual colors
#  07/12/18        Michael Nunez           Use data without cutoffs

## To do: Generate plotly

## References:
# http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3

## Necessary packages
library(ggplot2)
library(R.matlab)

## Code
fontsize = 18;

# Read the tables and samples
trialdata = read.csv(
  '../Data/N200_rt_window_150_275.csv')
colnames(trialdata)
for (i in seq(1,dim(trialdata)[1])) {
	if (trialdata$Experiment[i] == 1) {
		expstr = 'Exp. 1'
	} else {
		expstr = 'Exp. 2'
	}
	if (trialdata$SNR.condition[i] == 0) {
		condstr = 'High'
	} else if (trialdata$SNR.condition[i] == 1) {
		condstr = 'Med'
	} else {
		condstr = 'Low'
	}
	trialdata$NoiseCondition[i] = sprintf('%s %s',expstr,condstr)
}
trialdata$NoiseCondition = factor(trialdata$NoiseCondition, levels = c("Exp. 1 Low", "Exp. 1 Med", "Exp. 1 High","Exp. 2 Low", "Exp. 2 Med", "Exp. 2 High"))

sesdata = read.csv(
  '../Data/N1deflec2_allSNR_window_150_275.csv')
colnames(sesdata)
for (i in seq(1,dim(sesdata)[1])) {
	if (sesdata$Experiment[i] == 1) {
		expstr = 'Exp. 1'
	} else {
		expstr = 'Exp. 2'
	}
	if (sesdata$SNR.condition[i] == 0) {
		condstr = 'High'
	} else if (sesdata$SNR.condition[i] == 1) {
		condstr = 'Med'
	} else {
		condstr = 'Low'
	}
	sesdata$NoiseCondition[i] = sprintf('%s %s',expstr,condstr)
}
sesdata$NoiseCondition = factor(sesdata$NoiseCondition, levels = c("Exp. 1 Low", "Exp. 1 Med", "Exp. 1 High","Exp. 2 Low", "Exp. 2 Med", "Exp. 2 High"))

samples = readMat('../Models/jagsmodel_all_n1lat_random_lapseJul_11_18_10_26.mat')
sesdata$EEG.Session.Counter[sesdata$EEG.Session.Counter < 16]  = sesdata$EEG.Session.Counter[sesdata$EEG.Session.Counter < 16] +1
sesdata$EEG.Session.Counter[sesdata$EEG.Session.Counter > 23]  = sesdata$EEG.Session.Counter[sesdata$EEG.Session.Counter > 23] -1
sesdata$EEG.Session.Counter[sesdata$EEG.Session.Counter > 35]  = sesdata$EEG.Session.Counter[sesdata$EEG.Session.Counter > 35] -1
for (n in 1:length(sesdata$EEG.Session.Counter)) {
	sesdata$NDT[n] = median(samples$tersub[(sesdata$SNR.condition[n]+1),(sesdata$EEG.Session.Counter[n]),,])
	}
#Convert to milliseconds
sesdata$NDT = sesdata$NDT*1000

# Get linear model equation and R^2_adj as a string
# SOURCE: http://goo.gl/K4yh
lm_eqn <- function(lmresults){
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2[adj]~"="~r2adj,
         list(a = format(coef(lmresults)[1], digits=0), 
              b = format(coef(lmresults)[2], digits = 4), 
             r2adj = format(summary(lmresults)$adj.r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

# http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=7
cbbPalette <- c("#7570B3", "#1B9E77", "#D95F02", "#E7298A", "#66A61E", "#E6AB02")
shapePalette <- c(21, 24, 22, 25, 23, 21)


## Plots
lm1 = lm(X10th.RT.percentiles ~ X..N1.latencies,sesdata)
png('n1_rt10per.png',units="in",width=10,height=10,res=300)
plot1 = ggplot(sesdata,aes(x=X..N1.latencies, y=X10th.RT.percentiles)) +
geom_point(aes(shape=NoiseCondition, fill=NoiseCondition),size=3) +
geom_smooth(method="lm", colour='black',size=1.5) +
xlim(150,275) + 
coord_cartesian(ylim=c(400,900)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_text(size=fontsize,face='bold')) +
labs(x='Trial-averaged N200 peak-latency (ms)',y='10th Reaction time percentiles (ms)',color='Noise\ncondition') +
scale_fill_manual(values = cbbPalette) + 
scale_shape_manual(values = shapePalette) + 
geom_abline(slope=1,intercept=lm1$coefficients[1]+30,colour="#A6761D",size=2, linetype=2)
# plot1 + geom_text(x = 200, y = 1200, label = lm_eqn(lm1), parse = TRUE,size=8)
plot(plot1)
dev.off()

#Note: use coord_cartesian(ylim=c(400,900)) instead of ylim(400,900) to 
lm2 = lm(RT ~ X..Single.trial.N200.latencies,trialdata)
png('n200_rt.png',units="in",width=10,height=10,res=300)
plot2 = ggplot(trialdata,aes(x=X..Single.trial.N200.latencies, y=RT)) +
geom_point(aes(shape=NoiseCondition, fill=NoiseCondition)) +
geom_smooth(method="lm", colour='black',size=1.5) +
xlim(150,275) + 
coord_cartesian(ylim=c(400,1400)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_text(size=fontsize,face='bold')) +
labs(x='Single-trial N200 peak-latency (ms)',y='Reaction time (ms)',color='Noise\ncondition') +
scale_fill_manual(values = cbbPalette) + 
scale_shape_manual(values = shapePalette) + 
geom_abline(slope=1,intercept=lm2$coefficients[1]+10,colour="#A6761D",size=2,linetype=2)
# plot2 + geom_text(x = 200, y = 1200, label = lm_eqn(lm2), parse = TRUE,size=8)
plot(plot2)
dev.off()

lm3 = lm(X10th.RT.percentiles ~ N1.deflection,sesdata)
png('n1deflec_rt10per.png',units="in",width=10,height=10,res=300)
plot3 = ggplot(sesdata,aes(x=N1.deflection, y=X10th.RT.percentiles)) +
geom_point(aes(shape=NoiseCondition, fill=NoiseCondition),size=3) +
geom_smooth(method="lm", colour='black',size=1.5) +
xlim(50,150) + 
coord_cartesian(ylim=c(400,900)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_text(size=fontsize,face='bold')) +
labs(x='Trial-averaged N200 deflection time (ms)',y='10th Reaction time percentiles (ms)',color='Noise\ncondition') +
scale_fill_manual(values = cbbPalette) + 
scale_shape_manual(values = shapePalette) + 
geom_abline(slope=1,intercept=lm3$coefficients[1]-100,colour="#A6761D",size=2, linetype=2)
# plot3 + geom_text(x = 200, y = 1200, label = lm_eqn(lm3), parse = TRUE,size=8)
plot(plot3)
dev.off()

lm4 = lm(NDT ~ X..N1.latencies,sesdata)
png('n1_NDT.png',units="in",width=10,height=10,res=300)
plot4 = ggplot(sesdata,aes(x=X..N1.latencies, y=NDT)) +
geom_point(aes(shape=NoiseCondition, fill=NoiseCondition),size=3) +
geom_smooth(method="lm", colour='black',size=1.5) +
xlim(150,275) + 
coord_cartesian(ylim=c(300,700)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_text(size=fontsize,face='bold')) +
labs(x='Trial-averaged N200 peak-latency (ms)',y='Non-decision time estimates (ms)',color='Noise\ncondition') +
scale_fill_manual(values = cbbPalette) + 
scale_shape_manual(values = shapePalette) + 
geom_abline(slope=1,intercept=lm4$coefficients[1]+110,colour="#A6761D",size=2, linetype=2)
# plot1 + geom_text(x = 200, y = 1200, label = lm_eqn(lm4), parse = TRUE,size=8)
plot(plot4)
dev.off()