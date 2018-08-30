# pdm5b_trainingplots.R - Creates scatter plots of N200 latencies and RT percentiles over training sessions
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
#  08/30/18        Michael Nunez              Converted from pdm5b_scatterplots.R

## To do: Generate plotly

## References:
# http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12

## Necessary packages
library(ggplot2)
library(R.matlab)

## Code
fontsize = 18;

# Read the tables and samples

sesdata = read.csv(
  '../Data/N1deflec2_allSNR_window_150_275.csv')
colnames(sesdata)

exp2data = sesdata[(sesdata$Experiment == 2),]


for (i in seq(1,dim(exp2data)[1])) {
	if (exp2data$True.Subject.Index[i] == 1) {
		expstr = 'Part. 1'
	} else if (exp2data$True.Subject.Index[i] == 11) {
		expstr = 'Part. 2'
	} else if (exp2data$True.Subject.Index[i] == 13) {
		expstr = 'Part. 3'
	} else if (exp2data$True.Subject.Index[i] == 14) {
		expstr = 'Part. 4'
	}
	if (exp2data$SNR.condition[i] == 0) {
		condstr = 'High'
	} else if (exp2data$SNR.condition[i] == 1) {
		condstr = 'Med'
	} else {
		condstr = 'Low'
	}
	exp2data$NoiseCondition[i] = sprintf('%s %s',expstr,condstr)
}
exp2data$NoiseCondition = factor(exp2data$NoiseCondition, levels = c("Part. 1 Low", "Part. 1 Med", "Part. 1 High","Part. 2 Low", "Part. 2 Med", "Part. 2 High", "Part. 3 Low", "Part. 3 Med", "Part. 3 High","Part. 4 Low", "Part. 4 Med", "Part. 4 High"))

samples = readMat('../Models/jagsmodel_all_n1lat_random_lapseJul_11_18_10_26.mat')
exp2data$EEG.Session.Counter[exp2data$EEG.Session.Counter < 16]  = exp2data$EEG.Session.Counter[exp2data$EEG.Session.Counter < 16] +1
exp2data$EEG.Session.Counter[exp2data$EEG.Session.Counter > 23]  = exp2data$EEG.Session.Counter[exp2data$EEG.Session.Counter > 23] -1
exp2data$EEG.Session.Counter[exp2data$EEG.Session.Counter > 35]  = exp2data$EEG.Session.Counter[exp2data$EEG.Session.Counter > 35] -1
for (n in 1:length(exp2data$EEG.Session.Counter)) {
	exp2data$NDT[n] = median(samples$tersub[(exp2data$SNR.condition[n]+1),(exp2data$EEG.Session.Counter[n]),,])
	}
#Convert to milliseconds
exp2data$NDT = exp2data$NDT*1000

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
cbbPalette <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")
shapePalette <- c(21, 24, 22, 25, 23, 21, 21, 24, 22, 25, 23, 21)

## Plots
lm1 = lm(X10th.RT.percentiles ~ X..N1.latencies,exp2data)
png('n1_rt10per_exp2.png',units="in",width=10,height=10,res=300)
plot1 = ggplot(exp2data,aes(x=X..N1.latencies, y=X10th.RT.percentiles)) +
geom_point(aes(shape=NoiseCondition, fill=NoiseCondition),size=3) +
geom_smooth(method="lm", colour='black',size=1.5) +
xlim(150,275) + 
coord_cartesian(ylim=c(400,900)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_text(size=fontsize,face='bold')) +
labs(x='Trial-averaged N200 peak-latency (ms)',y='NDT estimates: 10th RT percentiles (ms)',color='Noise\ncondition') +
scale_fill_manual(values = cbbPalette) + 
scale_shape_manual(values = shapePalette) + 
geom_abline(slope=1,intercept=lm1$coefficients[1]+175,colour="#A6761D",size=2, linetype=2)
# plot1 + geom_text(x = 200, y = 1200, label = lm_eqn(lm1), parse = TRUE,size=8)
plot(plot1)
dev.off()


lm2 = lm(NDT ~ X..N1.latencies,exp2data)
png('n1_NDT_exp2.png',units="in",width=10,height=10,res=300)
plot2 = ggplot(exp2data,aes(x=X..N1.latencies, y=NDT)) +
geom_point(aes(shape=NoiseCondition, fill=NoiseCondition),size=3) +
geom_smooth(method="lm", colour='black',size=1.5) +
xlim(150,275) + 
coord_cartesian(ylim=c(300,700)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_text(size=fontsize,face='bold')) +
labs(x='Trial-averaged N200 peak-latency (ms)',y='Non-decision time estimates (ms)',color='Noise\ncondition') +
scale_fill_manual(values = cbbPalette) + 
scale_shape_manual(values = shapePalette) + 
geom_abline(slope=1,intercept=lm2$coefficients[1]+210,colour="#A6761D",size=2, linetype=2)
# plot1 + geom_text(x = 200, y = 1200, label = lm_eqn(lm2), parse = TRUE,size=8)
plot(plot2)
dev.off()

lm3 = lm( X..N1.latencies ~ Session,exp2data)
png('Session.png',units="in",width=10,height=10,res=300)
plot3 = ggplot(exp2data,aes(x=Session, y=X..N1.latencies, shape=NoiseCondition, color=NoiseCondition, fill=NoiseCondition)) +
geom_point(size=3) +
geom_smooth(method="lm", aes(fill=NoiseCondition),size=1.5,se=FALSE) +
coord_cartesian(ylim=c(150,275)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_blank()) +
labs(y='Trial-averaged N200 peak-latency (ms)',color='Noise\ncondition') +
scale_x_continuous(name='Training session',breaks=unique(exp2data$Session),lim=c(1,7),labels=unique(exp2data$Session)) +
scale_fill_manual(values = cbbPalette) + 
scale_shape_manual(values = shapePalette) + 
scale_color_manual(values = cbbPalette)
# plot1 + geom_text(x = 200, y = 1200, label = lm_eqn(lm3), parse = TRUE,size=8)
plot(plot3)
dev.off()
confint(lm3, 'Session', level=0.95)

lm4 = lm( NDT ~ Session,exp2data)
png('Session_NDT_exp2.png',units="in",width=10,height=10,res=300)
# plot4 = ggplot(exp2data,aes(x=Session, y=NDT, shape=NoiseCondition, color=NoiseCondition, fill=NoiseCondition)) +
# geom_point(size=3) +
# geom_smooth(method="lm", aes(fill=NoiseCondition),size=1.5,se=FALSE) +
plot4 = ggplot(exp2data,aes(x=Session, y=NDT)) +
geom_point(aes(shape=NoiseCondition, fill=NoiseCondition),size=3) +
geom_smooth(method="lm", colour='black',size=1.5) +
coord_cartesian(ylim=c(300,700)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_blank()) +
labs(y='Non-decision time estimates (ms)',color='Noise\ncondition') +
scale_x_continuous(name='Training session',breaks=unique(exp2data$Session),lim=c(1,7),labels=unique(exp2data$Session)) +
scale_fill_manual(values = cbbPalette) + 
scale_shape_manual(values = shapePalette) + 
scale_color_manual(values = cbbPalette)
# plot1 + geom_text(x = 200, y = 1200, label = lm_eqn(lm4), parse = TRUE,size=8)
plot(plot4)
dev.off()
confint(lm4, 'Session', level=0.95)

lm5 = lm( X10th.RT.percentiles ~ Session,exp2data)
png('Session_rt10per_exp2.png',units="in",width=10,height=10,res=300)
# plot5 = ggplot(exp2data,aes(x=Session, y=X10th.RT.percentiles, color=NoiseCondition, fill=NoiseCondition)) +
# geom_point(size=3) +
# geom_smooth(method="lm", aes(fill=NoiseCondition),size=1.5,se=FALSE) +
plot5 = ggplot(exp2data,aes(x=Session, y=X10th.RT.percentiles)) +
geom_point(aes(shape=NoiseCondition, fill=NoiseCondition),size=3) +
geom_smooth(method="lm", colour='black',size=1.5) +
coord_cartesian(ylim=c(400,900)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_blank()) +
labs(y='NDT estimates: 10th RT percentiles (ms)',color='Noise\ncondition') +
scale_x_continuous(name='Training session',breaks=unique(exp2data$Session),lim=c(1,7),labels=unique(exp2data$Session)) +
scale_fill_manual(values = cbbPalette) + 
scale_shape_manual(values = shapePalette) + 
scale_color_manual(values = cbbPalette)
# plot1 + geom_text(x = 200, y = 1200, label = lm_eqn(lm5), parse = TRUE,size=8)
plot(plot5)
dev.off()
confint(lm5, 'Session', level=0.95)