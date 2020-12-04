# pdm5b_scatterplots.R - Creates scatter plots of N200 latencies versus RT percentiles
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
#  01/30/17      Michael Nunez             Conversion from pdm5b_scatterplots.R
#  06/18/18      Michael Nunez             Addition of results with modeled lapse trials
#  06/20/18      Michael Nunez             Addition of results with no cutoffs (no change)

## References:
# http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
# h#ttp://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization

## Necessary packages
library(ggplot2)
library(ggjoy)
library(R.matlab)

## Code

fontsize = 18
cbbPalette <- c("#000000", "#E6AB02", "#A6761D", "#66A61E", "#E7298A", "#7570B3", "#D95F02", "#1B9E77")

# Read in the data points
recovery = readMat('recovery_changingbetas.mat')


## Plots
runindex = seq(1,30,length=30)
posteriorbetas = as.vector(recovery[7]$posteriorbetas)
percentilebetas = as.vector(recovery[6]$percentilebetas)
posteriorbetaslapse350 =  as.vector(recovery[5]$posteriorbetaslapse350)
posteriorbetaslapse =  as.vector(recovery[4]$posteriorbetaslapse)
percentilebetas.cutoffs = as.vector(recovery[3]$percentilebetas.cutoffs)
lm1 = lm(percentilebetas ~ runindex)
lm2 = lm(posteriorbetas ~ runindex)
lm3 = lm(posteriorbetaslapse ~ runindex)
lm4 = lm(posteriorbetaslapse350 ~ runindex)
lm5 = lm(percentilebetas.cutoffs ~ runindex)

plotdata1 = data.frame(matrix(ncol = 3, nrow = 90))
plotdata1$betas = c(percentilebetas,posteriorbetas,posteriorbetaslapse)
plotdata1$runindex = c(runindex,runindex,runindex)
plotdata1$Estimate = c(rep("10th Percentile",30),rep("Simple DDM Fit",30),rep("DDM Fit with Lapses",30))

# plotdata1 = data.frame(matrix(ncol = 5, nrow = 150))
# plotdata1$betas = c(percentilebetas,posteriorbetas,posteriorbetaslapse,posteriorbetaslapse350,percentilebetas.cutoffs)
# plotdata1$runindex = c(runindex,runindex,runindex,runindex,runindex)
# plotdata1$Estimate = c(rep("10th Percentile",30),rep("Simple DDM Fit",30),rep("DDM Fit with Lapses",30),rep("DDM Fit with Lapses & Cutoffs",30),rep("10th Percentile with Cutoffs",30))


png('contam_influence.png',units="in",width=10,height=10,res=300)
plot1 = ggplot(plotdata1,aes(x=runindex, y=betas, color=Estimate, shape=Estimate)) +
geom_point(size=3) + 
geom_smooth(method="lm", size=1.5, aes(fill=Estimate)) +
scale_x_continuous(name='Percentage of contaminant trials (%)',breaks=as.vector(recovery[2]$converttoindex),lim=c(1,30),labels=as.vector(recovery[1]$percentcontam)) +
coord_cartesian(ylim=c(0.5,1.25)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_blank(),legend.position="top") +
labs(y='Relationship between real NDT and NDT estimate (slope parameter)') +
geom_abline(slope=0,intercept=1,colour=cbbPalette[2],size=2,linetype=2)+
scale_color_manual(values = tail(cbbPalette,n=5))+
scale_fill_manual(values = tail(cbbPalette,n=5))
plot(plot1)
dev.off()