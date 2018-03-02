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
#  01/30/17        Michael Nunez             Conversion from pdm5b_scatterplots.R

## Necessary packages
library(ggplot2)
library(ggjoy)
library(R.matlab)

## Code

# Read in the data points
recovery = readMat('recovery_changingbetas.mat')

#Color-blind friendly palette (http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Plots
runindex = seq(1,30,length=30)
percentilebetas = as.vector(recovery[3]$percentilebetas)
posteriorbetas =  as.vector(recovery[4]$posteriorbetas)
lm1 = lm(percentilebetas ~ runindex)
lm2 = lm(posteriorbetas ~ runindex)
plotdata1 = data.frame(matrix(ncol = 3, nrow = 60))
plotdata1$betas = c(percentilebetas,posteriorbetas)
plotdata1$runindex = c(runindex,runindex)
plotdata1$Estimate = c(rep("10th Percentile",30),rep("Simple DDM Fit",30))


png('contam_influence.png',units="in",width=10,height=10,res=300)
plot1 = ggplot(plotdata1,aes(x=runindex, y=betas, color=Estimate, shape=Estimate)) +
geom_point() + 
geom_smooth(method="lm", size=1.5) +
scale_x_continuous(name='Percentage of contaminant trials (%)',breaks=as.vector(recovery[2]$converttoindex),lim=c(1,30),labels=as.vector(recovery[1]$percentcontam)) +
coord_cartesian(ylim=c(0.5,1.25)) +
theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face='bold'), legend.text = element_text(size=14),
    legend.title=element_text(size=14,face='bold')) +
labs(y='Relationship between real NDT and NDT estimate (slope parameter)') +
geom_abline(slope=0,intercept=1,colour=cbbPalette[2],size=2,linetype=2)+
scale_fill_manual(palette = tail(cbbPalette,n=2))
plot(plot1)
dev.off()