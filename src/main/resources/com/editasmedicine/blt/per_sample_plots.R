################################################################################
# Copyright (c) 2017-2018 Editas Medicine, Inc.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification,
# are permitted (subject to the limitations in the disclaimer below) provided
# that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, this
#   list of conditions and the following disclaimer in the documentation and/or
#   other materials provided with the distribution.
# * Neither the name of [Owner Organization] nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
# 
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
# LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
################################################################################

################################################################################
# Generates plots for a single sample from the detailed sample/target ouput
# file for a single sample.
#
# Arguments (positional):
# 1: The input file for a single sample (the .targets.txt.gz / SampleTargetMetrics)
# 2: The output PDF file to write
################################################################################

options(warn = -1)
library(ggplot2)

args    = commandArgs(trailingOnly=T)
input  = args[1]
output = args[2]

# Load and massage the data
data   = read.table(gzfile(input), sep="\t", header=T, stringsAsFactors=F)
sample = unique(data$sample)
data   = subset(data, data$indel_bases == 0)
data   = subset(data, data$mismatches <= 10)
zeroMm = subset(data, data$mismatches == 0)
zeroMmCutRate = sum(zeroMm$obs_cut) / (sum(zeroMm$obs_cut) + sum(zeroMm$obs_uncut))
data$normalized_cut_rate = data$cut_rate / zeroMmCutRate
genomic = subset(data, data$genomic_location != "")

pdf(output)


# Plot of cut rate by mismatch count
title  = paste("Distribution of Normalized Cut Rate vs. Mismatches for", sample)
ggplot(subset(data, data$mismatches <= 10)) + aes(x=factor(mismatches), y=normalized_cut_rate) +
  geom_boxplot() +
  geom_point(data=genomic, color="orange", size=3) + 
  labs(x="Number of Mismatches", y="Normalized Cut Rate", title=title) +
  theme(plot.title = element_text(hjust = 0.5))

# Plot of cut rate by average mismatch position and mismatch count
sub = subset(data, data$mismatches >= 1 & data$mismatches <= 4, select=c("mismatches", "mean_mismatch_position", "obs_cut", "obs_uncut"))
plotdata = aggregate(sub, by=list(sub$mismatches, sub$mean_mismatch_position), FUN=sum)
names(plotdata) = c("mismatches", "mean_mismatch_position", "XXX", "YYY", "obs_cut", "obs_uncut")
plotdata$cut_rate = plotdata$obs_cut / (plotdata$obs_cut + plotdata$obs_uncut)
plotdata$normalized_cut_rate = plotdata$cut_rate / zeroMmCutRate
plotdata$Mismatches = factor(plotdata$mismatches)
  
title = paste("Cut Rate by Mean Mismatch Position in", sample)
ggplot(plotdata) + aes(x=mean_mismatch_position, y=normalized_cut_rate, color=Mismatches) +
    geom_path() +
    labs(x="Mean Mismatch Position", y="Normalized Cut Rate", title=title) +
    theme(plot.title = element_text(hjust = 0.5))

dev.off()
