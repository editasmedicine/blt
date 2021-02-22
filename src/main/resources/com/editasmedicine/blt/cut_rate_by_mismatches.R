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
# Generates a pair of plots of the cut rate and the normalized cut rate
# by the number of mismatches, with a line for each sample.
#
# Arguments (positional):
# 1  : The output PDF file to create
# 2-n: Per sample input files (the .summary.txt / SampleMetric files)
################################################################################

options(warn = -1)
library(ggplot2)

args    = commandArgs(trailingOnly=T)
output  = args[1]
inputs  = args[2:length(args)]

# Load up all the inputs
first = TRUE

for (input in inputs) {
  tmp = read.table(input, sep="\t", header=T)
  if (first) {
    data = tmp
    first = FALSE
  }
  else {
    data = rbind(data, tmp)
  }
}

# Make the plots
xmax = 10

pdf(output, width=11.5, height=8)
ggplot(data) + aes(x=mismatches, y=cut_rate, color=sample) +
  geom_line() +
  scale_x_continuous(limits=c(0,xmax), breaks=0:xmax) +
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(x="Number of Mismatches", y="Cut Rate", title="Cut Rate By Number of Mismatches") +
  theme(plot.title = element_text(hjust = 0.5))
  
ggplot(data) + aes(x=mismatches, y=normalized_cut_rate, color=sample) +
  geom_line() +
  scale_x_continuous(limits=c(0,xmax), breaks=0:xmax) +
  scale_y_continuous(limits=c(0,max(1, data$normalized_cut_rate)), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(x="Number of Mismatches", y="Cut Rate", title="Normalized Cut Rate By Number of Mismatches") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
