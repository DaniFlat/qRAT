####
# qPCR - Relative Expression Analysis Tool
# version: 0.1.6
#
# PLEASE CITE
# Please cite the published manuscript in all studies using qRAT
# For authors and journal information, please refer to the qRAT website https://uibk.ac.at/microbiology/services/qrat
#
#
# MIT License
#
# Copyright (c) 2022 Daniel Flatschacher
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
####

####
# BESTKEEPER
####


## using ctrlGene package

housekeeperAnalysis <- function(x) {


  if (is.null(dt)) NULL

  data <- read.table("fileSingle.csv", sep = ",")
  data <- select(data, -c(Well, rp.num))
  print(data)
  df1_new<-dcast(data, Sample ~ Gene)

  print(df1_new)

}
