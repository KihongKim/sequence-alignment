

setwd("/Users/Potenza/Documents/01_Dissertation/Analysis2")

library(TraMineR)

load("wardTree.OM.cons.RData")
load("wardTree.OM.tran.RData")
load("wardTree.HAM.RData")
load("wardTree.DHD.RData")

seqtreedisplay(wardTree.OM.cons, type="d", cex.legend=0.6, border=NA, showdepth=TRUE)
seqtreedisplay(wardTree.OM.tran, type="d", cex.legend=0.6, border=NA, showdepth=TRUE)
seqtreedisplay(wardTree.HAM, type="d", cex.legend=0.6, border=NA, showdepth=TRUE)
seqtreedisplay(wardTree.DHD, type="d", cex.legend=0.6, border=NA, showdepth=TRUE)


