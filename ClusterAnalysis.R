setwd("/home/kihong/Dissertation/Analysis2") # on the sapporo server
setwd("/Users/Potenza/Documents/01_Dissertation/Analysis2") # on my mac air
options(scipen=999) # to remove any scientific notations

library(TraMineR)
library(WeightedCluster)
library(nnet)

### load data
# load("diary.RData")					# diaries in SPELL
# load("diary.sts.RData")				# diaries in STS
load("diary.sts.ivs.RData")				# diaries in STS with IVs
load("diary.seq288.RData")				# a sequence object of 5-min time intervals
load("diary.seq288.OM.cons.RData")		# a pairwise distance matrix by OM with a constant substitution cost
# load("diary.seq288.OM.tran.RData")	# a pairwise distance matrix by OM with transitioning substitution costs
# load("diary.seq288.DHD.RData")		# a pairwise distance matrix by DHD

####################
# Cluster Analysis #
####################

# >>> Hierachical Agglomerative Clustering
# 1. each of the observations is considered as a group
# 2. on each iteration, the two closet groups are grouped, until all the observations form a single group
# 3. the agglomeration schedule (i.e., the succession of groupings performed) can be represented in a dendrogram
# 4. with the dendrogram, the number of groups is selected
# 5. cut the grouping tree at the corresponding level

# create the agglomeration schedule
wardClust.OM.cons <- hclust(as.dist(diary.seq288.OM.cons), method = "ward")
wardClust.OM.tran <- hclust(as.dist(diary.seq288.OM.tran), method = "ward")
wardClust.HAM <- hclust(as.dist(diary.seq288.HAM), method = "ward")
wardClust.DHD <- hclust(as.dist(diary.seq288.DHD), method = "ward")

png("dendrogram.png")
plot(wardCluster)
dev.off()


# creat tree
wardTree.OM.cons <- as.seqtree(wardClust.OM.cons, seqdata = diary.seq288, diss = diary.seq288.OM.cons, ncluster = 10)
wardTree.OM.tran <- as.seqtree(wardClust.OM.tran, seqdata = diary.seq288, diss = diary.seq288.OM.tran, ncluster = 10)
wardTree.HAM <- as.seqtree(wardClust.HAM, seqdata = diary.seq288, diss = diary.seq288.HAM, ncluster = 10)
wardTree.DHD <- as.seqtree(wardClust.DHD, seqdata = diary.seq288, diss = diary.seq288.DHD, ncluster = 10)

save(wardTree.OM.cons, file="wardTree.OM.cons.RData")
save(wardTree.OM.tran, file="wardTree.OM.tran.RData")
save(wardTree.HAM, file="wardTree.HAM.RData")
save(wardTree.DHD, file="wardTree.DHD.RData")


# to create a tree, run 03_seqtreedisplay on my macbook air via GraphViz



# cluster quality
wardRange.OM.cons <- as.clustrange(wardClust.OM.cons, diss = diary.seq288.OM.cons, ncluster = 10)
wardRange.OM.tran <- as.clustrange(wardClust.OM.tran, diss = diary.seq288.OM.tran, ncluster = 10)
wardRange.HAM <- as.clustrange(wardClust.HAM, diss = diary.seq288.HAM, ncluster = 10)
wardRange.DHD <- as.clustrange(wardClust.DHD, diss = diary.seq288.DHD, ncluster = 10)

save(wardRange.OM.cons, file="wardRange.OM.cons.RData")
save(wardRange.OM.tran, file="wardRange.OM.tran.RData")
save(wardRange.HAM, file="wardRange.HAM.RData")
save(wardRange.DHD, file="wardRange.DHD.RData")

# select the number of clusters by statistics
pdf("wardRange.pdf")
par(mfrow=c(2,2))
plot(wardRange.OM.cons, stat=c("ASWw","HG","PBC","HC"), main="OM with CONSTANT")
plot(wardRange.OM.tran, stat=c("ASWw","HG","PBC","HC"), main="OM with TRATE")
plot(wardRange.HAM, stat=c("ASWw","HG","PBC","HC"), main="HAM")
plot(wardRange.DHD, stat=c("ASWw","HG","PBC","HC"), main="DHD")
dev.off()

pdf("wardRangeZ.pdf")
par(mfrow=c(2,2))
plot(wardRange.OM.cons, stat=c("ASWw","HG","PBC","HC"), main="OM with CONSTANT", norm="zscore")
plot(wardRange.OM.tran, stat=c("ASWw","HG","PBC","HC"), main="OM with TRATE", norm="zscore")
plot(wardRange.HAM, stat=c("ASWw","HG","PBC","HC"), main="HAM", norm="zscore")
plot(wardRange.DHD, stat=c("ASWw","HG","PBC","HC"), main="DHD", norm="zscore")
dev.off()

summary(wardRange.OM.cons, max.rank=3)
summary(wardRange.OM.tran, max.rank=3)
summary(wardRange.OM.HAM, max.rank=3)
summary(wardRange.OM.DHD, max.rank=3)

# show activity state distributions by the number of clusters
# 5 clusters
pdf("wardClust5.OM.cons.pdf")
seqdplot(diary.seq288, group = wardRange.OM.cons$clustering$cluster5, border = NA, title="OM with CONSTANT")
dev.off()

pdf("wardClust5.OM.tran.pdf")
seqdplot(diary.seq288, group = wardRange.OM.tran$clustering$cluster5, border = NA, title="OM with TRATE")
dev.off()

pdf("wardClust5.HAM.pdf")
seqdplot(diary.seq288, group = wardRange.HAM$clustering$cluster5, border = NA, title="HAM")
dev.off()

pdf("wardClust5.DHD.pdf")
seqdplot(diary.seq288, group = wardRange.DHD$clustering$cluster5, border = NA, title="DHD")
dev.off()

# 6 clusters
pdf("wardClust6.OM.cons.pdf")
seqdplot(diary.seq288, group = wardRange.OM.cons$clustering$cluster6, border = NA, title="OM with CONSTANT")
dev.off()

pdf("wardClust6.OM.tran.pdf")
seqdplot(diary.seq288, group = wardRange.OM.tran$clustering$cluster6, border = NA, title="OM with TRATE")
dev.off()

pdf("wardClust6.HAM.pdf")
seqdplot(diary.seq288, group = wardRange.HAM$clustering$cluster6, border = NA, title="HAM")
dev.off()

pdf("wardClust6.DHD.pdf")
seqdplot(diary.seq288, group = wardRange.DHD$clustering$cluster6, border = NA, title="DHD")
dev.off()

# 7 clusters
pdf("wardClust7.OM.cons.pdf")
seqdplot(diary.seq288, group = wardRange.OM.cons$clustering$cluster7, border = NA, title="OM with CONSTANT")
dev.off()

pdf("wardClust7.OM.tran.pdf")
seqdplot(diary.seq288, group = wardRange.OM.tran$clustering$cluster7, border = NA, title="OM with TRATE")
dev.off()

pdf("wardClust7.HAM.pdf")
seqdplot(diary.seq288, group = wardRange.HAM$clustering$cluster7, border = NA, title="HAM")
dev.off()

pdf("wardClust7.DHD.pdf")
seqdplot(diary.seq288, group = wardRange.DHD$clustering$cluster7, border = NA, title="DHD")
dev.off()


# 8 clusters
pdf("wardClust8.OM.cons.pdf")
seqdplot(diary.seq288, group = wardRange.OM.cons$clustering$cluster8, border = NA, title="OM with CONSTANT", withlegend=FALSE)
dev.off()

pdf("wardClust8.OM.tran.pdf")
seqdplot(diary.seq288, group = wardRange.OM.tran$clustering$cluster8, border = NA, title="OM with TRATE", withlegend=FALSE)
dev.off()

pdf("wardClust8.HAM.pdf")
seqdplot(diary.seq288, group = wardRange.HAM$clustering$cluster8, border = NA, title="HAM", withlegend=FALSE)
dev.off()

pdf("wardClust8.DHD.pdf")
seqdplot(diary.seq288, group = wardRange.DHD$clustering$cluster8, border = NA, title="DHD", withlegend=FALSE)
dev.off()

# 9 clusters
pdf("wardClust9.OM.cons.pdf")
seqdplot(diary.seq288, group = wardRange.OM.cons$clustering$cluster9, border = NA, title="OM with CONSTANT", withlegend=FALSE)
dev.off()

pdf("wardClust9.OM.tran.pdf")
seqdplot(diary.seq288, group = wardRange.OM.tran$clustering$cluster9, border = NA, title="OM with TRATE", withlegend=FALSE)
dev.off()

pdf("wardClust9.HAM.pdf")
seqdplot(diary.seq288, group = wardRange.HAM$clustering$cluster9, border = NA, title="HAM", withlegend=FALSE)
dev.off()

pdf("wardClust9.DHD.pdf")
seqdplot(diary.seq288, group = wardRange.DHD$clustering$cluster9, border = NA, title="DHD", withlegend=FALSE)
dev.off()

# 10 clusters
pdf("wardClust10.OM.cons.pdf")
seqdplot(diary.seq288, group = wardRange.OM.cons$clustering$cluster10, border = NA, title="OM with CONSTANT", withlegend=FALSE)
dev.off()

pdf("wardClust10.OM.tran.pdf")
seqdplot(diary.seq288, group = wardRange.OM.tran$clustering$cluster10, border = NA, title="OM with TRATE", withlegend=FALSE)
dev.off()

pdf("wardClust10.HAM.pdf")
seqdplot(diary.seq288, group = wardRange.HAM$clustering$cluster10, border = NA, title="HAM", withlegend=FALSE)
dev.off()

pdf("wardClust10.DHD.pdf")
seqdplot(diary.seq288, group = wardRange.DHD$clustering$cluster10, border = NA, title="DHD", withlegend=FALSE)
dev.off()




#
wardClust5.OM.cons <- cutree(wardClust.OM.cons, k=5)
wardClust5.OM.tran <- cutree(wardClust.OM.tran, k=5)
wardClust5.HAM <- cutree(wardClust.HAM, k=5)
wardClust5.DHD <- cutree(wardClust.DHD, k=5)

wardClust6.OM.cons <- cutree(wardClust.OM.cons, k=6)
wardClust6.OM.tran <- cutree(wardClust.OM.tran, k=6)
wardClust6.HAM <- cutree(wardClust.HAM, k=6)
wardClust6.DHD <- cutree(wardClust.DHD, k=6)

pdf("wardClust5.OM.cons.pdf"); seqdplot(diary.seq288, group = wardClust5.OM.cons, border = NA); dev.off()
pdf("wardClust5.OM.tran.pdf"); seqdplot(diary.seq288, group = wardClust5.OM.tran, border = NA); dev.off()
pdf("wardClust5.HAM.pdf"); seqdplot(diary.seq288, group = wardClust5.HAM, border = NA); dev.off()
pdf("wardClust5.DHD.pdf"); seqdplot(diary.seq288, group = wardClust5.DHD, border = NA); dev.off()

cq5.OM.cons <- wcClusterQuality(diary.seq288.OM.cons, wardClust5.OM.cons)
cq5.OM.cons$stats
cq5.OM.tran <- wcClusterQuality(diary.seq288.OM.tran, wardClust5.OM.tran)
cq5.OM.tran$stats
cq5.HAM <- wcClusterQuality(diary.seq288.HAM, wardClust5.HAM)
cq5.HAM$stats
cq5.DHD <- wcClusterQuality(diary.seq288.DHD, wardClust5.DHD)
cq5.DHD$stats






# Linking trajectory types and explanatory factors
tb5 <- with(diary.sts.ivs, table(h_size, wardClust5.OM.cons))
chisq.test(tb5)

tb6 <- with(diary.sts.ivs, table(h_size, wardClust6.OM.cons))
chisq.test(tb6)

ds <- cbind(diary.sts.ivs,
			wardClust5.OM.cons,
			wardClust5.OM.tran,
			wardClust5.HAM,
			wardClust5.DHD,
			wardClust6.OM.cons,
			wardClust6.OM.tran,
			wardClust6.HAM,
			wardClust6.DHD)

ds.wardClust6.DHD <- mlogit.data(ds, shape="wide", choice="wardClust6.DHD")

mnl.wardClust6.DHD <- mlogit(wardClust6.DHD ~ 0 | 
	p_gender+p_age+p_disabled+p_licensed+p_transitpass+p_biker+
	p_worker+p_retiree+p_homemaker+p_ftcollege+p_ptcollege+p_k12+p_g9to12+p_child0to4+p_adult65over,
	data=ds.wardClust6.DHD)
					



MNL <- multinom(
	wardClust6.DHD ~ 
		p_gender+p_age+p_disabled+p_licensed+p_transitpass+p_biker+
		p_worker+p_retiree+p_homemaker+p_ftcollege+p_ptcollege+p_k12+p_g9to12+p_child0to4+p_adult65over+
		h_size+h_vehicles+h_zerovehicle+h_lowincome+h_highincome+h_white+
		h_anykids0to4+h_anykids5to14+h_anykids15to17+h_numkids0to4+h_numkids5to14+h_numkids15to17+
		h_1a0k+h_1awk+h_2a0k+h_2awk,
	data=ds)


		s_multnomah+s_clackamas+s_washington+
		s_beaverton+s_portland+s_hillsboro+s_gresham+s_oregoncity+s_downtown,


MNL$prog2 <- relevel(ml$prog, ref = "academic")
MNL <- multinom(prog2 ~ ses + write, data = ml)
summary(test)


# >>> Partitioning Around Medoids



pamclust4 <- wcKMedoids(diary.seq288.OM.cons, k=4)
png("pamclust4.png")
seqdplot(diary.seq288, group = pamclust4$clustering, border=NA)
dev.off()








