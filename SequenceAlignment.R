# sequence alignment followed by ANODI (Analysis of Dissimilarities) and tree-structured analysis
# by Kihong Kim

# 1. load data and manipulation
# 2. transform from SPELL to STS (STate Sequence)
# 3. sequence alignment aka optimal matching
# 4. ANODI
# 5. tree-structured analysis

###########
# 0. set up
###########
setwd("/home/kihong/Dissertation/Analysis2") # on the sapporo server
options(scipen=999) # to remove any scientific notations

library(TraMineR)
library(sqldf)
library(reshape)

###############################
# 1. load data and manipualtion
###############################
# load data
activity <- read.csv("../Data/ACTIVITY.csv", header=T)
activity <- with(activity, activity[order(SAMPN,PERNO,PLANO),])
save(activity, file="activity.RData")
diary <- activity

# select households only in Metro (AREA=11; 4799), not RTC (AREA=21; 1650)
diary <- subset(diary, SAMPN>8000000)

# add 'arrive' and 'depart' on the 1440-minute time scale
diary$arrive <- with(diary, ARR_HR*60+ARR_MIN)
diary$depart <- with(diary, DEP_HR*60+DEP_MIN)
# correct 'arrive' and 'depart' for activities after midnight
diary$arrive <- with(diary, ifelse(arrive<180, arrive+1440, arrive))
diary$depart <- with(diary, ifelse(depart<180, depart+1440, depart))

# add 'actNo'
diary$actNo <- with(diary, ave(PLANO, SAMPN,PERNO, FUN=seq))
### tip: 'ave' is to do group averages over the combined levels of factors

# add 'maxAct'
diary$maxAct <- with(diary, ave(actNo, SAMPN,PERNO, FUN=length))

# add 'lastDepart'
diary$lastDepart <- c(NA, diary$depart[-length(diary$depart)])
diary$lastDepart[diary$actNo==1] <- NA
### tip: -length() is all elements except for the last observation

# add 'tripDur'
diary$tripDur <- with(diary, arrive-lastDepart)

# add 'actDur'
diary <- rename(diary, c(ACTDUR="ACTDUR_OLD"))
diary$actDur <- with(diary, depart-arrive)

# add 'uid'
diary$uid <- with(diary, SAMPN*10+PERNO) # maximum number of household persons is 8

# define the sample
# delete persons who have no trip during the day
NoTrip <- sqldf("select uid from diary where actNo==1 and maxAct==1")
diary <- subset(diary, !(diary$uid %in% NoTrip$uid))
# delete persons who do not start at home
NotAtHome180 <- sqldf("select uid from diary where actNo==1 and TPURP!=1 and TPURP!=2")
diary <- subset(diary, !(diary$uid %in% NotAtHome180$uid))
# delete persons who do not end at home
NotAtHome1619 <- sqldf("select uid from diary where actNo==maxAct and TPURP!=1 and TPURP!=2")
diary <- subset(diary, !(diary$uid %in% NotAtHome1619$uid))
# delete persons whose depart!=1619 in the end
NotDepart1619 <- sqldf("select uid from diary where actNo==maxAct and depart!=1619")
diary <- subset(diary, !(diary$uid %in% NotDepart1619$uid))

# aggregate activity types, partially following MIT's aggregation (thisActG1)
# load the table of activity aggregation
ActivityAggregation = read.csv("../Data/ActivityAggregation1.csv", stringsAsFactors=FALSE)
# aggregate activities
diary$thisActG = ActivityAggregation$ActG[match(diary$TPURP, ActivityAggregation$Code)]
diary$thisActG_Name = ActivityAggregation$ActG_Name[match(diary$TPURP, ActivityAggregation$Code)]
# disaggregate home activities into HB, HR, and HE
attach(diary)
diary$thisActG[thisActG=="HM" & actNo==1] <- "HB"
diary$thisActG[thisActG=="HM" & actNo!=1 & actNo!=maxAct] <- "HR"
diary$thisActG[thisActG=="HM" & actNo!=1 & actNo==maxAct] <- "HE"
diary$thisActG_Name[thisActG=="HB"] <- "HomeBegining"
diary$thisActG_Name[thisActG=="HR"] <- "HomeMiddle"
diary$thisActG_Name[thisActG=="HE"] <- "HomeEnd"
detach(diary)

# compute average duration by thisActG
sqldf("select thisActG,count(*),round(avg(actDur),1),round(stdev(actDur),1),median(actDur),min(actDur),max(actDur) from diary group by thisActG")

# delete pesons with negative activity duration (1 case)
diary <- subset(diary, uid!=80791657)

# delete pesons with too long duration for 'CM' (1 case)
diary <- subset(diary, uid!=80466342)

#
save(diary, file="diary.RData")

###############################
# 2. transform from SPEE to STS
###############################
# add BEGIN and END to transform diaries from SPELL to STS
# note: the activity type of individuals while traveling is equal to be that of their destination activity type
diary$BEGIN <- with(diary, ifelse(actNo==1, arrive, lastDepart)) 
diary$END <- with(diary, depart)

# 'diary.sts': transform from SPELL to STS
diary.sts <- seqformat(diary,
  id="uid", begin="BEGIN", end="END", status="thisActG",
  from="SPELL", to="STS", process=FALSE)

names(diary.sts) <- gsub("y","t",names(diary.sts)) # change column names from y... to t...
save(diary.sts, file="diary.sts.RData")

# 'diary.sts.ivs': attach 'household' and 'person' tables to diary.sts
person <- read.csv("../Data/PER.csv", header=T)
household <- read.csv("../Data/HH.csv", header=T)
ivs <- merge(person, household, by="SAMPN", all.x=T)
ivs$uid <- with(ivs, SAMPN*10+PERNO)
ivs <- subset(ivs, (ivs$uid %in% diary$uid))
ivs <- with(ivs, ivs[order(SAMPN,PERNO),])
diary.sts.ivs <- cbind(diary.sts, ivs)
save(diary.sts.ivs, file="diary.sts.ivs.RData")

# 'diary.seq1440': create a sequence object in 1-min time intervals
diary.alphabet <- c("HM","HR","WK","WR","SC","CM","ES","EO","SH","HE","PB","SR","OT")
diary.label <- c("home","home returning temporarily",
                 "work","work related","school",
                 "change mode","escort","eat out",
                 "shopping","household errands","personal business",
                 "social recreation","other")
diary.seq1440 <- seqdef(diary.sts.ivs, 1:1440, alphabet=diary.alphabet, labels=diary.label)
attr(diary.seq1440, "cpal") <- c("light grey","grey",
								 "blue","cyan","yellow",
								 "black","purple","orange",
								 "red","violet","golden rod",
								 "green","white")
summary(diary.seq1440)
save(diary.seq1440, file="diary.seq1440.RData")

# 'diary.seq288': create a sequence object in 5-min time intervals
diary.seq288 <- diary.seq1440[,seq(from=1, to=1440, by=5)]
save(diary.seq288, file="diary.seq288.RData")

# graphs
pdf("diary.seq288.seq10.pdf")
seqiplot(diary.seq288, border=NA, cex.legend=0.6, title="Activity Sequence Index Plot (first 10 sequences)")
dev.off()
pdf("diary.seq288.seqAll.pdf")
seqIplot(diary.seq288, border=NA, cex.legend=0.6, title="Activity Sequence Index Plot (all sequences)")
dev.off()
pdf("diary.seq288.seqD.pdf")
seqdplot(diary.seq288, border=NA, cex.legend=0.6, title="Activity State Distribution Plot")
dev.off()

print(diary.seq288[9:10,], format="SPS")

###############################################
# 3. sequence alignment a.k.a. optimal matching
###############################################

# produce a pairwise OM distance matrix with indels=1 & substitutions=2
diary.seq288.OM.cons <- seqdist(diary.seq288, method="OM", indel=1, sm="CONSTANT", full.matrix=FALSE)
save(diary.seq288.OM.cons, file="diary.seq288.OM.cons.RData")

# produce a pairwise OM distance matrix with indels=1 & transitioning substitutions
diary.seq288.OM.tran <- seqdist(diary.seq288, method="OM", indel=1, sm="TRATE", full.matrix=FALSE)
save(diary.seq288.OM.tran, file="diary.seq288.OM.tran.RData")

# produce a pairwise HAM matrix (Hamming distance)
diary.seq288.HAM <- seqdist(diary.seq288, method="HAM", full.matrix=FALSE)
save(diary.seq288.HAM, file="diary.seq288.HAM.RData")

# produce a pairwise DHD matrix (dynamic Hamming distance)
diary.seq288.DHD <- seqdist(diary.seq288, method="DHD", full.matrix=FALSE)
save(diary.seq288.DHD, file="diary.seq288.DHD.RData")



