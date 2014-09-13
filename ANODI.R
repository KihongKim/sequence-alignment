# ANODI (Anlysis of Dissimilarities) after sequence alignment
# by Kihong Kim


setwd("/home/kihong/Dissertation/Analysis2") # on the sapporo server

library(TraMineR)

### load data for ANODI
# load("diary.RData")					# diaries in SPELL (128122 rows / 42013 rows)
# load("diary.sts.RData")				# diaries in STS (15286 rows / 4971 rows)
load("diary.sts.ivs.RData")				# diaries in STS with IVs (15286 rows / 4968 rows)
# load("diary.seq1440.RData")			# a sequence objcet of 1-min interval diaries in STS with IVs (15286 rows / 4968 rows)
load("diary.seq288.RData")				# a sequence objcet of 5-min interval diaries in STS with IVs (15286 rows / 4968 rows)

load("diary.seq288.OM.cons.RData")		# a pairwise distance matrix by OM with constant substitution costs
load("diary.seq288.OM.tran.RData")		# a pairwise distance matrix by OM with transitioning substitution costs
load("diary.seq288.HAM.RData")			# a pairwise distance matrix by HAM (Hamming distance)
load("diary.seq288.DHD.RData")			# a pairwise distance matrix by DHD (dynamic Hamming distance)


###########################################
# to create independent variables for ANODI
###########################################
attach(diary.sts.ivs)

# individuals' person characteristics
diary.sts.ivs$p_gender[GEND==1] <- 0													# gender - male
diary.sts.ivs$p_gender[GEND==2] <- 1													# gender - female

diary.sts.ivs$p_age[AGE<5] <- 1															# age - under 5 years
diary.sts.ivs$p_age[AGE>=5&AGE<=14] <- 2												# age - 5 to 14 years
diary.sts.ivs$p_age[AGE>=15&AGE<=17] <- 3												# age - 15 to 17 years
diary.sts.ivs$p_age[AGE>=18&AGE<=29] <- 4												# age - 18 to 29 years
diary.sts.ivs$p_age[AGE>=30&AGE<=49] <- 5												# age - 30 to 49 years
diary.sts.ivs$p_age[AGE>=50&AGE<=64] <- 6												# age - 50 to 64 years
diary.sts.ivs$p_age[AGE>=65] <- 7														# age - over 65 years
diary.sts.ivs$p_age[AGE==99] <- NA  													# age - unknown

diary.sts.ivs$p_disabled <- 0							  								# disability - no
diary.sts.ivs$p_disabled[DISAB==1] <- 1													# disability - yes

diary.sts.ivs$p_licensed <- 0															# driver's license - no
diary.sts.ivs$p_licensed[LIC==1] <- 1													# driver's license - yes

diary.sts.ivs$p_transitpass <- 0						  								# transit pass - no
diary.sts.ivs$p_transitpass[TRANS==1] <- 1												# transit pass - yes

diary.sts.ivs$p_biker <- 0							      								# regular biker - no
diary.sts.ivs$p_biker[PBIKE==1] <- 1													# regular biker - yes

diary.sts.ivs$p_employed <- 0															# employed either full-time or part-time - no
diary.sts.ivs$p_employed[EMPLY==1] <- 1													# employed either full-time or part-time - yes

diary.sts.ivs$p_worker <- 0							    								# work either employed or volunteer - no
diary.sts.ivs$p_worker[WORKS==1] <- 1													# work either employed or volunteer - yes

diary.sts.ivs$p_retiree <- 0							    							# retiree - no
diary.sts.ivs$p_retiree[WKSTAT==1] <- 1													# retiree - yes

diary.sts.ivs$p_homemaker <- 0						    								# homemaker - no
diary.sts.ivs$p_homemaker[WKSTAT==3] <- 1												# homemaker - yes

diary.sts.ivs$p_ftcollege <- 0								                            # full-time college student - no
diary.sts.ivs$p_ftcollege[STUDE==1&(SCHOL==6|SCHOL==7|SCHOL==8)] <- 1					# full-time college student - yes

diary.sts.ivs$p_ptcollege <- 0															# part-time college student - no
diary.sts.ivs$p_ptcollege[STUDE==2&(SCHOL==6|SCHOL==7|SCHOL==8)] <- 1					# part-time college student - yes

diary.sts.ivs$p_k12 <- 0								      							# K-12 student - no
diary.sts.ivs$p_k12[SCHOL==3|SCHOL==4] <- 1												# K-12 student - yes

diary.sts.ivs$p_g9to12 <- 0							    								# grade 9 to 12 student - no
diary.sts.ivs$p_g9to12[SCHOL==4] <- 1													# grade 9 to 12 student - yes

diary.sts.ivs$p_child0to4 <- 0						    								# child under 5 years - no
diary.sts.ivs$p_child0to4[AGE<5] <- 1													# child under 5 years - yes

diary.sts.ivs$p_adult65over <- 0					      	  							# adult over 65 years - no
diary.sts.ivs$p_adult65over[AGE>=65&AGE<=95] <- 1										# adult over 65 years - yes

# individuals' household characteristics
diary.sts.ivs$h_size[HHSIZ==1] <- 1					# 1 person in household
diary.sts.ivs$h_size[HHSIZ==2] <- 2					# 2 persons in household
diary.sts.ivs$h_size[HHSIZ==3] <- 3					# 3 persons in household
diary.sts.ivs$h_size[HHSIZ==4] <- 4					# 4 persons in household
diary.sts.ivs$h_size[HHSIZ>=5] <- 5					# 5 or more persons in household

diary.sts.ivs$h_vehicles[HHVEH==0] <- 0			# 0 vehicles in household
diary.sts.ivs$h_vehicles[HHVEH==1] <- 1			# 1 vehicle in household
diary.sts.ivs$h_vehicles[HHVEH==2] <- 2			# 2 vehicles in household
diary.sts.ivs$h_vehicles[HHVEH==3] <- 3			# 3 vehicles in household
diary.sts.ivs$h_vehicles[HHVEH>=4] <- 4			# 4 or more vehicles in household

diary.sts.ivs$h_zerovehicle <- 0						  # zero vehicle in household - no
diary.sts.ivs$h_zerovehicle[HHVEH==0] <- 1		# zero vehicle in household - yes

diary.sts.ivs$h_lowincome <- 0							          # low-income (less than $25,000) household - no
diary.sts.ivs$h_lowincome[INCOME==1|INCOME==2] <- 1	# low-income (less than $25,000) household - yes

diary.sts.ivs$h_highincome <- 0											         	 # high-income (over $75,000) household - no
diary.sts.ivs$h_highincome[INCOME==6|INCOME==7|INCOME==8] <- 1	 # high-income (over $75,000) household - yes

diary.sts.ivs$h_white <- 0								    # white household - no
diary.sts.ivs$h_white[RACE==2] <- 1					# white household - yes

diary.sts.ivs$h_workers[HHWRK==0] <- 0				# 0 workers in household
diary.sts.ivs$h_workers[HHWRK==1] <- 1				# 1 worker in household
diary.sts.ivs$h_workers[HHWRK==2] <- 2				# 2 workers in household
diary.sts.ivs$h_workers[HHWRK>=3] <- 3				# 3 or more workers in household

diary.sts.ivs$h_anykids0to4 <- with(diary.sts.ivs, ave(AGE, SAMPN, FUN=function(x) any(x<=4)))				  # household with children <= 4 years
diary.sts.ivs$h_anykids5to14 <- with(diary.sts.ivs, ave(AGE, SAMPN, FUN=function(x) any(x>=5&x<=14)))		# household with children 5 to 14 years
diary.sts.ivs$h_anykids15to17 <- with(diary.sts.ivs, ave(AGE, SAMPN, FUN=function(x) any(x>=15&x<=17)))	# household with children 15 to 17 years
diary.sts.ivs$h_anykids0to17 <- with(diary.sts.ivs, ave(AGE, SAMPN, FUN=function(x) any(x<=17)))			  # household with children <= 17 years

diary.sts.ivs$h_numkids0to4 <- with(diary.sts.ivs, ave(AGE, SAMPN, FUN=function(x) length(which(x<=4))))			    # num of children <= 4 years
diary.sts.ivs$h_numkids5to14 <- with(diary.sts.ivs, ave(AGE, SAMPN, FUN=function(x) length(which(x>=5&x<=14))))		# num of children 5 to 14 years
diary.sts.ivs$h_numkids15to17 <- with(diary.sts.ivs, ave(AGE, SAMPN, FUN=function(x) length(which(x>=15&x<=17))))	# num of children 15 to 17 years
diary.sts.ivs$h_numkids0to17 <- with(diary.sts.ivs, ave(AGE, SAMPN, FUN=function(x) length(which(x<=17))))			  # num of children <= 17 years

diary.sts.ivs$h_1a0k <- 0									# single without children under 18 years - no
diary.sts.ivs$h_1a0k[HHSIZ==1] <- 1 						# single without children under 18 years - yes

diary.sts.ivs$h_1awk <- 0												# single with children under 18 years - no
diary.sts.ivs$h_1awk[(HHSIZ>1)&((HHSIZ-diary.sts.ivs$h_numkids0to17)==1)] <- 1		# single with children under 18 years - yes

diary.sts.ivs$h_2a0k <- 0									# two or more adults without children under 18 years - no
diary.sts.ivs$h_2a0k[(HHSIZ>1)&(diary.sts.ivs$h_numkids0to17==0)] <- 1		# two or more adults without children under 18 years - yes

diary.sts.ivs$h_2awk <- 0																# two or more adults with children under 18 years - no
diary.sts.ivs$h_2awk[(HHSIZ>1)&((HHSIZ-diary.sts.ivs$h_numkids0to17)>1)&(diary.sts.ivs$h_numkids0to17>0)] <- 1		# two or more adults with children under 18 years - yes

diary.sts.ivs$h_type[diary.sts.ivs$h_1a0k==1] <- 1
diary.sts.ivs$h_type[diary.sts.ivs$h_1awk==1] <- 2
diary.sts.ivs$h_type[diary.sts.ivs$h_2a0k==1] <- 3
diary.sts.ivs$h_type[diary.sts.ivs$h_2awk==1] <- 4

# individuals' activity-travel characteristics

# individuals' residential location characteristics
diary.sts.ivs$s_multnomah <- 0							# household residing in Multnomah County - no
diary.sts.ivs$s_multnomah[CTFIP==41051] <- 1				# household residing in Multnomah County - yes

diary.sts.ivs$s_clackamas <- 0							# household residing in Clackamas County - no
diary.sts.ivs$s_clackamas[CTFIP==41005] <- 1				# household residing in Clackamas County - yes

diary.sts.ivs$s_washington <- 0							# household residing in Washington County - no
diary.sts.ivs$s_washington[CTFIP==41067] <- 1				# household residing in Washington County - yes

diary.sts.ivs$s_beaverton <- 0							            # household residing in Beaverton City - no
diary.sts.ivs$s_beaverton[HCITY=="BEAVERTON"] <- 1			# household residing in Beaverton City - yes

diary.sts.ivs$s_portland <- 0								            # household residing in Portland City - no
diary.sts.ivs$s_portland[HCITY=="PORTLAND"] <- 1			  # household residing in Portland City - yes

diary.sts.ivs$s_hillsboro <- 0  							          # household residing in Hillsboro City - no
diary.sts.ivs$s_hillsboro[HCITY=="HILLSBORO"] <- 1	    # household residing in Hillsboro City - yes

diary.sts.ivs$s_gresham <- 0    						            # household residing in Gresham City - no
diary.sts.ivs$s_gresham[HCITY=="GRESHAM"] <- 1	        # household residing in Gresham City - yes

diary.sts.ivs$s_oregoncity <- 0      					          # household residing in Oregon City - no
diary.sts.ivs$s_oregoncity[HCITY=="OREGON CITY"] <- 1   # household residing in Oregon City - yes

diary.sts.ivs$s_downtown <- 0        				            # household residing in Portland Downtown - no
diary.sts.ivs$s_downtown[HZIP==97201|HZIP==97205|HZIP==97209] <- 1   # household residing in Portland Downtown - yes

detach(diary.sts.ivs)



#####################################
# ANODI (analysis of dissimilarities)
#####################################

system.time(diary.seq288.OM.cons.ANODI <- dissmfacw(
	diary.seq288.OM.cons ~
		p_gender+p_age+p_disabled+p_licensed+p_transitpass+p_biker+
		p_worker+p_retiree+p_homemaker+p_ftcollege+p_ptcollege+p_k12+p_g9to12+p_child0to4+p_adult65over+
		h_size+h_vehicles+h_zerovehicle+h_lowincome+h_highincome+h_white+
		h_anykids0to4+h_anykids5to14+h_anykids15to17+h_numkids0to4+h_numkids5to14+h_numkids15to17+
		h_1a0k+h_1awk+h_2a0k+h_2awk+
		s_multnomah+s_clackamas+s_washington+
		s_beaverton+s_portland+s_hillsboro+s_gresham+s_oregoncity+s_downtown,
	data=diary.sts.ivs, R=1000))
save(diary.seq288.OM.cons.ANODI, file="diary.seq288.OM.cons.ANODI.RData")

diary.seq288.OM.tran.ANODI <- dissmfacw(
	diary.seq288.OM.tran ~
		p_gender+p_age+p_disabled+p_licensed+p_transitpass+p_biker+
		p_worker+p_retiree+p_homemaker+p_ftcollege+p_ptcollege+p_k12+p_g9to12+p_child0to4+p_adult65over+
		h_size+h_vehicles+h_zerovehicle+h_lowincome+h_highincome+h_white+
		h_anykids0to4+h_anykids5to14+h_anykids15to17+h_numkids0to4+h_numkids5to14+h_numkids15to17+
		h_1a0k+h_1awk+h_2a0k+h_2awk+
		s_multnomah+s_clackamas+s_washington+
		s_beaverton+s_portland+s_hillsboro+s_gresham+s_oregoncity+s_downtown,
	data=diary.sts.ivs, R=1000)
save(diary.seq288.OM.tran.ANODI, file="diary.seq288.OM.tran.ANODI.RData")

diary.seq288.HAM.ANODI <- dissmfacw(
	diary.seq288.HAM ~
		p_gender+p_age+p_disabled+p_licensed+p_transitpass+p_biker+
		p_worker+p_retiree+p_homemaker+p_ftcollege+p_ptcollege+p_k12+p_g9to12+p_child0to4+p_adult65over+
		h_size+h_vehicles+h_zerovehicle+h_lowincome+h_highincome+h_white+
		h_anykids0to4+h_anykids5to14+h_anykids15to17+h_numkids0to4+h_numkids5to14+h_numkids15to17+
		h_1a0k+h_1awk+h_2a0k+h_2awk+
		s_multnomah+s_clackamas+s_washington+
		s_beaverton+s_portland+s_hillsboro+s_gresham+s_oregoncity+s_downtown,
	data=diary.sts.ivs, R=1000)
save(diary.seq288.HAM.ANODI, file="diary.seq288.HAM.ANODI.RData")

diary.seq288.DHD.ANODI <- dissmfacw(
	diary.seq288.DHD ~
		p_gender+p_age+p_disabled+p_licensed+p_transitpass+p_biker+
		p_worker+p_retiree+p_homemaker+p_ftcollege+p_ptcollege+p_k12+p_g9to12+p_child0to4+p_adult65over+
		h_size+h_vehicles+h_zerovehicle+h_lowincome+h_highincome+h_white+
		h_anykids0to4+h_anykids5to14+h_anykids15to17+h_numkids0to4+h_numkids5to14+h_numkids15to17+
		h_1a0k+h_1awk+h_2a0k+h_2awk+
		s_multnomah+s_clackamas+s_washington+
		s_beaverton+s_portland+s_hillsboro+s_gresham+s_oregoncity+s_downtown,
	data=diary.sts.ivs, R=1000)
save(diary.seq288.DHD.ANODI, file="diary.seq288.DHD.ANODI.RData")







