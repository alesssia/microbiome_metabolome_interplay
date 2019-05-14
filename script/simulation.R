#Code used as part of the following manuscript:
#  "Interplay between the human gut microbiome and host metabolism"
# doi: https://doi.org/10.1101/561787

# Author: Mario Falchi
# Developer: Mario Falchi, Alessia Visconti

library(qvalue)
library(car)

minObs <- 100 #We need at least this overlap between metagenome and metabolome
minMissing <- 30 #At least these missing when evaluating the correlation
ncycles <- 1000 #Number of times we need to run the simulation

#-----------------------------
#	LOADs (corrected) metabolites
#-----------------------------

#Metabolites were corrected by age, sex and family structure
FaecalBloodCorrectAgeSexFam <- read.table("FaecalBloodCorrectAgeSexFam.RData", header=T, sep="\t")

#-----------------------------
#	LOADs ASSOCIATION RESULTs
#-----------------------------

#This can be run also on the results of the association with species by loading the correspondant 
#results files
assocF <- read.table("pathways_vs_faecal.tsv", sep="\t", header=T)
assocB <- read.table("pathways_vs_blood.tsv", sep="\t", header=T)

#Inner joint by common pathways, joining results in blood and faeces (those passing FDR 5%)
assocFB <- merge(assocB[assocB$adj_P < 0.05, ], assocF[assocF$adj_P < 0.05, ], by="Pathway", all=F)

#Unique ID for the pair of faecal and blood metabolites (co-associated metabolites)
assocFB$key <- paste(assocFB$Faecal_metabolite, assocFB$Blood_metabolite, sep="-")


#-----------------------------
#	Generates NULL sample
#-----------------------------

#Divide faecal from blood metabolites and generates all their possible pairs
null_sample <- colnames(FaecalBloodCorrectAgeSexFam)
f <- setdiff(null_sample[grepl("^F", null_sample)], "FID")
p <- null_sample[grepl("^B", null_sample)]

null_sample <- expand.grid(f, p)
colnames(null_sample) <- c("code_F", "code_B")
null_sample$code_F <- as.character(null_sample$code_F)
null_sample$code_B <- as.character(null_sample$code_B)

#Evaluates the correlation (along with its pvalue) for all the pairs
corr <- NULL
P <- NULL
for (i in 1:nrow(null_sample))
{
	tmp <- na.omit(FaecalBloodCorrectAgeSexFam[, c(null_sample[i, 1], null_sample[i, 2])])
	if (nrow(tmp) >= minObs)
	{
		ctest <- cor.test(tmp[, 1], tmp[, 2])
		corr <- c(corr, ctest$estimate)
		P <- c(P, ctest$p.value)
	}	else
	{
		corr <- c(corr, NA)
		P <- c(P, NA)
	}
}
null_sample <- cbind(null_sample, corr, P)

#Removes pairs with no enough observations
null_sample <- na.omit(null_sample)

#Qvalues
null_sample$qv <- qvalue(null_sample$P)$qvalues

#Generates unique keys and remove co-associated metabolites
null_sample$key <- paste(null_sample$code_F, null_sample$code_B, sep="-")
null_sample <- null_sample[!null_sample$key %in% assocFB$key, ]

#-----------------------------
#	LOAD PATHWAYS 
#-----------------------------

#If species is to be evaluated, use here the corresponding file
pathways <- read.table("community_pathways.tsv", sep="\t", header=T)
pathways$sex <- as.factor(pathways$sex)
pathways$age <- scale(pathways$age)


#-----------------------------
#	EVALUATING ACTUAL CORRELATIONS
#-----------------------------

sink("simulation_correlation_actual.txt")
for(i in 1:nrow(assocFB))
{
	#Selecting the bacteria and the pathway to which they associate
	code_B <- assocFB$Blood_metabolite[i]
	code_F <- assocFB$Faecal_metabolite[i]
	response <- assocFB$Pathway[i]
	
	#Selecting non null data with both metabolites
	tmp <- na.omit(FaecalBloodCorrectAgeSexFam[,c("IID", code_B, code_F)])
        
	if(nrow(tmp)>minObs)
	{
		#Merging with the pathway
		tmpT <- merge(tmp, pathways[,c("IID", response)], by="IID")
		names(tmpT)[4]<-"Y"
		
		#Al least minMissing missing data
		if(sum(is.na(tmpT$Y)) > minMissing)
		{ 
			#Global correlation, then with the pathway and without
			all  <- cor.test(tmpT[,2], tmpT[,3])
			withY <- cor.test(tmpT[!is.na(tmpT$Y),2], tmpT[!is.na(tmpT$Y),3])
			withoutY <- cor.test(tmpT[is.na(tmpT$Y),2], tmpT[is.na(tmpT$Y),3])
            
			#Levene's test
			leveneF <- leveneTest(tmpT[,2] ~ as.factor(is.na(tmpT$Y)))$`Pr(>F)`[1]
			leveneP <- leveneTest(tmpT[,3] ~ as.factor(is.na(tmpT$Y)))$`Pr(>F)`[1]
            
			#Saving the results
			cat(i, "\t", response , "\t", assocFB$code_F[i], "\t", assocFB$code_B[i], "\t", all$estimate,"\t", all$p.value, "\t", leveneF, "\t", leveneP, "\t", withY$estimate,"\t", withY$p.value, "\t",  nrow(tmp[!is.na(tmpT$Y),]), "\t", withoutY$estimate,"\t", withoutY$p.value, "\t",  nrow(tmp[is.na(tmpT$Y),]), "\n")
		}
	}	
}
sink()

#Results just generated, Mario sank them
correlation_actual <- read.table("simulation_correlation_actual.txt", as.is=T, header=F, sep="\t")
names(correlation_actual)<-c("index", "bacteria", "faecal", "plasma", "corrALL", "PAll", "LeveneF","LeveneP","CorrY", "PY", "NY", "CorrNoY", "PNoY", "NnoY")

#Stats for doing the match (percentage of incomplete and abs of global correlation)
correlation_actual$percIncomplete<-correlation_actual$NnoY/correlation_actual$NY
correlation_actual$correlation  <-abs(round(correlation_actual$corrALL,2)) 


set.seed(42)
best=0
for(cycle in 1:ncycles)
{

	#This is used to sample for similar correlations and number of observation,
	completeness <- data.frame(correlation_actual$correlation, (correlation_actual$NY+correlation_actual$NnoY), correlation_actual$NnoY, correlation_actual$index)
	names(completeness)=c("correlation","n","nincomplete", "index")
    
	sink("simulation_correlation_null.txt")
    
	#Random sampling
	null_sample <- random <- null_sample[(sample(nrow(null_sample))),]

	#Count the matched
	matched = 0
    
	i=0
	while (nrow(completeness)>0 & matched<1000 & i < nrow(random)) 
	{
		i=i+1
        
		#Selecting the sampled metabolites
		pos1<-which(names(FaecalBloodCorrectAgeSexFam)==random$code_B[i])
		pos2<-which(names(FaecalBloodCorrectAgeSexFam)==random$code_F[i])
		tmp<-FaecalBloodCorrectAgeSexFam[,c(2,pos1,pos2)]
		tmp<-tmp[complete.cases(tmp[,2], tmp[,3]),]
            
		if(dim(tmp)[1] > minObs)
		{ 
            #Stats for matching   
			all <- cor.test(tmp[,2], tmp[,3])
			correlation <- round(abs(all$estimate),2)
			size <- nrow(tmp)
			y <- 0
                
			found <- completeness[ which(completeness$n == size & completeness$correlation == correlation),]
            
			#At least one actual value with the same number of observation and correlation 
			if(nrow(found)>0) 
			{
				matched <- matched+1
				RandomN <- sample(length(found$nincomplete),1)
                
				tmp <- found$nincomplete[RandomN]
				pres = rep(0, tmp)
				Y <- sample(c(pres, rep(1,nrow(tmp)-length(pres))))
				tmp$Y <- Y
                
				#Doing tests done before
				leveneF<-leveneTest(tmp[,2] ~ as.factor((tmp$Y)))$`Pr(>F)`[1]
				leveneP<-leveneTest(tmp[,3] ~ as.factor((tmp$Y)))$`Pr(>F)`[1]
				if(leveneF>0.05 & leveneP>0.05)
				{
					#Removing from sample pool
					completeness[which(completeness$index==found$index[RandomN]),]<-NA
                    
					#Correlation as before
					withY <- cor.test(tmp[tmp$Y>0,2],tmp[tmp$Y>0,3])
					withoutY <- cor.test(tmp[tmp$Y==0,2],tmp[tmp$Y==0,3])
                    
					#Saving the results
					cat(i, "\t", "-", "\t", random$code_F[i] , "\t", random$code_B[i], "\t", all$estimate,"\t", all$p.value, "\t", withY$estimate,"\t", withY$p.value, "\t", nrow(tmp[tmp$Y>0,]) , "\t", withoutY$estimate,"\t", withoutY$p.value, "\t", nrow(tmp[tmp$Y==0,]), "\n")
				}
			}
		}
	}
	sink();sink();sink()
    

	
	#Results just generated, Mario sank them
	simulation_correlation_null <- read.table("simulation_correlation_null.txt", as.is=T, header=F, sep="\t")
	names(simulation_correlation_null) <- c("index", "bacteria", "faecal", "plasma", "corrALL", "PAll","CorrY", "PY", "NY","CorrNoY", "PNoY", "NnoY")

	#Summary stats
	used <- correlation_actual[ -which(correlation_actual$index %in% completeness[,4]),]
	actual <- sum(abs(used $CorrY)> abs(used $CorrNoY))/nrow(used)
	null <- sum(abs(simulation_correlation_null$CorrY)>abs(simulation_correlation_null$CorrNoY))/nrow(simulation_correlation_null)

	#Saving the results
	if (cycle>1)
	{
		cat(round(sum(prevUsed$index %in% used$index)/nrow(used),2),"\t",round(sum(simulation_correlation_null.tmp$index %in% simulation_correlation_null$index)/nrow(simulation_correlation_null),2),"\t")
	}
	else {
		cat("","\t","","\t")
	}
	cat(actual,"\t")
	cat(null,"\t")
	prevUsed<-used
	
	#Stats
	simulation_correlation_null.tmp <- simulation_correlation_null
	cat(wilcox.test(abs(simulation_correlation_null$corrALL), abs(used$corrALL))$p.value,"\t")
	cat(wilcox.test(abs(simulation_correlation_null$NY), abs(used$NY))$p.value,"\t")
	cat(wilcox.test(abs(simulation_correlation_null$NnoY), abs(used$NnoY))$p.value,"\t")
	cat(wilcox.test(-log10(simulation_correlation_null$PAll), -log10(used$PAll))$p.value,"\n")

	sink()
	if(null >= actual) best= best+1
}

#Final results (probability)
best/cycle

