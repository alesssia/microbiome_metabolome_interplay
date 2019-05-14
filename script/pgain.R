#Code used as part of the following manuscript:
#  "Interplay between the human gut microbiome and host metabolism"
# doi: https://doi.org/10.1101/561787

# Author: Mario Falchi
# Developer: Mario Falchi, Alessia Visconti

options(stringsAsFactors = FALSE)

library(lmerTest)
library(parallel)

#This function allows selecting the proper model accordingly to the 
#prensence of both males and females (sex is included as covariate) 
#and of family structures (family ID is included as random effect and 
#a LMM is used)
select.models <- function(tmp)
{
	if (length(unique(tmp$FID)) == nrow(tmp) )
	{
		mymodel <- lm
		if (length(unique(tmp$sex)) == 2)
		{
			mycovariates <- " age + sex"
		} else
		{
			mycovariates <- " age"
		}
	} else
	{
		mymodel <- lmer
		if (length(unique(tmp$sex)) == 2)
		{
			mycovariates <- " age + sex + (1|FID)"
		} else
		{
			mycovariates <- " age + (1|FID)"
		}
	}

	list(mymodel=mymodel, mycovariates=mycovariates)
}


minObs <- 100 #We need at least this overlap between metagenome and metabolome
todo <- 1000 #Number of observation in each null
ncycles <- 1000 #Number of times we need to run the simulation
 

#-----------------------------
#	LOADs METABOLITES
#-----------------------------

fmet <- read.csv("Faecal_metabolites.csv")
bmet <- read.csv("Blood_metabolites.csv")
met <- merge(fmet, bmet, by=c("FID", "IID"), all=F)

#-----------------------------
#	LOADs ASSOCIATION RESULTS
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
#	LOAD PATHWAYS 
#-----------------------------

#If species is to be evaluated, use here the corresponding file
pathways <- read.table("community_pathways.tsv", sep="\t", header=T)
pathways$sex <- as.factor(pathways$sex)
pathways$age <- scale(pathways$age)

#-----------------------------
#	EVALUATING THE PGAIN
#-----------------------------

#This is the results of the real model, to be compared to the results of the null 
#model generated below
pgain <- mclapply(1:nrow(assocFB), function(i, df)
{ 
	codeB <- assocFB[i, "Blood_metabolite"]
	codeF <- assocFB[i, "Faecal_metabolite"]
	mypathway <- assocFB[i, "Pathway"]

	tmp <- na.omit(met[, c("FID", "IID", codeB, codeF)])
	tmp <- na.omit(merge(tmp, pathways[, c("FID", "IID", mypathway, "age", "sex")], by=c("FID", "IID"),  all=F))
	#Evaluating the ratio, that since the data are log transformed, can be done using simply
	#log(A/B) == log(A)-log(B)
	tmp$ratio <- tmp[, 4]-tmp[, 3]
	
	#Cheching if there are enough observation
	if (nrow(tmp) < minObs)
	{
		return (c(i, mypathway, codeF, codeB, nrow(tmp), rep(NA, 8)))
	}

	#Selecting the correct model
	a <- select.models(tmp)
	mymodel <- a$mymodel
	mycovariates <- a$mycovariates
	
	#Taking the residual of the corrected metabolites and evaluating their correlation
	#This will be used in the matching, since the metabolites correlation influence the
	#pgain
	rF <- residuals(mymodel(paste0(codeF, " ~ ", mycovariates), data=tmp))
	rB <- residuals(mymodel(paste0(codeB, " ~ ", mycovariates), data=tmp))
	ctest <- cor.test(rF, rB)

	#assiciating the metabolites and their ratio
	Fm <- summary(mymodel(paste0(codeF, " ~ ", mypathway, "+", mycovariates), data=tmp))$coefficient
	Bm <- summary(mymodel(paste0(codeB, " ~ ", mypathway, "+", mycovariates), data=tmp))$coefficient
	Gm <- summary(mymodel(paste0("ratio ~ ", mypathway, "+", mycovariates), data=tmp))$coefficient 
          
	c(i, mypathway, codeF, codeB, nrow(tmp), Fm[2, c(1,5)], Bm[2, c(1,5)], Gm[2, c(1,5)], ctest$estimate, ctest$p.value)
}, mc.cores=4)
pgain <- do.call(rbind, pgain)
colnames(pgain)<-c("index", "Response", "Fmet", "Pmet", "Num", "FaecalBeta", "FaecalP", "BloodB","BloodP", "GainB", "GainP", "corr", "corrP")


pgain <- as.data.frame(pgain)
pgain$Num <- as.numeric(pgain$Num)

pgain$FaecalP <- as.numeric(pgain$FaecalP)
pgain$BloodP <- as.numeric(pgain$BloodP)
pgain$GainP <- as.numeric(pgain$GainP)
pgain$corr <- as.numeric(pgain$corr)

#Evaluating the pgain as min(P_faecal, P_blood)/P_ratio
pgain$pgain <- apply(pgain[,c("FaecalP","BloodP","GainP")], 1, function(v) {min(v[1],v[2])/v[3]})
#Removing those with no enough observation
pgain <- na.omit(pgain)

#-----------------------------
#	GENERATING THE RANDOM SAMPLES AND COMPARING WITH ACTUAL RESULTS
#-----------------------------

#Set of metabolites I can sample from
assocFFDR <- assocF[assocF$adj_P < 0.05, ]
assocBFDR <- assocB[assocB$adj_P < 0.05, ]

#This is used to sample for similar correlations and number of observation,
#both influencing the pgain
completeness <- data.frame(round(pgain$corr,2), pgain$Num, pgain$index)
names(completeness) <- c("corr", "n", "index")

set.seed(42)	
for (cycle in 1:ncycles)
{	
	#Samples two metabolites (actually two results of the association tests)
	assocFFDR<-assocFFDR[sample(1:nrow(assocFFDR)),]
	assocBFDR<-assocBFDR[sample(1:nrow(assocBFDR)),]
    	 
	sink(paste0("pgainNull_pathways_cycle_", cycle, ".txt"))
    
	#Counts how many have been done
	done = 0
	
	#TODO: To be optimised using a function, and discarding the FOR
	for (i in 1:nrow(assocFFDR))
	{
		#Checks if I have still someting to sample, or sample to generate
		if (sum(!is.na(completeness)) == 0 | done >= todo)  { break }
		
		for (j in 1:nrow(assocBFDR))
		{
			#Checks if I have still someting to sample, or sample to generate
			if (sum(!is.na(completeness)) == 0 | done >= todo) { break }
            
			#Selects the metabolites and the two pathways they are associated
			codeF <- assocFFDR$Faecal_metabolite[i]
			codeB <- assocBFDR$Blood_metabolite[j]
			mypathwayF <- assocFFDR$Pathway[i]
			mypathwayB <- assocBFDR$Pathway[j]
			
			key<-paste(codeF, codeB, sep="-")
			
			#They should not be a real pair!
			if (!key %in% assocFB$key)
			{	
				#Selects the metabolites and evaluates their ratio
				tmp1 <- tmp2 <- na.omit(met[,c("FID","IID", codeF, codeB)])
				tmp1$ratio <- tmp2$ratio <- tmp1[, 4]-tmp1[, 3]
  
  			    #Merges with the two pathways generating two datasets (one is the pathways
				#associate withe the randomly selected faecal metabolites, the other with 
				#the blood metabolites). Remember that the metabolites are actually selected
				#from the results file
				tmp1 <- na.omit(merge(tmp1, pathways[, c("FID", "IID", mypathwayF, "age", "sex")], by=c("FID", "IID"), all=F))
				tmp2 <- na.omit(merge(tmp2, pathways[, c("FID", "IID", mypathwayB, "age", "sex")], by=c("FID", "IID"), all=F))
				
				#As done for the real data, count the number of overlapping omics, 
				#selectes the model and evaluates their correlation
				dim1 <- nrow(tmp1)
				dim2 <- nrow(tmp2)
				
				a1 <- select.models(tmp1)
				mymodel1 <- a1$mymodel
				mycovariates1 <- a1$mycovariates

				a2 <- select.models(tmp2)
				mymodel2 <- a2$mymodel
				mycovariates2 <- a2$mycovariates
  
				if (dim1 > minObs)
				{
					rF <- residuals(mymodel1(paste0(codeF, " ~ ", mycovariates1), data=tmp1))
					rP <- residuals(mymodel1(paste0(codeB, " ~ ", mycovariates1), data=tmp1))
					ctest1 <- cor.test(rF, rP)
					ctest1 <- c(round(ctest1$estimate, 2), ctest1$p.value)
				} else 
				{
					ctest1 <- c(NA, 1)
				}
    
				if (nrow(tmp2) > minObs)
				{
					rF <- residuals(mymodel2(paste0(codeF, " ~ ", mycovariates2), data=tmp2))
					rP <- residuals(mymodel2(paste0(codeB, " ~ ", mycovariates2), data=tmp2))
					ctest2 <- cor.test(rF, rP)
					ctest2 <- c(round(ctest2$estimate, 2), ctest2$p.value)
				} else 
				{
					ctest2 <- c(NA, 1)
				}
                
				#Checks if one of the real value has the same number of observations of the selected
				#pair of metabolites and the same correlation
				if(nrow(completeness[ which(completeness$n==dim1 & completeness$corr==ctest1[1]),])>0) 
				{
                    #If yes, I take it out from the list
					found<-completeness[ which(completeness$n==dim1 & completeness$corr==ctest1[1]), ]
					RandomN<-sample(nrow(found),1)
					completeness[ which(completeness$index==found$index[RandomN]),]<-NA
                    
					#Do the tests
					Fm <- summary(mymodel1(paste0(codeF, " ~ ", mypathwayF, " + ", mycovariates1), data=tmp1))$coefficient
					Pm <- summary(mymodel1(paste0(codeB, " ~ ", mypathwayF, " + ", mycovariates1), data=tmp1))$coefficient
					Gm <- summary(mymodel1(paste0("ratio ~ ", mypathwayF, " + ", mycovariates1), data=tmp1))$coefficient 
  
  				    #Store the data
					r1 <- c(responseF, codeF, codeB, nrow(tmp1), Fm[2, c(1,5)], Pm[2, c(1,5)], Gm[2, c(1,5)], ctest1)
					cat(r1, "\n", sep="\t")
					
					#Count a test done                    
					done <- done + 1
				}
	        
				#If the number of test is enough I don't need to test the second pathway
				if (sum(!is.na(completeness))==0 | done >= todo) { break }
            
				#As before
				if(nrow(completeness[ which(completeness$n==dim2 & completeness$corr==ctest2[1]),])>0) {
                    
					found<-completeness[ which(completeness$n==dim2 & completeness$corr==ctest2[1]),]
					RandomN<-sample(nrow(found),1)
					completeness[ which(completeness$index==found$index[RandomN]),]<-NA
                    
					Fm <- summary(mymodel2(paste0(codeF, " ~ ", mypathwayB, " + ", mycovariates2), data=tmp2))$coefficient
					Pm <- summary(mymodel2(paste0(codeB, " ~ ", mypathwayB, " + ", mycovariates2), data=tmp2))$coefficient
					Gm <- summary(mymodel2(paste0("ratio ~ ", mypathwayB, " + ", mycovariates2), data=tmp2))$coefficient 
  
					r2 <- c(responseB, codeF, codeB, nrow(tmp2), Fm[2, c(1,5)], Pm[2, c(1,5)], Gm[2, c(1,5)], ctest2)
                    
					cat(r2, "\n", sep="\t")
					
					done <-done + 1                   
				}
			}
		}
	}
    
	sink(); sink(); sink()
	
	#Those matched 
	used <- pgain[-which(pgain$index %in% completeness$index),]
    write.table(used, paste0("pgainNull_pathways_used_cycle_", cycle, ".txt"), col.names=F, row.names=F, quote=F, sep="\t")
	
	#Results just generated, Mario sank them
	pgainNull <- read.table(paste0("pgainNull_pathways_cycle_", cycle, ".txt"), as.is=T, header=F, sep="\t")[, 1:12]
	names(pgainNull)<-c("Response", "faecal", "blood","Num","FB", "FP", "PB","PP", "GB", "GP", "corr", "Pcorr")
	
	#Generates pgain
	pgainNull$pgain<-apply(pgainNull[,c("FP","FP","GP")], 1, function(x){min(x[1],x[2])/x[3]})
	pgainNull$key<-paste(pgainNull$faeca,"-", pgainNull$blood, sep="")
    
	#For each cycle, compare real results with those randomly sampled
	if(cycle==1)
	{
			r <- c("NA", "NA", index, cycle, wilcox.test(used$corr, pgainNull$corr)$p.value, wilcox.test(used$Num, pgainNull$Num)$p.value, wilcox.test(used$pgain, pgainNull$pgain, alternative = "greater")$p.value)
		} else 
		{
			r <- c(sum(used$index %in% usedPrev$index), sum(pgainNull$key %in% pgainNullPrev$key), index, cycle, wilcox.test(used$corr, pgainNull$corr)$p.value, wilcox.test(used$Num, pgainNull$Num)$p.value, wilcox.test(used$pgain, pgainNull$pgain, alternative = "greater")$p.value)
		}
	
		#Saves the cycle results
		print(r)
		write.table(paste(r, collapse="\t"), file="pgainNull_results_pathways.txt", col.names=F, row.names=F, append=T, quote=F)
		
		#Reset for next cycle
		pgainNullPrev <- pgainNull
		usedPrev <- used
}

