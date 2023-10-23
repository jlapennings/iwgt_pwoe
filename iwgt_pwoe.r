### Bayesian modelling to calculate predictive weight of evidence for a combination of (gentox) test results
### Jeroen Pennings and Tom Aldenberg, July 2021
### version 202310
###
### Input is csv file with the (in vivo) outcome in the first column and the various input assays  in other columns
### The file should only contain -1 (negative), 1 (positive) or NA (missing)
###


### set parameters !
myfolder <- "J:/myfolder"   # or whatever
myfile <- "gtxdb309set.csv"   # or whatever
nperm <- 4000     # number of permutations
myprior <- c(0.5,0.5)     # prior neg resp pos, sum should be 1, so 0.5-0.5 means odds ratio of 1
set.seed(1729)


### load file and do some preprocessing
library(MuMIn)
setwd(myfolder)
getwd()
mydata <- read.csv(myfile, header=TRUE, row.names=1)
head(mydata)
#
dim(mydata)
mydata <- mydata[which(apply(is.na(mydata),1,sum)==0),]   # remove NAs
dim(mydata)
unique(unlist(mydata))   # should be only -1 and 1
apply(mydata ==1,2,sum, na.rm=TRUE); apply(mydata ==-1,2,sum, na.rm=TRUE)
mydata[which(mydata[,1] == -1),1] <- 0   # outcome values should be 0 or 1 for logistic regression


### prepare
ntests <- ncol(mydata)-1   # npar (for AIC and BIC) = ntests + 1 
myformula <- as.formula(paste0(colnames(mydata)[1], " ~ ", paste(colnames(mydata)[2:ncol(mydata)], collapse=" + ")))
#
myoutputA <- matrix(NA, 2^ntests, ntests + 1)
rownames(myoutputA) <- paste0("row", 1:2^ntests)
colnames(myoutputA) <- colnames(mydata) ; colnames(myoutputA)[1] <- "Intercept"
for(i in 1:ncol(myoutputA)) myoutputA[,i] <- (-1)^(trunc(((1:(2^ntests))-1)/2^(ntests+1-i)))
myoutputA[,-1] <- -1*myoutputA[,-1]   # negatives first
#
myoutputB <- matrix(NA, 2^ntests, 18)
rownames(myoutputB) <- paste0("row", 1:2^ntests)
colnames(myoutputB) <- c("abbrev","dataNeg","dataPos","dataSum", "binprNeg","binprPos", 
                          "modelNeg","modelPos","modelNegc","modelPosc",
                          "testprNeg","testprPos", "postprNeg","postprPos", "LR","WoE","WoE.LL","WoE.UL")
for(i in 1:nrow(myoutputB)) myoutputB[i,"abbrev"] <- rawToChar(as.raw(79+myoutputA[i,-1]))
myoutputB <- as.data.frame(myoutputB)
#
myoutputW <- matrix(NA, 2^ntests, nperm)
rownames(myoutputW) <- paste0("row", 1:2^ntests)


### calculate and bootstrap
myprogress <- winProgressBar(title = "progress bar", min = 0, max = nperm + 1, width = 300)
for(jj in nperm:0) {   # jj to 1 = bootstrap, 0 = full set
   setWinProgressBar(myprogress, (nperm + 1 - jj), title=paste(trunc((nperm + 1 - jj)/(nperm + 1)*100, 0), "% done"))
   if(jj == 0) mydata2 <- mydata
   if(jj >  0) mydata2 <- mydata[ceiling(runif(nrow(mydata))*nrow(mydata)),]
   myoutputB[,-1] <- NA   # reset all except abbrev
   #
   mydata2abbrev <- NULL 
   for(i in 1:nrow(mydata2)) mydata2abbrev <- c(mydata2abbrev, rawToChar(as.raw(79+mydata2[i,-1])))
   mydata2abbrev <- cbind(outcome=mydata2[,1], mydata2abbrev)
   for(i in 1:nrow(myoutputB)) {
      myoutputB[i,"dataNeg"] <- nrow(subset(mydata2abbrev, mydata2abbrev[,1] == 0 & mydata2abbrev[,2] == myoutputB[i,"abbrev"]))
      myoutputB[i,"dataPos"] <- nrow(subset(mydata2abbrev, mydata2abbrev[,1] == 1 & mydata2abbrev[,2] == myoutputB[i,"abbrev"]))
      myoutputB[i,"dataSum"] <- nrow(subset(mydata2abbrev, mydata2abbrev[,2] == myoutputB[i,"abbrev"]))
      }
   myglm <- glm(myformula, family="binomial", data = mydata2)   # NB bootstrap model too
   myoutputC <- t(summary(myglm)$coef)
   for(i in 1:nrow(myoutputB)) {
      mybetas <- sum(myoutputA[i,] * myoutputC[1,])
      myoutputB[i,"binprNeg"] <- exp(-mybetas) / (1 + exp(-mybetas))
      myoutputB[i,"binprPos"] <- 1 / (1 + exp(-mybetas))
      }
      myoutputB[,"modelNeg"] <- myoutputB[,"binprNeg"] * myoutputB[,"dataSum"]
      myoutputB[,"modelPos"] <- myoutputB[,"binprPos"] * myoutputB[,"dataSum"]
      #myoutputB[,"modelNegc"] <- myoutputB[,"modelNeg"]   # correct for zeros or low values?
      #myoutputB[,"modelPosc"] <- myoutputB[,"modelPos"]   # correct for zeros or low values?
      myoutputB[,"testprNeg"] <-  (myoutputB[,"modelNeg"] / sum(myoutputB[,"modelNeg"])) + 1E-9
      myoutputB[,"testprPos"] <-  (myoutputB[,"modelPos"] / sum(myoutputB[,"modelPos"])) + 1E-9
   for(i in 1:nrow(myoutputB)) {
      myoutputB[i,"postprNeg"] <- (myprior[1]*myoutputB[i,"testprNeg"])/
         sum((myprior[1]*myoutputB[i,"testprNeg"]) + (myprior[2]*myoutputB[i,"testprPos"]))
      myoutputB[i,"postprPos"] <- (myprior[2]*myoutputB[i,"testprPos"])/
         sum((myprior[1]*myoutputB[i,"testprNeg"]) + (myprior[2]*myoutputB[i,"testprPos"]))
      myoutputB[i,"LR"]  <- myoutputB[i,"postprPos"]/myoutputB[i,"postprNeg"]
      myoutputB[i,"WoE"] <- 10*log10(myoutputB[i,"postprPos"]/myoutputB[i,"postprNeg"])
      }
   if(jj >  0) myoutputW[,jj] <- myoutputB[,"WoE"]
   } # end jj
close(myprogress)
myoutputW <- t(apply(myoutputW, 1, sort))
myoutputB[,"WoE.LL"]  <- myoutputW[,trunc(nperm*0.025)+1]
myoutputB[,"WoE.UL"] <- myoutputW[,trunc(nperm*0.975)]


### combine output
myoutputsums <- c(rep("",ncol(myoutputA)), "", apply(myoutputB[,2:4],2,sum), "","", 
                apply(myoutputB[,7:12],2,sum),  "","","","","","")
myoutputAB <- rbind(cbind(myoutputA,myoutputB), myoutputsums, spacer = rep("",ncol(myoutputA)+ncol(myoutputB)))
myoutputAB[,1] <- ""
rownames(myoutputAB)[nrow(myoutputAB)-1] <- "sums"
myoutputD <- rep("", 4 * ncol(myoutputB)); dim(myoutputD) <- c(4, ncol(myoutputB))
myoutputD[1,3] <- "n records";    myoutputD[1,4] <- nrow(mydata)
myoutputD[2,3] <- "n parameters"; myoutputD[2,4] <- ncol(myoutputA)   # yes, include intercept
myoutputD[3,3] <- "prior ratio";  myoutputD[3,4] <- model_priorratio <- myprior[2]/myprior[1]   # here: pos/neg
myoutputD[4,3] <- "AICc";         myoutputD[4,4] <- AICc(myglm)
myoutputCD <- cbind(myoutputC,myoutputD); colnames(myoutputCD) <- colnames(myoutputAB)
myoutputABCD <- rbind(myoutputAB, myoutputCD)
# save, but leave out column modelNegc and modelPosc
outputfile <- strsplit(myfile,"")[[1]]
outputfile <- paste0(outputfile[1:(length(outputfile)-4)], collapse="")
outputfile <- paste0(outputfile, "_output.txt", sep="")
write.table(myoutputABCD[,-which(colnames(myoutputABCD) %in% c("modelNegc","modelPosc"))], 
            outputfile, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

### end



