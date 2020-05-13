# fold of SSEIH data

library(data.table)
library(plotly)
library(fields)
library(Hmisc)
library(purrr)

prepath <- "."
filenames <- list.files(paste0(prepath,"/runs"))
pngmode <- FALSE
pdfmode <- TRUE

BSsel <- .00 # throw out bootstraps that does not fit the initial phase

ResDK <- as.numeric(as.Date("2020-03-11"))-as.numeric(as.Date("2020-01-01"))
Res2DK <-  as.numeric(as.Date("2020-04-15"))-as.numeric(as.Date("2020-01-01"))
Res3DK <-  as.numeric(as.Date("2020-05-18"))-as.numeric(as.Date("2020-01-01"))

endTimes <- as.numeric(as.Date("2020-08-01"))-as.numeric(as.Date("2020-01-01"))
times <- seq(ResDK,endTimes,1)

xt <- c(ResDK,as.numeric(as.Date(format(as.Date((4:24)*30,origin="2020-01-01"),
                                        "%Y-%m-01")))-
          as.numeric(as.Date("2020-01-01")))
timelabels <- format(as.Date(xt,origin="2020-01-01"),"%d %b")
timelabels <-gsub("Jul","jul",gsub("Jun","jun",gsub("Mar","mar",
                                       gsub("May","maj",gsub("Apr",replacement = "apr",timelabels)))))



## The data that is loaded as 'mdt' cannot be shared at this point.
## We have asked if can be shared and are awaiting the answer.
prepath <- "."
mdt <- fread(file = paste0(prepath,"/data_for_kaare2.csv"))
mdt[Event=="Hospitalised",Event:="HospitalisedCase"]
mdt$Date <- as.Date(mdt$Date)
all <- mdt[,.(New=sum(New),Cumulative=sum(Cumulative),Instant=sum(Instant),RegionName="all"),
           by=.(AgeGroup,Event,Date)]
amdt <- rbind(mdt,all)
rm(list=c("all","mdt"))

regions=c("Hovedstaden","Zealand","Syddanmark",
          "Midtjylland","Nordjylland","all")
regionsTitle=c("Region Hovedstaden","Region Sjælland","Region Syddanmark",
          "Region Midtjylland","Region Nordjylland","Hele Danmark")


mainTitles <- c("Uden yderligere åbninger",
"Grundblok (G)",
"G + skole 6.-10 kl.",
"G + skole 6.-10 kl. + efterskole",
"G + skole 6.-10 kl. + restauration (ude)",
"G + skole 6.-10 kl. + restauration (inde- og ude)",
"G + skole 6.-10 kl. + private",
"G + skole 6.-10 kl. + efterskole + restauration (ude) + private",
"G + skole 6.-10 kl. + efterskole + restauration (inde og ude) + private",
"G + skole 9.-10 kl.",
"G + skole 9.-10 kl. + efterskole + restauration (ude) + private",
"G + skole 9.-10 kl. + efterskole + restauration (inde og ude) + private",
"Grundblok (G) + halv fysisk afstand",
"Grundblok (G) + normal fysisk afstand",
"G + skole 6.-10 kl. + efterskole + rest.(i & u) + private + halv fysisk afstand",
"G + skole 6.-10 kl. + efterskole + rest.(i & u) + private + normal fysisk afstand",
"Uden yderligere åbninger + halv fysisk afstand")


testmode <- !TRUE

colorTable<- c("lightgray","darkgray","lightgray")

prepath <- "."
allProbs.old.early <- read.csv(paste0(prepath,"/currentStateProbs/stateProbs_old_early.csv"))
allProbs.young.early <- read.csv(paste0(prepath,"/currentStateProbs/stateProbs_young_early.csv"))
allProbs.old.weighted <- read.csv(paste0(prepath,"/currentStateProbs/stateProbs_old_weighted.csv"))
allProbs.young.weighted <- read.csv(paste0(prepath,"/currentStateProbs/stateProbs_young_weighted.csv"))


stateAll.fun <- function(t.val, probsData, N=1)
{
  state1.fun <- approxfun(x = probsData$time, y = probsData$pstate1, yright = 0)
  state2.fun <- approxfun(x = probsData$time, y = probsData$pstate2, yright = 0)
  state3.fun <- approxfun(x = probsData$time, y = probsData$pstate3, rule =2)
  state4.fun <- approxfun(x = probsData$time, y = probsData$pstate4, rule =2)
  temp <- N*cbind(inHosp=state1.fun(t.val),inICU=state2.fun(t.val),outsideHosp=state3.fun(t.val),dead=state4.fun(t.val))
  row.names(temp) <- NULL
  return(temp)
}

HosFoldTheis2 <- function(cumIndHosU60,cumIndHosO60,initIndHosU60,initIndHosO60)
{
  newHosU60 <- c(initIndHosU60,diff(cumIndHosU60 ))
  newHosO60 <- c(initIndHosO60,diff(cumIndHosO60 ))
  
  nHorPred <- 100
  
  allHosU60 <- matrix(0,nrow=length(newHosU60)+nHorPred,ncol=4)
  allHosO60 <- matrix(0,nrow=length(newHosO60)+nHorPred,ncol=4)
  
  oldProb <- stateAll.fun(0:nHorPred, probsData = allProbs.old.early, N=1)
  youngProb <- stateAll.fun(0:nHorPred, probsData = allProbs.young.early, N=1)
  
  
  for (i in 1:21) 
  {
    allHosU60[i:(i+nHorPred),] <- allHosU60[i:(i+nHorPred),] + youngProb*newHosU60[i]
    allHosO60[i:(i+nHorPred),] <- allHosO60[i:(i+nHorPred),] + oldProb*newHosO60[i]
  }
  
  oldProb <- stateAll.fun(0:nHorPred, probsData = allProbs.old.weighted, N=1)
  youngProb <- stateAll.fun(0:nHorPred, probsData = allProbs.young.weighted, N=1)
  
  
  for (i in 22:length(newHosO60)) 
  {
    allHosU60[i:(i+nHorPred),] <- allHosU60[i:(i+nHorPred),] + youngProb*newHosU60[i]
    allHosO60[i:(i+nHorPred),] <- allHosO60[i:(i+nHorPred),] + oldProb*newHosO60[i]
  }
  
  allHosU60 <- allHosU60[1:length(newHosU60),]
  allHosO60 <- allHosO60[1:length(newHosO60),]
  
  return( rowSums(allHosU60[,1:2]) + rowSums(allHosO60[,1:2]) )
}


HosFoldTheis2ICU <- function(cumIndHosU60,cumIndHosO60,initIndHosU60,initIndHosO60)
{
  newHosU60 <- c(initIndHosU60,diff(cumIndHosU60 ))
  newHosO60 <- c(initIndHosO60,diff(cumIndHosO60 ))
  
  nHorPred <- 100
  
  allHosU60 <- matrix(0,nrow=length(newHosU60)+nHorPred,ncol=4)
  allHosO60 <- matrix(0,nrow=length(newHosO60)+nHorPred,ncol=4)
  
  oldProb <- stateAll.fun(0:nHorPred, probsData = allProbs.old.early, N=1)
  youngProb <- stateAll.fun(0:nHorPred, probsData = allProbs.young.early, N=1)
  
  
  for (i in 1:21) 
  {
    allHosU60[i:(i+nHorPred),] <- allHosU60[i:(i+nHorPred),] + youngProb*newHosU60[i]
    allHosO60[i:(i+nHorPred),] <- allHosO60[i:(i+nHorPred),] + oldProb*newHosO60[i]
  }
  
  oldProb <- stateAll.fun(0:nHorPred, probsData = allProbs.old.weighted, N=1)
  youngProb <- stateAll.fun(0:nHorPred, probsData = allProbs.young.weighted, N=1)
  
  
  for (i in 22:length(newHosO60)) 
  {
    allHosU60[i:(i+nHorPred),] <- allHosU60[i:(i+nHorPred),] + youngProb*newHosU60[i]
    allHosO60[i:(i+nHorPred),] <- allHosO60[i:(i+nHorPred),] + oldProb*newHosO60[i]
  }
  
  allHosU60 <- allHosU60[1:length(newHosU60),]
  allHosO60 <- allHosO60[1:length(newHosO60),]

  return( allHosU60[,2] + allHosO60[,2] )
}



# for kk loop

prepath <- "."

feRetSamlet <- list()
for(kk in filenames){
  load(paste0(prepath,"/runs/",kk))
  for (i in 1:length(idScenarier)){
    feRetSamlet[[idScenarier[i]]] <- feRet[[i]]
  }
}
rm(feRet)

library(doParallel)

numCores<- detectCores()
registerDoParallel(cores=numCores)


tic <- Sys.time()
listInfHos<-foreach(id =1:17, .packages = "data.table") %dopar%{
  ldt <- feRetSamlet[[id]]$ldt
  pdt <- feRetSamlet[[id]]$pdt
  
  scenario <- pdt$scenario[id]
  
  InfHos <- array(NA,dim=c(uniqueN(ldt$Date),
                           uniqueN(ldt$rep)-1,
                           uniqueN(ldt$RegionName)))
  
  for (j in 1:(uniqueN(ldt$rep)-1))
  {
    for (i in 1:uniqueN(ldt$RegionName))
    {
      InfHos[,j,i] <- HosFoldTheis2(ldt[RegionName==regions[i] & AgeGroup=="0-59" & rep==j &
                                          Event=="HospitalisedCase",Cumulative],
                                    ldt[RegionName==regions[i] & AgeGroup=="60+" & rep==j &
                                          Event=="HospitalisedCase",Cumulative],
                                    amdt[RegionName==regions[i] & AgeGroup=="0-59" & Date>="2020-03-11" &
                                           Event=="HospitalisedCase",Instant][1],
                                    amdt[RegionName==regions[i] & AgeGroup=="60+" & Date>="2020-03-11" &
                                           Event=="HospitalisedCase",Instant][1])
      
    }
  }
  InfHos
}
Sys.time()-tic
#save(listInfHos,file="listInfHos.RData")
#listInfHos[[id]] <- InfHos

##ICU
tic <- Sys.time()
listInfICU <- foreach(id =1:17, .packages = "data.table") %dopar%{
  ldt <- feRetSamlet[[id]]$ldt
  pdt <- feRetSamlet[[id]]$pdt
  
  scenario <- pdt$scenario[id]
  
  InfICU <- array(NA,dim=c(uniqueN(ldt$Date),
                           uniqueN(ldt$rep)-1,
                           uniqueN(ldt$RegionName)))
  
  for (j in 1:(uniqueN(ldt$rep)-1))
  {
    for (i in 1:uniqueN(ldt$RegionName))
    {
      
      InfICU[,j,i] <- HosFoldTheis2ICU(ldt[RegionName==regions[i] & AgeGroup=="0-59" & rep==j &
                                             Event=="HospitalisedCase",Cumulative],
                                       ldt[RegionName==regions[i] & AgeGroup=="60+" & rep==j &
                                             Event=="HospitalisedCase",Cumulative],
                                       amdt[RegionName==regions[i] & AgeGroup=="0-59" & Date>="2020-03-11" &
                                              Event=="HospitalisedCase",Instant][1],
                                       amdt[RegionName==regions[i] & AgeGroup=="60+" & Date>="2020-03-11" &
                                              Event=="HospitalisedCase",Instant][1])
      
    }
  }
  InfICU
}
Sys.time()-tic
#listInfICU[[id]] <- InfICU

#save(listInfICU,file="listInfICU.RData")

## ICU done HERE

## Plotting ####

# Pick a scenario and run the inner part of loops
sce <- 3
InfHos <- listInfHos[[3]]

## Order of scenarios
mapSce <- c(1:3,10,4:8,11,9,12:16,17)
sceId <- 17

png(paste0("fase1_indlagte.png"), width = 1600, height=2100, pointsize = 36)
par(mar=c(2.,3.3,4.3,1)+.1, mgp=c(2.1,0.7,0))

par(mfcol=c(3,2))


if(pdfmode) {
  pdf(file = "B4_Indlagte_alle.pdf", width=6, height=8)
    par(mfrow=c(3,2))
  par(mar=c(2.,3.3,4.3,1)+.1, mgp=c(2.1,0.7,0))
}
for (sceId in 1:16){
  sce <- mapSce[sceId]
  InfHos <- listInfHos[[sce]]
for (j in 1:6)
{
solBS <- InfHos[,,j]
weights <- rep(1,500)
weights[weights < quantile(weights,BSsel)] <- 0

BSsolQuan <-  apply(solBS, 1, function(x) 
  wtd.quantile(x, probs = c(.05,.25,.75,.95), weights = weights ) )

# Observational noise
PredMed <-  apply(solBS, 1, function(x) 
  wtd.quantile(x, probs = c(.5), weights = weights ) )
ObsNoise <- apply(t(PredMed),2,function(x) qpois(c(.05,.25,.75,.95),x)-x)

BSsolQuan <-  BSsolQuan + ObsNoise
BSsolQuan[BSsolQuan<0] <- 0


if (pngmode) {
  png(paste0("./Figures/",substr(filenames[kk],1,nchar(filenames[kk])-6),"_hospital.png"))
par(mfrow=c(1,1))
par(mar=c(5,2,4,2)+.1,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))
}

ymax <- ifelse(j==6,2000,700)
plot(NULL, xlim=c(times[1],times[115]), ylim=c(0,ymax), ylab="Antal indlagte på hospital", # var dag 130
     xlab=" ",xaxt="n", frame.plot = F)
axis(1,at=xt,labels=timelabels)

itimes <- (as.numeric(as.Date("2020-05-02"))-as.numeric(as.Date("2020-03-11"))):length(times)
ntimes <- times[itimes]
for (i in 4:2) { 
  polygon(c(ntimes, times[length(times)], rev(ntimes)), c(BSsolQuan[i,itimes], 0, 
                                                          rev(BSsolQuan[i-1,itimes])),col = adjustcolor(colorTable[5-i],alpha.f=0.75), border = NA)
}

grid()
legend("top",c( "SI:","5 - 95%","25 - 75%"),
       col=c(NA,colorTable[1:2] ),
       lty = c(NA,1,1),lwd=5) 

IndlagteHos <- amdt[RegionName==regions[j]  & Date>="2020-03-11" &
       Event=="HospitalisedCase",sum(Instant),by=Date]$V1

points(times[1:length(IndlagteHos)],IndlagteHos,pch=19,col=4)
title(regionsTitle[j])  #LAEC
if (pngmode) dev.off()
}
mtext(paste0("B4.",sceId-1,": ",mainTitles[sce]), side=3, line=-1.5, outer=TRUE, cex=0.9)
}


if(pdfmode) dev.off()



################################### #
## Plot ICU ####
png(paste0("fase1_ICU.png"), width = 1600, height=2100, pointsize = 36)
par(mar=c(2.,3.3,4.3,1)+.1, mgp=c(2.1,0.7,0))

par(mfcol=c(3,2))

if(pdfmode) {
  pdf(file = "B5_ICU_alle.pdf", width=6, height=8)
  par(mfrow=c(3,2))
  par(mar=c(2.,3.3,4.3,1)+.1, mgp=c(2.1,0.7,0))#,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))
}
for (sceId in 1:16){
  sce <- mapSce[sceId]
  InfICU <- listInfICU[[sce]]
  for (j in 1:6)
  {
    solBS <- InfICU[,,j]
    weights <- rep(1,500)
    weights[weights < quantile(weights,BSsel)] <- 0
    
    BSsolQuan <-  apply(solBS, 1, function(x) 
      wtd.quantile(x, probs = c(.05,.25,.75,.95), weights = weights ) )

    # Observational noise
    PredMed <-  apply(solBS, 1, function(x) 
      wtd.quantile(x, probs = c(.5), weights = weights ) )
    ObsNoise <- apply(t(PredMed),2,function(x) qpois(c(.05,.25,.75,.95),x)-x)
    
    BSsolQuan <-  BSsolQuan + ObsNoise
    BSsolQuan[BSsolQuan<0] <- 0
    
    
    if (pngmode) {
      png(paste0("./Figures/",substr(filenames[kk],1,nchar(filenames[kk])-6),"_hospital.png"))
      par(mfrow=c(1,1))
      par(mar=c(5,2,4,2)+.1,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))
    }
    
    ymax <- ifelse(j==6,400,400)
    plot(NULL, xlim=c(times[1],times[115]), ylim=c(0,ymax), ylab="Antal indlagte på intensiv", # var dag 130
         xlab=" ",xaxt="n", frame.plot = F)
    axis(1,at=xt,labels=timelabels)
    
    itimes <- (as.numeric(as.Date("2020-05-02"))-as.numeric(as.Date("2020-03-11"))):length(times)
    ntimes <- times[itimes]
    for (i in 4:2) { 
      polygon(c(ntimes, times[length(times)], rev(ntimes)), c(BSsolQuan[i,itimes], 0, 
                                                              rev(BSsolQuan[i-1,itimes])),col = adjustcolor(colorTable[5-i],alpha.f=0.75), border = NA)
    }
    
    grid()
    legend("top",c( "SI:","5 - 95%","25 - 75%"),
           col=c(NA,colorTable[1:2] ),
           lty = c(NA,1,1),lwd=5) 
    
    IndlagteHos <- amdt[RegionName==regions[j]  & Date>="2020-03-11" &
                          Event=="ICU",sum(Instant),by=Date]$V1
    
    points(times[1:length(IndlagteHos)],IndlagteHos,pch=19,col=4)
    title(regionsTitle[j])
    if (pngmode) dev.off()
  }
  mtext(paste0("B5.",sceId-1,": ",mainTitles[sce]), side=3, line=-1.5, outer=TRUE, cex=0.9)
}

if(pdfmode) dev.off()


#### png ####
## Selected plots for the report

forFig <- c(2,13,14, 9,15,16)
mainTitlesSub <- c("Grundblok (G)",
                "Grundblok (G) + \n halv fysisk afstand",
                "Grundblok (G) + \n normal fysisk afstand",
                "G + skole 6.-10 kl. + efterskole + restauration (i & u) + \n private",
                "G + skole 6.-10 kl. + efterskole + rest.(i & u) +\n private + halv fysisk afstand",
                "G + skole 6.-10 kl. + efterskole + rest.(i & u) +\n private + normal fysisk afstand")

if (pngmode) {
  png(paste0("./Figures/",substr(filenames[kk],1,nchar(filenames[kk])-6),"_hospital.png"))
  par(mfrow=c(1,1))
  par(mar=c(5,2,4,2)+.1,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))
}


## følsomhed indlagte ####

png(paste0("følsomhed_indlagte.png"), width = 1600, height=2100, pointsize = 36)
par(mar=c(2.,3.3,4.3,1)+.1, mgp=c(2.1,0.7,0))#,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))

par(mfcol=c(3,2))
for (sceId in 1:length(forFig)){
  sce <- forFig[sceId]
  InfHos <- listInfHos[[sce]]
  j<- 6
  
    solBS <- InfHos[,,j]
    weights <- rep(1,500)
    weights[weights < quantile(weights,BSsel)] <- 0
    
    BSsolQuan <-  apply(solBS, 1, function(x) 
      wtd.quantile(x, probs = c(.05,.25,.75,.95), weights = weights ) )

    # Observational noise
    PredMed <-  apply(solBS, 1, function(x) 
      wtd.quantile(x, probs = c(.5), weights = weights ) )
    ObsNoise <- apply(t(PredMed),2,function(x) qpois(c(.05,.25,.75,.95),x)-x)
    
    BSsolQuan <-  BSsolQuan + ObsNoise
    BSsolQuan[BSsolQuan<0] <- 0

    ymax <- ifelse(j==6,2000,700)
    plot(NULL, xlim=c(times[1],times[115]), ylim=c(0,ymax), ylab="Antal indlagte på hospital", # var dag 130
         xlab=" ",xaxt="n", frame.plot = F)
    axis(1,at=xt,labels=timelabels)
    
    itimes <- (as.numeric(as.Date("2020-05-02"))-as.numeric(as.Date("2020-03-11"))):length(times)
    ntimes <- times[itimes]

    for (i in 4:2) { 
      polygon(c(ntimes, times[length(times)], rev(ntimes)), c(BSsolQuan[i,itimes], 0, 
                                                              rev(BSsolQuan[i-1,itimes])),col = adjustcolor(colorTable[5-i],alpha.f=0.75), border = NA)
    }
    
    grid()
    legend("top",c( "SI:","5 - 95%","25 - 75%"),
           col=c(NA,colorTable[1:2] ),
           lty = c(NA,1,1),lwd=5) 
    
    IndlagteHos <- amdt[RegionName==regions[j]  & Date>="2020-03-11" &
                          Event=="HospitalisedCase",sum(Instant),by=Date]$V1
    
    points(times[1:length(IndlagteHos)],IndlagteHos,pch=19,col=4)
    title(mainTitlesSub[sceId])
    if (pngmode) dev.off()
} 
dev.off()



## følsomhed intensiv ####

png(paste0("følsomhed_intensiv.png"), width = 1600, height=2100, pointsize = 36)
par(mar=c(2.,3.3,4.3,1)+.1, mgp=c(2.1,0.7,0))#,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))

par(mfcol=c(3,2))
for (sceId in 1:length(forFig)){
  sce <- forFig[sceId]
  
  InfICU <- listInfICU[[sce]]
  j<- 6
  
  solBS <- InfICU[,,j]
  weights <- rep(1,500)
  weights[weights < quantile(weights,BSsel)] <- 0
  
  BSsolQuan <-  apply(solBS, 1, function(x) 
    wtd.quantile(x, probs = c(.05,.25,.75,.95), weights = weights ) )
  
  # Observational noise
  PredMed <-  apply(solBS, 1, function(x) 
    wtd.quantile(x, probs = c(.5), weights = weights ) )
  ObsNoise <- apply(t(PredMed),2,function(x) qpois(c(.05,.25,.75,.95),x)-x)
  
  BSsolQuan <-  BSsolQuan + ObsNoise
  BSsolQuan[BSsolQuan<0] <- 0
  
  ymax <- ifelse(j==6,400,400)
  plot(NULL, xlim=c(times[1],times[115]), ylim=c(0,ymax), ylab="Antal indlagte på intensiv", # var dag 130
       xlab=" ",xaxt="n", frame.plot = F)
  axis(1,at=xt,labels=timelabels)
  
  itimes <- (as.numeric(as.Date("2020-05-02"))-as.numeric(as.Date("2020-03-11"))):length(times)
  ntimes <- times[itimes]
  
  for (i in 4:2) { 
    polygon(c(ntimes, times[length(times)], rev(ntimes)), c(BSsolQuan[i,itimes], 0, 
                                                            rev(BSsolQuan[i-1,itimes])),col = adjustcolor(colorTable[5-i],alpha.f=0.75), border = NA)
  }
  
  grid()
  legend("top",c( "SI:","5 - 95%","25 - 75%"),
         col=c(NA,colorTable[1:2] ),
         lty = c(NA,1,1),lwd=5) 
  
  IndlagteHos <- amdt[RegionName==regions[j]  & Date>="2020-03-11" &
                        Event=="ICU",sum(Instant),by=Date]$V1
  
  points(times[1:length(IndlagteHos)],IndlagteHos,pch=19,col=4)
  title(mainTitlesSub[sceId])
  if (pngmode) dev.off()
} 
dev.off()




#### png fase1 fortsat  ####
## sce <- 1 eller 17 
sce <- 17
png(paste0("fase1_fortsat_sce",sce,".png"), width = 1600, height=700, pointsize = 36)
par(mar=c(2.,3.3,1.3,1)+.1, mgp=c(2.1,0.7,0))#,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))

par(mfcol=c(1,2))

  InfHos <- listInfHos[[sce]]
  j<- 6
  
  solBS <- InfHos[,,j]
  weights <- rep(1,500)
  weights[weights < quantile(weights,BSsel)] <- 0
  
  BSsolQuan <-  apply(solBS, 1, function(x) 
    wtd.quantile(x, probs = c(.05,.25,.75,.95), weights = weights ) )

  # Observational noise
  PredMed <-  apply(solBS, 1, function(x) 
    wtd.quantile(x, probs = c(.5), weights = weights ) )
  ObsNoise <- apply(t(PredMed),2,function(x) qpois(c(.05,.25,.75,.95),x)-x)
  
  BSsolQuan <-  BSsolQuan + ObsNoise
  BSsolQuan[BSsolQuan<0] <- 0
  
  ymax <- 700#ifelse(j==6,2000,700)
  plot(NULL, xlim=c(times[1],times[115]), ylim=c(0,ymax), ylab="Antal indlagte på hospital", # var dag 130
       xlab=" ",xaxt="n", frame.plot = F)
  axis(1,at=xt,labels=timelabels)
  
  itimes <- (as.numeric(as.Date("2020-05-02"))-as.numeric(as.Date("2020-03-11"))):length(times)
  ntimes <- times[itimes]
  
  for (i in 4:2) { 
    polygon(c(ntimes, times[length(times)], rev(ntimes)), c(BSsolQuan[i,itimes], 0, 
                                                            rev(BSsolQuan[i-1,itimes])),col = adjustcolor(colorTable[5-i],alpha.f=0.75), border = NA)
  }
  
  grid()
  legend("topright",c( "SI:","5 - 95%","25 - 75%"),
         col=c(NA,colorTable[1:2] ),
         lty = c(NA,1,1),lwd=5) 
  
  IndlagteHos <- amdt[RegionName==regions[j]  & Date>="2020-03-11" &
                        Event=="HospitalisedCase",sum(Instant),by=Date]$V1
  
  points(times[1:length(IndlagteHos)],IndlagteHos,pch=19,col=4)

  ## ICU part
  
  InfICU <- listInfICU[[sce]]
  j<- 6
    solBS <- InfICU[,,j]
    weights <- rep(1,500)
    weights[weights < quantile(weights,BSsel)] <- 0
    
    BSsolQuan <-  apply(solBS, 1, function(x) 
      wtd.quantile(x, probs = c(.05,.25,.75,.95), weights = weights ) )
    
    # Observational noise
    PredMed <-  apply(solBS, 1, function(x) 
      wtd.quantile(x, probs = c(.5), weights = weights ) )
    ObsNoise <- apply(t(PredMed),2,function(x) qpois(c(.05,.25,.75,.95),x)-x)
    
    BSsolQuan <-  BSsolQuan + ObsNoise
    BSsolQuan[BSsolQuan<0] <- 0
    
    
    ymax <- 200 #ifelse(j==6,400,400)
    plot(NULL, xlim=c(times[1],times[115]), ylim=c(0,ymax), ylab="Antal indlagte på intensiv", # var dag 130
         xlab=" ",xaxt="n", frame.plot = F)
    axis(1,at=xt,labels=timelabels)
    
    itimes <- (as.numeric(as.Date("2020-05-02"))-as.numeric(as.Date("2020-03-11"))):length(times)
    ntimes <- times[itimes]
    for (i in 4:2) { 
      polygon(c(ntimes, times[length(times)], rev(ntimes)), c(BSsolQuan[i,itimes], 0, 
                                                              rev(BSsolQuan[i-1,itimes])),col = adjustcolor(colorTable[5-i],alpha.f=0.75), border = NA)
    }
    
    grid()
    legend("topright",c( "SI:","5 - 95%","25 - 75%"),
           col=c(NA,colorTable[1:2] ),
           lty = c(NA,1,1),lwd=5) 
    
    IndlagteHos <- amdt[RegionName==regions[j]  & Date>="2020-03-11" &
                          Event=="ICU",sum(Instant),by=Date]$V1
    
    points(times[1:length(IndlagteHos)],IndlagteHos,pch=19,col=4)
    
dev.off()
