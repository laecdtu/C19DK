# fold of SSEIH data

library(data.table)
library(plotly)
library(fields)
library(Hmisc)
library(purrr)


prepath <-  "."
filenames <- list.files(paste0(prepath,"/runs"))
pngmode <- FALSE
pdfmode <- TRUE

BSsel <- 0.20 # throw out bootstraps that does not fit the initial phase

ResDK <- as.numeric(as.Date("2020-03-11"))-as.numeric(as.Date("2020-01-01"))
Res2DK <-  as.numeric(as.Date("2020-04-15"))-as.numeric(as.Date("2020-01-01"))
Res3DK <-  as.numeric(as.Date("2020-05-18"))-as.numeric(as.Date("2020-01-01"))
Res4DK <-  as.numeric(as.Date("2020-06-01"))-as.numeric(as.Date("2020-01-01"))
ResSummerOn <-  as.numeric(as.Date("2020-06-27"))-as.numeric(as.Date("2020-01-01"))
ResSummerOff <-  as.numeric(as.Date("2020-08-09"))-as.numeric(as.Date("2020-01-01"))



endTimes <- as.numeric(as.Date("2020-10-01"))-as.numeric(as.Date("2020-01-01"))
times <- seq(ResDK,endTimes,1)

xt <- c(ResDK,as.numeric(as.Date(format(as.Date((4:24)*30,origin="2020-01-01"),
                                        "%Y-%m-01")))-
          as.numeric(as.Date("2020-01-01")))
timelabels <- format(as.Date(xt,origin="2020-01-01"),"%d %b")
# Check 'timelabels' depends on locale on system
#timelabels <-gsub("Jul","jul",gsub("Jun","jun",gsub("Mar","mar",
#                                       gsub("May","maj",gsub("Apr",replacement = "apr",timelabels)))))


## The data that is loaded as 'mdt' cannot be shared at this point.
## We have asked if can be shared and are awaiting the answer.
prepath <- "."
#mdt <- fread(file = paste0(prepath,"/data_for_kaare.csv"))
mdt <- fread(file = paste0(prepath,"/linelist_snapshot.csv"))
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


mainTitles <- c("Fortsat anden fase af genåbning",
"Fortsat anden fase af genåbning + halvdelen af voksne holder fysisk afstand",
"Fortsat anden fase af genåbning + fysisk afstand som før corona",
"Fase 3 åbner d. 1/6",
"Fase 3 åbner d. 1/6 + halvdelen af voksne holder fysisk afstand",
"Fase 3 åbner d. 1/6 + fysisk afstand som før corona",
"Udvidet fase 3 åbner d. 1/6",
"Udvidet fase 3 åbner d. 1/6 + halvdelen af voksne holder fysisk afstand",
"Udvidet fase 3 åbner d. 1/6 + fysisk afstand som før corona",
"Fase 3 & 4 åbner d. 1/6",
"Fase 3 & 4 åbner d. 1/6 + halvdelen af voksne holder fysisk afstand",
"Fase 3 & 4 åbner d. 1/6 + fysisk afstand som før corona"
)

mainTitles <- c("Fase 2",
                "Fase 2 + halvdelen af voksne holder fysisk afstand",
                "Fase 2 + fysisk afstand som før corona",
                "Udvidelse med fase 3",
                "Udvidelse med fase 3 + halvdelen af voksne holder fysisk afstand",
                "Udvidelse med fase 3 + fysisk afstand som før corona",
                "Udvidelse med fase 3 samt udvidet fase 3",
                "Udvidelse med fase 3 samt udvidet fase 3 + halvdelen af voksne holder fysisk afstand",
                "Udvidelse med fase 3 samt udvidet fase 3 + fysisk afstand som før corona",
                "Udvidelse med fase 3, udvidet fase 3 samt fase 4",
                "Udvidelse med fase 3, udvidet fase 3 samt fase 4 + halvdelen af voksne holder fysisk afstand",
                "Udvidelse med fase 3, udvidet fase 3 samt fase 4 + fysisk afstand som før corona"
)



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

nHorPred <- 100
oldProb.early <- stateAll.fun(0:nHorPred, probsData = allProbs.old.early, N=1)
youngProb.early <- stateAll.fun(0:nHorPred, probsData = allProbs.young.early, N=1)
oldProb.weighted <- stateAll.fun(0:nHorPred, probsData = allProbs.old.weighted, N=1)
youngProb.weighted <- stateAll.fun(0:nHorPred, probsData = allProbs.young.weighted, N=1)



HosFoldTheis2 <- function(cumIndHosU60,cumIndHosO60,initIndHosU60,initIndHosO60)
{
  newHosU60 <- c(initIndHosU60,diff(cumIndHosU60 ))
  newHosO60 <- c(initIndHosO60,diff(cumIndHosO60 ))
  
  nHorPred <- 100
  
  allHosU60 <- matrix(0,nrow=length(newHosU60)+nHorPred,ncol=4)
  allHosO60 <- matrix(0,nrow=length(newHosO60)+nHorPred,ncol=4)
  
  for (i in 1:21) 
  {
    allHosU60[i:(i+nHorPred),] <- allHosU60[i:(i+nHorPred),] + youngProb.early*newHosU60[i]
    allHosO60[i:(i+nHorPred),] <- allHosO60[i:(i+nHorPred),] + oldProb.early*newHosO60[i]
  }
  
  for (i in 22:length(newHosO60)) 
  {
    allHosU60[i:(i+nHorPred),] <- allHosU60[i:(i+nHorPred),] + youngProb.weighted * newHosU60[i]
    allHosO60[i:(i+nHorPred),] <- allHosO60[i:(i+nHorPred),] + oldProb.weighted * newHosO60[i]
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
  
  for (i in 1:21) 
  {
    allHosU60[i:(i+nHorPred),] <- allHosU60[i:(i+nHorPred),] + youngProb.early * newHosU60[i]
    allHosO60[i:(i+nHorPred),] <- allHosO60[i:(i+nHorPred),] + oldProb.early * newHosO60[i]
  }
  
  for (i in 22:length(newHosO60)) 
  {
    allHosU60[i:(i+nHorPred),] <- allHosU60[i:(i+nHorPred),] + youngProb.weighted * newHosU60[i]
    allHosO60[i:(i+nHorPred),] <- allHosO60[i:(i+nHorPred),] + oldProb.weighted * newHosO60[i]
  }
  
  allHosU60 <- allHosU60[1:length(newHosU60),]
  allHosO60 <- allHosO60[1:length(newHosO60),]

  return( allHosU60[,2] + allHosO60[,2] )
}


## Loading saved data:
feRetSamlet <- list()
for(kk in filenames){
  load(paste0(prepath,"/runs/",kk))
  for (i in 1:length(idScenarier)){
    feRetSamlet[[idScenarier[i]]] <- feRet[[i]]
  }
}
rm(feRet)

# Simple check
sapply(feRetSamlet,function(x)nrow(x$ldt))
sapply(feRetSamlet,function(x)length(x$pdt))


library(doParallel)

numCores<- detectCores()
registerDoParallel(cores=numCores-1)

tic <- Sys.time()
listInfHos<-foreach(id =1:12, .packages = "data.table") %dopar%{
  ldt <- feRetSamlet[[id]]$ldt

  InfHos <- array(NA,dim=c(uniqueN(ldt$Date),
                           uniqueN(ldt$rep)-1,
                           uniqueN(ldt$RegionName)))
  
  for (j in 1:(uniqueN(ldt$rep)-1))
  {
    for (i in 1:uniqueN(ldt$RegionName))
    {
      InfHos[,j,i] <-  HosFoldTheis2(ldt[RegionName==regions[i] & AgeGroup=="0-59" & rep==j &
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

##ICU
tic <- Sys.time()
listInfICU <- foreach(id =1:12, .packages = "data.table") %dopar%{
  ldt <- feRetSamlet[[id]]$ldt

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

## ICU done HERE


sce <- 1
InfHos <- listInfHos[[3]]

mapSce <- c(1:12)
sceId <- 1

## saving quantiles to plot later
listBSQHos <- list()
for (sceId in 1:12){
  sce <- mapSce[sceId]
  InfHos <- listInfHos[[sceId]] 
  pdt <- as.data.table(t(sapply(feRetSamlet[[sce]]$pdt, function(x)c(region=x$region, rep=x$rep, nll=x$nll))))
  print(c(sceId, sce))
  listTmpReg <- list()
  for (j in 1:6)
  {
    solBS <- InfHos[,,j]
    if(j<6){
      weights <-   pdt[region==j,-nll]
    } else {
      weights <- pdt[,-mean(nll),by=rep]$V1
    }
    weights <- (1-0.9*(weights-min(weights))/(max(weights)-min(weights)))
    weights[weights < quantile(weights,BSsel)] <- 0
    
    BSsolQuan <-  apply(solBS, 1, function(x) 
      wtd.quantile(x, probs = c(.05,.25,.75,.95), weights = weights ) )

    # Observational noise
    PredMed <-  apply(solBS, 1, function(x) 
      wtd.quantile(x, probs = c(.5), weights = weights ) )
    ObsNoise <- apply(t(PredMed),2,function(x) qpois(c(.05,.25,.75,.95),x)-x)
    
    BSsolQuan <-  BSsolQuan + ObsNoise
    BSsolQuan[BSsolQuan<0] <- 0
    listTmpReg[[j]] <- BSsolQuan
  }
  listBSQHos[[sceId]] <- listTmpReg
}

listBSQICU <- list()
for (sceId in 1:12){
  sce <- mapSce[sceId]
  InfICU <- listInfICU[[sceId]] 
  pdt <- as.data.table(t(sapply(feRetSamlet[[sce]]$pdt, function(x)c(region=x$region, rep=x$rep, nll=x$nll))))
  print(c(sceId, sce))
  listTmpReg <- list()
  for (j in 1:6)
  {
    solBS <- InfICU[,,j]
    if(j<6){
      weights <-   pdt[region==j,-nll]
    } else {
      weights <-pdt[,-mean(nll),by=rep]$V1
    }
    weights <- (1-0.9*(weights-min(weights))/(max(weights)-min(weights)))
    weights[weights < quantile(weights,BSsel)] <- 0
    
    BSsolQuan <-  apply(solBS, 1, function(x) 
      wtd.quantile(x, probs = c(.05,.25,.75,.95), weights = weights ) )
    
    # Observational noise
    PredMed <-  apply(solBS, 1, function(x) 
      wtd.quantile(x, probs = c(.5), weights = weights ) )
    ObsNoise <- apply(t(PredMed),2,function(x) qpois(c(.05,.25,.75,.95),x)-x)
    
    BSsolQuan <-  BSsolQuan + ObsNoise
    BSsolQuan[BSsolQuan<0] <- 0
    listTmpReg[[j]] <- BSsolQuan
  }
  listBSQICU[[sceId]] <- listTmpReg
}

#save.image()

## Plotting ####
## Credit: https://stackoverflow.com/questions/30765866/get-margin-line-locations-in-log-space/30835971#30835971
line2user <- function(line, side){
  lh <- par('cin')[2] * par('cex')*par('lheight')
  x_off <- diff(grconvertX(0:1, 'inches','user'))
  y_off <- diff(grconvertY(0:1, 'inches','user'))
  switch(side,
         '1' = par('usr')[3] - line * y_off * lh,
         '2' = par('usr')[1] - line * x_off * lh,
         '3' = par('usr')[4] + line * y_off * lh,
         '4' = par('usr')[2] + line * x_off * lh,
          stop("side must be in 1:4", call. = FALSE))
}


png(paste0("indlagteC.png"), width = 1600, height=2100, pointsize = 36)
par(mfrow=c(4,2))
par(mar=c(2.,3.3,3.3,1)+.1, mgp=c(2.1,0.7,0))

for (sceId in (1:4)*3-0){
  for(reg in (2:1)){
      
    BSsolQuan<- listBSQHos[[sceId]][[reg]]
    ymax <- ifelse(reg==6,2000,1500)
    plot(NULL, xlim=c(times[1],times[205]), ylim=c(0,ymax), ylab="Antal indlagte på hospital", 
         xlab=" ",xaxt="n", frame.plot = F)
    axis(1,at=xt,labels=timelabels)
    
    itimes <- (as.numeric(as.Date("2020-05-15"))-as.numeric(as.Date("2020-03-11"))):length(times)
    ntimes <- times[itimes]
    for (i in 4:2) { 
      polygon(c(ntimes, times[length(times)], rev(ntimes)), c(BSsolQuan[i,itimes], 0, 
                                                              rev(BSsolQuan[i-1,itimes])),col = adjustcolor(colorTable[5-i],alpha.f=0.75), border = NA)
    }
    
    grid()
    legend("topleft",c( "SI:","5 - 95%","25 - 75%"),
           col=c(NA,colorTable[1:2] ),
           lty = c(NA,1,1),lwd=5) 
    
    IndlagteHos <- amdt[RegionName==regions[reg]  & Date>="2020-03-11" &
                          Event=="HospitalisedCase",sum(Instant),by=Date]$V1
    
    points(times[1:length(IndlagteHos)],IndlagteHos,pch=19,col=4)
    title(regionsTitle[reg], cex.main=1, line=0.7)  
  }
  text(line2user(line=mean(par('mar')[c(2,4)]), side = 2),
       line2user(line=2.3, side=3), mainTitles[sceId], xpd=NA, cex=1, font=2)
}
dev.off()

## Plot ICU ####
png(paste0("ICU_A.png"), width = 1600, height=2100, pointsize = 36)
par(mfrow=c(4,2))
par(mar=c(2.,3.3,3.3,1)+.1, mgp=c(2.1,0.7,0))

for (sceId in (1:4)*3-2){  #Update index
  for(reg in (2:1)){
    BSsolQuan<- listBSQICU[[sceId]][[reg]]
    ymax <- ifelse(reg==6,500,500)
    plot(NULL, xlim=c(times[1],times[205]), ylim=c(0,ymax), ylab="Antal indlagte på intensiv", 
         xlab=" ",xaxt="n", frame.plot = F)
    axis(1,at=xt,labels=timelabels)
    
    itimes <- (as.numeric(as.Date("2020-05-15"))-as.numeric(as.Date("2020-03-11"))):length(times)
    ntimes <- times[itimes]
    for (i in 4:2) { 
      polygon(c(ntimes, times[length(times)], rev(ntimes)), c(BSsolQuan[i,itimes], 0, 
                                                              rev(BSsolQuan[i-1,itimes])),col = adjustcolor(colorTable[5-i],alpha.f=0.75), border = NA)
    }
    
    grid()
    legend("topleft",c( "SI:","5 - 95%","25 - 75%"),
           col=c(NA,colorTable[1:2] ),
           lty = c(NA,1,1),lwd=5) 
    
    IndlagteHos <- amdt[RegionName==regions[reg]  & Date>="2020-03-11" &
                          Event=="ICU",sum(Instant),by=Date]$V1
    
    points(times[1:length(IndlagteHos)],IndlagteHos,pch=19,col=4)
    title(regionsTitle[reg], cex.main=1, line=0.7)  
  }
  text(line2user(line=mean(par('mar')[c(2,4)]), side = 2),
       line2user(line=2.3, side=3), mainTitles[sceId], xpd=NA, cex=1, font=2)
}
dev.off()


## Status udvalgte datoer ####

idLab <- c(5:7)
rowID <- xt[idLab]-xt[1]
dagLab <- timelabels[idLab]

matSum <- matrix(NA, nrow=72, ncol=12)

for (fa in 1:3){
  for (sce in (1:4)){  #Update index
    for(reg in (2:1)){
      for(dag in 1:3){
      tmp <- (reg-1)*36 + (fa-1)*12 + (dag-1)*4 + sce
      print(paste(paste(tmp, collapse = " "), mainTitles[(sce-1)*3 + fa]))
      matSum[tmp, 5:8] <- round(t(listBSQHos[[ (sce-1)*3 + fa ]][[reg]][,rowID[dag]]))
      matSum[tmp, 9:12] <- round(t(listBSQICU[[ (sce-1)*3 + fa ]][[reg]][,rowID[dag]]))
      matSum[tmp, 1:4] <- c(reg, c("FA","halv FA","normal FA")[fa],dagLab[dag], mainTitles[(sce-1)*3+1])
      }
    }
  }
}
write.csv(matSum, file="matSum_combined_quantiles_20200518.csv")

## Estimation af Reff på udvalgte datoer ####

rapDates <- c("2020-07-01", "2020-08-01","2020-09-01")
timesReff <- (as.numeric(as.Date(rapDates))-as.numeric(as.Date("2020-03-11")))

sce <- 1
idDag<- 1
a <- matrix(NA, nrow=36,ncol=12)
for (region in c("Hovedstaden","Zealand")){
for (idDag in 1:3) {
  print(c(region, idDag))
Reff <- sapply(1:12,function(sce){
  
  tid <- 0:5
  tmp <- feRetSamlet[[sce]]$ldt[-1,]
  tmp  <- tmp[RegionName==region,]
  gr <- sapply(1:max(tmp$rep), function(i){
    tmp2 <- tmp[ rep==i,sum(New),keyby=Date]$V1[timesReff[idDag]+(-5:0)]
    coef(lm(log(tmp2)~tid))[2]
  } )
  
  pdt <- feRetSamlet[[sce]]$pdt
  weights <-pdt[,-mean(as.numeric(nll)),by=rep]$V1
  weights <- (1-0.9*(weights-min(weights))/(max(weights)-min(weights)))
  #weights <- rep(1,500)
  weights[weights < quantile(weights,BSsel)] <- 0
  
  1+wtd.quantile(gr, probs = c(.05,.25,.75,.95), weights = weights ) *6.5
})

a[(idDag-1)*12+1:12, 5:8+ifelse(region=="Zealand",0,4)] <- t(Reff)
}
}
a[,4] <- 1:12

a[,3] <- rep(1:3, each=12) # Dato
a[,2] <- 1:3 # Fysisk afstand
  
round(a,2)
#write.csv(a, file="Reff_20200519.csv")

## Eksempel på scenarier ####
sceId <- 3
mainTitles[sceId]
reg<- 2
BSsolQuan<- listBSQHos[[sceId]][[reg]]
ymax <- ifelse(reg==6,2000, 2000)

png("eksempel_scenarier.png", width=1400, height = 1000, pointsize = 30)
plot(NULL, xlim=c(times[1],times[205]), ylim=c(0,ymax), ylab="Antal indlagte på hospital", 
     xlab=" ",xaxt="n", frame.plot = F)
axis(1,at=xt,labels=timelabels)

itimes <- (as.numeric(as.Date("2020-05-15"))-as.numeric(as.Date("2020-03-11"))):length(times)
ntimes <- times[itimes]
for (i in 4:2) { 
  polygon(c(ntimes, times[length(times)], rev(ntimes)), c(BSsolQuan[i,itimes], 0, 
                                                          rev(BSsolQuan[i-1,itimes])),col = adjustcolor(colorTable[5-i],alpha.f=0.75), border = NA)
}

grid()
legend("topleft",c("5 - 95%","25 - 75%",paste("Simulation",1:6)),title = "SI:",
       col=c(colorTable[1:2], 1:3,5:7 ),
       lty = 1,lwd=c(5,5,rep(2,6))*1.5, bg="white") 


IndlagteHos <- amdt[RegionName==regions[reg]  & Date>="2020-03-11" &
                      Event=="HospitalisedCase",sum(Instant),by=Date]$V1

points(times[1:length(IndlagteHos)],IndlagteHos,pch=19,col=4)
title(paste(mainTitles[sceId],"for",regionsTitle[reg]), cex.main=1, line=0.7)  



InfHos <- listInfHos[[sceId]] 
pdt <- feRetSamlet[[sce]]$pdt
j <- reg
  solBS <- InfHos[,,j]
    weights <-   pdt[region==j,-as.numeric(nll)]
  weights <- (1-0.9*(weights-min(weights))/(max(weights)-min(weights)))
  weights[weights < quantile(weights,BSsel)] <- 0
selBS <- solBS[,weights!=0]

itimes <- (as.numeric(as.Date("2020-05-15"))-as.numeric(as.Date("2020-03-11"))):length(times)
ntimes <- times[itimes]
subset <- c(1,5,6,8:10)
for (i in 1:length(subset)) { 
  lines(ntimes, selBS[itimes,subset[i]*10],col=i+(i>3), lwd=3)
}

dev.off()


#################################### #
## Old plots  (prior to May 20, 2020)####

if(pdfmode) {
  pdf(file = "Indlagte_alle.pdf", width=6, height=8)
    par(mfrow=c(3,2))
  par(mar=c(2.,3.3,4.3,1)+.1, mgp=c(2.1,0.7,0))#,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))
}
for (sceId in 1:12){
  sce <- mapSce[sceId]
  InfHos <- listInfHos[[sceId]] 
  pdt <- feRetSamlet[[sce]]$pdt
  print(c(sceId, sce))
for (j in 1:6)
{
solBS <- InfHos[,,j]
if(j<6){
  weights <-   pdt[region==j,-as.numeric(nll)]
} else {
  weights <-pdt[,-mean(as.numeric(nll)),by=rep]$V1
}
weights <- (1-0.9*(weights-min(weights))/(max(weights)-min(weights)))
#weights <- rep(1,500)
weights[weights < quantile(weights,BSsel)] <- 0

BSsolQuan <-  apply(solBS, 1, function(x) 
  wtd.quantile(x, probs = c(.05,.25,.75,.95), weights = weights ) )
#  wtd.quantile(x, probs = c(seq(0,.9,0.1),.95), weights = weights ) )
#quantile(x, probs = c(.05,seq(0.1,1,0.1))) )

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

ymax <- ifelse(j==6,2000,1000)
plot(NULL, xlim=c(times[1],times[205]), ylim=c(0,ymax), ylab="Antal indlagte på hospital", # var dag 130
     xlab=" ",xaxt="n", frame.plot = F)
axis(1,at=xt,labels=timelabels)

#    polygon(c(times, times[length(times)]), c(BSsolQuan[1,], 0),col = adjustcolor(colorTable[10],alpha.f=0.75), border = NA)
itimes <- (as.numeric(as.Date("2020-05-02"))-as.numeric(as.Date("2020-03-11"))):length(times)
ntimes <- times[itimes]
#    polygon(c(times, times[length(times)]), c(BSsolQuan[1,], 0),col = adjustcolor(colorTable[10],alpha.f=0.75), border = NA)
for (i in 4:2) { 
  polygon(c(ntimes, times[length(times)], rev(ntimes)), c(BSsolQuan[i,itimes], 0, 
                                                          rev(BSsolQuan[i-1,itimes])),col = adjustcolor(colorTable[5-i],alpha.f=0.75), border = NA)
}

#lines(times, BSsolQuan[6,],type="l",lwd=3)
#abline(h=Hoscap,col=adjustcolor("black",alpha.f = 0.75),lwd=3,lty=2)
grid()
legend("top",c( "SI:","5 - 95%","25 - 75%"),
       col=c(NA,colorTable[1:2] ),
       lty = c(NA,1,1),lwd=5) 

IndlagteHos <- amdt[RegionName==regions[j]  & Date>="2020-03-11" &
       Event=="HospitalisedCase",sum(Instant),by=Date]$V1

#IndlagteHos <- c(10,NA,23,NA,28,62,82,129,153,183,206,232,254,301,350,386,430,459,499,533,529,535,525,517,507,504)
points(times[1:length(IndlagteHos)],IndlagteHos,pch=19,col=4)
#abline(v=as.numeric(as.Date("2020-03-28"))-as.numeric(as.Date("2020-01-01"))) # 
title(regionsTitle[j])  
if (pngmode) dev.off()
}
#mtext(paste0("B4.",sceId-1,": ",mainTitles[sce]), side=3, line=-1.5, outer=TRUE, cex=0.9)
mtext(paste0(mainTitles[sce]), side=3, line=-1.5, outer=TRUE, cex=0.9)
}

if(pdfmode) dev.off()



################################### #
## Plot ICU ####
if(pdfmode) {
  pdf(file = "ICU_alle.pdf", width=6, height=8)
  par(mfrow=c(3,2))
  par(mar=c(2.,3.3,4.3,1)+.1, mgp=c(2.1,0.7,0))#,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))
}
for (sceId in 1:12){
  sce <- mapSce[sceId]
  InfICU <- listInfICU[[sce]]
  pdt <- feRetSamlet[[sce]]$pdt
  
  for (j in 1:6)
  {
    solBS <- InfICU[,,j]
    if(j<6){
      weights <-   pdt[region==j,-as.numeric(nll)]
    } else {
      weights <-pdt[,-mean(as.numeric(nll)),by=rep]$V1
    }
    weights <- (1-0.9*(weights-min(weights))/(max(weights)-min(weights)))
    #weights <- rep(1,500)
    weights[weights < quantile(weights,BSsel)] <- 0
    
    BSsolQuan <-  apply(solBS, 1, function(x) 
      wtd.quantile(x, probs = c(.05,.25,.75,.95), weights = weights ) )
    #  wtd.quantile(x, probs = c(seq(0,.9,0.1),.95), weights = weights ) )
    #quantile(x, probs = c(.05,seq(0.1,1,0.1))) )
    
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
    
    ymax <- ifelse(j==6,500,500)
    plot(NULL, xlim=c(times[1],times[205]), ylim=c(0,ymax), ylab="Antal indlagte på intensiv", # var dag 130
         xlab=" ",xaxt="n", frame.plot = F)
    axis(1,at=xt,labels=timelabels)
    
    #    polygon(c(times, times[length(times)]), c(BSsolQuan[1,], 0),col = adjustcolor(colorTable[10],alpha.f=0.75), border = NA)
    itimes <- (as.numeric(as.Date("2020-05-02"))-as.numeric(as.Date("2020-03-11"))):length(times)
    ntimes <- times[itimes]
    #    polygon(c(times, times[length(times)]), c(BSsolQuan[1,], 0),col = adjustcolor(colorTable[10],alpha.f=0.75), border = NA)
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
    #  wtd.quantile(x, probs = c(seq(0,.9,0.1),.95), weights = weights ) )
    #quantile(x, probs = c(.05,seq(0.1,1,0.1))) )
    
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
  #  wtd.quantile(x, probs = c(seq(0,.9,0.1),.95), weights = weights ) )
  #quantile(x, probs = c(.05,seq(0.1,1,0.1))) )
  
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
  #title(mainTitlesSub[sceId])
  
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
