# SSEIH model of COVID 19 
# only model until arrival at hospital
# optimize on initial conditions and infectionrate per region

rm(list=ls())
library(deSolve)
library(data.table)
library(doParallel)

#### set params ####
numCores<- detectCores()
registerDoParallel(cores=numCores-1)

#### init par ####
#setwd(<where the files are>)
source("./inputsSSEIH.R")
source("./betasSSEIH.R")

ResDK <- as.numeric(as.Date("2020-03-11"))-as.numeric(as.Date("2020-01-01"))
Res2DK <-  as.numeric(as.Date("2020-04-15"))-as.numeric(as.Date("2020-01-01"))
Res3DK <-  as.numeric(as.Date("2020-05-18"))-as.numeric(as.Date("2020-01-01"))
Res4DK <-  as.numeric(as.Date("2020-06-01"))-as.numeric(as.Date("2020-01-01"))
ResSummerOn <-  as.numeric(as.Date("2020-06-27"))-as.numeric(as.Date("2020-01-01"))
ResSummerOff <-  as.numeric(as.Date("2020-08-09"))-as.numeric(as.Date("2020-01-01"))

endTimes <- as.numeric(as.Date("2020-10-01"))-as.numeric(as.Date("2020-01-01"))
times <- seq(ResDK,endTimes,1)

## The data that is loaded as 'mdt' cannot be shared at this point.
## We have asked if it can be shared and are awaiting the answer.

# load hospital data
#mdt <- fread(file = "data_for_kaare.csv")
mdt <- fread(file = "linelist_snapshot.csv")
mdt[Event=="Hospitalised",Event:="HospitalisedCase"]

mdt$Date <- as.Date(mdt$Date)
all <- mdt[,.(New=sum(New),Cumulative=sum(Cumulative),Instant=sum(Instant),RegionName="all"), by=.(AgeGroup,Event,Date)]

amdt <- rbind(mdt,all)
rm(list=c("all","mdt"))

endOptim1 <- as.numeric(max(amdt$Date))-as.numeric(as.Date("2020-01-01"))
timesOptim1 <- seq(ResDK,endOptim1,1)

#### init ldt ####

ldt <- data.table(RegionName = "all",
                  AgeGroup = "0-59",
                  Event = "HospitalisedCase",
                  Date = as.Date(ResDK,origin="2020-01-01"),
                  New = NA_real_,
                  Cumulative = NA_real_,
                  rep = NA_integer_,
                  nll = NA_real_)


#### functions ####

## SSEIH
SSEIH <- function(t,x,pl)
{
  ## Unpack state by hand 
  n   <- length(x) / 12
  S   <- x[1:n]
  E1  <- x[n + (1:n)]
  E2  <- x[2*n + (1:n)]
  E3  <- x[3*n + (1:n)]
  I1R <- x[4*n + (1:n)] # infected at home - recovering
  I2R <- x[5*n + (1:n)]
  I3R <- x[6*n + (1:n)]
  I1M <- x[7*n + (1:n)] # infected at home - will go to hospital
  I2M <- x[8*n + (1:n)]
  I3M <- x[9*n + (1:n)]
  HC  <- x[10*n + (1:n)] # cumulative hospital
  U   <- x[11*n + (1:n)] # time
  
  with(pl,
       {
         # contact reduction
         if (U[1]<Res2DK) {beta <- betaNed} else {beta <- betaFase1}
         if (U[1]>Res3DK) {beta <- betaFase2}
         if (U[1]>Res4DK) {beta <- betaFase3}
         if (U[1]>ResSummerOn & U[1]<ResSummerOff) {beta <- betaFaseSummer}
         
         # Scaling of Risk of transmission by contact
         RR <- ER*0.01
         
         
         # nyeSmittede
         nyeSmittede <- numeric(n)
         for(i in 1:n){
           nyeSmittede[i] <-S[i] * beta[i,] %*% (I1R+I1M+I2R+I2M+I3R+I3M) * RR
         }
         
         dX <- c( 
           dS = - nyeSmittede,
           
           dE1 =   nyeSmittede - 3*gammaEI1 * E1,
           dE2  = 3*gammaEI1*E1 - 3*gammaEI1*E2,
           dE3  = 3*gammaEI1*E2 - 3*gammaEI1*E3,
           
           dI1R = 3*gammaEI1 * E3*pI1R       - 3*recI1 * I1R,
           dI2R = 3*recI1*I1R - 3*recI1*I2R,
           dI3R = 3*recI1*I2R - 3*recI1*I3R,
           
           dI1M = 3*gammaEI1 * E3*(1-pI1R)   - 3*gammaI12 * I1M,
           dI2M = 3*gammaI12 * I1M - 3*gammaI12 * I2M,
           dI3M = 3*gammaI12 * I2M - 3*gammaI12 * I3M,
           
           dHC  = 3*gammaI12 * I3M,
           
           dU   = rep(1,n)
         )
         
         return(list(dX))
       }
  )
  
}

## optimizer function

nloglike3 <- function(theta,pl,j)
{
  pl['ER'] <- exp(theta[3])
  
  region <- regions[j]
  ## Initial states
  pop <- as.numeric(popDK[j,1:2])
  Ntot <- sum(pop)
  
  Inf0 <- exp(theta[1:2])
  ### init conditions ###
  S0 <- (pop-Inf0-c(5,5))/Ntot # Five recovered in each age group
  U0 <- c(ResDK,ResDK)
  
  HC0 <- c(amdt[Date=="2020-03-11" &
                  AgeGroup=="0-59" &
                  Event=="HospitalisedCase" &
                  RegionName==regions[j], Cumulative ],
           amdt[Date=="2020-03-11" &
                  AgeGroup=="60+" &
                  Event=="HospitalisedCase" &
                  RegionName==regions[j], Cumulative ])/Ntot

  E10  <-c(10,19)/24*Inf0/Ntot
  E20  <-c(9,2)/24*Inf0/Ntot
  E30  <-c(0,0)/24*Inf0/Ntot
  
  I1R0 <-c(1,2)/24*Inf0/Ntot*pI1R
  I1M0 <-c(1,2)/24*Inf0/Ntot*(1-pI1R)
  
  I2R0 <-c(2,1)/24*Inf0/Ntot*pI1R
  I2M0 <-c(2,1)/24*Inf0/Ntot*(1-pI1R)
  
  I3R0 <-c(2,1)/24*Inf0/Ntot*pI1R
  I3M0 <-c(2,1)/24*Inf0/Ntot*(1-pI1R)
  
  init0 <- c(S=S0,
             
             E1=E10,
             E2=E20, 
             E3=E30,
             
             I1R=I1R0,
             I2R=I2R0,
             I3R=I3R0,
             
             I1M=I1M0,
             I2M=I2M0,
             I3M=I3M0,
             
             HC=HC0,
             
             U=U0)
  
  sol <- ode(init0, timesOptim1, SSEIH, pl)
  
  tmp <- sol[,22:23]*Ntot
  yp <- cbind(c(NA,diff(tmp[,1])),c(NA,diff(tmp[,2])))
  
  ydU <- amdt[Date>="2020-03-11" & 
                AgeGroup=="0-59" & 
                Event=="HospitalisedCase" &
                RegionName==regions[j], New ] 
  ydO <- amdt[Date>="2020-03-11" & 
                AgeGroup=="60+" & 
                Event=="HospitalisedCase" &
                RegionName==regions[j], New ] 
  
  res <- sum(dpois(ydU,yp[,1], log=TRUE) , na.rm=TRUE) + sum(dpois(ydO,yp[,2], log=TRUE) , na.rm=TRUE)
  return(-res)
}


#### bootstrap ####

sceNames <- c("fase2", "fase2f50", "fase2f100", "fase3","fase3f50",  "fase3f100", 
              "fase3+", "fase3+f50","fase3+f100", "fase4", "fase4f50","fase4f100" )
#sceNames <- c("fase4", "fase2")

nrep<- 500*1.25
tic <- Sys.time()
sceIndex<-12
feRet <- foreach(sceIndex=1:length(sceNames), .packages = c("data.table", "deSolve"), .verbose = TRUE) %dopar% {
  sceName <- sceNames[sceIndex]
  print(paste(Sys.time(),sceName))
  #### init ldt ####
  
  ldt <- data.table(RegionName = "all",
                    AgeGroup = "0-59",
                    Event = "HospitalisedCase",
                    Date = as.Date(ResDK,origin="2020-01-01"),
                    New = NA_real_,
                    Cumulative = NA_real_,
                    rep = NA_integer_,
                    nll = NA_real_)
  plAll <- list()

# region
for (j in 1:5)
{
  
  set.seed(1234)
  
  for (i in 1:nrep)
  {
    # efficiency of restrictions
    ER <- NA
    
    ## Two daily interacting groups 

    n <- 2
    boffdia <- crunif(input$r.betaBR1)
    betaNed  <- array(c(crunif(input$r.betaWIltR1),boffdia,
                        boffdia,crunif(input$r.betaWIgeR1)),c(2,2))
    
    boffdia <- crunif(input$rd.betaBR2)
    betaRes2 <- array(c(crunif(input$rd.betaWIltR2),boffdia,
                        boffdia,crunif(input$rd.betaWIgeR2)),c(2,2))
    betaFase1 <- betaRes2 + betaNed
    
    boffdia <- runif(3,0,2)
    mscal <- array( boffdia[c(1,3,3,2)], c(2,2))
    betaFase2 <- array((betasPast[[3]]-betasPast[[2]])[c(1,3,3,2)], c(2,2)) * mscal + betaFase1

    boffdia <- runif(3,0,2)
    mscal <- array( boffdia[c(1,3,3,2)], c(2,2))
    betaFase3 <- array((betas[[sceName]]-betasPast[[3]])[c(1,3,3,2)], c(2,2)) * mscal + betaFase2
    
    boffdia <- runif(3,0,2)
    mscal <- array( boffdia[c(1,3,3,2)], c(2,2))
    betaFaseSummer <- array((betas[[paste0(sceName,"S")]] - betas[[sceName]])[c(1,3,3,2)], c(2,2)) * mscal + betaFase3
        
    betaFaseSummer <- array(pmax(0.0, betaFaseSummer),c(2,2))
    # recovery rates for different stages
    recI1 <- c(1/crunif(input$r.recI11),1/crunif(input$r.recI12)) # 1/dage før rask udenfor hospital
    
    # rates going from one state to the next
    gammaEI1 <- c(1/crunif(input$r.kE1),1/crunif(input$r.kE2)) # rate going from E to I
    gammaI12 <- c(1/crunif(input$r.k1),1/crunif(input$r.k2)) # rate going from IxM to HC
    
    # percentage Recover or Moving on
    pI1R <- 1-c(crunif(input$r.pI1R1),crunif(input$r.pI1R2))*0.01
    
    pl <- list(betaNed=betaNed, betaFase1=betaFase1, betaFase2=betaFase2, betaFase3=betaFase3, betaFaseSummer=betaFaseSummer,
           gammaEI1=gammaEI1, gammaI12=gammaI12, recI1=recI1, pI1R=pI1R,ER=ER) #  variables

    nll3w <- function(theta)
    {
      nloglike3(theta,pl,j)
    }
    
    # optimize on initial conditions
    theta3 <- log(c(mean(input$a.InfInit1[j,]),mean(input$a.InfInit2[j,]),mean(input$a.ER[j,])))
    rll <- nlminb(theta3, nll3w, lower=c(6,6,2), upper = c(12,12,5.3))
    
    pl['ER'] <- exp(rll$par)[3]
    pl$rep <- i
    pl$region <- j
    pl$sceName <- sceName
    pl$nll <-  rll$objective
    pl$Inf01 <- exp(rll$par)[1]
    pl$Inf02 <- exp(rll$par)[2]

    ## Saving parameters
    plAll[[i + nrep*(j-1)]] <- pl

    theta <- rll$par
    
    # run full time
    
    region <- regions[j]
    ## Initial states
    pop <- as.numeric(popDK[j,1:2])
    Ntot <- sum(pop)
    
    Inf0 <- exp(theta[1:2])
    ### init conditions ###
    S0 <- (pop-Inf0-c(5,5))/Ntot
    U0 <- c(ResDK,ResDK)
    
    HC0 <- c(amdt[Date=="2020-03-11" &
                    AgeGroup=="0-59" &
                    Event=="HospitalisedCase" &
                    RegionName==regions[j], Cumulative ],
             amdt[Date=="2020-03-11" &
                    AgeGroup=="60+" &
                    Event=="HospitalisedCase" &
                    RegionName==regions[j], Cumulative ])/Ntot

    E10  <-c(10,19)/24*Inf0/Ntot
    E20  <-c(9,1)/24*Inf0/Ntot
    E30  <-c(0,0)/24*Inf0/Ntot
    
    I1R0 <-c(1,2)/24*Inf0/Ntot*pI1R
    I1M0 <-c(1,2)/24*Inf0/Ntot*(1-pI1R)
    
    I2R0 <-c(2,1)/24*Inf0/Ntot*pI1R
    I2M0 <-c(2,1)/24*Inf0/Ntot*(1-pI1R)
    
    I3R0 <-c(2,1)/24*Inf0/Ntot*pI1R
    I3M0 <-c(2,1)/24*Inf0/Ntot*(1-pI1R)
    
    init0 <- c(S=S0,
               
               E1=E10,
               E2=E20,
               E3=E30,
               
               I1R=I1R0,
               I2R=I2R0,
               I3R=I3R0,
               
               I1M=I1M0,
               I2M=I2M0,
               I3M=I3M0,
               
               HC=HC0,
               
               U=U0)
    
    sol <- ode(init0,times,SSEIH,pl)

    tmp <- sol[,22:23]*Ntot
    yp <- cbind(c(NA,diff(tmp[,1])),c(NA,diff(tmp[,2])))

    tmp2 <- data.table(RegionName = region,
                       AgeGroup = "0-59",
                       Event = "HospitalisedCase",
                       Date = as.Date(sol[,1],origin="2020-01-01"),
                       New = yp[,1],
                       Cumulative = tmp[,1],
                       rep = i,
                       nll = rll$objective)

    ldt <- rbind(ldt,tmp2)

    tmp2 <- data.table(RegionName = region,
                       AgeGroup = "60+",
                       Event = "HospitalisedCase",
                       Date = as.Date(sol[,1],origin="2020-01-01"),
                       New = yp[,2],
                       Cumulative = tmp[,2],
                       rep = i,
                       nll = rll$objective)

    ldt <- rbind(ldt,tmp2)
    
  }
}

all <- ldt[,.(New=sum(New), Cumulative=sum(Cumulative),
              RegionName="all", nll = sum(nll)),
           by=.(AgeGroup,Event,Date,rep)]
ldt <- rbind(ldt,all)


 list(ldt=ldt, pdt=plAll)

}

save(feRet,sceNames,nrep, file=paste0("runAll_pl_",length(sceNames),sceNames[1],".RData"))


toc <- Sys.time()
toc-tic


