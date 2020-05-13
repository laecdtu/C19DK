
library(shiny)
library(shinydashboard)
#library(shinyalert)
#library(plotly)
library(deSolve)
library(epiR)
library(fields)
library(Hmisc)
library(purrr)

testmode <- !TRUE

# colors for plot
colorTable<- rev(designer.colors(10, c( "#1FD082","#F6D04D", "#E83F48") ))

#### parameters ####
ResDK <- as.numeric(as.Date("2020-03-11"))-as.numeric(as.Date("2020-01-01"))
Res2DK <-  as.numeric(as.Date("2020-04-14"))-as.numeric(as.Date("2020-01-01"))
InitPer <- 1:35
xt <- c(ResDK,as.numeric(as.Date(format(as.Date((4:24)*30,origin="2020-01-01"),
                                       "%Y-%m-01")))-
         as.numeric(as.Date("2020-01-01")))
timelabels <- format(as.Date(xt,origin="2020-01-01"),"%d %b")
itl <- format(as.Date(ResDK+InitPer-1,origin="2020-01-01"),"%d %b")
endOfTime <- as.numeric(as.Date("2020-09-01"))-as.numeric(as.Date("2020-01-01"))   #730
x0 <- seq(ResDK,endOfTime,1)
times <- seq(ResDK,endOfTime,1)

cols <- c("#990000","#FC7634","#F6D04D")
popDK <- c(4.3e6,1.8e6)
NtotDK <- sum(popDK)
popHSJ <- c(2e6,6e5) # H+SJ
NtotHSJ <- sum(popHSJ)
pop <- popDK
Ntot <- sum(pop)

IndlagteHos <- c(10,NA,23,NA,28,62,82,129,153,183,206,232,254,301,350,386,430,459,499,533,529,535,525,517,507,504,503,472,453, 433,401, 408, 396, 388)  # 13/4
IndlagteInt <- c(0,0,4,NA,2,10,18,24,30,37,42,46,55,69,87,94,109,121,131,139,145,146,153,143,142,144,139,127,127, 120, 113,  106, 104, 100)
IndlagteHos <- IndlagteHos - IndlagteInt

nDaysData <- length(IndlagteHos)
if (length(IndlagteInt)!=nDaysData) stop("Input data of different lengths")

# function to easier call runif
crunif <- function(invec)
{
  if (length(invec)!=2) stop("wrong size")
  runif(1,invec[1],invec[2])
}

#### input dates for interactive mitigation ####
x <- c(ResDK,Res2DK,c(as.numeric(as.Date(format(as.Date((5:11)*30,origin="2020-01-01"), "%Y-%m-01"))),
                      as.numeric(as.Date("2022-01-01",origin="2020-01-01"))) - as.numeric(as.Date("2020-01-01")))
#y <- c(1,1,1,.8,.6,.4,.1,.1,rep(0,2))
y <- rep(1,10)
pltimelabels <- format(as.Date(x,origin="2020-01-01"),"%d %b %Y")
xx <- x

# Define server logic 
shinyServer(function(input, output, session) {
  
  #### dataBS #### 
  dataBS <- reactive({
    
    set.seed(1234)

    
    SSEIRD <- function(t,x,p)
    {
      ## Unpack state by hand 
      n   <- length(x) / 19
      S   <- x[1:n]
      E1   <- x[n + (1:n)]
      I1R <- x[2*n + (1:n)] # infected at home - recovering
      I1M <- x[3*n + (1:n)] # infected at home - will go to hospital
      HR <- x[4*n + (1:n)] # at hospital - recovering
      HM <- x[5*n + (1:n)] # at hospital - will go to ICU
      CR <- x[6*n + (1:n)] # infected intensive care - recovering
      CD <- x[7*n + (1:n)] # infected intensive care - dying
      R   <- x[8*n + (1:n)]
      D   <- x[9*n + (1:n)]
      U   <- x[10*n + (1:n)]
      HC <- x[11*n + (1:n)] # cumulative hospital
      CC <- x[12*n + (1:n)] # cumulative hospital
      E2 <- x[13*n + (1:n)]
      E3 <- x[14*n + (1:n)]
      I2R <- x[15*n + (1:n)]
      I3R <- x[16*n + (1:n)]
      I2M <- x[17*n + (1:n)]
      I3M <- x[18*n + (1:n)]
      
      
      dS <- numeric(n)
      dE1 <- numeric(n)
      dI1R <- numeric(n)
      dI1M <- numeric(n)
      dHR <- numeric(n)
      dHM <- numeric(n)
      dCR <- numeric(n)
      dCD <- numeric(n)
      dR <- numeric(n)
      dD <- numeric(n)
      dU <- numeric(n)
      dHC <- numeric(n)
      dCC <- numeric(n)
      dE2 <- numeric(n)
      dE3 <- numeric(n)
      dI2R <-numeric(n)
      dI3R <-numeric(n)
      dI2M <-numeric(n)
      dI3M <-numeric(n)
      
      
      
      with(as.list(p),
           {
             # contact reduction
             if (U[1]<Res2DK) {beta <- betaNow} else {beta <- betaNew}
             
             # interactive adjust
             RR <- ER*0.01 

             for(i in 1:n)
             {
               dS[i]   <- - S[i] * beta[i,] %*% (I1R+I1M+I2R+I2M+I3R+I3M) * RR
               dE1[i]  <-   S[i] * beta[i,] %*% (I1R+I1M+I2R+I2M+I3R+I3M) * RR - 3*gammaEI1[i] * E1[i]
               dI1R[i] <- 3*gammaEI1[i] * E3[i]*pI1R[i]       - 3*recI1[i] * I1R[i]
               dI1M[i] <- 3*gammaEI1[i] * E3[i]*(1-pI1R[i])   - 3*gammaI12[i] * I1M[i]
               dHR[i]  <- 3*gammaI12[i] * I3M[i]*pI2R[i]      + recI3[i]*CR[i] - recI2[i]* HR[i]
               dHM[i]  <- 3*gammaI12[i] * I3M[i]*(1-pI2R[i])  - gammaI23[i] * HM[i]
               dCR[i]  <- gammaI23[i] * HM[i]*pI3R[i]         - recI3[i] * CR[i]
               dCD[i]  <- gammaI23[i] * HM[i]*(1-pI3R[i])     - muI3[i] * CD[i]
               dR[i]   <- 3*recI1[i]*I3R[i] + recI2[i]*HR[i] 
               dD[i]   <- muI3[i]*CD[i] 
               dU[i]   <- 1
               dHC[i]  <- 3*gammaI12[i] * I3M[i]
               dCC[i]  <- gammaI23[i] * HM[i]
               dE2[i]  <- 3*gammaEI1[i]*E1[i] - 3*gammaEI1[i]*E2[i]
               dE3[i]  <- 3*gammaEI1[i]*E2[i] - 3*gammaEI1[i]*E3[i]
               dI2R[i] <- 3*recI1[i]*I1R[i] - 3*recI1[i]*I2R[i]
               dI3R[i] <- 3*recI1[i]*I2R[i] - 3*recI1[i]*I3R[i]
               dI2M[i] <- 3*gammaI12[i] * I1M[i] - 3*gammaI12[i] * I2M[i]
               dI3M[i] <- 3*gammaI12[i] * I2M[i] - 3*gammaI12[i] * I3M[i]
             }
             dX <- c(dS,dE1,dI1R,dI1M,dHR,dHM,dCR,dCD,dR,dD,dU,dHC,dCC,
                     dE2,dE3,dI2R,dI3R,dI2M,dI3M)
             #if(abs(sum(dX)-2)>1e-8) stop("Non-conservation!")
             return(list(dX))
           }
      )
      
    }
    
    withProgress(message = 'Simulerer -', value = 0, {
      
      nrep <- as.numeric(input$nrep)
      
      for (i in 1:nrep)
      {
        ## Two daily interacting groups 
        n <- 2
        betaNow <- array(c(crunif(input$r.betaWIlt),crunif(input$r.betaB),
                        crunif(input$r.betaB),crunif(input$r.betaWIge)),c(2,2))
        beta1 <- array(c(crunif(input$r.betaWIlt1),crunif(input$r.betaB1),
                           crunif(input$r.betaB1),crunif(input$r.betaWIge1)),c(2,2))
        beta4 <- array(c(crunif(input$r.betaWIlt4),crunif(input$r.betaB4),
                         crunif(input$r.betaB4),crunif(input$r.betaWIge4)),c(2,2))
        beta2 <- array(c(crunif(input$r.betaWIlt2),crunif(input$r.betaB2),
                         crunif(input$r.betaB2),crunif(input$r.betaWIge2)),c(2,2))
        beta3 <- array(c(crunif(input$r.betaWIlt3),crunif(input$r.betaB3),
                         crunif(input$r.betaB3),crunif(input$r.betaWIge3)),c(2,2))
        betaNew <- switch(as.numeric(input$scenarie),betaNow,beta1,beta2,beta3,beta4)
        
        # efficiency of restrictions
        ER <- crunif(input$r.ER)
        
        # mortality rates for different stages
        muI1 <- c(0,0)
        muI2 <- c(0,0)
        muI3 <- c(1/crunif(input$r.mu31),1/crunif(input$r.mu32)) #c(0.14,0.14) 1/Dage til død
        
        # recovery rates for different stages
        recI1 <- c(1/crunif(input$r.recI11),1/crunif(input$r.recI12)) # 1/dage før rask udenfor hospital
        recI2 <- c(1/crunif(input$r.recI21),1/crunif(input$r.recI22)) # 1/dage på hospital før udskrivelse
        recI3 <- c(1/crunif(input$r.recI31),1/crunif(input$r.recI32)) # 1/dage på ICU før udskrivelse
        
        # rates going from one state to the next
        gammaEI1 <- c(1/crunif(input$r.kE1),1/crunif(input$r.kE2)) # rate going from E to I1
        gammaI12 <- c(1/crunif(input$r.k1),1/crunif(input$r.k2)) # rate going from I1 to I2
        gammaI23 <- c(1/crunif(input$r.gI231),1/crunif(input$r.gI232)) # rate going from I2 to I3 - move to ICU
        
        # percentage Recover or Moving on
        pI1R <- 1-c(crunif(input$r.pI1R1),crunif(input$r.pI1R2))*0.01
        pI2R <- c(crunif(input$r.pI2R1),crunif(input$r.pI2R2)) #1-c(.05,.5)
        pI3R <- c(crunif(input$r.pI3R1),crunif(input$r.pI3R2))
        
        ## Initial states
        Inf0 <- c(crunif(input$r.InfInit1),crunif(input$r.InfInit2))
        
        
        S0 <- (pop-Inf0-c(11,8))/Ntot
        HR0 <- c(11,8)/Ntot
        HM0 <- 0*S0
        CR0 <- 0*S0
        CD0 <- 0*S0
        R0 <- 0*S0
        D0 <- 0*S0
        U0 <- c(ResDK,ResDK)
        HC0 <- c(29,30)/Ntot
        CC0 <- c(0,3)/Ntot
        
        E10  <-c(10,19)/24*Inf0/Ntot
        E20  <-c(9,1)/24*Inf0/Ntot
        E30  <-c(0,0)/24*Inf0/Ntot
        
        I1R0 <-c(1,2)/24*Inf0/Ntot*pI1R
        I1M0 <-c(1,2)/24*Inf0/Ntot*(1-pI1R)
        
        I2R0 <-c(2,1)/24*Inf0/Ntot*pI1R
        I2M0 <-c(2,1)/24*Inf0/Ntot*(1-pI1R)
        
        I3R0 <-c(2,1)/24*Inf0/Ntot*pI1R
        I3M0 <-c(2,1)/24*Inf0/Ntot*(1-pI1R)
        
        init0 <- c(S=S0,E1=E10,I1R=I1R0,I1M=I1M0,HR=HR0,HM=HM0,
                   CR=CR0,CD=CD0,R=R0,D=D0,U=U0,HC=HC0,CC=CC0,
                   E2=E20, 
                   E3=E30 ,
                   I2R=I2R0,
                   I3R=I3R0,
                   I2M=I2M0,
                   I3M=I3M0)
        
        p <- c(betaNow=betaNow,betaNew=betaNew,gammaEI1=gammaEI1,gammaI12=gammaI12,gammaI23=gammaI23, 
               muI1=muI1, muI2=muI2, muI3=muI3,
               recI1=recI1, recI2=recI2, recI3=recI3,
               pI1R=pI1R,pI2R=pI2R,pI3R=pI3R,ER=ER) # 15 variables
        
        sol <- ode(init0,times,SSEIRD,p)
        
        if (i==1)
        {
          InfTot <- matrix(NA,nrow=NROW(sol),ncol = nrep)
          InfHos <- matrix(NA,nrow=NROW(sol),ncol = nrep)
          InfInt <- matrix(NA,nrow=NROW(sol),ncol = nrep)
        }
        
        InfTot[,i] <- rowSums(sol[,c(4:17,28:39)])*Ntot
        InfHos[,i] <- rowSums(cbind(sol[,seq(10,12,2)]*Ntot,sol[,seq(11,13,2)]*Ntot))
        InfInt[,i] <- rowSums(cbind(sol[,seq(14,16,2)]*Ntot,sol[,seq(15,17,2)]*Ntot))
        
        if (i==1)
        {
          Sens <- c(p,Inf01=Inf0[1],Inf02=Inf0[2],maxInfInt=max(InfInt[,i]),wmaxInfInt=which.max(InfInt[,i]))
        } else {
          Sens <- rbind(Sens,c(p,Inf01=Inf0[1],Inf02=Inf0[2],maxInfInt=max(InfInt[,i]),wmaxInfInt=which.max(InfInt[,i])))
        }
        
        incProgress(1/nrep, detail = "Vent venligst")
      }
      
    })
 
    resid <- (InfInt[6:nDaysData,] - matrix(rep(IndlagteInt[6:nDaysData],nrep),ncol=nrep))^2
    weights <- colSums(resid,na.rm=TRUE)
    weights <- (1-0.9*(weights-min(weights))/(max(weights)-min(weights)))
    
    list(InfTot=InfTot,InfHos=InfHos,InfInt=InfInt,Sens=Sens,weights=weights)
    
  }) # dataBS 
 
    #### plot1BS ####
    output$plot1BS <- renderPlot({
      
      solBS <- dataBS()$InfTot
      weights <- dataBS()$weights
      weights[weights < quantile(weights,input$BSsel)] <- 0
      
      BSsolQuan <-  apply(solBS, 1, function(x) 
        wtd.quantile(x, probs = c(seq(0,.9,0.1),.95), weights = weights ) )
      
      par(mfrow=c(1,1))
      par(mar=c(5,2,4,2)+.1,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))
      
      plot(NULL, xlim=c(times[1],times[length(times)]), ylim=c(0,ceiling(max(max(BSsolQuan))) ), ylab=" ", xlab=" ",xaxt="n", frame.plot = F)
      axis(1,at=xt,labels=timelabels)
      
      for (i in 11:2) { 
        polygon(c(times, times[length(times)], rev(times)), c(BSsolQuan[i,], 0, rev(BSsolQuan[i-1,])),col = adjustcolor(colorTable[12-i],alpha.f=0.75), border = NA)
      }
      lines(times, BSsolQuan[6,],type="l",lwd=3)
      legend("topright",c( "Risiko:","50%","0-10 %","10-20 %","20-30 %","30-40 %","40-50 %","50-60 %","60-70 %","70-80 %","80-90 %","90-100 %"),
             col=c( NA,"black",colorTable ),lty = c(NA,1,rep(1,11)),lwd=3)
      if (testmode) mtext("Test Resultat! Data mangler!",side=3,col=adjustcolor("red",alpha.f = 0.75),padj=5,cex=2)
      
      
    })
    
    #### plot2BS ####
    output$plot2BS <- renderPlot({
      
      
      solBS <- dataBS()$InfHos
      weights <- dataBS()$weights
      weights[weights < quantile(weights,input$BSsel)] <- 0
      
      BSsolQuan <-  apply(solBS, 1, function(x) 
        wtd.quantile(x, probs = c(seq(0,.9,0.1),.95), weights = weights ) )
      
      par(mfrow=c(1,1))
      par(mar=c(5,2,4,2)+.1,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))
      
      plot(NULL, xlim=c(times[1],times[length(times)]), ylim=c(0,ceiling(max(max(BSsolQuan))) ), ylab=" ", xlab=" ",xaxt="n", frame.plot = F)
      axis(1,at=xt,labels=timelabels)
      
      for (i in 11:2) { 
        polygon(c(times, times[length(times)], rev(times)), c(BSsolQuan[i,], 0, rev(BSsolQuan[i-1,])),col = adjustcolor(colorTable[12-i],alpha.f=0.75), border = NA)
      }
      lines(times, BSsolQuan[6,],type="l",lwd=3)
      abline(h=input$Hoscap,col=adjustcolor("black",alpha.f = 0.75),lwd=3,lty=2)
      legend("topright",c("Kapacitet", "Risiko:","50%","5-10 %","10-20 %","20-30 %","30-40 %","40-50 %","50-60 %","60-70 %","70-80 %","80-90 %","90-100 %"),
             col=c(adjustcolor("black",alpha.f = 0.75), NA,"black",colorTable ),
             lty = c(2,NA,1,rep(1,11)),lwd=3)
      if (testmode) mtext("Test Resultat! Data mangler!",side=3,col=adjustcolor("red",alpha.f = 0.75),padj=5,cex=2)
      
    })
    
    #### plot3BS ####
    output$plot3BS <- renderPlot({
      
      solBS <- dataBS()$InfInt
      weights <- dataBS()$weights
      weights[weights < quantile(weights,input$BSsel)] <- 0
      
      BSsolQuan <-  apply(solBS, 1, function(x) 
        wtd.quantile(x, probs = c(seq(0,.9,0.1),.95), weights = weights ) )
      #quantile(x, probs = c(.05,seq(0.1,1,0.1))) )
      
      par(mfrow=c(1,1))
      par(mar=c(5,2,4,2)+.1,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))
      
      plot(NULL, xlim=c(times[1],times[length(times)]), ylim=c(0,ceiling(max(max(BSsolQuan))) ), ylab=" ", xlab=" ",xaxt="n", frame.plot = F)
      axis(1,at=xt,labels=timelabels)
      
      for (i in 11:2) { 
        polygon(c(times, times[length(times)], rev(times)), c(BSsolQuan[i,], 0, rev(BSsolQuan[i-1,])),col = adjustcolor(colorTable[12-i],alpha.f=0.75), border = NA)
      }
      lines(times, BSsolQuan[6,],type="l",lwd=3)
      abline(h=input$ICUcap,col=adjustcolor("black",alpha.f = 0.5),lwd=3,lty=2)
      legend("topright",c("Kapacitet", "Risiko:","50%","5-10 %","10-20 %","20-30 %","30-40 %","40-50 %","50-60 %","60-70 %","70-80 %","80-90 %","90-100 %"),
             col=c(adjustcolor("black",alpha.f = 0.75), NA,"black",colorTable ),
             lty = c(2,NA,1,rep(1,11)),lwd=3) 
      if (testmode) mtext("Test Resultat! Data mangler!",side=3,col=adjustcolor("red",alpha.f = 0.75),padj=5,cex=2)
    })
    
    #### plot Start fasen ####
    output$plotInit <- renderPlot({
      
      par(mfrow=c(1,2))
      # Hos
      solBS <- dataBS()$InfHos
      times2 <- times[InitPer]
      weights <- dataBS()$weights
      weights[weights < quantile(weights,input$BSsel)] <- 0
      
      BSsolQuan <-  apply(solBS[InitPer,], 1, function(x) 
        wtd.quantile(x, probs =c(seq(0,.9,0.1),.95), weights = weights ) )

      par(mar=c(5,2,4,2)+.1,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))
      
      plot(NULL, xlim=c(times2[1],times2[length(times2)]), ylim=c(0,ceiling(max(max(BSsolQuan))) ), ylab=" ", xlab=" ",xaxt="n", frame.plot = F)
      axis(1,at=ResDK+InitPer-1,labels=itl)
      
      for (i in 11:2) { 
        polygon(c(times2, times2[length(times2)], rev(times2)), c(BSsolQuan[i,], 0, rev(BSsolQuan[i-1,])),col = adjustcolor(colorTable[12-i],alpha.f=0.75), border = NA)
      }
      lines(times2, BSsolQuan[6,],type="l",lwd=3)
      abline(h=input$Hoscap,col=adjustcolor("black",alpha.f = 0.75),lwd=3,lty=2)
      legend("topleft",c("Kapacitet", "Risiko:","50%","5-10 %","10-20 %","20-30 %","30-40 %","40-50 %","50-60 %","60-70 %","70-80 %","80-90 %","90-100 %"),
             col=c(adjustcolor("black",alpha.f = 0.75), NA,"black",colorTable ),
             lty = c(2,NA,1,rep(1,11)),lwd=3)
      points(0:(length(IndlagteHos)-1)+ResDK,IndlagteHos,col=4,pch=19)
      title("indlagte Hospital ikke medregnet intensiv")
      
      # ICU
      solBS <- dataBS()$InfInt
      
      BSsolQuan <-  apply(solBS[InitPer,], 1, function(x)
        wtd.quantile(x, probs = c(seq(0,.9,0.1),.95), weights = weights ) )
      
      
      par(mar=c(5,2,4,2)+.1,mai=c(1.02,0.82,0.82,2.02) + c(-0.6,0,-0.5,-1.5))
      
      plot(NULL, xlim=c(times2[1],times2[length(times2)]), ylim=c(0,ceiling(max(max(BSsolQuan))) ), ylab=" ", xlab=" ",xaxt="n", frame.plot = F)
      axis(1,at=ResDK+InitPer-1,labels=itl)
      
      for (i in 11:2) { 
        polygon(c(times2, times2[length(times2)], rev(times2)), c(BSsolQuan[i,], 0, rev(BSsolQuan[i-1,])),col = adjustcolor(colorTable[12-i],alpha.f=0.75), border = NA)
      }
      lines(times2, BSsolQuan[6,],type="l",lwd=3)
      abline(h=input$ICUcap,col=adjustcolor("black",alpha.f = 0.5),lwd=3,lty=2)
      legend("topleft",c("Kapacitet", "Risiko:","50%","5-10 %","10-20 %","20-30 %","30-40 %","40-50 %","50-60 %","60-70 %","70-80 %","80-90 %","90-100 %"),
             col=c(adjustcolor("black",alpha.f = 0.75), NA,"black",colorTable ),
             lty = c(2,NA,1,rep(1,11)),lwd=3) 
      points(0:(length(IndlagteInt)-1)+ResDK,IndlagteInt,col=4,pch=19)
      title("indlagte Intensiv")
      
      par(mfrow=c(1,1))
    })
    
    #### information tab ####
    output$teori_text <- renderText({
      'Antagelser og yderligere informationer om arbejdet forud for første genåbning den 14. April 2020
      er tilgængeligt fra SSIs hjemmeside: <a href="https://files.ssi.dk/Ekspertrapport%20Matematisk%20modellering%20COVID19%20smittespredning%20og%20sygehusbelastning%201604202"> fuldt notat </a><br><br>
      Hvis man ønsker at køre flere bootstraps er man velkommen til at hente <a href="http://people.compute.dtu.dk/laec/C19REST030.zip" >kildekoden</a>. <br><br>
      <B>Kort beskrivelse af modellen:</B><br>
      Modellen er baseret på en udvidet aldersstruktureret SSEIR model, hvor rater kan reduceres i givne tidsrum.
    Oprindeligt udviklet i simpel form i starten af forrige århundrede (Se eks. William Ogilvy Kermack, A. G. McKendrick 
    and Gilbert Thomas Walker 1997A contribution to the mathematical theory of
    epidemics. Proc. R. Soc. Lond. A115700–721.).
    En let forståelig gennemgang kan ses på Wikipedia 
    (https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology).
    Modellen er udvidet så der  anvendes en dynamisk komponent
    således at man kan ændrer på effektiviteten af restriktioner i et givent tidsrum fra 11 Marts. 
    Desuden modelleres
    antal indlagte på hospitaler. Aldersstrukturering betyder at alle parametre der beskriver sygdommen
    er afhængig af alder. 
    Bemærk at disse dynamiske ændringer 
    medfører at en del teoretiske størrelser fra simplere modeller ikke bare lige kan overføres.
    Den dynamiske simulering bør således være mere præcis end et teoretisk resultat.
      '
    })
    
    output$flow_text <- renderText({
      'Her ses et flow diagram af underliggende model, <B>S</B> beskriver modtagelige individer, <B>E</B> er eksponerede individer, 
    <B>I<sub>R</sub></B> er inficerede, som vil blive raske mens <B>I<sub>H</sub></B> er inficerede, som vil blive indlagt på 
    hospitalet, <B>H<sub>R</sub></B> og <B>H<sub>C</sub></B> dækker over individer på hospitalet, hvor R og C angiver, om de vil blive raske
    eller udvikle en kritisk tilstand, på samme måde bekriver <B>C<sub>R</sub></B> og <B>C<sub>D</sub></B> om individer i kritisk vil overleve (recover)
    eller dø henholdsvis, til slut er <B>R</B> raske (immune) og <B>D</B> er døde. <br> <br>

    <B>S</B> svarer til den modtagelige del af populationen, som kan blive smittet gennem kontakt med inficerede individer. Inficerede individer (<B>I</B>)
    bliver inddelt i to kategorier, inficerede med milde symptomer, som ikke kræver hospitalsindlæggelse og derfor efterfølgende bliver raske 
    (<B>I<sub>R</sub></B>) og inficerede som kræver hospitalsindlæggelse efterfølgende (<B>I<sub>H</sub></B>). De hospitalsindlagte er også delt ind 
    i to, dem som bliver raske (<B>H<sub>R</sub></B>) samt dem, hvor symptomerne udvikler sig til en kritisk tilstand (<B>H<sub>C</sub></B>). 
    Individer i kritisk tilstand (<B>C</B>) vil enten få forbedret tilstanden og gå fra kritisk tilstand  (<B>C<sub>R</sub></B>) til hospitalsindlæggelse med 
    udsigt til at blive rask (<B>H<sub>R</sub></B>), eller 
    være i kritisk tilstand (<B>C<sub>D</sub></B>) med udgang død (<B>D</B>). Hvis et individ er rask og immun betegnes dette med <B>R</B> og <B>D</B> betegner 
    individer afgået ved døden.'
    })
    
    #### Sens Analysis ####
    output$plotSensP <- renderTable({
      
      Sens <- dataBS()$Sens[,-c(15:18)]
      prcc <- epi.prcc(Sens[,-ncol(Sens)])
      prcc <- cbind(prcc,colnames(Sens)[1:(ncol(Sens)-2)])
      prcc <- prcc[order(abs(prcc$gamma),decreasing = TRUE),]
      names(prcc)[5] <- "param"
      prcc

    })
    

})