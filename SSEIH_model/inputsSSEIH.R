# input data

popDK <- data.frame(lt60=c(1435432,590495,880891,991793,423790),
                    ge60=c(410591,246864,342214,334547,166146),
                    region=c("Hovedstaden","Zealand","Syddanmark",
                             "Midtjylland","Nordjylland")) 
regions=c("Hovedstaden","Zealand","Syddanmark",
          "Midtjylland","Nordjylland")

input <- list()

# Boot:
input$r.ER <- c(50, 150) # Skalering

input$a.ER <- outer(c(85,70,65,60,75),c(.75,1.25))

input$r.betaWIlt <- c(.12, .22) # Smitterate indenfor aldersgrupper <60
input$r.betaWIge <- c(.14, .24) # Smitterate indenfor aldersgrupper >=60
input$r.betaB <- c(.02, .12) # Smitterate mellem aldersgrupper

# in restriction periods
# Period 1: Lock down
input$r.betaWIltR1 <- c(.13, .23) # Smitterate indenfor aldersgrupper <60
input$r.betaWIgeR1 <- c(.13, .23) # Smitterate indenfor aldersgrupper >=60
input$r.betaBR1 <- c(.03, .12) # Smitterate mellem aldersgrupper

# Period 2: First reopening
input$r.betaWIltR2 <- c(.20, .30) # Smitterate indenfor aldersgrupper <60
input$r.betaWIgeR2 <- c(.15, .25) # Smitterate indenfor aldersgrupper >=60
input$r.betaBR2 <- c(.04, .14) # Smitterate mellem aldersgrupper

#diff of first two periods (The code samples from this rather than from the second period)
input$rd.betaWIltR2 <- c(0, 0.144) # Smitterate indenfor aldersgrupper <60
input$rd.betaWIgeR2 <- c(0, 0.044) # Smitterate indenfor aldersgrupper >=60
input$rd.betaBR2 <- c(0, 0.021) # Smitterate mellem aldersgrupper


# Initially infected: c("Hovedstaden","Zealand","Syddanmark","Midtjylland","Nordjylland")
input$a.InfInit1 <- outer(c(40000,15000,15000,15000,4000),c(.75,1.25))
input$a.InfInit2 <- outer(c(4000,1500,1000,800,700),c(.75,1.25))

# Under 60:
input$r.kE1 <- c(4,6) # Latensperiode
input$r.k1 <- c(6,10) # Antal dage efter smitte før indlæggelse
input$r.gI231 <- c(.5,2.5) # Antal dage til ICU efter indlæggelse
input$r.pI1R1 <- c(0.05, 0.5) # Risiko for indlæggelse [%]
input$r.pI2R1 <- c(0.77, 0.97) # Sandsynlighed for rask efter indlæggelse
input$rc.pI2R1 <- c(0.5, 1) # Sandsynlighed for rask efter indlæggelse -change 1.april
input$r.pI3R1 <- c(.7, .95) # Sandsynlighed for at overleve ICU
input$r.mu31 <-c(14,28) # Dage på ICU før død
input$r.recI11 <- c(4,6) #Antal dage før rask udenfor hospital
input$r.recI21 <- c(5,9) # Antal dage på hospital før udskrivelse
input$r.recI31 <- c(14,28) # Antal dage på ICU før udskrivelse
# Over 60:
input$r.kE2 <- c(4,6) # Latensperiode
input$r.k2 <- c(5,9) # Antal dage efter smitte før indlæggelse
input$r.gI232 <- c(0.5,1.5) # Antal dage til ICU efter indlæggelse
input$r.pI1R2 <- c(5, 6.2) # Risiko for indlæggelse [%]
input$r.pI2R2 <- c(0.7, 0.9) # Sandsynlighed for rask efter indlæggelse
input$rc.pI2R2 <- c(0.2, 1) # Sandsynlighed for rask efter indlæggelse
input$r.pI3R2 <- c(.45, .55) # Sandsynlighed for at overleve ICU
input$r.mu32 <-c(14,28) # Dage på ICU før død
input$r.recI12 <- c(4,6) # Antal dage før rask udenfor hospital
input$r.recI22 <- c(5,15) # Antal dage på hospital før udskrivelse
input$r.recI32 <- c(14,28) # Antal dage på ICU før udskrivelse

# Sampling uniform based on vector with [min ; max]
crunif <- function(invec)
{
  if (length(invec)!=2) stop("wrong size")
  if (nrep>1) {runif(1,invec[1],invec[2])} else {mean(invec)}
}