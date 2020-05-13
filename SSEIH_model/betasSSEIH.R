# make betas for use
# indtast mode værdi, der ligger automatisk pm 0.05 (bpm, værdien kan ændres i script)
# 
# betas <- list()
# scenarie <- list()
# betasPast <- list()

# #                      betaWIlt, betaWIge, betaB
# i<-1; betas[[i]]  <- c(.288    , .206    , .095  ); scenarie[[i]] <- "Example"
# 
# #                            betaWIlt,   betaWIge, betaB
# i<-1; betasPast[[i]]  <- c(0.17599952, 0.17951076, 0.07467067) # ned
# i<-2; betasPast[[i]]  <- c(0.25376615, 0.18864727, 0.08404281) # fase 1

load("beta_til_Kaare_20200506.RData")

names(beta)

scenarie <- as.list(names(beta)[-(1:2)])
betas <- lapply(beta[-(1:2)], function(x)c(x[c(1,4,2)])) 
betasPast <- lapply(beta[(2:3)], function(x)c(x[c(1,4,2)])) 
## 'betasPast' is used in this implementation - and should manually by given as inputs in 'inputsSSEIH.R'


