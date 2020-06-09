## Gruppering og konvertering af kontaktrater
#setwd(<location of file>)
library(data.table)
# Population data is downloaded from: mortality.org To have comparable data for UK and DK
popDK <- fread("Pop_DK_2017.txt", skip=2)
popDK
popUK <- fread("Pop_UK_2017.txt", skip=2)


library(openxlsx)
# Credit for this data from the BBC Pandemic:
# https://doi.org/10.1101/2020.02.16.20023754
CTypes <- getSheetNames("BBC_matrices_reciprocal_filled.xlsx")[1:8]
tmp <- read.xlsx("BBC_matrices_reciprocal_filled.xlsx", rowNames = TRUE, sheet = CTypes[1])

listCT <- list()
for (type in CTypes){
  listCT[[type]] <- read.xlsx("BBC_matrices_reciprocal_filled.xlsx", rowNames = TRUE, sheet = type)
}

rn <- rownames(tmp)
lower.breaks <- as.numeric(sapply(strsplit(rn,"[-+]"), function(x)x[[1]]))
popDK[,Age:=as.numeric(substr(Age,1,3))]
popUK[,Age:=as.numeric(substr(Age,1,3))]

pop <- rbind(popDK[,Country:="DK"], popUK[,Country:="UK"])

pop[,AG1:=lower.breaks[rowSums(outer(Age, lower.breaks,">="))]]

popCont <- pop[,.(Total=sum(Total)), keyby=.(Country,AG1)][,Popu:=sum(Total),by=Country]
popCont[,AG1p:=Total/Popu]

propUK <- popCont[Country=="UK",AG1p]
propDK <- popCont[Country=="DK",AG1p]

######################### ##
merge.contacts <- function(groups, lower.breaks, prop, m, make.symmetric=TRUE){ 
  # Groups: liste med vektorer for hver ønsket gruppe
  # lower.breaks: Vektor med nedre ende af alle intervaller - det sidste er åbent til uendelig alder ;-)
  # prop: fordeling af population på aldersgrupper
  # m: antal kontakter per individ per dag (afsender i kolonner og modtager i rækker)
  # make.symmertric: Hvis TRUE laves T symmetrisk inden der grupperes og transformeres tilbage
  
  # Først omregnes til antal kontakter i original skala:
  cont <- as.matrix(m)*rep(prop, each=length(prop))
  if (make.symmetric)  
    cont <- (cont + t(cont))/2
  
  ## Nye gruppenavne  
  rn <- sapply(groups,function(x)lower.breaks[x[1]])  
  dimNames <- paste0(rn, "-",c(rn[-1]-1,""))
  
  ng <- length(groups)
  gc <- matrix(NA, ng, ng)
  gp <- numeric(ng)
  for (i in 1:ng){
    gp[i] <- sum(prop[groups[[i]]])
    for (j in 1:ng){
      gc[i,j] <- sum(cont[groups[[i]], groups[[j]]])
    }
  }
  mgc <- gc/rep(gp, each=length(gp))
  dimnames(mgc) <- list(dimNames, dimNames)
  names(gp) <- dimNames
  return(list(m=mgc, prop=gp)) # Kontakter og andel af population i hver gruppe
}


## Version 2 med vægte på polymod ####
# use polymod data for 0-4x5-0
polymod_all <- matrix(c(1.92, 0.65,0.95,6.64), nrow = 2, byrow = TRUE) *0.25   ## Ganget med 0.25 fordi PolyMod i 10-14 og 15-19 er på samme niveau som 5-9
polymod_phys <- matrix(c(1.49, 0.59, 0.74, 3.82), nrow = 2, byrow = TRUE) *0.25
polymod_con <-  polymod_all-polymod_phys

## Først samtaler: Der normeres, så 5-9 årige har samme antal samtaler som 10-19 årige
#ch
(ch <- merge.contacts(groups = list(1,2,3:4,5:6),lower.breaks = c(0,5,10,13,15,18), prop = propUK[1:6], m = listCT[[1]][1:6,1:6], make.symmetric = TRUE))
listCT[[1]][1:2,1:2] <- polymod_con/polymod_con[2,2]*mean(diag(ch$m)[3:4])
#cs
(cs <- merge.contacts(groups = list(1,2,3:4,5:6),lower.breaks = c(0,5,10,13,15,18), prop = propUK[1:6], m = listCT[[3]][1:6,1:6], make.symmetric = TRUE))
listCT[[3]][1:2,1:2] <- polymod_con/polymod_con[2,2]*mean(diag(cs$m)[3:4])

## Dernæst de fysisk kontakter
#ph
(ph <- merge.contacts(groups = list(1,2,3:4,5:6),lower.breaks = c(0,5,10,13,15,18), prop = propUK[1:6], m = listCT[[5]][1:6,1:6], make.symmetric = TRUE))
listCT[[5]][1:2,1:2] <- polymod_phys/polymod_phys[2,2]*mean(diag(ph$m)[3:4])
#ps
(ps <- merge.contacts(groups = list(1,2,3:4,5:6),lower.breaks = c(0,5,10,13,15,18), prop = propUK[1:6], m = listCT[[7]][1:6,1:6], make.symmetric = TRUE))
listCT[[7]][1:2,1:2] <- polymod_phys/polymod_phys[2,2]*mean(diag(ps$m)[3:4])



########################### ##
## Kigger på aldersmønstre i kontakttyper
# groups <- list(1,2, 3:4,5:6 ,7:9, 10:11, 12:13, 14:15, 16:17, 18:19) # Angiv grupper baseret på lower breaks
# Split hvert fem år
groups5yr <- list(1,2,3:4,5:6 ,7:8,9, 10,11, 12,13, 14,15, 16,17, 18,19) # Angiv grupper baseret på lower breaks

listTypes <- list()
for (i in names(listCT)){
  listTypes[[i]] <- merge.contacts(groups = groups5yr, lower.breaks = lower.breaks, prop = propDK, m = listCT[[i]] )
}
ageMin <- c((0:15)*5)

## Plot kontakt matricer ####
# outfolder <- "~/proj/Corona/kontaktrater/"
# png(file=paste0(outfolder,"BBC_kontaktmatricer_5yr_0.2conv.png"), width=2000, height=1000, pointsize = 30)
# ageMin <- (0:15)*5
par(mfrow=c(2,4), mar=c(3,3,2,1), mgp=c(2, 0.7,0))
zmax <- max(sapply(listTypes,function(x)max(x$m))*rep(c(0.2,1), each=4))
for (i in 1:length(listTypes)){
  image(ageMin,ageMin,listTypes[[i]]$m*ifelse(i<=4,0.2,1), main=names(listTypes)[i], zlim=c(0,zmax))
}
# dev.off()



## Plan: Skalere til symmetri og derefter reducere  og re-skalere ###

## Counts for all og samlet i fire ####
#cont <- as.matrix(m)*rep(prop, each=length(prop))
nAgeGr <- length(listTypes[[1]]$prop)
listCounts<- list()
for (i in names(listTypes)){
  listCounts[[i]] <- listTypes[[i]]$m*rep(listTypes[[i]]$prop, each=nAgeGr)
}

weightConversational <- 0.2
names(listCounts)
counts4<- list()
counts4$home <- listCounts$conversational_home * weightConversational + listCounts$physical_home
counts4$work <- listCounts$conversational_work * weightConversational + listCounts$physical_work
counts4$school <- listCounts$conversational_school * weightConversational + listCounts$physical_school
counts4$other <- listCounts$conversational_other * weightConversational + listCounts$physical_other
prop <- listTypes[[1]]$prop


############################## #
#### End of initializing ##

## Generate symmetric weight matrix based on diagonal. 
## Ensures that you can add so that when the diagonal is one then the rest is as well
##   it is done by adding contacts to and from a given age group to those that are younger (In Danish this could be labelled as "sildebensparket")
genWeight <- function(dw){
  n <- length(dw)
  if (n==1) return(dw)
  w <- diag(dw)
  for (i in 1:n){
    w[i,(1:i)] <- w[(1:i),i] <- dw[i]
  }
  return(w)
}

############################## #
## Structure for activities ####
# Further details can be found in the reports.
# All activites are handled in 5 year age groups

#Template: listAkt[[]] <- list(AktName = "", home = , work = ,school = , other=)
listAkt <- list()
arbStyrk <- 3348000 # Opdateret 16/5, var 5180000

################################## #
# Tiden fra d. 11/3 til 14/4:
listAkt[["ned"]] <- list(AktName = "Nedlukning d. 11/3", home = 0.45, work = genWeight( rep(c(0.45, 1), c(13,3)) ) ,
                         school = 0, other=genWeight( rep(c(0.2,0.1), c(14,2)) ) )

################################## #
## Tiltag i fase 1: fra d. 15/4 ####
# 5800 ansatte ud af 17000 er mødt fysisk på gym+EUD -- KILDE?
listAkt[["UngUdd.f1"]] <- list(AktName = "UngUdd til eksamen", 
                               home = 0, 
                               work = genWeight( rep(c(0,5800/arbStyrk,0),c(5,8,3)) ), # Kilde: BUVM
                               school = genWeight(rep(c(0,0.16,0.04,0),c(3, 1,9,3))) , # 80% af årgang
                               other=genWeight(rep(c(0,0.1,0),c(3, 1,12))))  

# FM forventer at 13500 går i arbejde. Kunder håndteres som 0.025 i other
listAkt[["LibErv"]] <- list( AktName = "Liberale erhverv", 
                             home = 0, 
                             work = genWeight( rep(c(13500/arbStyrk,0),c(13,3)) ), # Kilde: FM
                             school = 0, 
                             other = 0.025)

listAkt[["Praksis"]] <- list( AktName = "Praksis sektoren", 
                              home = 0,
                              work = genWeight( rep(c(33000*0.8/arbStyrk,0),c(13,3)) ), 
                              school = 0, 
                              other = 0.025)

# Domstole 5000 personer + 7900 på udvalgte videregående uddannelser jf UFM. (7900 tilføjet 2020-05-12)
listAkt[["Domstole"]] <- list(AktName = "Domstole og videregående udd.", 
                              home = 0, 
                              work = genWeight( rep(c((5000+7900)/arbStyrk,0),c(13,3)) ),
                              school = 0, 
                              other=0)

listAkt[["PrivArb.f1"]] <- list(AktName = "Privat arbejdsmarked efter påske", 
                                home =0 ,
                                work = genWeight( rep(c(0.12,0),c(13,3)) ) ,  # Trafik Google 60% i arbejde OK med +0.15 Kilde: FM
                                school = 0, 
                                other=0)

# Alle børn til og med 11 år. +0.05 i other for 10-14-årige
# 71341 lærer på arbejde (ud af 102502) 
# Antager, at det giver +0.05 home kontakter
listAkt[["DagInstSkole.f1"]] <- list(AktName = "Daginstitutioner og skole 0.-5. kl.", 
                                     home = 0.05, 
                                     work = genWeight( rep(c(0,71341/arbStyrk,0),c(5,8,3)) ),
                                     school =  genWeight(c(1,1,0.4,0,rep(0.6,9),0,0,0)), 
                                     other=genWeight(rep(c(0,0.05,0),c(2,1,13))))

################################## #
## Tiltag, som undersøges i fase 2 ####
# Normalt 450000 nu 150000 kunder per dag.
# 60000 (Kilde: FM)
listAkt[["Detail"]] <- list(AktName = "Detailhandel: Storcentre og stormagasiner", 
                            home = 0, 
                            work = genWeight( rep(c(60000/arbStyrk,0),c(13,3)) ),
                            school =0 , 
                            other= 0.05)

# Professionel idræt (3000), biblioteker (kun lån) (BIB7: 3500 ansatte på folkebiblioteker), offentlig sektor AT mf 5000 ÅV i alt
# Zoo 2000 ansatte og 4500000 årlige besøgene ~ 3 timer per dansker per år ~ 0.5 min/dag = 0.002 i other (Dette blev brugt d. 6/5, men i genberegningen var det kun kørende)
# Opdateret 2020-05-12 til kun adgang med biler så 340 ansatte (KM) og ingen i other:
listAkt[["Anden.f2"]] <- list(AktName = "Prof idræt + biblioteker + 5000 offentligt ansatte",
                              home = 0, 
                              work = genWeight( rep(c((3000+3500+5000+340)/arbStyrk,0),c(13,3)) ) ,
                              school = 0, 
                              other=0)
# Skoler
# 9.-10. er ca 1,5 årgang, hvoraf 0,5 årgang er på efterskole
#  lærer på arbejde, må være resten: 102502-71341 = 31161
# Home + 0.05
listAkt[["6.-10.kl"]] <- list(AktName = "6.-10. klasse samt klubtilbud", 
                              home = 0.05, 
                              work = genWeight( rep(c(0,31161/arbStyrk,0),c(5,8,3)) ), # OBS
                              school = genWeight(c(0,0,0.6,0.2,rep(0.212,9),0,0,0)), 
                              other= genWeight(rep(c(0,0.1,0.03,0), c(2,1,1,12))) )
# 1/4 af lærerne ift 6.-10. klasse
listAkt[["9.-10.kl"]] <- list(AktName = "9.-10. klasse samt klubtilbud", 
                              home = 0.05, 
                              work = genWeight( rep(c(0,31161/4/arbStyrk,0),c(5,8,3)) ), 
                              school = genWeight(c(0,0,0,0.2,rep(0.05,9),0,0,0)), 
                              other=genWeight(rep(c(0,0.03,0), c(3,1,12))) )

listAkt[["Efterskoler"]] <- list(AktName = "Efterskoler - skoledel", 
                                 home = 0, 
                                 work = genWeight(rep(c(5000/arbStyrk,0),c(13,3))),
                                 school = genWeight(c(0,0,0,0.1,rep(0.028,9),0,0,0)), 
                                 other=genWeight(rep(c(0,0.03,0), c(3,1,12))) )

# Restaurant og cafe
# Ref. HORESTA publikation "danskernes udespisevaner i 2019"
c(109, 136, 101, 98, 62,33)* (c(0.27, 0.26, 0.26, 0.31, 0.4, 0.55)+0.1)  # Antal * <andel cafe+restaurant> for [12-19, 20-29, 30-39, 40-49, 50-64, 65-]
# UPDATE 20200512: +0.1 svarer til værtshuse + cocktailbarer + lidt fastfood.
tmp <- c(10, 20, 30, 40,  49, 49, 36,36, 40,40, 31,31,31 , 21,21,21) # Antaget for de unge, som ikke fremgår af publikation fra HORESTA.
# Gennemsnit 1/2 besøg per uge ~ 30min per uge ~ 5min per dag svarende til 2% af tiden i 'other' *2 pga højere antal kontakter/tid
tmp2 <- tmp/mean(tmp)*0.04*0.7  #0.7 pga 70% aktivitet
# Ved ~70% aktivitet
# 40.000 forventes at gå i arbejde Kilde:EM
listAkt[["RestCafe"]] <- list(AktName = "Restauranter og caféer", 
                              home = 0, 
                              work = genWeight( rep(c(0,40000/arbStyrk,0),c(3,10,3)) ),
                              school = 0, 
                              other=genWeight(tmp2))
## Udeservering er 10% af fuld omsætning (Kilde: FM)
listAkt[["RestCafeUde"]] <- list(AktName = "Restauranter og caféer - udeservering", 
                                 home = 0, 
                                 work = genWeight( rep(c(0,40000/7/arbStyrk,0),c(3,10,3)) ),
                                 school = 0, 
                                 other=genWeight(tmp2/7))

## For det private arbejdsmarked Kilde: FM 
# 2020-05-17: Justeret til 0.05 Kilde: FM 
listAkt[["PrivArb.f2"]] <- list(AktName = "Privat arbejdesmarked efter anden åbning", 
                                home =0 , 
                                work = genWeight( rep(c(0.05,0),c(13,3)) ),
                                school = 0, 
                                other=0)

listAkt[["divUdd"]] <- list(AktName = "STU + EUD + FGU + skolehjem", 
                            home = 0, 
                            work = 0,
                            school = genWeight(rep(c(0,0.15,0.04,0),c(3, 1,9,3))),  
                            other=genWeight(rep(c(0,0.1,0),c(3, 1,12))))

listAkt[["fysEksam"]] <- list(AktName = "Videregående udd med fysisk eksamen", 
                              home = 0, 
                              work = genWeight( rep(c(0,13400/765000,0),c(4,2,10)) ),  #765000 personer i 20-29-årige
                              school = 0,
                            other=0)

# 650.000 personer som begynder udendørsaktiviteter. Vi antager 15 min per dag svarende til 5% af 'other'
listAkt[["Idræt"]] <- list(AktName = "Idræt - udendørs", home = 0, work = 0,school = 0, other = 0.05)

listAkt[["home.225"]] <- list(AktName = "Home + 22.5%", home = 0.225, work =0 ,school = 0, other=0)
listAkt[["home.25"]] <- list(AktName = "Home + 25%", home = 0.25, work =0 ,school = 0, other=0)
listAkt[["home.275"]] <- list(AktName = "Home + 27.5%", home = 0.275, work =0 ,school = 0, other=0)


###################################### #
## Fase 3 ####
# Se notat på SSIs hjemmeside for dokumentation, hvor der ikke er nævnt detaljer
# Sommeraktiviter håndteres for sig selv
# 35.8 mio besøgende per år ~ 3min i gennemsnit per dag ~ 0.01 i other, men da biograf, teater og zoo har højere tæthed benyttes 0.02
listAkt[["Kultur"]] <- list(AktName = "Museer, teater, biograf, film og TV, zoo", 
                            home = 0, 
                            work = genWeight( rep(c(0,(6660+6400+2000+2500+960)/arbStyrk,0),c(3,10,3)) ),
                            school = 0, 
                            other = 0.02)

listAkt[["OffSag"]] <- list(AktName = "Offentlig sektor med sagspukler", 
                            home = 0, 
                            work = genWeight( rep(c(0,(33000)/arbStyrk,0),c(3,10,3)) ), # Kilde: FM
                            school = 0, 
                            other = 0) # Ca 90.000 borgere til møder, men gennemsnitligt tidsforbrug er meget lavt

listAkt[["OffForsk"]] <- list(AktName = "Offentlig forskning som kræver tilstædeværelse", 
                              home = 0, 
                              work = genWeight( rep(c(0,(14500+300)/arbStyrk,0),c(3,10,3)) ), # Kilde: UFM+KEFM+KUM+VIVE
                              school = 0, 
                              other = 0)

listAkt[["VoksUdd"]] <- list(AktName = "Voksenuddannelse for ledige", 
                    home = 0, 
                    work = genWeight( rep(c(0,(21800+1500+3720+117)/arbStyrk,0),c(3,10,3)) ),
                    school = 0, 
                    other = 0)

listAkt[["Hojskoler"]] <- list(AktName = "Højskoler mv", 
                    home = 0, 
                    work = genWeight( rep(c(0,(1500)/arbStyrk,0),c(3,10,3)) ) + genWeight( rep(c(0,(5000)/765000,0),c(4,2,10)) )+
                      + genWeight( rep(c(0,(1000)/640000,0),c(13,2,1)) ), # ansatte + unge på lange + ældre(65-75-årige) på korte
                    school = 0, 
                    other = 0.005) 

# Svarer til 1,37 min per danske per dag ~ 0.005 i 'other' med anslået firdobbelt risiko per tid: 0.02
listAkt[["IndIdræt"]] <- list(AktName = "Indendørs idræt og foreningsliv", 
                    home = 0, 
                    work = genWeight( rep(c(0,(2000+7000+3500)/arbStyrk,0),c(3,10,3)) ),
                    school = 0, 
                    other = 0.02)

################################# #
## Fase 3+: Udvidet 3. fase ####

listAkt[["DRTV2"]] <- list(AktName = "DR og TV2 åbnes fuldt", 
                    home = 0, 
                    work = genWeight( rep(c(0,(4695)/arbStyrk,0),c(3,10,3)) ),
                    school = 0, 
                    other = 0)

listAkt[["OffAnsat.f3+"]] <- list(AktName = "Offentlige ansatte i fase 3+", 
                    home = 0.05, # Normalisering af samfundet
                    work = genWeight( rep(c(0,(190000)/arbStyrk,0),c(3,10,3)) ), # Kilde: FM 
                    school = 0 , 
                    other = 0.05) # Kombination af 5% mere i arbejde, som giver transport og borgerkontakt

################################# #
## Fase 4:  ####

listAkt[["AllUdd"]] <- list(AktName = "Alle øvrige uddannelser inkl 1.g og 2.g", 
                    home = 0.1, 
                    work = genWeight( rep(c(0,(320000)/arbStyrk,0),c(3,10,3)) ),
                    school = genWeight(rep(c(0,0.39,0.08,0),c(3, 1,9,3))),  #OBS anslået to årgange 15-19-årige
                    other = 0.1)

# Kilde: Dansk Svømmebadsteknisk Forening. 
# Hver dansker bruger 4min/dag ~ 0.015 pp i other. Der anslås firdobbelt risiko, så 0.06
listAkt[["Fitness"]] <- list(AktName = "Fitness + svømmebade + legeland", 
                    home = 0, 
                    work = genWeight( rep(c(0,(1900)/arbStyrk,0),c(3,10,3)) ),
                    school = 0, 
                    other = 0.06)



funSumAkt <- function(included, activityList, faWork){
  # included: activities to include (names of elements in activityList)
  # activityList: list of activities. Must include elements named 'home', 'work', 'school' and 'other'
  # faWork: Physical distance factor at work. 'work' gives how many are at work and faWork defines how effective control measures at work are.
  if(any(!included %in% names(activityList))){
    warning(paste(included[!included %in% names(activityList)],"not in activityList"))
  }
  relCont <- list(home =0 , work =0 ,school =0 , other=0)
  for (akt in included){
    for (nm in names(relCont)){
      fysAfstandFaktor <- ifelse(nm=="work", faWork, 1)
      relCont[[nm]] <- relCont[[nm]] + activityList[[akt]][[nm]] * fysAfstandFaktor
    }
  }
  relCont[["activities"]] <- included
  return(relCont)
}

dummy <- c("ned", "LibErv")
funSumAkt(dummy, listAkt, faWork = 0.45)

names(listAkt)
scenarier <- list(basis=list(1,1,1,1,"Basis"))
scenarier[["ned"]] <- funSumAkt("ned", listAkt, faWork=0.45)
sc.fase1 <- c("ned","UngUdd.f1", "LibErv","Praksis", "Domstole", "PrivArb.f1", "DagInstSkole.f1")
scenarier[["fase1"]] <- funSumAkt(sc.fase1, listAkt, faWork=0.45)
# The following lines relates to the report on May 6th 2020
# scenarier[["fase2.01"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2"), listAkt, faWork=0.5)
# scenarier[["fase2.02"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "6.-10.kl"), listAkt, faWork=0.5)
# scenarier[["fase2.03"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "6.-10.kl", "Efterskoler"), listAkt, faWork=0.5)
# scenarier[["fase2.04"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "6.-10.kl", "RestCafeUde"), listAkt, faWork=0.5)
# scenarier[["fase2.05"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "6.-10.kl", "RestCafe"), listAkt, faWork=0.5)
# scenarier[["fase2.06"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "6.-10.kl", "PrivArb.f2"), listAkt, faWork=0.55)
# scenarier[["fase2.07"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "6.-10.kl", "Efterskoler", "RestCafeUde","PrivArb.f2"), listAkt, faWork=0.55)
# scenarier[["fase2.08"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "6.-10.kl", "Efterskoler", "RestCafe","PrivArb.f2"), listAkt, faWork=0.55)
# ## 9.-10. kl i stedet for 6.10. kl.
# scenarier[["fase2.09"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "9.-10.kl"), listAkt, faWork=0.45)
# scenarier[["fase2.10"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "9.-10.kl", "Efterskoler", "RestCafeUde","PrivArb.f2"), listAkt, faWork=0.55)
# scenarier[["fase2.11"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "9.-10.kl", "Efterskoler", "RestCafe","PrivArb.f2"), listAkt, faWork=0.55)

## Følsomheds scenarier
#"fase2.01" og "fase2.08" home skal 1/2 og 1/1 til 100%
## Først: Hvad har vi egentlig:
# sapply(scenarier, function(x)sapply(x[1:4],max))
# ## Skal bruge home+0.225
# scenarier[["fase2.01f50"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "home.25"), listAkt, faWork=0.5+0.25)
# scenarier[["fase2.01f100"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "home.25", "home.25"), listAkt, faWork=1)
# 
# scenarier[["fase2.08f50"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "6.-10.kl", "Efterskoler", "RestCafe","PrivArb.f2","home.225"), listAkt, faWork=0.55+.225)
# scenarier[["fase2.08f100"]] <- funSumAkt(c(sc.fase1, "Detail", "Anden.f2", "6.-10.kl", "Efterskoler", "RestCafe","PrivArb.f2","home.225","home.225"), listAkt, faWork=1)
# 
# scenarier[["fase1f50"]] <- funSumAkt(c(sc.fase1, "home.275"), listAkt, faWork=0.45)

## Hvad blev vedtaget
sc.fase2 <- c(sc.fase1, "Detail", "Anden.f2", "6.-10.kl", "Efterskoler", "RestCafe","PrivArb.f2", "divUdd","Idræt","fysEksam")
scenarier[["fase2"]] <- funSumAkt(sc.fase2, listAkt, faWork=0.55)
#scenarier[["fase2.1.1f50"]] <- funSumAkt(c(sc.fase2,"home.225"), listAkt, faWork=0.55)
#scenarier[["fase2.1.1f100"]] <- funSumAkt(c(sc.fase2,"home.225","home.225"), listAkt, faWork=0.55)
add.fase3 <- c("Kultur", "OffSag", "OffForsk", "VoksUdd", "Hojskoler", "IndIdræt")
scenarier[["fase3"]] <- funSumAkt(c(sc.fase2, add.fase3 ), listAkt, faWork=0.60) # + 0.05 faWork

add.fase3p <- c("DRTV2", "OffAnsat.f3+")
scenarier[["fase3+"]] <- funSumAkt(c(sc.fase2, add.fase3, add.fase3p ), listAkt, faWork=0.65) # + 0.05 faWork

add.fase4 <- c("AllUdd", "Fitness" )
scenarier[["fase4"]] <- funSumAkt(c(sc.fase2, add.fase3, add.fase3p, add.fase4 ), listAkt, faWork=0.65) # + 0.05 faWork

# Følsomhedsscenarier
tmp <- lapply(scenarier[["fase4"]][1:4],function(xx)diag(as.matrix(xx)))
full <- lapply(tmp, function(x)1-x)
half <- lapply(tmp, function(x)(1-x)*0.5)

addScenarios <- function(x,y){
  # x is a scenario
  # y are diagonal elements (vectors) to be added
  z <- x
  for (i in 1:4){
    z[[i]] <- z[[i]] + genWeight( y[[i]])
  }
  return(z)
}

scenarier[["fase2f100"]] <- addScenarios( scenarier[["fase2"]], full)
scenarier[["fase2f50"]] <- addScenarios( scenarier[["fase2"]], half)
scenarier[["fase3f100"]] <- addScenarios( scenarier[["fase3"]], full)
scenarier[["fase3f50"]] <- addScenarios( scenarier[["fase3"]], half)
scenarier[["fase3+f100"]] <- addScenarios( scenarier[["fase3+"]], full)
scenarier[["fase3+f50"]] <- addScenarios( scenarier[["fase3+"]], half)
scenarier[["fase4f100"]] <- addScenarios( scenarier[["fase4"]], full)
scenarier[["fase4f50"]] <- addScenarios( scenarier[["fase4"]], half)

makeSummer <- function(x, wSchool=0.5, wWork=0.75){
  # x is a scenario
  # wSchool, wWork weights for those two
    z <- x
    z$school <- z$school * wSchool
    z$work <- z$work * wWork
}

scenarier[["fase2S"]] <- makeSummer(scenarier[["fase2"]])
scenarier[["fase2f100S"]] <- makeSummer(scenarier[["fase2f100"]])
scenarier[["fase2f50S"]] <- makeSummer(scenarier[["fase2f50"]])

scenarier[["fase3S"]] <- makeSummer(scenarier[["fase3"]])
scenarier[["fase3f100S"]] <- makeSummer(scenarier[["fase3f100"]])
scenarier[["fase3f50S"]] <- makeSummer(scenarier[["fase3f50"]])

scenarier[["fase3+S"]] <- makeSummer(scenarier[["fase3+"]])
scenarier[["fase3+f100S"]] <- makeSummer(scenarier[["fase3+f100"]])
scenarier[["fase3+f50S"]] <- makeSummer(scenarier[["fase3+f50"]])

scenarier[["fase4S"]] <- makeSummer(scenarier[["fase4"]])
scenarier[["fase4f100S"]] <- makeSummer(scenarier[["fase4f100"]])
scenarier[["fase4f50S"]] <- makeSummer(scenarier[["fase4f50"]])




## Først: Hvad har vi egentlig:
sapply(scenarier, function(x)sapply(x[1:4],max))

lapply(scenarier, function(x)sapply(x[1:4],function(xx)round(diag(as.matrix(xx)),3)))

## Eigenvalues
max.eigen.4w <- function(vaekt=rep(1,4), listW, prop){
  m2 <- counts4[[1]]*listW[[1]]/rep(prop,each=length(prop))*vaekt[1]
  for (i in 2:4){
    m2 <- m2 + counts4[[i]]*listW[[i]]/rep(prop,each=length(prop)) * vaekt[i]
  }
  
  return(list(me = eigen(m2)$values, eig=eigen(m2)$vec))
}


## Eigenvalues
(tmp <- max.eigen.4w(vaekt = rep(1, each=4), listW = scenarier[["fase4"]],prop=prop)$eig[,1])

full <-max.eigen.4w(vaekt = rep(1, each=4), listW = list(1,1,1,1),prop=prop)$me[1]

(max.eigen.4w(vaekt = rep(1, each=4), listW = scenarier[["ned"]],prop=prop)$me[1]/full*3.3)

## Beregning af R0 for alle scenarier:
sapply(scenarier, function(x){max.eigen.4w(vaekt = rep(1, each=4), listW = x,prop=prop)$me[1]/full*3.3})



######################### ####
merge.counts2contacts <- function(groups, lower.breaks, prop, cont, make.symmetric=TRUE){ 
  # Groups: liste med vektorer for hver ønsket gruppe
  # lower.breaks: Vektor med nedre ende af alle intervaller - det sidste er åbent til uendelig alder ;-)
  # prop: fordeling af population på aldersgrupper
  # cont: antal kontakter (evt skaleret) mellem aldersgrupper antages at være symmetrisk (afsender i kolonner og modtager i rækker)
  # make.symmertric: Skal T laves symmetrisk inden der grupperes og transformeres tilbage?
  
  ## Nye gruppenavne  
  rn <- sapply(groups,function(x)lower.breaks[x[1]])  
  dimNames <- paste0(rn, "-",c(rn[-1]-1,""))
  
  ng <- length(groups)
  gc <- matrix(NA, ng, ng)
  gp <- numeric(ng)
  for (i in 1:ng){
    gp[i] <- sum(prop[groups[[i]]])
    for (j in 1:ng){
      gc[i,j] <- sum(cont[groups[[i]], groups[[j]]])
    }
  }
  mgc <- gc/rep(gp, each=length(gp))
  dimnames(mgc) <- list(dimNames, dimNames)
  names(gp) <- dimNames
  return(list(m=mgc, prop=gp)) # Kontakter og andel af population i hver gruppe
}


groupsKG <- list(1:12 , 13:16) # Angiv grupper baseret på lower breaks

listEffCont <- list()
for(red in names(scenarier)){
  EffCont <- matrix(0, nrow=length(groupsKG), ncol=length(groupsKG))
  for (i in 1:4){
    tmp <- merge.counts2contacts(groups = groupsKG, lower.breaks = ageMin, prop = prop, cont = counts4[[i]]*scenarier[[red]][[i]] )
    EffCont <- EffCont + tmp$m
  }
  listEffCont[[red]] <- EffCont / tmp$prop*0.15*0.8
}
## Faktor er justeret til samme egenværdi som tidligere svarende til 20% vækst per dag
beta <- listEffCont

save(beta, scenarier, file="beta_til_Kaare_20200517.RData")


################################# #
## Plot til teksnisk appendiks ####

library(viridis)

png(filename = "Counts_per_type.png", width=1600, height=1600, pointsize = 30)
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(1.8,0.7,0))
zmax <- 1
for (i in 1:4){
  tmp<-counts4[[i]]/rep(prop, each=length(prop))
  print(range(tmp))
  image(ageMin,ageMin,tmp , main=names(counts4)[i], zlim=c(0,zmax), xlab="Fra aldersgruppe", ylab=" Til aldersgruppe", col = viridis(20,direction = -1))
  box()
}
dev.off()

png(filename = "Contacts_fase.png", width=1600, height=1600, pointsize = 30)
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(1.8,0.7,0))
zmax <- 0.1
sce <- c("basis","ned","fase1", "fase3")
sceNames <- c("Normale kontakter", "Nedlukket", "Efter første genåbning", "Efter udvidet anden genåbning")
for (red in 1:4){
  EffCont <- matrix(0, nrow=length(prop), ncol=length(prop))
  for (i in 1:4){
    EffCont <- EffCont + counts4[[i]]*scenarier[[ sce[red] ]][[i]]
  }
  print(range(EffCont))
  image(ageMin,ageMin,EffCont , main=sceNames[red], zlim=c(0,zmax), xlab="Fra aldersgruppe", ylab=" Til aldersgruppe", col = viridis(20,direction = -1))
  box()
 
}
dev.off()
