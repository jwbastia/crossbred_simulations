## Simulation of Xbred pedigree with purebred and crossbred performance
## 
setwd("/home//basti015/JohnB/current/projects/B4F/rg_accuracy")
rm(list=ls())
require(pedigree)

##
## DESIGN

h2pb <- 0.40
h2xb <- h2pb

ngenpb <- 10
sexratio <- 8  # number of females per male
littersize <- 10 # number of offspring per mating

# population size for purebred lines
      Ne  <- c( 300 , 300 , 300 , 300)
names(Ne) <- c("A","B","C","D")

damxb  <- TRUE    # TRUE : 2 dam lines will be simulated and crossbred dams will be simulated to produce the terminal product
                  # FALSE: 1 dam line will be simulated and pureline dams will be simulated to produce the terminal product

sirexb <- FALSE   # TRUE : 2 sire lines will be simulated and crossbred sires will be simulated to produce the terminal product
                  # FALSE: 1 sire line will be simulated and pureline sires will be simulated to produce the terminal product
# crossbred parents are only made in one direction, i.e. AxB but NOT BxA
# terminal product is made from A or AB "males" mated to C or CD "females" (sex assignment is arbitrary, so interpretation could be reversed)

# number of animals per generation
  # pb phenotyped animals are from line A
  # xb phenotyped animals' line composition depends on damxb and sirexb
  # (F,F)->AC, (F,T)->A(CD), (T,F)->(AB)C, (T,T)->(AB)(CD)

      N  <- c( NA,   100,  NA,   160,  NA,   NA,    NA,    NA,   NA  )  # Am and Cf must be provided. Others can be calculated
                                                                            #Cf should be at least Am * (sexratio / littersize) * 2
names(N) <- c( "Af", "Am", "Bf", "Cf", "Dm", "ABm", "CDf", "pb", "xb")  

# line A has males (Am) and females (Af) simulated because purebred line A animals are needed for phenos
# Bf, Cf, and CDf are females.  ABm, and Dm are males. pb, and xb are both males and females. 
# line letter refer to position in the (AB)*(CD) mating structure

##################################################
# For "damxb == sirexb == TRUE" the numbers of   #
# animals in the mating structure are as follows #  
##################################################
#                                           #
#   gen -1    NeA     NeB     NeC    NeD    # 
#             / \      |       |      |     #
#   gen 0   Af * Am * B(f)    C(f) * D(m)   #
#              |    |              |        #
#   gen 1      pb   AB(m)    *    CD(f)     # 
#                            |              #
#   gen 2                    xb             #
#                                           #
###########################################################################################
# For sirexb == F : Am (gen=0) is mated to CD(gen=1) (or  C(gen=0) if  damxb = F) and NeB, B and AB do not exist #
# for  damxb == F :  C (gen=0) is mated to AB(gen=1) (or Am(gen=0) if sirexb = F) and NeD, D and CD do not exist #
###########################################################################################


# if numbers are not specified they are calculated based on sexratio and littersize
# calculate number of dams "Af", if not provided
if(is.na(N[["Af"]])) N[["Af"]] <- N[["Am"]]*sexratio
# calculate number of line B females for production of AB males, if required and not provided
if(!sirexb) N[["Bf"]] <- 0 else if(is.na(N[["Bf"]])) N[["Bf"]] <- N[["Am"]]*sexratio
# calculate number of line D males for production of CD females, if required and not provided
if(!damxb)  N[["Dm"]] <- 0 else if(is.na(N[["Dm"]])) N[["Dm"]] <- ceiling(N[["Cf"]]/sexratio)
# calculate number of AB animals (all males), if required and not provided
if(!sirexb) N[["ABm"]] <- 0 else if(is.na(N[["ABm"]])) N[["ABm"]] <- N[["Bf"]]*littersize/2 # /2 because AB are only male offspr.
# calculate number of CD animals (all females), if required and not provided
if(!damxb)  N[["CDf"]] <- 0 else if(is.na(N[["CDf"]])) N[["CDf"]] <- N[["Cf"]]*littersize/2 # /2 because CD are only female offspr.
# calculate number of pb offspring, if not provided (= always required)
if(is.na(N[["pb"]])) N[["pb"]] <- N[["Af"]]*littersize
# calculate number of xb offspring, if not provided (= always required)
if(is.na(N[["xb"]])) N[["xb"]] <- (  damxb  * N[["CDf"]] + 
                                   (!damxb) * N[["Cf"]]  ) * littersize


## correlations
            # Gpb    Gxb    Epb    Exb
R <- matrix(c(1.00 , 0.70 , 0.00 , 0.00 ,  # Gpb
              NA   , 1.00 , 0.00 , 0.00 ,  # Gxb
              NA   , NA   , 1.00 , 0.00 ,  # Epb
              NA   , NA   , NA   , 1.00 ), # Exb
            nrow = 4, byrow = TRUE)
cholR <- chol(R)

#collect <- data.frame(rG = NA, rP = NA, vGpb = NA , vGxb= NA, vPpb = NA, vPxb = NA)


#i <- 1
for(i in c(1))   
{  
  pedigree <- data.frame(line = NA , gen = NA ,
                         id   = 0  , sid = NA , did = NA , sex = NA ,
                         Gpb  = NA , Gxb = NA , Epb = NA , Exb = NA )
  pednames <- names(pedigree)
  
  ## purebred pedigree
  for(pp in names(Ne[c(TRUE,sirexb,TRUE,damxb)]))    # tst #  pp <- names(Ne)[1]
  {
    #founders
    # random deviates for G and E
    rd <- matrix( c(rnorm(Ne[[pp]], sd = sqrt(    h2pb) ) , 
                    rnorm(Ne[[pp]], sd = sqrt(    h2xb) ) ,
                    rnorm(Ne[[pp]], sd = sqrt(1 - h2pb) ) ,
                    rnorm(Ne[[pp]], sd = sqrt(1 - h2xb) ) ) , 
                  byrow = FALSE, ncol = 4 )
    # add correlations of G and E
    rd <- rd %*% cholR 
    
    # append founders to ped
    ped <- data.frame(line = rep(pp, Ne[[pp]])                 ,
                      gen  = -ngenpb                           ,
                      id   = seq(max(pedigree$id) + 1 , 
                                 max(pedigree$id) + Ne[[pp]] ) ,
                      sid  = NA                                , 
                      did  = NA                                ,
                      sex  = sample(c("M","F")      , 
                                    size = Ne[[pp]] , 
                                    replace = T)              
                      )

    ped <- cbind(ped,rd)  
    names(ped) <- pednames
      
    # later generations
    for(gg in -(ngenpb-1):-1)
    {  ## tst :  gg <- -9
       curr_sires <- ped[sample(which(ped$gen == gg-1 & 
                                       ped$sex =="M"), 
                               size = Ne[[pp]] , 
                               replace = T    
                               ), 
                        match(c("id","Gpb","Gxb"), pednames)]
      
      curr_dams  <- ped[sample(which(ped$gen == gg-1 & 
                                       ped$sex =="F"), 
                               size = Ne[[pp]] , 
                               replace = T
                               ), 
                        match(c("id","Gpb","Gxb"), pednames)]
      
      curr_gen   <- data.frame(line = rep(pp, Ne[[pp]])            ,
                               gen  = gg                           ,
                               id   = seq(max(ped$id) + 1          , 
                                          max(ped$id) + Ne[[pp]] ) ,
                               sid  = curr_sires$id                , 
                               did  = curr_dams$id                 ,
                               sex  = sample(c("M","F")      ,
                                             size = Ne[[pp]] ,
                                             replace = T     )     )
      
      rd <- matrix( c(rnorm(Ne[[pp]], sd = sqrt(    h2pb/2) ) , 
                      rnorm(Ne[[pp]], sd = sqrt(    h2xb/2) ) ,
                      rnorm(Ne[[pp]], sd = sqrt(1 - h2pb)   ) ,
                      rnorm(Ne[[pp]], sd = sqrt(1 - h2xb)   ) ) , 
                    byrow = FALSE                               , 
                    ncol = 4                                    )
      rd <- rd %*% cholR
      curr_gen$Gpb <- 0.5*curr_sires$Gpb + 0.5*curr_dams$Gpb + rd[,1]
      curr_gen$Gxb <- 0.5*curr_sires$Gxb + 0.5*curr_dams$Gxb + rd[,2]
      curr_gen <- cbind(curr_gen,rd[,3:4])
      names(curr_gen) <- pednames
      ped <- rbind(ped,curr_gen)        
    } # end later generations 
    
    # add purebred line to pedigree
    pedigree <- rbind(pedigree, ped)

  } # end of purebred lines 

  # remove placeholder row 1 of pedigree
  pedigree <- pedigree[-1,]
  
  ## GENERATION 0
  ##  
  gen0 <- names(N)[1:5]
  gen0 <- gen0[N[1:5]>0]
  for(pp in gen0 ) # pp <- "Af"
  {
    
    curr_line <- substr(pp,1,1)
    curr_sex  <- substr(pp,2,2)
    curr_sires <- pedigree[sample(which(pedigree$gen  == -1        &
                                        pedigree$line == curr_line &
                                        pedigree$sex  == "M"       ) , 
                                  size = N[[pp]]                     , 
                                  replace = T                        ) ,
                           match(c("id","Gpb","Gxb"), pednames)]
    
    curr_dams  <- pedigree[sample(which(pedigree$gen  == -1        & 
                                        pedigree$line == curr_line &
                                        pedigree$sex  == "F"       ) , 
                                  size = N[[pp]]                     , 
                                  replace = T                        ) ,
                           match(c("id","Gpb","Gxb"), pednames)]
    
    ped  <- data.frame(line = rep(pp, N[[pp]])                ,
                       gen  = 0                               ,
                       id   = seq(max(pedigree$id) + 1 , 
                                  max(pedigree$id) + N[[pp]]) ,
                       sid  = curr_sires$id                   , 
                       did  = curr_dams$id                    ,
                       sex  = toupper(curr_sex)               )
    
    rd <- matrix( c(rnorm(N[[pp]], sd = sqrt(    h2pb/2) ) , 
                    rnorm(N[[pp]], sd = sqrt(    h2xb/2) ) ,
                    rnorm(N[[pp]], sd = sqrt(1 - h2pb)   ) ,
                    rnorm(N[[pp]], sd = sqrt(1 - h2xb)   ) ) ,
                  byrow = FALSE, ncol = 4 )
    rd <- rd %*% cholR
    ped$Gpb <- 0.5*curr_sires$Gpb + 0.5*curr_dams$Gpb + rd[,1]
    ped$Gxb <- 0.5*curr_sires$Gxb + 0.5*curr_dams$Gxb + rd[,2]
    ped     <- cbind(ped,rd[,3:4])
    names(ped) <- pednames
    pedigree <- rbind(pedigree , ped)  
  }
  
  ## GENERATION 1
  ## 
  gen12 <- list(ABm = c("ABm","Am","Bf"), 
                CDf = c("CDf","Dm","Cf"),
                pb  = c("pb","Am","Af"),
                xb  = c("xb", ifelse(sirexb,"ABm","Am"),
                              ifelse(damxb ,"CDf","Cf")))
  gen12 <- gen12[c(sirexb,damxb,TRUE,TRUE)]
  
  for(pp in gen12) # pp <- gen12[[1]]
  {
    curr_sires <- pedigree[sample(which(pedigree$line == pp[2] & pedigree$sex  == "M") , 
                                  size = ceiling(N[[pp[1]]]/(sexratio*littersize/ifelse(pp[1]%in%c("pb","xb"),1,2))) , 
                                  replace = F                                          ) ,  # FALSE to get exact sex ratio 
                           match(c("id","Gpb","Gxb"), pednames)]                              # gives error with too few sires
    
    curr_dams  <- pedigree[sample(which(pedigree$line == pp[3] & pedigree$sex  == "F") , 
                                  size = ceiling(N[[pp[1]]]/(littersize/ifelse(pp[1]%in%c("pb","xb"),1,2)))            , 
                                  replace = F                                          ) ,  # FALSE to get exact sex ratio
                           match(c("id","Gpb","Gxb"), pednames)]
    
    ped  <- data.frame(line = rep(pp[1], N[[pp[1]]])                    ,
                       gen  = 1                                      ,
                       id   = seq(max(pedigree$id) + 1 , 
                                  max(pedigree$id) + N[[pp[1]]])     ,
                       sid  = curr_sires$id                          , 
                       did  = curr_dams$id                           ,
                       sex  = sample(c(ifelse(pp[1]=="ABm","M","F")  ,
                                       ifelse(pp[1]=="CDf","F","M")) ,
                                     size = N[[pp[1]]], replace = T) )
    
    rd <- matrix( c(rnorm(N[[pp[1]]], sd = sqrt(    h2pb/2) ) , 
                    rnorm(N[[pp[1]]], sd = sqrt(    h2xb/2) ) ,
                    rnorm(N[[pp[1]]], sd = sqrt(1 - h2pb)   ) ,
                    rnorm(N[[pp[1]]], sd = sqrt(1 - h2xb)   ) ) ,
                  byrow = FALSE, ncol = 4 )
    rd <- rd %*% cholR
    ped$Gpb <- 0.5*curr_sires$Gpb + 0.5*curr_dams$Gpb + rd[,1]
    ped$Gxb <- 0.5*curr_sires$Gxb + 0.5*curr_dams$Gxb + rd[,2]
    ped     <- cbind(ped,rd[,3:4])
    names(ped) <- pednames
    pedigree <- rbind(pedigree , ped)  
  }
  
  ### ADD phenotypes
  pedigree$Ppb[pedigree$line =="pb"] <-  round(pedigree$Gpb[pedigree$line =="pb"] +  pedigree$Epb[pedigree$line =="pb"],6)
  pedigree$Pxb[pedigree$line =="xb"] <-  round(pedigree$Gxb[pedigree$line =="xb"] +  pedigree$Exb[pedigree$line =="xb"],6)
  
  ### set missing parents to 0
  pedigree$sid[is.na(pedigree$sid)] <- 0
  pedigree$did[is.na(pedigree$did)] <- 0
  
  ## TST ## 
  #  with(pedigree[pedigree$line%in%c("pb","xb"),], tapply(did,sid,function(x)length(table(x))))
  print(c( cor(with(pedigree[pedigree$line=="xb",], tapply(Gxb,sid,mean)),with(pedigree[pedigree$line=="pb",], tapply(Gpb,sid,mean))) ,
           with(pedigree[pedigree$line=="Am",], cor(Gxb,Gpb))) )
    
  ## write files for ASreml
  file.create("repPbXbRg.as")
  asfile <- c("replicate of simulations for PBXB rg accuracy" ,
              " animal !P",
              " Ppb",
              " Pxb",
              "repPbXbRg.ped !MAKE !SORT !SKIP 1" ,
              "repPbXbRg.dat !SKIP 1 !ASUV" ,
              "Ppb Pxb ~ Trait !r Tr.animal !f mv",
              "1 2 1",
              "0 !S2==1",
              "Trait 0 DIAG",
              "1 1",
              "Tr.animal 2",
              "Trait 0 US",
              "1",
              "0.5 1",
              "animal"
              )
  writeLines(asfile, con="repPbXbRg.as" )
  
  file.create("repPbXbRg.pin")
  asfile <- c("F PvarPb 1+3",  ### component  6
              "F PvarXb 2+5",  ### component  7
              "H h2Pb 3 6"  ,  ### component  8
              "H h2Xb 5 7"  ,  ### component  9
              "R Gcor 3:5"  )  ### component 10
  writeLines(asfile, con="repPbXbRg.pin" )
  
  write.table(pedigree[pedigree$line %in% c("pb","xb"),match(c("id" , "Ppb" , "Pxb"), names(pedigree))], 
              file = "repPbXbRg.dat", na = "NA", sep = "\t" , row.names = FALSE , quote = FALSE )

  write.table(pedigree[,match(c("id" , "sid" , "did"), names(pedigree))], 
              file = "repPbXbRg.ped", na = "NA", sep = "\t" , row.names = FALSE , quote = FALSE )
  
  ## call ASreml
  system("asreml -nw1000 repPbXbRg.as")
  system("asreml -nP repPbXbRg.pin")
  
  ## read ASreml estimates?
  curr_estimates <- readLines("repPbXbRg.pvc")
  curr_h2Pb <- unlist(strsplit( curr_estimates[grep("h2Pb", curr_estimates)], split = " "))
  curr_h2Pb <- curr_h2Pb[curr_h2Pb != ""]
  curr_h2Pb <- as.numeric(curr_h2Pb[(length(curr_h2Pb)-1):length(curr_h2Pb)])
 
  curr_h2Xb <- unlist(strsplit( curr_estimates[grep("h2Xb", curr_estimates)], split = " "))
  curr_h2Xb <- curr_h2Xb[curr_h2Xb != ""]
  curr_h2Xb <- as.numeric(curr_h2Xb[(length(curr_h2Xb)-1):length(curr_h2Xb)])
  
  curr_Gcor <- unlist(strsplit( curr_estimates[grep("Gcor", curr_estimates)], split = " "))
  curr_Gcor <- curr_Gcor[curr_Gcor != ""]
  curr_Gcor <- as.numeric(curr_Gcor[(length(curr_Gcor)-1):length(curr_Gcor)])
  
  curr_estimates <- c(curr_h2Pb, curr_h2Xb, curr_Gcor)
  ## 

  write.table()

}  # end replicate i


# ## inspect data
# if(FALSE){
#   par(mfrow=c(2,2))
#   hist(by(pedigree,list(pedigree$line,pedigree$gen), function(x) cor(x$Gpb,x$Gxb)),
#        main = "correlation", xlab = "Gpb-Gxb")
#   hist(by(pedigree,list(pedigree$line,pedigree$gen), function(x) cor(x$Epb,x$Exb)),
#        main = "correlation", xlab = "Epb-Exb")
#   hist(by(pedigree,list(pedigree$line,pedigree$gen), function(x) cor(x$Gpb,x$Exb)),
#        main = "correlation", xlab = "Gpb-Exb")
#   hist(by(pedigree,list(pedigree$line,pedigree$gen), function(x) cor(x$Gxb,x$Epb)),
#        main = "correlation", xlab = "Gxb-Epb")
#   
#   par(mfrow=c(2,2))
#   hist(by(pedigree,list(pedigree$line,pedigree$gen), function(x) var(x$Gpb)),
#        main = "variance", xlab = "Gpb")
#   hist(by(pedigree,list(pedigree$line,pedigree$gen), function(x) var(x$Gxb)),
#        main = "variance", xlab = "Gxb")
#   hist(by(pedigree,list(pedigree$line,pedigree$gen), function(x) var(x$Epb)),
#        main = "variance", xlab = "Epb")
#   hist(by(pedigree,list(pedigree$line,pedigree$gen), function(x) var(x$Exb)),
#        main = "variance", xlab = "Exb")
# }
#   
#   
