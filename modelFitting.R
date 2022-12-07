library(tidyverse)
library(lme4)
library(lme4qtl)

run_models <- function(models = c('bivar', 'slope_add_12', 'slope_add_21'), dat, G, ldak_loc){
  #Convenience function that fits models and wrap up the results
  stopifnot(packageVersion("lme4qtl") >= "0.2.2")
  
  lapply(models, function(model) {
    tab <- switch(model,
                  "bivar" = {
                    bdat <- gather(dat, tname, tvalue, trait1, trait2)
                    start <- Sys.time()
                    mod <- relmatLmer(tvalue ~ -1 + tname + (0 + tname | id) + (0 + tname | rid), 
                                      data = bdat, relmat = list(id = G), calc.derivs = FALSE)
                    end <- Sys.time()
                    time <- end - start
                    
                    #Estimated variance components
                    vf <- VarProp(mod) %>% filter(grp == "id")
                    sigma.t1 <- vf$vcov[vf$var1 == 'tnametrait1' & is.na(vf$var2)]
                    sigma.t2 <- vf$vcov[vf$var1 == 'tnametrait2' & is.na(vf$var2)]
                    sigma.t1t2 <- vf$vcov[which(vf$var1 == 'tnametrait1' & vf$var2 == 'tnametrait2')]
                    
                    #Estimated h2s and rhog
                    rhog <- sigma.t1t2 / sqrt(sigma.t1 * sigma.t2)
                    h2.add.t1 <- vf$prop[vf$var1 == 'tnametrait1' & is.na(vf$var2)]
                    h2.add.t2 <- vf$prop[vf$var1 == 'tnametrait2' & is.na(vf$var2)]
                    
                    data_frame(rhog = rhog, h2.add.t1 = h2.add.t1, h2.add.t2 = h2.add.t2, runtime = as.numeric(time, unit = 'mins'))
                  },
                  "slopeG_add_12" = {
                    start <- Sys.time()
                    mod <- relmatLmer(trait2 ~ (1|rid) + trait1 + (0 + trait1|id), data = dat, 
                                      relmat = list(id = G, rid = G))
                    end <- Sys.time()
                    time <- end - start
                    
                    #Estimates of h2
                    vf <- VarProp(mod)
                    h2.add.t2 <- vf$prop[vf$grp == 'rid']
                    h2.cov_12 <- vf$prop[vf$grp == 'id']
                    
                    data_frame(h2.cov_12 = h2.cov_12, h2.add.t2 = h2.add.t2, runtime = as.numeric(time, unit = 'mins'))
                  },
                  "slopeG_add_21" = {
                    start <- Sys.time()
                    mod <- relmatLmer(trait1 ~ (1|rid) + trait2 + (0 + trait2|id), data = dat, 
                                      relmat = list(id = G, rid = G))
                    end <- Sys.time()
                    time <- end - start
                    
                    #Estimates of h2
                    vf <- VarProp(mod)
                    h2.add.t1 <- vf$prop[vf$grp == 'rid']
                    h2.cov_21 <- vf$prop[vf$grp == 'id']
                    
                    data_frame(h2.cov_21 = h2.cov_21, h2.add.t1 = h2.add.t1, runtime = as.numeric(time, unit = 'mins'))
                  },
                  "slopeGandE_add_12" = {
                    dat$eid <- dat$id #Grouping variable for the iid random slopes

                    start <- Sys.time()
                    mod <- relmatLmer(trait2 ~ (1|rid) + trait1 + (0 + trait1|id) + (0 + trait1|eid), data = dat, 
                                      relmat = list(id = G, rid = G))
                    end <- Sys.time()
                    time <- end - start
                    
                    #Estimates of h2
                    vf <- VarProp(mod)
                    h2.add.t2 <- vf$prop[vf$grp == 'rid']
                    h2.cov_12 <- vf$prop[vf$grp == 'id']
                    e.cov_12 <- vf$prop[vf$grp == 'eid']
                    
                    data_frame(h2.cov_12 = h2.cov_12, h2.add.t2 = h2.add.t2, e.cov_12 = e.cov_12, runtime = as.numeric(time, unit = 'mins'))
                  },
                  "slopeE_add_12" = {
                    start <- Sys.time()
                    mod <- relmatLmer(trait2 ~ (1|rid) + trait1 + (0 + trait1|id), data = dat)
                    end <- Sys.time()
                    time <- end - start
                    
                    #Estimates of h2
                    vf <- VarProp(mod)
                    icc.t2 <- vf$prop[vf$grp == 'rid']
                    e.cov_12 <- vf$prop[vf$grp == 'id']

                    data_frame(icc.t2 = icc.t2, e.cov_12 = e.cov_12, runtime = as.numeric(time, unit = 'mins'))
                  },
                  "gxemm_iid" = {
                    if(length(unique(dat$rep)) > 1) #If repeated observations, use just one per individual
                      dat.gxemm <- dat[!duplicated(dat$id),]
                    else
                      dat.gxemm <- dat
                    
                    start <- Sys.time()
                    mod		<- GxEMM(y = dat.gxemm$trait2, 
                                  K = G, 
                                  X = dat.gxemm$trait1,
                                  Z = dat.gxemm$trait1, gtype='iid', ldak_loc=ldak_loc)
                    end <- Sys.time()
                    time <- end - start
                    
                    data.frame(h2.hom_gxemm = mod$h2[1], h2.het_gxemm = mod$h2[2], runtime = as.numeric(time, unit = 'mins'))
                  },
                  "gxemm_giid_eiid" = {
                    if(length(unique(dat$rep)) > 1) #If repeated observations, use just one per individual
                      dat.gxemm <- dat[!duplicated(dat$id),]
                    else
                      dat.gxemm <- dat
                    
                    start <- Sys.time()
                    mod		<- GxEMM(y = dat.gxemm$trait2, 
                                  K = G, 
                                  X = dat.gxemm$trait1,
                                  Z = dat.gxemm$trait1, gtype='iid', etype = 'iid', ldak_loc=ldak_loc)
                    end <- Sys.time()
                    time <- end - start
                    
                    data.frame(h2.hom_gxemm = mod$h2[1], h2.het_gxemm = mod$h2[2], runtime = as.numeric(time, unit = 'mins'))
                  },
                  stop("error in switch: unknown `model`"))
    
    mutate(tab, model = model)
  }) %>% bind_rows
}
