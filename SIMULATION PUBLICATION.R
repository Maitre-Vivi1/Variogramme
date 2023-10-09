# Library -----------------------------------------------------------------

library(MASS)
library(dplyr)
library(ggplot2)
library(nlme)
library(reshape2)
library(AICcmodavg)
library(parallel)
library(corpcor)
require(reshape2)
require(np)


# Simulation --------------------------------------------------------------

simu <- function(type_corr, nb_mesure, n_patient, temps=1:nb_mesure, effet_temps=rnorm(1)*temps) {
  
  require(corpcor)
  require(MASS)
  
  rho <- runif(1,0.5,0.9)
  vari <- runif(1, 15,30)
  pos_def <- F
  
  if (type_corr == "AR1") {
     while (!pos_def) {
      corr <-  matrix(0, nrow=nb_mesure, ncol=nb_mesure)
      for (i in temps) {
        for (j in temps) {
          
          corr[which(temps == i),which(temps == j)] <- vari * rho**(abs(i-j))
        }
      }
      pos_def <- is.positive.definite(corr)
    }
  }
  
  else if (type_corr == "AR1H") {
    corr <-  matrix(0, nrow=nb_mesure, ncol=nb_mesure)
    while (!pos_def) {
      
      sigma <- runif(n = nb_mesure, max = 1.5, min = 0.5)
      
      for (i in temps) {
        for (j in temps) {
          corr[which(temps == i),which(temps == j)] <- vari*sigma[i]*sigma[j]*rho**(abs(i-j))
        }
      }
      pos_def <- is.positive.definite(corr)
    }
  }
  
  else if (type_corr == "ARMA11") {
    corr <-  matrix(0, nrow=nb_mesure, ncol=nb_mesure)
    while (!pos_def) {
      autocorr <- runif(1)
      for (i in temps) {
        for (j in temps) {
          if (i == j) {
            corr[i,j] <- vari
          } else {
            corr[i,j] <- vari * autocorr**abs(i-j)* rho
          }
        }
      }
      pos_def <- is.positive.definite(corr)
    }
  }
  
  else if (type_corr == "UN") {
    while (!pos_def) {
      m <- matrix(runif(nb_mesure^2)*2-1, ncol=nb_mesure) 
      corr <- t(m) %*% m
      pos_def <- is.positive.definite(corr)
    }
  }
  
  else if (type_corr == "CS") {
    while (!pos_def) {
      
      corr <- matrix(rho*vari*0.3, nrow = nb_mesure, ncol = nb_mesure)
      diag(corr) <- diag(corr) + vari*runif(1)
      pos_def <- is.positive.definite(corr)
    }
  }
  
  else if (type_corr == "CSH") {
    corr <-  matrix(0, nrow=nb_mesure, ncol=nb_mesure)
    while (!pos_def) {
      
      sigma <- runif(n = nb_mesure, max = 2)
      autocorr <- c(1,rep(runif(n = 1, max = 2),nb_mesure))
      
      for (i in temps) {
        for (j in temps) {
          
          corr[which(temps == i),which(temps == j)] <- vari*sigma[i]*sigma[j]*autocorr[abs(i-j)+1]
        }
      }
      pos_def <- is.positive.definite(corr)
    }
  }
  
  else if (type_corr == "VC") {
    corr <- matrix(0, nrow = nb_mesure, ncol = nb_mesure)
    diag(corr) <- vari*runif(1)
  }
  
  else if (type_corr == "VCH") {
    while (!pos_def) {
      corr <- matrix(0, nrow = nb_mesure, ncol = nb_mesure)
      diag(corr) <- vari*runif(nb_mesure)
      pos_def <- is.positive.definite(corr)
    }
  }
  
  
  else if (type_corr == "TOEP") {  
    while (!pos_def) {
      corr <- toeplitz(x = sample(c(-1,1), nb_mesure, replace = T)* runif(n = nb_mesure)) *vari
      pos_def <- is.positive.definite(corr)
    }
  }
  
  
  else if (type_corr == "ToepH") {
    corr <-  matrix(0, nrow=nb_mesure, ncol=nb_mesure)
    while (!pos_def) {
      
      sigma <- runif(n = nb_mesure, max = 2)
      autocorr <- c(1,runif(n = nb_mesure-1, max = 2))
      
      for (i in temps) {
        for (j in temps) {
          
          corr[which(temps == i),which(temps == j)] <- sigma[i]*sigma[j]*autocorr[abs(i-j)+1]*vari
        }
      }
      pos_def <- is.positive.definite(corr)
    }
  }
  
  
  else if (type_corr == "Puis") {
    while (!pos_def) {
      corr <- matrix(0, nrow = nb_mesure, ncol = nb_mesure)
      
      for (i in temps) {
        for (j in temps) {
          
          corr[which(temps == i),which(temps == j)] <- vari*rho**(abs(temps[i]-temps[j]))
        }
      }
      pos_def <- is.positive.definite(corr)
    }
  }
  
  else if (type_corr == "Gaus") {
    while (!pos_def) {
      corr <- matrix(0, nrow = nb_mesure, ncol = nb_mesure)
      
      for (i in temps) {
        for (j in temps) {
          
          corr[which(temps == i),which(temps == j)] <-vari* exp(-abs(temps[i]-temps[j])**2/rho**2)
          
          
          
        }
      }
      pos_def <- is.positive.definite(corr)
    }
  }
  
  else if (type_corr == "GAUSH") {
    while (!pos_def) {
      corr <- matrix(0, nrow = nb_mesure, ncol = nb_mesure)
      sigma <- runif(n = nb_mesure, max = 2)
      for (i in temps) {
        for (j in temps) {
          
          corr[which(temps == i),which(temps == j)] <-vari* exp(-abs(temps[i]-temps[j])**2/rho**2)*sigma[i]*sigma[j]
        }
      }
      pos_def <- is.positive.definite(corr)
    }
  }
  
  info_struct <- list(Type = type_corr, Matrice = corr, rho = rho, vari = vari)
  
  data <- mvrnorm(n_patient, mu=effet_temps, Sigma=corr)
  
  trt <- sample(1:2, size = n_patient, replace = T)
  
  effet_trt <- rnorm(2, mean = 3, sd = 2)
  
  for (i in temps) {
    for (j in 1:2) {
      
      data[,i][which(trt == j)] <- data[,i][which(trt == j)] + effet_trt[j]
    }
  }
  
  wideData <- data.frame(id=1:n_patient, data)
  
  longData <- reshape(wideData, varying=names(wideData)[-1],
                      direction="long", sep="", idvar="id")
  longData <- longData[order(longData$id,longData$time),]
  
  longData$id <- as.factor(longData$id)
  
  names(longData)[length(names(longData))] <- "Y"
  
  info_struct <- c(info_struct, list(n_mesure = nb_mesure, n_patient = n_patient, effet_temps = effet_temps, effet_trt=effet_trt))
  
  longData$center <- rep(sample(1:3, size = n_patient, replace = T), each = nb_mesure, size = 1)
  
  for (i in levels(longData$center)) {
    longData$Y[which(longData$center == i)] <- longData$Y[which(longData$center == i)] + rnorm(1, sd = 3)
  }
  
  longData$trt <- rep(trt, each=nb_mesure, 1)
  
  vec_int <- rep(1, length(nb_mesure))
  
  info_struct$ZGZ <- matrix(nrow = length(temps), ncol = 1, data = c(vec_int)) %*% matrix(data = 3) %*% t(matrix(nrow = length(temps), ncol = 1, data = c(vec_int)))
  
  return(list(longData,info_struct))  
}


# Comparaison -------------------------------------------------------------

reussite_variog <- function(data.sim, formula = Y ~ time + trt) {
  
  require(reshape2)
  require(dplyr)
  require(np)
  require(AICcmodavg)
  require(nlme)
  
  df <- data.sim[[1]]
  info_struct <- data.sim[[2]]
  rho <- info_struct$rho
  mat <- info_struct$Matrice
  n.mesure <- info_struct$n_mesure
  
  seuil <- max(diag(mat))*0.10
  
  res.lm <- lme(data = df, fixed = formula, random = ~ 1 | center)
  
  df_res <- data.frame(id = df$id, res = res.lm$residuals[,2])
  df_dist <- NULL
  df_resid <- NULL
  
  for (i in levels(df$id)) {
    df_dist <- rbind(df_dist, melt(as.matrix(dist(df$time[which(df$id == i)]))))
    df_resid <- rbind(df_resid, melt(as.matrix(dist(df_res$res[which(df$id == i)]))))
  }
  
  df_resid$value <- (df_resid$value**2)/2
  
  df_res_temps <- data.frame(residus = df_resid$value, ecart_temp = df_dist$value, temps1 = df_dist$Var1, temps2 = df_dist$Var2)
  df_res_temps <- df_res_temps[-which(df_res_temps$ecart_temp == 0),]
  df_res_temps <- df_res_temps[-which(duplicated(df_res_temps$residus)),]
  
  moy <- df_res_temps %>%
    group_by(ecart_temp) %>%
    summarize(moy = mean(residus)) %>%
    as.data.frame()
  
  nb_repet <- sort(c(dist(x = unique(levels(as.factor(df$time))))))
  
  df4 <- data.frame(table(df_res_temps$ecart_temp), moy$moy)
  df4.2 <- NULL
  
  moy2 <- NULL
  for (i in levels(as.factor(moy$ecart_temp))) {
    
    for (j in 1:length(nb_repet[nb_repet==i])) {
      moy2 <- rbind(moy2, moy[moy$ecart_temp==i,])
      df4.2 <- rbind(df4.2, df4[df4$Var1 == i,])
    }
  }
  
  moy <- moy2
  df4 <- df4.2
  
  temps_unique <- unique(nb_repet)
  
  df_res2 <- data.frame(
    res = unname(resid(res.lm)),
    id = df$id,
    time = df$time
  )
  
  total_var_sigma2_tau2_nu2 <- function(df_res2) {
    
    dfff2 <- c()
    
    for (i in unique(df_res2$time)) {
      dfff2 <- c(dfff2, c(dist(df_res2$res[which(df_res2$time == i)], upper = F, method = "euclidean")))
    }
    
    sum(dfff2**2)/2/length(dfff2)
    
  }
  
  bw1 <- np::npregbw(xdat = moy$ecart_temp, ydat = moy$moy, regtype = "ll", ckertype="epanechnikov", bwscaling = F, bandwidth.compute = T, bwmethod = "cv.ls")
  kre1 <- np::npreg(bws = bw1, exdat =  temps_unique)
  
  kre1$mean <- c(0, kre1$mean)
  
  
  getvarcov2 <- function(model) {
    
    S <- corMatrix(model$modelStruct$corStruct)[[1]]
    
    vw <- rep(1, nrow(S))
    vars <- (model$sigma * vw)^2
    result <- t(S * sqrt(vars)) * sqrt(vars)
    class(result) <- c("marginal", "VarCov")
    attr(result, "group.levels") <- names(model$groups)
    as.matrix(result)
    
  }
  
  res.lme.toep <- lme(fixed = formula, data = df, correlation = corARMA(form = ~ 1| center/id, p = n.mesure - 1), random = ~ 1 | center, control = 
                        lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
  res.lme.ar1 <- lme(fixed = formula, data = df, correlation = corCAR1(form = ~ 1 | center/id), random = ~ 1 | center, control = 
                       lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
  res.lme.cs <- lme(fixed = formula, data = df, correlation = corCompSymm(form = ~ 1 | center/id), random = ~ 1 | center, control = 
                      lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
  res.lme.gaus <- lme(fixed = formula, data = df, correlation = corGaus(form = ~ 1 | center/id), random = ~ 1 | center, control = 
                        lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
  res.lme.vc <- lme(fixed = formula, data = df, correlation = corIdent(form = ~ 1 | center/id), random = ~ 1 | center, control = 
                      lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
  
  vec_aic <- sapply(list(res.lme.toep, res.lme.ar1, res.lme.cs, res.lme.gaus, res.lme.vc), AIC)
  vec_bic <- sapply(list(res.lme.toep, res.lme.ar1, res.lme.cs, res.lme.gaus, res.lme.vc), BIC)
  vec_aicc <- sapply(list(res.lme.toep, res.lme.ar1, res.lme.cs, res.lme.gaus, res.lme.vc), AICc)
  vec_aic_plus_bic <- vec_aic + vec_bic
  
  noms <- c("TOEP", "AR1", "CS", "Gaus", "VC")
  pred_aic <- noms[which.min(vec_aic)] ; pred_bic <- noms[which.min(vec_bic)]
  pred_aicc <- noms[which.min(vec_aicc)] ; pred_aic_plus_bic <- noms[which.min(vec_aic_plus_bic)]
  
  
  m.ar1 <- getvarcov2(res.lme.ar1)
  m.toep <- getvarcov2(res.lme.toep)
  m.cs <- getvarcov2(res.lme.cs)
  m.gaus <- getvarcov2(res.lme.gaus)
  m.vc <- diag(sigma(res.lme.vc)**2, nrow = dim(m.ar1)[1])
  
  
  m.ar1 <- matrix(data = m.ar1[1,1], ncol = n.mesure, nrow = n.mesure) - m.ar1
  m.toep <- matrix(data = m.toep[1,1], ncol = n.mesure, nrow = n.mesure) - m.toep
  m.cs <- matrix(data = m.cs[1,1], ncol = n.mesure, nrow = n.mesure) - m.cs
  m.gaus <- matrix(data = m.gaus[1,1], ncol = n.mesure, nrow = n.mesure) - m.gaus
  m.vc <- matrix(data = m.vc[1,1], ncol = n.mesure, nrow = n.mesure) - m.vc
  
  
  df_temps_loess <- data.frame(
    ecart_temp = c(0,temps_unique),
    ar1 = m.ar1[1,],
    toep = m.toep[1,],
    cs = m.cs[1,],
    vc = m.vc[1,],
    gaus = m.gaus[1,],
    var_np = kre1$mean
  )
  
  eqm <- function(fitted, np=df_temps_loess$var_np) {
    sum((fitted - np)**2)/length(fitted)
  }
  
  list_eqm <- as.list(df_temps_loess[,-c(1,7)])
  
  vec.eqm <- sapply(list_eqm, eqm)
  
  names(vec.eqm) <- c("AR1", "TOEP", "CS", "VC", "Gaus")
  pred_variog <- names(vec.eqm)[which.min(vec.eqm)]
  
  if (sqrt(vec.eqm[which(names(vec.eqm) == info_struct$Type)]) < sqrt(vec.eqm[which.min(vec.eqm)]) + seuil &
      sqrt(vec.eqm[which(names(vec.eqm) == info_struct$Type)]) > sqrt(vec.eqm[which.min(vec.eqm)]) - seuil) {
    pred_variog <- info_struct$Type
  }
  
  
  print(paste(info_struct$Type, "rho=", rho, "pat=", info_struct$n_patient, "ni=", info_struct$n_mesure, pred_variog))
  
  return(list(pred = pred_variog, pat = info_struct$n_patient, n_mesure = info_struct$n_mesure, type = info_struct$Type,
              pred_aic = pred_aic, pred_bic = pred_bic, pred_aicc = pred_aicc, pred_aic_plus_bic = pred_aic_plus_bic))
  
}

reussite_variog_total <- function(data.sim, formula.moyenne= Y ~ time + trt, seuil_iso = 0.05) {
  
  df <- data.sim[[1]]
  info_struct <- data.sim[[2]]
  rho <- info_struct$rho
  mat <- info_struct$Matrice
  n.mesure <- info_struct$n_mesure
  
  seuil <- max(diag(mat))*0.10
  
  require(ggpubr)
  require(MASS)
  require(dplyr)
  require(ggplot2)
  require(nlme)
  require(reshape2)
  require(AICcmodavg)
  require(plotly)
  
  id <- df$id
  time <- df$time
  
  temps <- sort(unique(df$time))
  
  # -------------------------------------------------------------------------
  
  
  z <- unname(lme(data = df, fixed =  formula.moyenne, random = ~ 1 | center)$residuals[,2])
  z_2 <- z**2
  
  
  df_res_time <- data.frame(
    res = z_2,
    time = time
  )
  
  pval.iso <- min(p.adjust(unname(summary(lm(data = df_res_time, formula = res ~ as.factor(time)))$coefficients[,4])[-1], method = "BH"))
  
  
  # -------------------------------------------------------------------------
  
  
  getvarcov2 <- function(model) {
    
    S <- corMatrix(model$modelStruct$corStruct)[[1]]
    
    if (!is.null(model$modelStruct$varStruct)) {
      ind <- model$groups == 1
      vw <- unique(unname(1/varWeights(model$modelStruct$varStruct)[ind]))
    } else { vw <- rep(1, nrow(S)) }
    
    vars <- (model$sigma * vw)^2
    result <- t(S * sqrt(vars)) * sqrt(vars)
    class(result) <- c("marginal", "VarCov")
    attr(result, "group.levels") <- names(model$groups)
    as.matrix(result)
    
  }
  
  
  
  if (info_struct$Type %in% c("AR1", "TOEP", "CS", "VC", "Gaus")) {
    
    
    res.lm <- lme(data = df, fixed = formula.moyenne, random = ~ 1 | center)
    
    df_res <- data.frame(id = df$id, res = res.lm$residuals[,2])
    df_dist <- NULL
    df_resid <- NULL
    
    for (i in levels(id)) {
      df_dist <- rbind(df_dist, melt(as.matrix(dist(time[which(id == i)]))))
      df_resid <- rbind(df_resid, melt(as.matrix(dist(df_res$res[which(id == i)]))))
    }
    
    
    df_resid$value <- (df_resid$value**2)/2
    
    df_res_temps <- data.frame(residus = df_resid$value, ecart_temp = df_dist$value, temps1 = df_dist$Var1, temps2 = df_dist$Var2)
    df_res_temps <- df_res_temps[-which(df_res_temps$ecart_temp == 0),]
    df_res_temps <- df_res_temps[-which(duplicated(df_res_temps$residus)),]
    
    moy <- df_res_temps %>%
      group_by(ecart_temp) %>%
      summarize(moy = mean(residus)) %>%
      as.data.frame()
    
    
    nb_repet <- sort(c(dist(x = unique(levels(as.factor(time))))))
    
    df4 <- data.frame(table(df_res_temps$ecart_temp), moy$moy)
    df4.2 <- NULL
    
    moy2 <- NULL
    for (i in levels(as.factor(moy$ecart_temp))) {
      
      for (j in 1:length(nb_repet[nb_repet==i])) {
        moy2 <- rbind(moy2, moy[moy$ecart_temp==i,])
        df4.2 <- rbind(df4.2, df4[df4$Var1 == i,])
      }
    }
    
    moy <- moy2
    df4 <- df4.2
    
    temps_unique <- unique(nb_repet)
    
    df_res2 <- data.frame(
      res = unname(resid(res.lm)),
      id =id,
      time = time
    )
    
    total_var_sigma2_tau2_nu2 <- function(df_res2) {
      
      dfff2 <- c()
      
      for (i in unique(df_res2$time)) {
        dfff2 <- c(dfff2, c(dist(df_res2$res[which(df_res2$time == i)], upper = F, method = "euclidean")))
      }
      
      sum(dfff2**2)/2/length(dfff2)
      
    }
    
    bw1 <- np::npregbw(xdat = moy$ecart_temp, ydat = moy$moy, regtype = "ll", ckertype="epanechnikov", bwscaling = F, bandwidth.compute = T, bwmethod = "cv.ls")
    kre1 <- np::npreg(bws = bw1, exdat =  as.numeric(as.character(temps))-as.numeric(as.character(temps[1])))
    
    kre1$mean[1] <- 0
    
    res.lme.toep <- lme(fixed = formula.moyenne, data = df, correlation = corARMA(form = ~ 1| center/id, p = n.mesure - 1), random = ~ 1 | center, control = 
                          lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
    res.lme.ar1 <- lme(fixed = formula.moyenne, data = df, correlation = corCAR1(form = ~ 1 | center/id), random = ~ 1 | center, control = 
                         lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
    res.lme.cs <- lme(fixed = formula.moyenne, data = df, correlation = corCompSymm(form = ~ 1 | center/id), random = ~ 1 | center, control = 
                        lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
    res.lme.gaus <- lme(fixed = formula.moyenne, data = df, correlation = corGaus(form = ~ 1 | center/id), random = ~ 1 | center, control = 
                          lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
    res.lme.vc <- lme(fixed = formula.moyenne, data = df, correlation = corIdent(form = ~ 1 | center/id), random = ~ 1 | center, control = 
                        lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
    
    vec_aic <- sapply(list(res.lme.toep, res.lme.ar1, res.lme.cs, res.lme.gaus, res.lme.vc), AIC)
    vec_bic <- sapply(list(res.lme.toep, res.lme.ar1, res.lme.cs, res.lme.gaus, res.lme.vc), BIC)
    vec_aicc <- sapply(list(res.lme.toep, res.lme.ar1, res.lme.cs, res.lme.gaus, res.lme.vc), AICc)
    vec_aic_plus_bic <- vec_aic + vec_bic
    
    noms <- c("TOEP", "AR1", "CS", "Gaus", "VC")
    pred_aic <- noms[which.min(vec_aic)] ; pred_bic <- noms[which.min(vec_bic)]
    pred_aicc <- noms[which.min(vec_aicc)] ; pred_aic_plus_bic <- noms[which.min(vec_aic_plus_bic)]
    
    
    m.ar1 <- getvarcov2(res.lme.ar1)
    m.toep <- getvarcov2(res.lme.toep)
    m.cs <- getvarcov2(res.lme.cs)
    m.gaus <- getvarcov2(res.lme.gaus)
    m.vc <- diag(sigma(res.lme.vc)**2, nrow = dim(m.ar1)[1])
    
    
    m.ar1 <- matrix(data = m.ar1[1,1], ncol = n.mesure, nrow = n.mesure) - m.ar1
    m.toep <- matrix(data = m.toep[1,1], ncol = n.mesure, nrow = n.mesure) - m.toep
    m.cs <- matrix(data = m.cs[1,1], ncol = n.mesure, nrow = n.mesure) - m.cs
    m.gaus <- matrix(data = m.gaus[1,1], ncol = n.mesure, nrow = n.mesure) - m.gaus
    m.vc <- matrix(data = m.vc[1,1], ncol = n.mesure, nrow = n.mesure) - m.vc
    
    df_temps_loess <- data.frame(
      ecart_temp = as.numeric(as.character(temps))-as.numeric(as.character(temps))[1],
      ar1 = m.ar1[1,],
      toep = m.toep[1,],
      cs = m.cs[1,],
      vc = m.vc[1,],
      gaus = m.gaus[1,],
      var_np = kre1$mean
    )
    
    
    eqm <- function(fitted, np=df_temps_loess$var_np) {
      sum((fitted - np)**2)/length(fitted)
    }
    
    vec.eqm <- unlist(lapply(as.list(df_temps_loess[,-c(1,7)]), eqm))
    
    names(vec.eqm) <- c("AR1", "TOEP", "CS", "VC", "Gaus")
    pred_variog <- names(vec.eqm)[which.min(vec.eqm)]
    
    if (sqrt(vec.eqm[which(names(vec.eqm) == info_struct$Type)]) < sqrt(vec.eqm[which.min(vec.eqm)]) + seuil &
        sqrt(vec.eqm[which(names(vec.eqm) == info_struct$Type)]) > sqrt(vec.eqm[which.min(vec.eqm)]) - seuil) {
      pred_variog <- info_struct$Type
    }
    
    
    print(paste(info_struct$Type, "rho=", rho, "pat=", info_struct$n_patient, "ni=", info_struct$n_mesure, pred_variog))
    
    return(list(pred = pred_variog, pat = info_struct$n_patient, n_mesure = info_struct$n_mesure, type = info_struct$Type,
                pred_aic = pred_aic, pred_bic = pred_bic, pred_aicc = pred_aicc, pred_aic_plus_bic = pred_aic_plus_bic, pval = pval.iso))
    
  } else {
    
    
    res.lm <- lme(data = df, fixed = formula.moyenne, random = ~ 1 | center)
    
    df_res <- data.frame(id = id, res = res.lm$residuals[,2])
    df_dist <- NULL
    df_resid <- NULL
    
    for (i in levels(id)) {
      df_dist <- rbind(df_dist, melt(as.matrix(dist(time[which(id == i)]))))
      df_resid <- rbind(df_resid, melt(as.matrix(dist(df_res$res[which(id == i)]))))
    }
    
    
    df_resid$value <- (df_resid$value**2)/2
    
    df_res_temps <- data.frame(residus = df_resid$value, ecart_temp = df_dist$value, temps1 = df_dist$Var1, temps2 = df_dist$Var2)
    df_res_temps <- df_res_temps[-which(df_res_temps$ecart_temp == 0),]
    
    
    moy <- df_res_temps %>%
      group_by(temps1, temps2, ecart_temp) %>%
      summarize(moy = mean(residus)) %>%
      as.data.frame()
    
    moy$temps1 <- as.factor(moy$temps1)
    moy$temps2 <- as.factor(moy$temps2)
    
    levels(moy$temps1) <- temps
    levels(moy$temps2) <- temps
    
    moy$temps1 <- as.numeric(as.character(moy$temps1))
    moy$temps2 <- as.numeric(as.character(moy$temps2))
    
    x_grid <- temps
    y_grid <- temps
    xy_grid <- as.matrix(expand.grid(x_grid, y_grid))
    
    bw1 <- np::npregbw(xdat = moy[,1:2], ydat = moy[,4], regtype = "ll", ckertype="epanechnikov", bwscaling = T, bandwidth.compute = T, bwmethod = "cv.aic")
    
    if (mean(bw1$bw) > 1) {
      bw1$bw <- c(0.1, 0.1)
    }
    
    kre1 <- np::npreg(bws = bw1, exdat = xy_grid)
    
    kre1$mean[which(kre1$eval[,1] == kre1$eval[,2])] <- 0
    
    mat2 <- matrix(kre1$mean, nrow = length(x_grid), ncol = length(x_grid))
    
    
    
    res_toeph <- lme(fixed = formula.moyenne,
                     data = df,
                     correlation = corARMA(form = ~ 1 | center/id, p = n.mesure-1),
                     weights = varIdent(form = ~ 1 | time),
                     random = ~ 1 | center,
                     control = lmeControl(niterEM = 500, msMaxEval = 5000, maxIter = 500, msMaxIter = 500))
    res_toeph_mat <- getvarcov2(res_toeph)
    
    
    res_ar1h <- lme(fixed = formula.moyenne,
                    data = df,
                    correlation = corCAR1(form = ~ 1 |  center/id),
                    weights = varIdent(form = ~ 1 | time),
                    random = ~ 1 | center,
                    control = lmeControl(niterEM = 500, msMaxEval = 5000, maxIter = 500, msMaxIter = 500))
    res_ar1h_mat <- getvarcov2(res_ar1h)
    
    res_csh <- lme(fixed = formula.moyenne,
                   data = df,
                   correlation = corCompSymm(form = ~ 1 |  center/id),
                   weights = varIdent(form = ~ 1 | time),
                   random = ~ 1 | center,
                   control = lmeControl(niterEM = 500, msMaxEval = 5000, maxIter = 500, msMaxIter = 500))
    res_csh_mat <- getvarcov2(res_csh)
    
    res_gaush <- lme(fixed = formula.moyenne,
                     data = df,
                     correlation = corGaus(form = ~ 1 |  center/id),
                     weights = varIdent(form = ~ 1 | time),
                     random = ~ 1 | center,
                     control = lmeControl(niterEM = 500, msMaxEval = 5000, maxIter = 500, msMaxIter = 500, opt = "optim"))
    res_gaush_mat <- getvarcov2(res_gaush)
    
    
    res_un <- lme(fixed = formula.moyenne,
                  data = df,
                  correlation = corSymm(form = ~ 1 |  center/id),
                  weights = varIdent(form = ~ 1 | time),
                  random = ~ 1 | center,
                  control = lmeControl(niterEM = 500, msMaxEval = 5000, maxIter = 500, msMaxIter = 500, opt = "optim"))
    res_un_mat <- getvarcov2(res_un)
    
    vec_aic <- sapply(list(res_toeph, res_ar1h, res_csh, res_gaush, res_un), AIC)
    vec_bic <- sapply(list(res_toeph, res_ar1h, res_csh, res_gaush, res_un), BIC)
    vec_aicc <- sapply(list(res_toeph, res_ar1h, res_csh, res_gaush, res_un), AICc)
    vec_aic_plus_bic <- vec_aic + vec_bic
    
    noms <- c("ToepH", "AR1H", "CSH", "GAUSH", "UN")
    pred_aic <- noms[which.min(vec_aic)] ; pred_bic <- noms[which.min(vec_bic)]
    pred_aicc <- noms[which.min(vec_aicc)] ; pred_aic_plus_bic <- noms[which.min(vec_aic_plus_bic)]
    
    
    for (i in 1:n.mesure) {
      for (j in 1:n.mesure) {
        if (i != j) {
          res_toeph_mat[i,j] <- (res_toeph_mat[i,i] + res_toeph_mat[j,j])/2 - res_toeph_mat[i,j]
          res_ar1h_mat[i,j] <- (res_ar1h_mat[i,i] + res_ar1h_mat[j,j])/2 - res_ar1h_mat[i,j]
          res_csh_mat[i,j] <- (res_csh_mat[i,i] + res_csh_mat[j,j])/2 - res_csh_mat[i,j]
          res_gaush_mat[i,j] <- (res_gaush_mat[i,i] + res_gaush_mat[j,j])/2 - res_gaush_mat[i,j]
          res_un_mat[i,j] <- (res_un_mat[i,i] + res_un_mat[j,j])/2 - res_un_mat[i,j]
        }
      }
    }
    diag(res_toeph_mat) <- 0
    diag(res_ar1h_mat) <- 0
    diag(res_csh_mat) <- 0
    diag(res_gaush_mat) <- 0
    diag(res_un_mat) <- 0
    
    
    eqm <- function(var_param, var_np=mat2){
      sum((var_param - var_np)**2) / (dim(var_np)[1] * dim(var_np)[2])
    }
    
    
    list_eqm <- lapply(list(un = res_un_mat, toeph = res_toeph_mat, ar1h = res_ar1h_mat, csh = res_csh_mat, gaush = res_gaush_mat) , eqm)
    
    nb_param_var <- function(model){
      length(model$modelStruct$corStruct)+length(model$modelStruct$varStruct)+1
    }
    
    list_nbparam <- lapply(list(un = res_un, toeph = res_toeph, ar1h = res_ar1h, csh = res_csh, gaush = res_gaush) , nb_param_var)
    
    score <- unlist(list_eqm) * unlist(list_nbparam)
    
    
    
    names(score) <- c("UN", "ToepH" , "AR1H", "CSH", "GAUSH")
    
    
    pred_variog <- names(score)[which.min(score)]
    
    if (sqrt(score[which(names(score) == info_struct$Type)]) < sqrt(score[which.min(score)]) + seuil &
        sqrt(score[which(names(score) == info_struct$Type)]) > sqrt(score[which.min(score)]) - seuil) {
      pred_variog <- info_struct$Type
    }
    
    
    print(paste(info_struct$Type, "rho=", rho, "pat=", info_struct$n_patient, "ni=", info_struct$n_mesure, pred_variog))
    
    return(list(pred = pred_variog, pat = info_struct$n_patient, n_mesure = info_struct$n_mesure, type = info_struct$Type,
                pred_aic = pred_aic, pred_bic = pred_bic, pred_aicc = pred_aicc, pred_aic_plus_bic = pred_aic_plus_bic, pval = pval.iso))
  }
  
}



# Replications ------------------------------------------------------------

concate <- function(dfs=list(), n_rep=1, vec_mesure, vec_patient, vec_mat) {
  
  for (i_patient in vec_patient) {
    for (i_mesure in vec_mesure) {
      for (i_mat in vec_mat) {
        
        dfs <- c(dfs, replicate(n_rep, simu(type_corr = i_mat, nb_mesure = i_mesure, n_patient = i_patient), simplify = F))
        
        # dfs <- c(dfs, parSapply(cl, 1:n_rep, FUN = simu, type_corr = i_mat, nb_mesure = i_mesure, n_patient = i_patient, simplify = F))
        
        # dfs <- c(dfs, sapply(1:n_rep, FUN = simu, type_corr = i_mat, nb_mesure = i_mesure, n_patient = i_patient))
        
        # parallel version of replicate
      }
    }
  }
  return(dfs)
}

# n.iter <- 1000
# 
# dfs <- concate(n_rep = n.iter, vec_mesure = c(3,5,7), vec_patient = c(7,25,50),  vec_mat = c("AR1", "Gaus", "CS", "VC", "TOEP", "ToepH", "GAUSH", "CSH", "UN", "AR1H"))



# Résultats ---------------------------------------------------------------

no_cores <- detectCores() # Calculate the number of cores

cl <- makeCluster(no_cores) # je prends tous les coeurs disponibles 

clusterSetRNGStream(cl) # Set a different seed on each member of the cluster (just in case)

clusterExport(cl, c("reussite_variog_total")) # exporter la fonction de simulation

resultat <- function(dfs) {
  
  res <- parLapply(cl = cl, X = dfs, fun = function(x){try(reussite_variog_total(x),T)})
  
  # res <- lapply(X = dfs, FUN = function(x) try(reussite_variog(x),T))

  for (j in length(res):1) {
    if (typeof(res[[j]])=="character") {
      res <- res[-j]
    }
  }

  return(res)
}

# 
# res <- resultat(dfs)

dfs_pred <- unlist(lapply(res, function(x){x$pred}))
dfs_nb_pat <- unlist(lapply(res, function(x){x$pat}))
dfs_nb_mesure <- unlist(lapply(res, function(x){x$n_mesure}))
dfs_simu <- unlist(lapply(res, function(x){x$type}))
dfs_pred_aic <- unlist(lapply(res, function(x){x$pred_aic}))
dfs_pred_bic <- unlist(lapply(res, function(x){x$pred_bic}))
dfs_pred_aicc <- unlist(lapply(res, function(x){x$pred_aicc}))
dfs_pred_aic_plus_bic <- unlist(lapply(res, function(x){x$pred_aic_plus_bic}))
dfs_test_iso <- unlist(lapply(res, function(x){x$pval})) <= 0.05

df_final <- data.frame(
  pred = dfs_pred,
  n_mesure = dfs_nb_mesure,
  pat = dfs_nb_pat,
  type = dfs_simu,
  bonne_pred = dfs_pred == dfs_simu,
  bonne_pred_aic = dfs_pred_aic == dfs_simu,
  bonne_pred_bic = dfs_pred_bic == dfs_simu,
  bonne_pred_aicc = dfs_pred_aicc == dfs_simu,
  bonne_pred_aic_plus_bic = dfs_pred_aic_plus_bic == dfs_simu,
  bonne_decision_iso = (dfs_simu %in% c("AR1", "Gaus", "CS", "VC", "TOEP") == !dfs_test_iso)
)


df_final <- df_final %>% group_by(pat, n_mesure, type) %>%
  summarize(taux_de_reussite = mean(bonne_pred), reu_aic = mean(bonne_pred_aic),
            reu_bic = mean(bonne_pred_bic), reu_aicc = mean(bonne_pred_aicc),
            reu_aic_plus_bic = mean(bonne_pred_aic_plus_bic),
            reu_test_iso = mean(bonne_decision_iso)) %>%
  as.data.frame()

colMeans(df_final[,-c(1:3)])



# Sauvegarde --------------------------------------------------------------

for (i in 1:100) {
  n.iter <- 10
  dfs <- concate(n_rep = n.iter, vec_mesure = c(3,5,7), vec_patient = c(10,25,50),  vec_mat = c("AR1", "Gaus", "CS", "VC", "TOEP", "ToepH", "GAUSH", "CSH", "UN", "AR1H"))
  

  res <- resultat(dfs)
  
  dfs_pred <- unlist(lapply(res, function(x){x$pred}))
  dfs_nb_pat <- unlist(lapply(res, function(x){x$pat}))
  dfs_nb_mesure <- unlist(lapply(res, function(x){x$n_mesure}))
  dfs_simu <- unlist(lapply(res, function(x){x$type}))
  dfs_pred_aic <- unlist(lapply(res, function(x){x$pred_aic}))
  dfs_pred_bic <- unlist(lapply(res, function(x){x$pred_bic}))
  dfs_pred_aicc <- unlist(lapply(res, function(x){x$pred_aicc}))
  dfs_pred_aic_plus_bic <- unlist(lapply(res, function(x){x$pred_aic_plus_bic}))
  dfs_test_iso <- unlist(lapply(res, function(x){x$pval})) <= 0.05
  
  df_final <- data.frame(
    pred = dfs_pred,
    n_mesure = dfs_nb_mesure,
    pat = dfs_nb_pat,
    type = dfs_simu,
    bonne_pred = dfs_pred == dfs_simu,
    bonne_pred_aic = dfs_pred_aic == dfs_simu,
    bonne_pred_bic = dfs_pred_bic == dfs_simu,
    bonne_pred_aicc = dfs_pred_aicc == dfs_simu,
    bonne_pred_aic_plus_bic = dfs_pred_aic_plus_bic == dfs_simu,
    bonne_decision_iso = (dfs_simu %in% c("AR1", "Gaus", "CS", "VC", "TOEP") == !dfs_test_iso)
  )
  
  write.csv(df_final, file = paste0("RESULTATS/comp_variog_", format(Sys.time(), "%d%m%Y"), "_", i, ".csv"), row.names = F)
}

 
# df_final <- read.csv("RESULTATS/comp_variog_08122022.csv")


# Ploting -----------------------------------------------------------------

df_ggplot <- data.frame(
  Valeur = c(df_final$taux_de_reussite, df_final$reu_aic, df_final$reu_bic, df_final$reu_aicc, df_final$reu_aic_plus_bic),
  crit = c(rep("Variogram", dim(df_final)[1]), rep("AIC", dim(df_final)[1]), rep("BIC", dim(df_final)[1]), rep("AICc", dim(df_final)[1]), rep("AIC+BIC", dim(df_final)[1])),
  pat = rep(df_final$pat, 5),
  n_mesure = rep(df_final$n_mesure, 5),
  type = rep(df_final$type, 5)
)

df_ggplot$type <- toupper(df_ggplot$type)
df_ggplot <- df_ggplot[-which(df_ggplot$crit == "AIC+BIC"),]
df_ggplot <- df_ggplot[-which(df_ggplot$crit == "AICc"),]

n.iter <- 1000

df_ggplot$type <- as.factor(df_ggplot$type)
levels(df_ggplot$type) <- c("AR(1)", "ARH(1)", "CS", "CSH", "GAUS", "GAUSH", "TOEP", "TOEPH", "UN", "VC")

ggplot(data = df_ggplot[which(df_ggplot$type %in% c("AR(1)", "CS", "GAUS", "TOEP", "VC")),], aes(x = pat, y = Valeur, group = crit, col = crit)) +
  # geom_errorbar(aes(ymin = Valeur - qnorm(0.975)*sqrt(Valeur*(1-Valeur)/(n.iter)), ymax = Valeur + qnorm(0.975)*sqrt(Valeur*(1-Valeur)/(n.iter))), position = position_dodge(width = 0.9)) +
  geom_point(position=position_dodge(width = 0.9), size = 3, aes(shape = crit)) +
  facet_grid(rows = vars(n_mesure), cols = vars(type)) +
  theme_minimal() + 
  scale_x_continuous(breaks = c(10,25,50)) +
  xlab("Number of patients") +
  ylab("Criterion variance-covariance \nrecovery rate") +
  labs(col = 'Criterion', shape = "Criterion") +
  ggtitle("Criterion comparison", subtitle =  "Lines : Number of repeated measurement\nColumns : Type of simulated covariance") +
  theme(legend.position = "bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black"))

ggplot(data = df_ggplot[which(!(df_ggplot$type %in% c("AR(1)", "CS", "GAUS", "TOEP", "VC"))),], aes(x = pat, y = Valeur, group = crit, col = crit)) +
  # geom_errorbar(aes(ymin = Valeur - qnorm(0.975)*sqrt(Valeur*(1-Valeur)/(n.iter)), ymax = Valeur + qnorm(0.975)*sqrt(Valeur*(1-Valeur)/(n.iter))), position = position_dodge(width = 0.9)) +
  geom_point(position=position_dodge(width = 0.9), size = 3, aes(shape = crit)) +
  facet_grid(rows = vars(n_mesure), cols = vars(type), ) +
  theme_minimal() + 
  scale_x_continuous(breaks = c(10,25,50)) +
  xlab("Number of patients") +
  ylab("Criterion variance-covariance \nrecovery rate") +
  labs(col = 'Criterion', shape = "Criterion") +
  ggtitle("Criterion comparison", subtitle =  "Lines : Number of repeated measurement\nColumns : Type of simulated covariance") +
  theme(legend.position = "bottom") +
  theme(panel.background = element_rect(fill = NA, color = "black"))

mean(df_final$reu_test_iso)
mean(df_final$taux_de_reussite)
mean(df_final$reu_aic)
mean(df_final$reu_bic)


# -------------------------------------------------------------------------

df <- df_final[,which(names(df_final) %in%
                        c("pat", "n_mesure", "type", "taux_de_reussite", "reu_aic", "reu_bic"))][,c(3,1,2,4,5,6)]

df <- df[order(df_final$type, df_final$pat, df_final$n_mesure),]

df$type <- as.factor(df$type)

levels(df$type) <- c("AR(1)", "ARH(1)", "CS", "CSH", "GAUS", "GAUSH", "TOEP", "TOEPH", "UN", "VC")

df$convergence <- convergence_type_pat_mesure$convergence

obj <- xtable::xtable(df, type = "latex", file = "filename2.tex", digits = 1, auto = T)


