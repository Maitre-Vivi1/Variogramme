# Graphique exemple du variogram ------------------------------------------


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
  
  longData$trt <- rep(trt, each=nb_mesure, 1)
  
  vec_int <- rep(1, length(nb_mesure))
  
  return(longData)  
}


AnalyseVariogramme::analyse_variog()
set.seed(3)

  df<- simu(type_corr = "AR1", nb_mesure = 9, n_patient = 100, temps = 1:9)
  id <- df$id
  time <- df$time
  n.mesure = length(unique(time))
  n.patient = length(unique(id))
  
  temps <- sort(unique(time))
  
  # -------------------------------------------------------------------------
  formula.moyenne <- Y ~ trt + time 
  

p <- get_vari()
p

p <- p +
  geom_point() +
  xlab("Lag") +
  ylab("Variance of the difference") +
  theme_minimal() +
  guides(linetype = guide_legend(title = ""), col = guide_legend(title = "")) +
  theme_minimal() +
  scale_x_continuous(breaks = 1:9) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))
  

p



# Graphique exemple fonction de variance ----------------------------------

data.sim <- simu(type_corr = "AR1", nb_mesure = 9, n_patient = 50)
df <- data.sim
id <- df$id
time <- df$time
n.mesure = length(unique(time))
n.patient = length(unique(id))

levels(id) <- 1:n.patient
require(ggpubr)
require(MASS)
require(dplyr)
require(ggplot2)
require(nlme)
require(reshape2)
require(plotly)


if ("time" %in% names(df)) {
  df <- df[,-which(names(df)=="time")]
}
if ("id" %in% names(df)) {
  df <- df[,-which(names(df)=="id")]
}

df <- cbind(df, id, time)

temps <- sort(unique(time))

  # -------------------------------------------------------------------------
formula.moyenne <- Y ~ trt + time 
  
  z <- unname(resid(lm(data = df, formula =  formula.moyenne)))
  z_2 <- z**2
  
  
  df_res_time <- data.frame(
    res = z_2,
    time = time
  )
  
  fct_variance <- lm(data = df_res_time, formula = res ~ as.factor(time) - 1)
  
  summary(lm(data = df_res_time, formula = res ~ as.factor(time)))
  anova(lm(data = df_res_time, formula = res ~ as.factor(time)))
  
  
  IC <- as.data.frame(confint(fct_variance))
  
  IC$moy <- coef(fct_variance)
  IC$temps <- unique(df_res_time$time)
  
ggplot(data = NULL, aes(x = as.numeric(as.character(IC$temps)))) +
  geom_line(aes( y= IC$moy)) +
  geom_ribbon(aes(ymin = IC$`2.5 %`, ymax = IC$`97.5 %`), alpha = 0.2) +
  geom_point(aes( x = as.numeric(as.character(df_res_time$time)), y = df_res_time$res)) +
  xlab("Timepoints") +
  ylab("Variance") +
  scale_x_continuous(breaks = 1:9) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme_minimal()




# Graphique variog iso exemple --------------------------------------------

get_vari <- function() {
  
       res.lm <- lm(data = df, formula = formula.moyenne)

       df_res <- data.frame(id = id, res = resid(res.lm))
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


       bw1 <- np::npregbw(xdat = moy$ecart_temp, ydat = moy$moy, regtype = "ll", ckertype="epanechnikov", bwscaling = T, bandwidth.compute = T, bwmethod = "cv.aic")
       kre1 <- np::npreg(bws = bw1, exdat =  as.numeric(as.character(temps))-as.numeric(as.character(temps[1])))

       kre1$mean[1] <- 0

       
       df_ggplot <- data.frame(
         vari = kre1$mean,
         lagf = 0:8
       )
       
    return(ggplot(data = df_ggplot, aes(x = lagf, y = vari)) +
      geom_line())
       
}

get_vari()

res.lme.toep <- gls(model = formula.moyenne, data = df, correlation = corARMA(form = ~ 1| id, p = n.mesure - 1))
res.lme.ar1 <- gls(model = formula.moyenne, data = df, correlation = corCAR1(form = ~ time | id))
res.lme.cs <- gls(model = formula.moyenne, data = df, correlation = corCompSymm(form = ~ 1 | id))
res.lme.gaus <- gls(model = formula.moyenne, data = df, correlation = corGaus(form = ~ 1 | id))
res.lme.vc <- gls(model = formula.moyenne, data = df)


m.ar1 <- as.matrix(getVarCov(res.lme.ar1, type = "conditional"))
m.toep <- as.matrix(getVarCov(res.lme.toep, type = "conditional"))
m.cs <- as.matrix(getVarCov(res.lme.cs, type = "conditional"))
m.gaus <- as.matrix(getVarCov(res.lme.gaus, type = "conditional"))
m.vc <- diag(total_var_sigma2_tau2_nu2(df_res2), nrow = n.mesure)


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

df_ggplot <- data.frame(
  semivariogramme = c(df_temps_loess$var_np, ar1 = df_temps_loess$ar1, toep = df_temps_loess$toep, cs = df_temps_loess$cs,
                      vc = df_temps_loess$vc, gaus = df_temps_loess$gaus),
  ecart_temp = rep(df_temps_loess$ecart_temp ,6),
  type = as.factor(rep(x = c("variog_np", "ar1", "toep", "cs", "vc", "gaus"), each = length(df_temps_loess$ecart_temp)))
)

p.variogs.iso <- ggplot(df_ggplot, aes(x = ecart_temp, y = semivariogramme, col = type, group = type, linetype = type)) +
  guides(linetype = guide_legend(title = "", nrow = 1), col = guide_legend(title = "", nrow = 1)) +
  labs(linetype = "") +
  geom_line(linewidth = 1.05) +
  theme_minimal() +
  theme(legend.position="bottom") +
  xlab("Lag") +
  ylab("Variogram") +
  scale_x_continuous(breaks = df_temps_loess$ecart_temp) +
  scale_colour_manual(labels = c("AR(1)", "CS", "GAUS", "Nonparametric", "TOEP", "VC"),
                      values = 1:6) +
  scale_linetype_manual(labels = c("AR(1)", "CS", "GAUS", "Nonparametric", "TOEP", "VC"),
                        values = c(1:6)) +
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14),
        legend.text = element_text(size=14)) +
  guides()
  
p.variogs.iso



geom_line() +
  scale_x_continuous(breaks = df_temps_loess$ecart_temp) +
  guides(linetype = guide_legend(title = ""), col = guide_legend(title = "")) +
  labs(linetype = "") +
  theme_minimal() +
  theme(legend.position="bottom") +
  xlab("Lag") +
  ylab("Variogram") +
  scale_colour_manual(labels = c("AR(1)", "CS", "GAUS", "Nonparametric", "TOEP", "VC"),
                      values = 1:6) +
  scale_linetype_manual(labels = c("AR(1)", "CS", "GAUS", "Nonparametric", "TOEP", "VC"),
                        values = c(1:6))




eqm <- function(fitted, np=df_temps_loess$var_np) {
  sum((fitted - np)**2)/length(fitted)
}

list_eqm <- lapply(as.list(df_temps_loess[,-c(1,7)]), eqm)

nb_param_var <- function(model){
  length(model$modelStruct$corStruct)+length(model$modelStruct$varStruct)+1
}


df_eqm_nbpar <- data.frame(
  type = names(unlist(list_eqm)),
  eqm = unlist(list_eqm),
  nb_par = unlist(lapply(list(ar1 = res.lme.ar1, toep = res.lme.toep, cs = res.lme.cs, vc=res.lme.vc, gaus=res.lme.gaus) , nb_param_var))
)

p.nb_param.eqm <- ggplot(data = df_eqm_nbpar, aes(x = nb_par, y = eqm, group = type, col = type, linetype = type)) +
  geom_point() +
  xlab("lag") +
  ylab("Variograms") +
  scale_x_continuous(breaks = 1:(n.mesure*(n.mesure+1)/2))





# Graphique heat map variog aniso -----------------------------------------

set.seed(4)

df <- simu(type_corr = "AR1H", n_patient = 50, nb_mesure = 7)
temps <- 1:7
n.mesure <- 7


formula.moyenne <- Y ~ time + trt

res.lm <- lme(data = df, fixed = formula.moyenne, random = ~ 1 | id)$residuals[,2]

id <- df$id
time <- df$time

df_res <- data.frame(id = id, res = res.lm)
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
if (bw1$bw[1] > 1 | bw1$bw[2] > 5) {
  bw1$bw <- c(1,1)
}
kre1 <- np::npreg(bws = bw1, exdat = xy_grid)

kre1$mean[which(kre1$eval[,1] == kre1$eval[,2])] <- 0

mat2 <- matrix(kre1$mean, nrow = length(x_grid), ncol = length(x_grid))


res_toeph <- gls(model = formula.moyenne,
                 data = df,
                 correlation = corARMA(form = ~ 1 | id, p = n.mesure-1),
                 weights = varIdent(form = ~ 1 | time),
                 control = lmeControl(niterEM = 500, msMaxEval = 5000, maxIter = 500, msMaxIter = 500))
res_toeph_mat <- getVarCov(res_toeph)


res_ar1h <- gls(model = formula.moyenne,
                data = df,
                correlation = corCAR1(form = ~ 1 | id),
                weights = varIdent(form = ~ 1 | time),
                control = lmeControl(niterEM = 500, msMaxEval = 5000, maxIter = 500, msMaxIter = 500))
res_ar1h_mat <- getVarCov(res_ar1h)

res_csh <- gls(model = formula.moyenne,
               data = df,
               correlation = corCompSymm(form = ~ 1 | id),
               weights = varIdent(form = ~ 1 | time),
               control = lmeControl(niterEM = 500, msMaxEval = 5000, maxIter = 500, msMaxIter = 500))
res_csh_mat <- getVarCov(res_csh)

res_gaush <- gls(model = formula.moyenne,
                 data = df,
                 correlation = corGaus(form = ~ 1 | id),
                 weights = varIdent(form = ~ 1 | time),
                 control = lmeControl(niterEM = 500, msMaxEval = 5000, maxIter = 500, msMaxIter = 500, opt = "optim"))
res_gaush_mat <- getVarCov(res_gaush)


res_un <- gls(model = formula.moyenne,
              data = df,
              correlation = corSymm(form = ~ 1 | id),
              weights = varIdent(form = ~ 1 | time),
              control = lmeControl(niterEM = 500, msMaxEval = 5000, maxIter = 500, msMaxIter = 500, opt = "optim"))
res_un_mat <- getVarCov(res_un)

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


vec_aic <- sapply(list(res_toeph, res_ar1h, res_csh, res_gaush, res_un), AIC)
vec_bic <- sapply(list(res_toeph, res_ar1h, res_csh, res_gaush, res_un), BIC)

noms <- c("TOEP", "ARH(1)", "CSG", "GAUSH", "UN")
pred_aic <- noms[which.min(vec_aic)] ; pred_bic <- noms[which.min(vec_bic)]




df_eqm_nbpar <- data.frame(
  type = names(unlist(list_eqm)),
  eqm = unlist(list_eqm),
  nb_par = unlist(list_nbparam)
)

ggplot(data = df_eqm_nbpar, aes(x = nb_par, y = eqm, group = type, col = type)) +
  geom_point() +
  scale_x_continuous(breaks = 1:(n.mesure*(n.mesure+1)/2))


df <- data.frame(
  variog = c(res_ar1h_mat, mat2),
  type = rep(c("ARH(1)", "Nonparametric"), each = 49),
  x = as.factor(rep(temps,n.mesure)),
  y = as.factor(rep(temps,rep(n.mesure,n.mesure)))
)


p3 <- ggplot(data = df, aes(x = x, y = y, fill = variog)) +
  geom_tile() +
  facet_wrap(~type) +
  xlab("Timepoints") +
  ylab("Timepoints") +
  scale_fill_viridis_c(option = 'B',  breaks = c(0,5,10,15,20,25,30,35,40,45,50)) +
  theme(legend.position="bottom") +
  guides(fill = guide_legend(title = "", nrow = 1))
  
p3


# -------------------------------------------------------------------------

