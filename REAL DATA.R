# Library -----------------------------------------------------------------

library(nlme)
library(lme4)
library(AnalyseVariogramme)


# Chargement des données --------------------------------------------------

data(sleepstudy)
data(BodyWeight)

BodyWeight <- BodyWeight[which(!BodyWeight$Diet == 2),]

df <- BodyWeight
# Variogramme -------------------------------------------------------------

analyse_variog(df = sleepstudy, id = sleepstudy$Subject, time = sleepstudy$Days, formula.moyenne = Reaction ~ Days, seuil_iso = 0.01)
analyse_variog(df = BodyWeight, id = BodyWeight$Rat, time = BodyWeight$Time, formula.moyenne = weight ~ Diet+Time, seuil_iso = 0.000005)
analyse_variog(df = paquid, id = paquid$ident, time = paquid$time, formula.moyenne = Isaacs ~ Age + Age_2 + Age_3 + niv.etu + sexe + v + compt_isaacs)
analyse_variog(df = paquid, id = paquid$ident, time = paquid$time, formula.moyenne = Isaacs ~ Age + Age_2 + Age_3 + niv.etu + sexe + v + compt_isaacs, seuil_iso = 0.00001)


# Mixed models ------------------------------------------------------------

summary(gls(model = Reaction ~ Days, data = sleepstudy, correlation = corCAR1(form =  ~ 1 | Subject)))
summary(gls(model = weight ~ Diet+Time, data = BodyWeight, correlation = corCAR1(form =  ~ Time | Rat)))
summary(gls(model = Isaacs ~ Age + I(Age**2) + I(Age**3) + niv.etu + sexe + v + compt_isaacs, data = paquid, correlation = corCompSymm(form =  ~ 1 | ident)))

# Publication -------------------------------------------------------------

ggplot(data = BodyWeight, aes(x= Time, y = weight, group = Rat, linetype = Diet)) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(linetype = "none") +
  xlab("Time") +
  ylab("Weight (g)") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))
  

id <- df$id
time <- df$time

temps <- sort(unique(df$time))

df <- BodyWeight

n.mesure <- length(unique(BodyWeight$Time))


formula.moyenne= Y ~ time + trt

seuil_iso = 0.05

names(df) <- c("Y", "time", "ID", "trt")

z <- unname(residuals(lm(data = df, formula =  formula.moyenne)))
z_2 <- z**2

time = df$time


df_res_time <- data.frame(
  res = z_2,
  time = time,
  trt = df$trt,
  name = df$ID
)

id  = df$ID
temps <- unique(df$time)
pval.iso <- car::Anova(lm(data = df_res_time, formula = res ~ as.factor(time)), type = "III")




fct_variance <- lm(data = df_res_time, formula = res ~ as.factor(time) -1  + trt)

IC <- as.data.frame(confint(fct_variance)[-12,])

IC$moy <- coef(fct_variance)[-12]
IC$temps <- unique(df_res_time$time)

p.fct_variance <- ggplot(data = NULL, aes(x = as.numeric(as.character(IC$temps)))) +
  geom_line(aes( y= IC$moy)) +
  geom_ribbon(aes(ymin = IC$`2.5 %`, ymax = IC$`97.5 %`), alpha = 0.2) +
  geom_point(aes( x = as.numeric(as.character(df_res_time$time)), y = df_res_time$res)) +
  xlab("Timepoints") +
  ylab("Squared residuals") +
  scale_x_continuous(breaks = as.numeric(as.character(temps))) +
  theme_minimal() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

p.fct_variance

res.lm <- lm(data = df, formula = formula.moyenne)

df_res <- data.frame(id = df$ID, res = resid(res.lm))
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

res.lme.toep <- gls(model = formula.moyenne, data = df, correlation = corARMA(form = ~ 1| ID, p = n.mesure - 1), control = 
                      lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
res.lme.ar1 <- gls(model = formula.moyenne, data = df, correlation = corCAR1(form = ~ time | ID), control = 
                     lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
res.lme.cs <- gls(model = formula.moyenne, data = df, correlation = corCompSymm(form = ~ 1 | ID), control = 
                    lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
res.lme.gaus <- gls(model = formula.moyenne, data = df, correlation = corGaus(form = ~ 1 | ID), control = 
                      lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))
res.lme.vc <- gls(model = formula.moyenne, data = df, correlation = corIdent(form = ~ 1 | ID), control = 
                    lmeControl(maxIter = 500,msMaxIter = 500,niterEM = 50))

vec_aic <- sapply(list(res.lme.toep, res.lme.ar1, res.lme.cs, res.lme.gaus, res.lme.vc), AIC)
vec_bic <- sapply(list(res.lme.toep, res.lme.ar1, res.lme.cs, res.lme.gaus, res.lme.vc), BIC)

noms <- c("TOEP", "AR1", "CS", "Gaus", "VC")
pred_aic <- noms[which.min(vec_aic)] ; pred_bic <- noms[which.min(vec_bic)]

m.ar1 <-  as.matrix(getVarCov(res.lme.ar1))
m.toep <- as.matrix(getVarCov(res.lme.toep))
m.cs <- as.matrix(getVarCov(res.lme.cs))
m.gaus <- as.matrix(getVarCov(res.lme.gaus))
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

names(vec.eqm) <- c("AR1", "TOEP", "CS", "VC", "GAUS")
pred_variog <- names(vec.eqm)[which.min(vec.eqm)]


df_ggplot <- data.frame(
  eqm = vec.eqm,
  type = names(vec.eqm),
  nb_param = c(2,n.mesure, 2, 1, 2)
)


ggplot(data = df_ggplot, aes(x= nb_param, y = log(eqm), group = type, shape = type, col = type)) +
  geom_point(size = 8) +
  labs(shape = "") +
  theme_minimal() +
  theme(legend.position="bottom") +
  guides(shape = guide_legend(title = ""), col = guide_legend(title = "")) +
  ylab("Mean square error (log)") +
  xlab("Number of parameters") +
  scale_x_continuous(breaks = c(1:11)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=14))


df_ggplot <- data.frame(
  semivariogramme = c(df_temps_loess$var_np, ar1 = df_temps_loess$ar1, toep = df_temps_loess$toep, cs = df_temps_loess$cs,
                      vc = df_temps_loess$vc, gaus = df_temps_loess$gaus),
  ecart_temp = rep(df_temps_loess$ecart_temp ,6),
  type = as.factor(rep(x = c("Nonparametric", "AR(1)", "TOEP", "CS", "VC", "GAUS"), each = length(df_temps_loess$ecart_temp)))
)




ggplot(df_ggplot, aes(x = ecart_temp, y = semivariogramme, group = type, linetype = type, col = type)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(breaks = df_temps_loess$ecart_temp) +
  guides(linetype = guide_legend(title = "", nrow = 1), col = guide_legend(title = "", nrow = 1)) +
  labs(linetype = "") +
  theme_minimal() +
  theme(legend.position="bottom") +
  xlab("Lag") +
  ylab("Variogram") +
  scale_colour_manual(labels = c("AR(1)", "CS", "GAUS", "Nonparametric", "TOEP", "VC"),
                      values = 1:6) +
  scale_linetype_manual(labels = c("AR(1)", "CS", "GAUS", "Nonparametric", "TOEP", "VC"),
                        values = c(1:6)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=14))

