files <- list.files(pattern = "*.csv", path = "~/PUBLICATION VARIOGRAMME/RESULTATS/11012023", recursive = TRUE, full.names = TRUE, all.files = TRUE, include.dirs = TRUE)

lst_of_frames <- lapply(files, read.csv)

oneframe <- do.call("rbind.data.frame", lst_of_frames)
df_final <- oneframe


convergence_type <- df_final %>% group_by(type) %>%
  summarize(convergence = round(n() / 9000 * 100,3))

convergence_type_pat_mesure <- df_final %>% group_by(type, pat, n_mesure) %>%
  summarize(convergence = round(n() / 10,3))




df_final <- df_final %>% group_by(pat, n_mesure, type) %>%
  summarize(taux_de_reussite = mean(bonne_pred), reu_aic = mean(bonne_pred_aic),
            reu_bic = mean(bonne_pred_bic), reu_aicc = mean(bonne_pred_aicc),
            reu_aic_plus_bic = mean(bonne_pred_aic_plus_bic),
            reu_test_iso = mean(bonne_decision_iso)) %>%
  as.data.frame()


colMeans(df_final[,-c(1:3)])

