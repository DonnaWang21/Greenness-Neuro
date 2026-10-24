# NDVI and incident neurological outcomes: Cox proportional hazards models

library(survival)

dat <- read.csv("data/analysis.csv", stringsAsFactors = FALSE)
outcomes <- c("dementia", "parkinson", "epilepsy", "migraine")

time_vars <- setNames(paste0("followup_", outcomes), outcomes)
exposure <- "followup_ndvi_mean"
categorical <- c("sex", "edu", "income", "ever_smo", "drink", "PA")
model1 <- c("age", "sex", "edu", "income")
model2 <- c(model1, "BMI", "ever_smo", "drink", "PA")

pollution <- c(dementia = "PM2.5_2006_2cat", parkinson = "PM2.5_2006_2cat",
               epilepsy = "PM2.5_2006_2cat", migraine = "PM2.5_mean")

for (y in outcomes) {
  stopifnot(all(na.omit(dat[[y]]) %in% c(0, 1)),
            all(na.omit(dat[[time_vars[[y]]]]) > 0))
}
dat[categorical] <- lapply(dat[categorical], factor)
pm_cut <- median(dat$PM2.5_2006, na.rm = TRUE)
dat$PM2.5_2006_2cat <- factor(ifelse(dat$PM2.5_2006 <= pm_cut, "Low", "High"),
                                 levels = c("Low", "High"))

ndvi_iqr <- IQR(dat[[exposure]], na.rm = TRUE)
dat$ndvi_iqr <- dat[[exposure]] / ndvi_iqr

results <- list()
for (y in outcomes) {
  adjustments <- list(Model1 = model1, Model2 = model2,
                      Model3 = c(model2, pollution[[y]]))
  for (m in names(adjustments)) {
    vars <- c(y, time_vars[[y]], "ndvi_iqr", adjustments[[m]])
    d <- droplevels(dat[complete.cases(dat[vars]), , drop = FALSE])
    if (!nrow(d) || sum(d[[y]]) == 0) stop("No usable events: ", y, "/", m)
    f <- reformulate(c("ndvi_iqr", adjustments[[m]]),
                     response = paste0("Surv(", time_vars[[y]], ", ", y, ")"))
    fit <- coxph(f, data = d, ties = "efron", na.action = na.fail,
                 singular.ok = FALSE)
    s <- coef(summary(fit))["ndvi_iqr", ]
    ci <- exp(confint(fit, parm = "ndvi_iqr"))
    stopifnot(all(is.finite(c(s, ci))))
    results[[paste(y, m)]] <- data.frame(
      outcome = y, model = m, n = fit$n, events = fit$nevent,
      ndvi_IQR = ndvi_iqr, HR = exp(s["coef"]),
      CI_lower = ci[1], CI_upper = ci[2], P = s["Pr(>|z|)"],
      row.names = NULL)
  }
}
results <- do.call(rbind, results)
results$P_FDR <- NA_real_
idx <- results$model == "Model3"

results$P_FDR[idx] <- p.adjust(results$P[idx], method = "BH")
dir.create("results", showWarnings = FALSE)
write.csv(results, "results/ndvi-neuro.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "results/ndvi-neuro-session.txt")
