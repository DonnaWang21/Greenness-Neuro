# NDVI and standardized brain volumes: Gaussian identity-link models

dat <- read.csv("data/analysis.csv", stringsAsFactors = FALSE)
brain <- read.csv("data/brain_structure.csv", stringsAsFactors = FALSE)
exposure <- "followup_ndvi_mean"
regions <- grep("_mean_z$", names(brain), value = TRUE)
categorical <- c("sex", "edu", "income", "ever_smo", "drink", "PA")
covariates <- c("age", "sex", "edu", "income", "BMI", "ever_smo", "drink", "PA")

d <- merge(dat[required], brain[c("eid", regions)], by = "eid")
d[categorical] <- lapply(d[categorical], factor)

ndvi_iqr <- IQR(d[[exposure]], na.rm = TRUE)
d$ndvi_iqr <- d[[exposure]] / ndvi_iqr

results <- lapply(regions, function(region) {
  vars <- c(region, "ndvi_iqr", covariates)
  dd <- droplevels(d[complete.cases(d[vars]), , drop = FALSE])
  if (!nrow(dd)) stop("No complete observations: ", region)
  fit <- glm(reformulate(c("ndvi_iqr", covariates), response = region),
             data = dd, family = gaussian(link = "identity"),
             na.action = na.fail, singular.ok = FALSE)
  stopifnot(fit$converged, fit$df.residual > 0)
  s <- coef(summary(fit))["ndvi_iqr", ]
  ci <- suppressMessages(confint(fit, parm = "ndvi_iqr", level = 0.95))
  stopifnot(all(is.finite(c(s, ci))))
  data.frame(region = region, n = nobs(fit), ndvi_IQR = ndvi_iqr,
             beta = s["Estimate"], CI_lower = ci[1], CI_upper = ci[2],
             P = s["Pr(>|t|)"], row.names = NULL)
})
results <- do.call(rbind, results)

results$P_FDR <- p.adjust(results$P, method = "BH")
dir.create("results", showWarnings = FALSE)
write.csv(results, "results/ndvi-brain.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "results/ndvi-brain-session.txt")
