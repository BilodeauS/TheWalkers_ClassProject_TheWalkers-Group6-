# =====================================================================
# BMEN 600 - Group 6 (The Walkers)
# Full analysis following the Midterm Research Plan (Sections 4.1, 4.3, 4.4)
#
#   Data: clean_data/sleepiness_eeg_clean.csv  (GitHub repo)
#         56 participants x 3 sessions = 168 eyes-closed recordings
#         (4 participants with missing KSS excluded - complete-case analysis)
#
#   Steps (numbered as in the plan):
#     1. Sample characteristics + descriptive statistics by session   (4.3)
#     2. STATE:  repeated-measures correlation, KSS vs log EEG          (4.1)
#     3. TRAIT:  Spearman, ESS vs each participant's mean log EEG       (4.1)
#     4. GROUPS: Mann-Whitney U, ESS <= 10 vs ESS > 10                  (4.3, 4.4)
#     5. Supplementary linear mixed-effects model (KSS + ESS + session) (4.1)
#     6. Comparison: pooled Pearson vs repeated-measures correlation    (4.4)
#     7. Figures
#   Primary outcome = log T/A ratio; secondary = log theta, log alpha
#   Bonferroni across the 3 EEG measures: alpha = 0.05 / 3 = 0.017
#
#   Run in RStudio: Code > Run Region > Run All (Cmd+Option+R on Mac)
#   Outputs: CSV tables + PNG figures in the folder printed at the end.
# =====================================================================

# ---- 0. Setup -----------------------------------------------------------
needed <- c("lme4", "lmerTest", "rmcorr")
for (p in needed) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(rmcorr)
})

ALPHA_BONF <- 0.05 / 3
EEG <- c(log_ratio = "log T/A ratio (primary)",
         log_theta = "log theta power",
         log_alpha = "log alpha power")

out_dir <- file.path(getwd(), "midterm_plan_results")
dir.create(out_dir, showWarnings = FALSE)

url <- paste0("https://raw.githubusercontent.com/BilodeauS/TheWalkers_ClassProject_TheWalkers-Group6-/",
              "refs/heads/main/clean_data/sleepiness_eeg_clean.csv")
d <- read.csv(url)
d$participant_id <- factor(d$participant_id)
d$session_f      <- factor(d$session)
d$log_ratio      <- log10(d$theta_alpha_ratio)      # log_theta, log_alpha already in file
stopifnot(nrow(d) == 168, nlevels(d$participant_id) == 56)

# One row per participant (used for trait analyses)
pl <- aggregate(cbind(log_ratio, log_theta, log_alpha, KSS) ~ participant_id, data = d, FUN = mean)
pl <- merge(pl, unique(d[, c("participant_id", "ESS_avg", "age", "sex")]), by = "participant_id")
pl$ESS_group <- factor(ifelse(pl$ESS_avg > 10, "ESS > 10", "ESS <= 10"),
                       levels = c("ESS <= 10", "ESS > 10"))

fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
header <- function(title) cat("\n", strrep("=", 70), "\n ", title, "\n", strrep("=", 70), "\n", sep = "")


# ---- 1. Sample characteristics + descriptives by session (4.3) ----------
header("1. SAMPLE AND DESCRIPTIVE STATISTICS")
cat(sprintf("Participants: %d (%d female, %d male); age %.2f +/- %.2f years\n",
            nrow(pl), sum(pl$sex == "f"), sum(pl$sex == "m"), mean(pl$age), sd(pl$age)))
cat(sprintf("ESS: mean %.2f +/- %.2f; ESS <= 10: n = %d, ESS > 10: n = %d\n",
            mean(pl$ESS_avg), sd(pl$ESS_avg), sum(pl$ESS_group == "ESS <= 10"), sum(pl$ESS_group == "ESS > 10")))

desc_vars <- c(KSS = "KSS", theta = "Theta power", alpha = "Alpha power",
               theta_alpha_ratio = "T/A ratio")
desc <- do.call(rbind, lapply(names(desc_vars), function(v) {
  do.call(rbind, lapply(1:3, function(s) {
    x <- d[d$session == s, v]
    data.frame(variable = desc_vars[[v]], session = s, n = sum(!is.na(x)),
               mean = mean(x), sd = sd(x), median = median(x), min = min(x), max = max(x))
  }))
}))
print(transform(desc, mean = round(mean, 3), sd = round(sd, 3), median = round(median, 3),
                min = round(min, 3), max = round(max, 3)), row.names = FALSE)
write.csv(desc, file.path(out_dir, "1_descriptives_by_session.csv"), row.names = FALSE)


# ---- 2. STATE: repeated-measures correlation (4.1) ---------------------
header("2. STATE SLEEPINESS - repeated-measures correlation (KSS vs EEG)")
state <- do.call(rbind, lapply(names(EEG), function(v) {
  tmp <- data.frame(pid = d$participant_id, x = d$KSS, y = d[[v]])
  r <- rmcorr(pid, x, y, tmp)
  data.frame(eeg = EEG[[v]], r_rm = r$r, ci_low = r$CI[1], ci_high = r$CI[2],
             df = r$df, p = r$p, significant_bonf = r$p < ALPHA_BONF)
}))
print(transform(state, r_rm = round(r_rm, 3), ci_low = round(ci_low, 3), ci_high = round(ci_high, 3),
                p = fmt_p(p)), row.names = FALSE)
write.csv(state, file.path(out_dir, "2_state_rmcorr.csv"), row.names = FALSE)


# ---- 3. TRAIT: Spearman correlation (4.1) -------------------------------
header("3. TRAIT SLEEPINESS - Spearman (ESS vs participant mean EEG)")
boot_spearman_ci <- function(x, y, n_boot = 2000, seed = 1) {
  set.seed(seed)
  vals <- replicate(n_boot, {
    i <- sample(seq_along(x), replace = TRUE)
    suppressWarnings(cor(x[i], y[i], method = "spearman"))
  })
  quantile(vals, c(0.025, 0.975), na.rm = TRUE)
}
trait <- do.call(rbind, lapply(names(EEG), function(v) {
  ct <- suppressWarnings(cor.test(pl$ESS_avg, pl[[v]], method = "spearman", exact = FALSE))
  ci <- boot_spearman_ci(pl$ESS_avg, pl[[v]])
  data.frame(eeg = EEG[[v]], rho = unname(ct$estimate), ci_low = ci[1], ci_high = ci[2],
             n = nrow(pl), p = ct$p.value, significant_bonf = ct$p.value < ALPHA_BONF)
}))
print(transform(trait, rho = round(rho, 3), ci_low = round(ci_low, 3), ci_high = round(ci_high, 3),
                p = fmt_p(p)), row.names = FALSE)
cat("(95% CI for rho from 2000 bootstrap resamples of participants)\n")
write.csv(trait, file.path(out_dir, "3_trait_spearman.csv"), row.names = FALSE)


# ---- 4. GROUPS: Mann-Whitney U, ESS <= 10 vs > 10 (4.3 / 4.4) -----------
header("4. GROUP COMPARISON - Mann-Whitney U (ESS <= 10 vs ESS > 10)")
groups <- do.call(rbind, lapply(names(EEG), function(v) {
  a <- pl[[v]][pl$ESS_group == "ESS <= 10"]; b <- pl[[v]][pl$ESS_group == "ESS > 10"]
  w <- wilcox.test(b, a, exact = FALSE)
  rank_biserial <- 2 * unname(w$statistic) / (length(a) * length(b)) - 1   # >0 means higher in ESS > 10
  data.frame(eeg = EEG[[v]], n_low = length(a), median_low = median(a),
             n_high = length(b), median_high = median(b),
             U = unname(w$statistic), p = w$p.value, rank_biserial = rank_biserial,
             significant_bonf = w$p.value < ALPHA_BONF)
}))
print(transform(groups, median_low = round(median_low, 3), median_high = round(median_high, 3),
                rank_biserial = round(rank_biserial, 3), p = fmt_p(p)), row.names = FALSE)
cat("(rank-biserial r: effect size; positive = higher in the ESS > 10 group)\n")
write.csv(groups, file.path(out_dir, "4_group_mannwhitney.csv"), row.names = FALSE)


# ---- 5. Supplementary linear mixed-effects model (4.1) ------------------
header("5. LINEAR MIXED-EFFECTS MODEL: log EEG ~ KSS + ESS + session + (1 | participant)")
lmm_rows <- list()
for (v in names(EEG)) {
  f <- as.formula(paste(v, "~ KSS + ESS_avg + session_f + (1 | participant_id)"))
  m <- lmer(f, data = d, REML = TRUE)
  co <- summary(m)$coefficients
  ci <- suppressMessages(confint(m, parm = rownames(co)[-1], method = "Wald"))
  cat("\n---", EEG[[v]], "---\n"); print(round(co, 4))
  vc <- as.data.frame(VarCorr(m))
  icc <- vc$vcov[1] / sum(vc$vcov)
  cat(sprintf("Intraclass correlation (share of variance between participants): %.2f\n", icc))
  for (term in rownames(co)[-1]) {
    lmm_rows[[length(lmm_rows) + 1]] <- data.frame(
      eeg = EEG[[v]], term = term, estimate = co[term, "Estimate"], std_err = co[term, "Std. Error"],
      df = co[term, "df"], t = co[term, "t value"], p = co[term, "Pr(>|t|)"],
      ci_low = ci[term, 1], ci_high = ci[term, 2],
      pct_change = (10^co[term, "Estimate"] - 1) * 100,
      pct_low = (10^ci[term, 1] - 1) * 100, pct_high = (10^ci[term, 2] - 1) * 100)
  }
}
lmm <- do.call(rbind, lmm_rows)
write.csv(lmm, file.path(out_dir, "5_mixed_model.csv"), row.names = FALSE)
cat("\nEstimates are in log10 units; pct_change = % change in the EEG measure per 1-point increase",
    "\n(per KSS point, per ESS point, or session 2/3 vs session 1).\n")


# ---- 6. Comparison: pooled Pearson vs repeated-measures (4.4) -----------
header("6. STATISTICAL APPROACH COMPARISON - pooled Pearson vs repeated-measures")
cmp <- do.call(rbind, lapply(names(EEG), function(v) {
  pc <- cor.test(d$KSS, d[[v]], method = "pearson")
  rr <- state[state$eeg == EEG[[v]], ]
  data.frame(eeg = EEG[[v]],
             pooled_r = unname(pc$estimate), pooled_p = pc$p.value, pooled_n = nrow(d),
             rm_r = rr$r_rm, rm_p = rr$p)
}))
print(transform(cmp, pooled_r = round(pooled_r, 3), rm_r = round(rm_r, 3),
                pooled_p = fmt_p(pooled_p), rm_p = fmt_p(rm_p)), row.names = FALSE)
cat("Pooled Pearson treats all 168 recordings as independent and mixes between- and\n",
    "within-person differences; the repeated-measures correlation uses only within-person change.\n")
write.csv(cmp, file.path(out_dir, "6_pooled_vs_rmcorr.csv"), row.names = FALSE)


# ---- 7. Figures ---------------------------------------------------------
BLUE <- "#2a78d6"; DARK <- "#104281"; ORANGE <- "#eb6834"; LIGHT <- "#cde2fb"; GREY <- "grey60"

# Fig A: each participant's KSS and T/A across sessions (spaghetti plots)
figA <- function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.8, 3, 1), las = 1)
  set.seed(2)
  for (v in c("KSS", "theta_alpha_ratio")) {
    yl <- if (v == "KSS") "KSS (1 = alert, 9 = very sleepy)" else "T/A ratio (log scale)"
    plot(NA, xlim = c(0.8, 3.2), ylim = range(d[[v]]), xaxt = "n", bty = "l",
         log = if (v == "KSS") "" else "y", xlab = "Session", ylab = yl,
         main = if (v == "KSS") "Reported sleepiness by session" else "Theta/alpha ratio by session")
    axis(1, at = 1:3, labels = c("1", "2\n(+90 min)", "3\n(+1 month)"), padj = 0.5)
    for (pid in levels(d$participant_id)) {
      g <- d[d$participant_id == pid, ]; g <- g[order(g$session), ]
      yy <- if (v == "KSS") g[[v]] + runif(3, -0.12, 0.12) else g[[v]]
      lines(g$session, yy, col = adjustcolor(BLUE, 0.25))
    }
    mu <- tapply(d[[v]], d$session, if (v == "KSS") mean else function(x) 10^mean(log10(x)))
    lines(1:3, mu, lwd = 3.5, col = DARK); points(1:3, mu, pch = 16, cex = 1.6, col = DARK)
  }
}

# Fig B: trait - ESS vs participant mean log T/A
figB <- function() {
  par(mar = c(4.5, 4.8, 3.5, 1), las = 1)
  plot(pl$ESS_avg, pl$log_ratio, pch = 16, cex = 1.2, bty = "l",
       col = ifelse(pl$ESS_group == "ESS > 10", ORANGE, BLUE),
       xlab = "Epworth Sleepiness Scale (average of 2 visits)",
       ylab = "Mean log T/A ratio (3 sessions)", main = "Trait sleepiness vs resting T/A ratio")
  abline(v = 10, lty = 2, col = GREY)
  abline(lm(log_ratio ~ ESS_avg, data = pl), lwd = 2, col = DARK)
  tr <- trait[1, ]
  mtext(sprintf("Spearman rho = %.2f, p = %s, n = %d", tr$rho, fmt_p(tr$p), tr$n),
        side = 3, line = 0.2, cex = 0.8, col = "grey30")
  legend("topleft", c("ESS <= 10", "ESS > 10"), pch = 16, col = c(BLUE, ORANGE), bty = "n", cex = 0.85)
}

# Fig C: groups - boxplots of mean log T/A, theta, alpha by ESS group
figC <- function() {
  par(mfrow = c(1, 3), mar = c(4, 4.8, 3.5, 1), las = 1)
  for (v in names(EEG)) {
    boxplot(pl[[v]] ~ pl$ESS_group, col = c(LIGHT, "#f9c9b3"), border = "grey30",
            xlab = "", ylab = sub(" \\(primary\\)", "", EEG[[v]]), main = sub(" \\(primary\\)", "", EEG[[v]]),
            outline = FALSE)
    stripchart(pl[[v]] ~ pl$ESS_group, vertical = TRUE, method = "jitter", add = TRUE,
               pch = 16, col = adjustcolor(c(BLUE, ORANGE), 0.7))
    g <- groups[groups$eeg == EEG[[v]], ]
    mtext(sprintf("Mann-Whitney p = %s", fmt_p(g$p)), side = 3, line = 0.2, cex = 0.7, col = "grey30")
  }
}

# Fig D: forest plot of the mixed model (KSS and ESS effects)
figD <- function() {
  fp <- lmm[lmm$term %in% c("KSS", "ESS_avg"), ]
  fp$measure <- sub(" \\(primary\\)", "", fp$eeg)
  fp <- fp[order(fp$term != "KSS", match(fp$eeg, EEG)), ]
  fp$y <- c(8, 7, 6, 3, 2, 1)
  sig <- fp$p < 0.05; col_sig <- ifelse(sig, DARK, GREY)
  par(mar = c(5.5, 9, 2.5, 2), las = 1)
  xr <- range(c(fp$pct_low, fp$pct_high, 0)); xr <- xr + c(-1, 1) * diff(xr) * 0.08
  plot(NA, xlim = xr, ylim = c(0.4, 9.6), yaxt = "n", bty = "n", ylab = "",
       xlab = "% change per 1-point higher score (95% CI)", main = "Mixed-effects model", mgp = c(2.6, 0.8, 0))
  abline(v = 0, lty = 2, col = "grey50"); abline(h = 4.6, col = "grey90")
  segments(fp$pct_low, fp$y, fp$pct_high, fp$y, lwd = 3, col = col_sig)
  points(fp$pct_change, fp$y, pch = 16, cex = 1.5, col = col_sig)
  text(fp$pct_change, fp$y + 0.35, sprintf("%+.1f%%%s", fp$pct_change, ifelse(sig, "*", "")),
       cex = 0.7, col = col_sig)
  axis(2, at = fp$y, labels = fp$measure, tick = FALSE, cex.axis = 0.85)
  text(xr[1], c(9.1, 4.1), c("KSS (state, per session)", "ESS (trait)"), adj = c(0, 0.5), font = 2, cex = 0.8)
  mtext("* p < 0.05. Adjusted for session; random intercept per participant (56 x 3).",
        side = 1, line = 4.3, cex = 0.65, col = "grey40")
}

save_fig <- function(fun, name, w, h) {
  png(file.path(out_dir, name), width = w, height = h, res = 250); fun(); dev.off()
}
save_fig(figA, "figA_sessions_spaghetti.png", 2800, 1300)
save_fig(figB, "figB_trait_ESS_vs_TA.png",    1800, 1500)
save_fig(figC, "figC_ESS_groups.png",         3000, 1300)
save_fig(figD, "figD_mixed_model_forest.png", 2400, 1600)

header("DONE")
cat("Tables (CSV) and figures (PNG) saved in:\n", out_dir, "\n")
cat("Significance: p < 0.05, and Bonferroni across the 3 EEG measures at p <", round(ALPHA_BONF, 3), "\n")
