# network_analysis.R
# ──────────────────────────────────────────────────────────────────────────────
# Crime Severity Network — Chile 2017
# Perceived severity of 15 crimes in a Chilean community sample (n ≈ 274)
#
# Produces all manuscript figures:
#   Figure 1  — Descriptive boxplots (ranking, years, seriousness)
#   Figure 2  — Primary network (EBICglasso, polychoric, gravedad)
#   Figure 3  — Sensitivity network (EBICglasso, Spearman, composite)
#   Figure 4  — Centrality plot (strength only; closeness/betweenness unstable)
#   Figure S1 — Bootstrap stability (supplementary)
#   Figure S2 — Edge accuracy bootstrap (supplementary)
#
# Authors: Herman E. Elgueta, Beatriz Pérez Sánchez
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(bootnet)
library(qgraph)
library(igraph)
library(psych)
library(scales)
library(gridExtra)

set.seed(415)   # RANDOM.ORG seed used in original data processing
dir.create("../outputs", showWarnings = FALSE, recursive = TRUE)

# ── Community palette ─────────────────────────────────────────────────────────
# 3-community solution: Walktrap on primary (gravedad) network
#   Contested    — Abortion, Euthanasia
#   Violence     — Rape, Child Sexual Abuse, Child Abuse, Partner Violence
#   Conventional — all remaining crimes (single cluster)
COMM_COLORS_3 <- c(
  "Contested"    = "#55A868",
  "Violence"     = "#DD8452",
  "Conventional" = "#4C72B0"
)

# 4-community solution: Walktrap on composite + Louvain/Spinglass on primary
#   Same first two, conventional cluster splits into:
#   Economic     — Fraud, Tax Evasion, Piracy, Bribery
#   ViolentPub   — Murder, Robbery, Terrorism, Drug Trafficking, Vandalism
COMM_COLORS_4 <- c(
  "Contested"   = "#55A868",   # green  — same as Figs 1 & 2
  "Violence"    = "#DD8452",   # orange — same as Figs 1 & 2
  "Economic"    = "#8172B3",   # purple — new (split from Conventional)
  "ViolentPub"  = "#64B5CD"    # teal   — new (split from Conventional)
)

# Keep COMM_COLORS pointing to 3-community palette (used for Fig 2 / Fig 1)
COMM_COLORS <- COMM_COLORS_3

community_of <- function(crime) {
  case_when(
    crime %in% c("Abortion", "Euthanasia")                                         ~ "Contested",
    crime %in% c("Rape", "Child Sexual Abuse", "Child Abuse", "Partner Violence")  ~ "Violence",
    TRUE                                                                             ~ "Conventional"
  )
}

community_of_4 <- function(crime) {
  case_when(
    crime %in% c("Abortion", "Euthanasia")                                         ~ "Contested",
    crime %in% c("Rape", "Child Sexual Abuse", "Child Abuse", "Partner Violence")  ~ "Violence",
    crime %in% c("Fraud", "Tax Evasion", "Piracy", "Bribery")                     ~ "Economic",
    TRUE                                                                             ~ "ViolentPub"
  )
}

# ── Load pre-computed correlation matrices (publicly available) ───────────────
# These are sufficient to reproduce Figures 2, 3, and 4 exactly.
# Individual-level data are not distributed (ethics restriction).
poly_mat   <- as.matrix(read.csv("../data/polychoric_gravedad.csv",  row.names = 1))
sp_mat     <- as.matrix(read.csv("../data/spearman_composite.csv",   row.names = 1))
NET_LABELS <- rownames(poly_mat)
N_SAMPLE   <- 274L   # analytic sample size used to estimate networks

cat(sprintf("Correlation matrices loaded: %d crimes\n", length(NET_LABELS)))

make_groups <- function(labels) {
  comm <- community_of(labels)
  list(
    Contested    = which(comm == "Contested"),
    Violence     = which(comm == "Violence"),
    Conventional = which(comm == "Conventional")
  )
}

# Layout helper: igraph Fruchterman-Reingold is truly stochastic (respects
# set.seed); qgraph's own spring layout is deterministic when fed a weight
# matrix, so we compute coords here and pass them in.
make_layout <- function(wt_mat, seed = 7) {
  g <- graph_from_adjacency_matrix(abs(wt_mat), mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  set.seed(seed)
  coords <- layout_with_fr(g)
  coords[, 1] <- (coords[, 1] - mean(coords[, 1])) /
                  max(abs(coords[, 1] - mean(coords[, 1])))
  coords[, 2] <- (coords[, 2] - mean(coords[, 2])) /
                  max(abs(coords[, 2] - mean(coords[, 2])))
  coords
}

# ── Figure 2: Primary network (polychoric EBICglasso) ─────────────────────────
net_g_wt  <- EBICglasso(poly_mat, n = N_SAMPLE, gamma = 0.5)
n_edges_g <- sum(net_g_wt != 0) / 2
groups_g  <- make_groups(NET_LABELS)
cat(sprintf("Primary network: %d edges retained\n", n_edges_g))

png("../outputs/fig2_network_gravedad.png",
    width = 1400, height = 1200, res = 160)
q_g <- qgraph(
  net_g_wt,
  layout      = make_layout(net_g_wt, seed = 42),
  labels      = NET_LABELS,
  label.scale = FALSE,
  label.cex   = 0.85,
  vsize       = 8,
  color       = unname(COMM_COLORS[c("Contested","Violence","Conventional")]),
  groups      = groups_g,
  legend      = FALSE,
  title       = "",
  posCol      = "#2166AC",
  negCol      = "#D6604D",
  edge.width  = 1.5,
  repulsion   = 0.9
)
dev.off()
cat("Saved: fig2_network_gravedad.png\n")

# ── Figure 3: Sensitivity network (Spearman EBICglasso on composite) ──────────
# Walktrap detects 4 communities: conventional cluster splits into economic
# (Fraud, Tax Evasion, Piracy, Bribery) and violent/public-order crime.
net_comp_wt <- EBICglasso(sp_mat, n = N_SAMPLE, gamma = 0.5)

g_comp_ig   <- graph_from_adjacency_matrix(
  pmax(net_comp_wt, 0), mode = "undirected", weighted = TRUE, diag = FALSE
)
n_comm_comp <- max(membership(cluster_walktrap(g_comp_ig)))

groups_comp_4 <- list(
  Contested  = which(community_of_4(NET_LABELS) == "Contested"),
  Violence   = which(community_of_4(NET_LABELS) == "Violence"),
  Economic   = which(community_of_4(NET_LABELS) == "Economic"),
  ViolentPub = which(community_of_4(NET_LABELS) == "ViolentPub")
)

cat(sprintf("Sensitivity network: %d edges, walktrap detects %d communities\n",
            sum(net_comp_wt != 0) / 2, n_comm_comp))

png("../outputs/fig3_network_composite.png",
    width = 1400, height = 1200, res = 160)
q_comp <- qgraph(
  net_comp_wt,
  layout      = make_layout(net_comp_wt, seed = 99),
  labels      = NET_LABELS,
  label.scale = FALSE,
  label.cex   = 0.85,
  vsize       = 8,
  color       = unname(COMM_COLORS_4[c("Contested","Violence","Economic","ViolentPub")]),
  groups      = groups_comp_4,
  legend      = FALSE,
  title       = "",
  posCol      = "#2166AC",
  negCol      = "#D6604D",
  edge.width  = 1.5,
  repulsion   = 0.9
)
dev.off()
cat("Saved: fig3_network_composite.png\n")

# ── Figure 4: Centrality ──────────────────────────────────────────────────────
cent     <- centrality(q_g)
cat("\n=== Strength centrality (gravedad network) ===\n")
str_vals <- sort(cent$OutDegree, decreasing = TRUE)
print(round(str_vals, 3))

png("../outputs/fig4_centrality.png",
    width = 900, height = 800, res = 160)
centralityPlot(q_g, include = "Strength", orderBy = "Strength", scale = "z-scores")
dev.off()
cat("Saved: fig4_centrality.png\n")

cat("\n=== Figures 2, 3, 4 complete (reproducible from correlation matrices) ===\n")

# ════════════════════════════════════════════════════════════════════════════════
# The following sections require individual-level data, which are not publicly
# distributed. Figure 1 (descriptive boxplots) and Figures S1–S2 (bootstrap)
# are skipped if the raw data file is not present locally.
# ════════════════════════════════════════════════════════════════════════════════
# Full cleaned dataset with _seriousness, _rank, _years columns for all crimes.
# Not publicly distributed — available from the corresponding author.
RAW_FILE <- "../data/crime_severity_chile2017.csv"
if (!file.exists(RAW_FILE)) {
  cat("\nFull cleaned dataset not found — skipping Figure 1 and bootstrap figures.\n")
  cat("Request data from the corresponding author (herman.elgueta@umag.cl).\n")
  quit(save = "no")
}

# ── Load raw data (full cleaned dataset with all three measures) ──────────────
# Figure 1 requires a dataset with columns named {Crime}_seriousness,
# {Crime}_rank, and {Crime}_years — the complete cleaned analytic file.
# Available from the corresponding author under a data use agreement.
df_clean <- read.csv(RAW_FILE)

# bootnet object for bootstrap figures
g_cols <- grep("_seriousness$", names(df_clean), value = TRUE)
dat_g  <- df_clean[, g_cols]
colnames(dat_g) <- sub("_seriousness$", "", colnames(dat_g)) |>
  gsub("\\.", " ", x = _)

net_g_boot <- estimateNetwork(dat_g, default = "EBICglasso",
                               corMethod = "cor_auto", tuning = 0.5, verbose = FALSE)

# ── Figure 1: Descriptive boxplots (ranking, sentencing, seriousness) ────────
make_long <- function(suffix) {
  cols <- grep(paste0("_", suffix, "$"), names(df_clean), value = TRUE)
  df_clean[, cols] |>
    setNames(sub(paste0("_", suffix, "$"), "", cols) |> gsub("\\.", " ", x = _)) |>
    pivot_longer(everything(), names_to = "crime", values_to = "score") |>
    mutate(community = community_of(crime))
}

rank_long  <- make_long("rank")
years_long <- make_long("years")
ser_long   <- make_long("seriousness")

crime_order <- rank_long |>
  group_by(crime) |>
  summarise(med = median(score, na.rm = TRUE), .groups = "drop") |>
  arrange(med) |>
  pull(crime)

pal <- setNames(COMM_COLORS[community_of(crime_order)], crime_order)

make_summary <- function(df) {
  df |>
    group_by(crime) |>
    summarise(m  = mean(score, na.rm = TRUE),
              se = sd(score, na.rm = TRUE) / sqrt(sum(!is.na(score))),
              .groups = "drop") |>
    mutate(crime = factor(crime, levels = crime_order))
}

rank_sum  <- make_summary(rank_long)
years_sum <- make_summary(years_long)
ser_sum   <- make_summary(ser_long)

bp_theme <- theme_linedraw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9),
        legend.position = "none")

make_bp <- function(long_df, sum_df, x_scale, x_lab, show_y = TRUE) {
  p <- long_df |>
    mutate(crime = factor(crime, levels = crime_order)) |>
    ggplot(aes(x = score, y = crime, fill = crime, color = crime)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.30, linewidth = 0.4) +
    geom_jitter(width = 0, height = 0.18, size = 0.45,
                alpha = 0.35, shape = 16, na.rm = TRUE) +
    geom_pointrange(
      data = sum_df,
      aes(x = m, xmin = m - 2*se, xmax = m + 2*se, y = crime, color = crime),
      size = 0.35, linewidth = 0.9, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    x_scale +
    labs(x = x_lab, y = NULL) +
    bp_theme
  if (!show_y) p <- p + theme(axis.text.y = element_blank())
  p
}

plot_rank <- make_bp(
  rank_long, rank_sum,
  scale_x_reverse(breaks = c(2, 4, 6, 8, 10, 12, 14)),
  "Ranking (1 = most severe)"
)

plot_years <- make_bp(
  years_long, years_sum,
  scale_x_continuous(
    breaks = c(0, 15, 30, 45, 60, 75, 90),
    labels = c("0", "15", "30", "45", "60", "75", "≥90"),
    limits = c(NA, 99)
  ),
  "Years in prison", show_y = FALSE
)

plot_ser <- make_bp(
  ser_long, ser_sum,
  scale_x_continuous(breaks = 0:4),
  "Seriousness (0–4)", show_y = FALSE
)

fig1 <- arrangeGrob(plot_rank, plot_years, plot_ser, ncol = 3,
                    widths = c(5, 4, 4))
ggsave("../outputs/fig1_descriptive.png", fig1,
       width = 12, height = 6, dpi = 200)
cat("Saved: fig1_descriptive.png\n")

# ── Figure S1: Bootstrap stability (CS-coefficients) ─────────────────────────
cat("\nRunning case-dropping bootstrap (n = 1000) — this takes a few minutes...\n")
boot_case <- bootnet(net_g_boot, nBoots = 1000, type = "case",
                     verbose = FALSE, nCores = 1)

safe_cs <- function(boot, stat) {
  r <- tryCatch(corStability(boot, statistics = stat, verbose = FALSE),
                error = function(e) NA_real_)
  if (length(r) == 0) NA_real_ else r[[1]]
}

cs_strength <- safe_cs(boot_case, "strength")
cat(sprintf("CS-coefficient  Strength: %.3f\n", cs_strength))

png("../outputs/figS1_bootstrap_stability.png", width = 1200, height = 800, res = 160)
tryCatch(
  plot(boot_case, statistics = "strength"),
  error = function(e) {
    plot.new(); text(0.5, 0.5, paste("Bootstrap plot error:\n", e$message), cex = 0.8)
  }
)
dev.off()
cat("Saved: figS1_bootstrap_stability.png\n")

# ── Figure S2: Edge accuracy bootstrap ────────────────────────────────────────
cat("\nRunning nonparametric edge bootstrap (n = 1000)...\n")
boot_np <- bootnet(net_g_boot, nBoots = 1000, type = "nonparametric",
                   verbose = FALSE, nCores = 1)
png("../outputs/figS2_edge_bootstrap.png", width = 1400, height = 1000, res = 160)
plot(boot_np, labels = FALSE, order = "sample")
dev.off()
cat("Saved: figS2_edge_bootstrap.png\n")

cat("\n=== Done. All figures saved to ../outputs/ ===\n")
