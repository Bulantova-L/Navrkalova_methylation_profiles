# Agreement of patient-vs-control methylation changes with the DLBCL direction
# from the literature (GSE54200), as shown on slide 15.
#
# Input: vysledky_pacienti_vs_kontrola.csv (output of DLBCL_metylace.R)
# Definition (same as in DLBCL_metylace.R):
#   among loci with FDR < 0.05 (patient vs. control), a locus "agrees" when
#   sign(delta) == -sign(delta_DLBCL_minus_CB)
#   (the -1 because GSE54200 values are HELP-assay HpaII/MspI ratios,
#    where a higher value means LESS methylation)

library(dplyr)
library(ggplot2)

# run from analyses/07_gse54200_dlbcl
res <- read.csv("output/vysledky_pacienti_vs_kontrola.csv")

# --- optional: drop duplicated probe rows -----------------------------------
# Some probes have several rows (liftOver fragments) and the full_join in
# DLBCL_metylace.R multiplies them. Better fix: aggregate per probe BEFORE the
# join, e.g. in agg_list:
#   df %>% group_by(probe_id) %>% summarise(Nmod = sum(Nmod), Ncov = sum(Ncov))
# Quick fix on the finished table:
dedup <- FALSE
if (dedup) res <- distinct(res, patient, probe_id, .keep_all = TRUE)

# --- agreement per patient ---------------------------------------------------
agree <- res %>%
  filter(fdr < 0.05, !is.na(delta_DLBCL_minus_CB)) %>%
  mutate(agrees = (delta * delta_DLBCL_minus_CB * -1) > 0) %>%
  group_by(patient) %>%
  summarise(n_sig = n(), agreement = 100 * mean(agrees), .groups = "drop") %>%
  mutate(sample = sub("^cfDNA_", "", patient),
         highlight = sample == "Pilot1")

print(agree)
# expected without dedup: BAR 47.8, MUD 55.2, Pilot1 64.9

# --- plot (styled like the slide) -------------------------------------------
p <- ggplot(agree, aes(x = sample, y = agreement, fill = highlight)) +
  geom_col(width = 0.62) +
  geom_hline(aes(yintercept = 50, linetype = "Chance (50 %)"),
             colour = "#C8452C", linewidth = 0.7) +
  geom_text(aes(y = 4, label = sprintf("%.0f%%", agreement)),
            vjust = 0, colour = "white", fontface = "bold", size = 5) +
  scale_fill_manual(values = c(`FALSE` = "#A7B3B8", `TRUE` = "#E07A1F"),
                    guide = "none") +
  scale_linetype_manual(name = NULL, values = c("Chance (50 %)" = "dashed")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Agreement with DLBCL direction",
       subtitle = "% of significant loci (FDR < 0.05) changed in the same\ndirection as in DLBCL (GSE54200)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", colour = "#0F2B33", size = 14),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "#5F6B76", size = 10),
        axis.text = element_text(colour = "#5F6B76"),
        legend.position = "bottom")

ggsave("plots/agreement_barplot.png", p, width = 5, height = 5, dpi = 300, bg = "white")
