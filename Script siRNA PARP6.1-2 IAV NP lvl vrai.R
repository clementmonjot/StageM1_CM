## --- SETUP ---

# Working directory
getwd()
setwd("./rawdata/")

# read arg



# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(car)

# Theme/style perso graphique
theme_unique_art <- function (base_size = 18, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(colour = "black"),
      title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.position = "none",
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# Options utilisation
control_sample <- "Control"
target_gene <- "NP"  # Modifier si besoin
#targets_to_exclude <- c("5h si6.1", "5h si6.2", "8h si6.1", "8h si6.2", "12h si6.1", "12h si6.2", "24h si6.1", "24h si6.2")
#targets_to_exclude <- c("Negative 8.1", "Negative 8.2", "RIG-I Control", "RIG-I 8.1", "RIG-I 8.2")
#targets_to_exclude <- c("Negative 8.1", "Negative 8.2", "cGAS Control", "cGAS 8.1", "cGAS 8.2")
targets_to_exclude <- c("")


## Importer données et nettoyage pour analyse

# lecture CSV (csv2 car format numbers français (,))
tableau <- read.csv2("Infection.csv")

# On garde GAPDH (ménage) et gene target (variable optionnelle)
tableau <- tableau %>% filter(X %in% c("Target Name", "GAPDH", target_gene))
colnames(tableau) <- tableau[1, ]
tableau <- tableau %>% filter(Well != "Well")
tableau <- tableau %>% filter(grepl("", `Sample Name`))

# Sélection de colonnes à analyser
tableauCT <- tableau %>% select("Sample Name", "Target Name", "Cт")
colnames(tableauCT) <- c("Sample", "Target", "CT")

# Variable tableauCT pour calcul Rq inside
tableauCT$CT <- gsub(",", ".", tableauCT$CT)
tableauCT$CT <- as.numeric(tableauCT$CT)

# Exclusion de taget_to_exclude (sinon peut avoir problème si pas assez de sample)
if (length(targets_to_exclude) > 0) {
  tableauCT <- tableauCT %>% filter(!(Sample %in% targets_to_exclude))
}

# Pivoter la table  et supprimer les NTC (s'il y en a)
tableauCTcast <- tableauCT %>% pivot_wider(names_from = Target, values_from = CT, values_fn = list)
tableauCTcast <- tableauCTcast %>% unnest(cols = everything()) %>% filter(Sample != "NTC")

# Calcul ∆CT, ∆∆CT, Rq
tableauCTmelt <- tableauCTcast %>% pivot_longer(cols = all_of(target_gene))
tableauCTmelt$DCT <- tableauCTmelt$value - tableauCTmelt$GAPDH
ref_target_gene <- as.numeric(tableauCTmelt[1, "DCT"])
tableauCTmelt$DDCT <- tableauCTmelt$DCT - ref_target_gene
tableauCTmelt$Rq <- 2^(-tableauCTmelt$DDCT)

# Résumé pour barplot (Regarder vidéo des ggPlots)
tableaueasy <- tableauCTmelt %>% select(Sample, name, Rq) %>% filter(!is.na(Sample), !is.na(Rq))
tableaueasy$Sample <- factor(tableaueasy$Sample, levels = unique(tableaueasy$Sample))
tableauPlot <- tableaueasy %>% group_by(Sample, name) %>% summarise(mean = mean(Rq), sd = sd(Rq), .groups = "drop")


## PARTIE STATS (Test si ANOVA possible)
anova_model <- aov(Rq ~ Sample, data = tableaueasy)
levene_res <- leveneTest(Rq ~ Sample, data = tableaueasy)
shapiro_res <- shapiro.test(residuals(anova_model))

text_signif <- ""

if (levene_res$`Pr(>F)`[1] > 0.05 && shapiro_res$p.value > 0.05) {
  message("ANOVA valide : ANOVA + Test de Tukey")
  
  tukey_res <- TukeyHSD(anova_model)
  pvals_df <- as.data.frame(tukey_res$Sample)
  pvals_df$comparison <- rownames(pvals_df)
  
  pvals_df <- pvals_df %>%
    mutate(group1 = sub("-.*", "", comparison),
           group2 = sub(".*-", "", comparison)) %>%
    filter(group1 == control_sample | group2 == control_sample) %>%
    filter(`p adj` < 0.05) %>%
    mutate(condition = ifelse(group1 == control_sample, group2, group1),
           p.signif = case_when(
             `p adj` <= 0.0001 ~ "****",
             `p adj` <= 0.001 ~ "***",
             `p adj` <= 0.01 ~ "**",
             `p adj` <= 0.05 ~ "*"
           ))
  
  text_signif <- ifelse(nrow(pvals_df) > 0, "Significatif", "Pas de différence significative")
  
  if (nrow(pvals_df) > 0) {
    pvals_df$condition <- factor(pvals_df$condition, levels = levels(tableaueasy$Sample))
    pvals_df$name <- unique(tableaueasy$name)[1]
    
    max_vals <- tableauPlot %>% group_by(Sample) %>%
      summarise(y.max = max(mean + sd), .groups = "drop")
    
    pvals_df <- pvals_df %>%
      left_join(max_vals, by = c("condition" = "Sample")) %>%
      mutate(y.position = y.max + 0.3)
  }
  
tableauPlot <- tableauPlot[!is.na(tableauPlot$mean), ]
  
  
  p <- ggplot(tableauPlot, aes(x = Sample, y = mean)) +
    geom_bar(stat = "identity", fill = "darkred", color = "black", width = 0.35) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    geom_point(data = tableaueasy, aes(x = Sample, y = Rq),
               inherit.aes = FALSE,
               position = position_jitter(width = 0.05),
               shape = 21, fill = "darkgray", color = "black", size = 2, stroke = 0.2) +
    facet_grid(~name) +
    (if (exists("pvals_df") && nrow(pvals_df) > 0)
      geom_text(data = pvals_df, aes(x = condition, y = y.position, label = p.signif),
                inherit.aes = FALSE, size = 6, fontface = "bold", vjust = 0)
     else NULL) +
    labs(y = "Relative expression ± SD", x = "Condition") +
    theme_unique_art() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # << C'est cette ligne qui change l'angle
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    )
  
  print(p)
  
} else {
  message("ANOVA invalide : Test de Wilcoxon")
  
  autres_conditions <- setdiff(levels(tableaueasy$Sample), control_sample)
  wilcox_results <- data.frame(condition = character(), p.value = numeric())
  
    for (cond in autres_conditions) {
      sub_data <- tableaueasy %>% filter(Sample %in% c(control_sample, cond))
      
      # Vérifie que deux niveaux sont bien présents
      if (length(unique(sub_data$Sample)) == 2 && all(table(sub_data$Sample) >= 1)) {
        test <- wilcox.test(Rq ~ Sample, data = sub_data)
        wilcox_results <- rbind(wilcox_results, data.frame(condition = cond, p.value = test$p.value))
      } else {
        message(paste("⚠️ Condition ignorée (pas 2 niveaux) :", cond))
      }
    }
    
  
  wilcox_results <- wilcox_results %>%
    mutate(signif = case_when(
      p.value <= 0.0001 ~ "****",
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  
  text_signif <- ifelse(any(wilcox_results$p.value < 0.05), "Significatif", "Pas de différence significative")
  
  max_vals <- tableauPlot %>%
    group_by(Sample) %>%
    summarise(y.max = max(mean + sd, na.rm = TRUE), .groups = "drop")
  
  annotations <- wilcox_results %>%
    filter(signif != "ns") %>%
    mutate(y.position = max_vals$y.max[match(condition, max_vals$Sample)] + 0.3,
           name = unique(tableaueasy$name)[1])
  
  p <- ggplot(tableauPlot, aes(x = Sample, y = mean)) +
    geom_bar(stat = "identity", fill = "darkred", color = "black", width = 0.35) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    geom_point(data = tableaueasy, aes(x = Sample, y = Rq),
               inherit.aes = FALSE,
               position = position_jitter(width = 0.05),
               shape = 21, fill = "darkgray", color = "black", size = 2, stroke = 0.2) +
    facet_grid(~name) +
    (if (exists("annotations") && nrow(annotations) > 0)
      geom_text(data = annotations, aes(x = condition, y = y.position, label = signif),
                inherit.aes = FALSE, size = 6, fontface = "bold", vjust = 0)
     else NULL) +
    labs(y = "Relative expression ± SD", x = "Condition") +
    theme_unique_art() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # << C'est cette ligne qui change l'angle
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    )
  print(p)
}


#ggsave("Timedep 12h PARP8.png", plot = p, width = 8, height = 6, dpi = 600)



