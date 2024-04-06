#libraries
library(tidyverse)
library(ggpubr)

#load data
patch_size <- read.table("patch_area.txt", header = T, sep = "\t")
site_data <-
  read.table("extracted_site_level_data.txt",
             header = T,
             sep = "\t")

monad_data <- subset(site_data, site_data$MDT == 1)

#combine data to 1 spreadsheet
patch_site_data <- bind_cols(monad_data, patch_size)

#plot plant richness

plant_richness_plot <- ggplot(data = patch_site_data) +
  geom_point(mapping = aes (x = HabPat_Ha, y = plant_rich, colour = Name)) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Habitat patch area (Ha)", y = "Plant species richness",
       colour = "Habitat") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )


#plot floral abundance
floral_abundance_plot <- ggplot(data = patch_site_data) +
  geom_point(mapping = aes (x = HabPat_Ha, y = openFU_abundance, colour = Name)) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(x = "Habitat patch area (Ha)", y = "Floral abundance",
       colour = "Habitat") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

#plot insect richness
insect_species_richness <- ggplot(data = patch_site_data) +
  geom_point(mapping = aes (x = HabPat_Ha, y = total_insect_sp_rich, colour = Name)) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Habitat patch area (Ha)", y = "Insect species richness",
       colour = "Habitat") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

#plot insect abundance

insect_abundance_plot <- ggplot(data = patch_site_data) +
  geom_point(mapping = aes (x = HabPat_Ha, y = total_insect_sp_abund, colour = Name)) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Habitat patch area (Ha)", y = "Insect abundance",
       colour = "Habitat") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )


#plot species evenness
species_evenness_plot <- ggplot(data = patch_site_data) +
  geom_point(mapping = aes (x = HabPat_Ha, y = species_even_insects, colour = Name)) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Habitat patch area (Ha)", y = "Insect species evenness",
       colour = "Habitat") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )


#plot interaction evenness
interaction_evenness_plot <- ggplot(data = patch_site_data) +
  geom_point(mapping = aes (x = HabPat_Ha, y = IE, colour = Name)) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Habitat patch area (Ha)", y = "Interaction evenness",
       colour = "Habitat") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

#all plots combined

pdf("Effect of surrounding habitat size.pdf")

ggarrange(
  plant_richness_plot,
  floral_abundance_plot,
  insect_species_richness,
  insect_abundance_plot,
  species_evenness_plot,
  interaction_evenness_plot,
  ncol = 2,
  nrow = 3,
  common.legend = T,
  legend = "right"
)
dev.off()
