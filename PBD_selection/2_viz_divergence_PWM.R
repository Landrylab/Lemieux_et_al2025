# Visualize the average AAD matrices for one domain group
library(tidyverse)
library(ggdendro)
library(tidyr)
library(grid)

################################
## Variables to modify ##
# directory of the repository with the output file of divergence_PWM.R
path = '~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/PBD_selection/'
output_dir = './'

# SH3, PDZ or WW
domain_group = 'SH3'
################################

aad_pwm <- read_rds(paste0(path, 'AAD_',domain_group, '.rds'))

# change format of the aad matrix
df_aad <- as_tibble(aad_pwm)
df_aad$PRM <- colnames(aad_pwm)

df_aad <- pivot_longer(df_aad, cols = 1:nrow(df_aad), values_to = 'mean_AAD')

ref_matrix <- aad_pwm

# Prepare the dendrogram describing the relations between all PWM
otter_matrix <- as.matrix(ref_matrix)[, -ncol(ref_matrix)]

rownames(otter_matrix)[1:(ncol(ref_matrix)-1)] <- colnames(ref_matrix)[-ncol(ref_matrix)]
otter_dendro <- as.dendrogram(hclust(d = dist(x = otter_matrix)))

otter_order <- order.dendrogram(otter_dendro)
df_aad$PRM <- factor(x = df_aad$PRM,
                      levels = row.names(otter_matrix)[otter_order], 
                      ordered = TRUE)

df_aad$name <- factor(x = df_aad$name,
                       levels = row.names(otter_matrix)[otter_order], 
                       ordered = TRUE)

# Create the dendro plot
dendro_plot <- ggdendrogram(data = otter_dendro, rotate = TRUE, labels = FALSE)+
  theme(axis.text.y = element_text(size = 0, colour = 'black'),
        axis.text.x = element_text(color = 'black', size = 8))

# Preview the plot
print(dendro_plot)

# Plot the minimal average aad between each comparison
df_aad %>%
  ggplot( aes(x = name, y = PRM)) +
  geom_tile(aes(fill = mean_AAD)) +
  scale_fill_viridis_c(limits = c(0,0.1)) +
  xlab(paste(domain_group, 'family PWM'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.5, hjust =1, color = 'black'),
        legend.position = "bottom", legend.title.position = 'top')+
  guides(fill = guide_colorbar(title = 'minimal average AAD',
                               theme = theme(
                                 legend.key.width  = unit(10, "lines"),
                                 legend.key.height = unit(1.5, "lines")
                               ))) -> heatmap_plot

# Put together the heatmap and dendrogram

library(cowplot)
library(svglite)
plot_grid(heatmap_plot, dendro_plot, rel_widths = c(1,0.3), 
          align = 'h')

ggsave(paste0('~/PL_projects/PL_papers/PPI_optimization_paper/figures/Fig_PWM_', domain_group, '.svg'), height = 10, width =18)
