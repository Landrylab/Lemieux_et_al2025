#import packages
library(ggplot2)
library(ggpubr)
library(magrittr)
library(stringr)
library(tidyverse)
library(gtools)
library(readxl)
library(cowplot)
library(svglite)

# Set theme for graphs
gt <- theme_classic2()+
  theme(axis.title = element_text(size = 10, color = 'black'), 
        axis.text = element_text(size = 8, color = 'black'), 
        axis.line = element_line(linewidth = 0.4), 
        axis.ticks = element_line(linewidth = 0.8),
        legend.title = element_text(size = 10, color = 'black'), 
        legend.text = element_text(size = 8, color = 'black'))

# Import data
df_growth <- read_csv('~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/growthcurves_analysis/gc_data/condensed_availability.csv')[, -1]

# Remove wells not inoculated
df_growth <- df_growth[df_growth$strain != 'control', ]

# Compute mean AUC and dgr for each strain/condition
df_growth %<>%
  dplyr::group_by(D12, D3, condition, strain) %>%
  dplyr::mutate(med_auc = mean(auc), 
         med_dgr = mean(dgr))

# Prepare data for visualization
df_viz <- df_growth

# Define protein binding domains and peptides
PBD <- c('PDZ',  'SH3', 'GBD','GRB2')
peptide <- c('pdz.1', 'pdz.2', 'pdz.3', 'sh3.1', 'sh3.2', 'sh3.3', 'sh3.4', 'gbd1', 'gab2')

# Classify orientation (TRUE if neither D12 nor D3 is PBD/peptide)
df_viz$orientation <- !(df_viz$D12 %in% peptide) & !(df_viz$D3 %in% PBD)

# Assign protein labels based on D12/D3
df_viz[ df_viz$D12 %in% PBD, 'protein'] <- df_viz[ df_viz$D12 %in% PBD, 'D12' ]
df_viz[ df_viz$D3 %in% PBD, 'protein'] <- df_viz[ df_viz$D3 %in% PBD, 'D3' ]
df_viz[ df_viz$D12 %in% peptide, 'protein'] <- df_viz[ df_viz$D12 %in% peptide, 'D12' ]
df_viz[ df_viz$D3 %in% peptide, 'protein'] <- df_viz[ df_viz$D3 %in% peptide, 'D3' ]

# Replace NA protein labels with 'empty'
df_viz$protein <- na.replace(df_viz$protein, 'empty')

# Set orientation factor labels
df_viz$orientation <- 
  factor(df_viz$orientation, 
         levels = c(TRUE, FALSE), 
         labels = c('PBD-F[1,2]_peptide-F[3]', 
                    'peptide-F[1,2]_PBD-F[3]'))

# Assign group labels for visualization
df_viz[df_viz$protein == 'SH3' | df_viz$protein %in% c('sh3.1', 'sh3.2', 'sh3.3', 'sh3.4'), 'group'] <- 'SH3'
df_viz[df_viz$protein == 'PDZ' | df_viz$protein %in% c('pdz.1', 'pdz.2', 'pdz.3'), 'group'] <- 'PDZ'
df_viz[df_viz$protein == 'GBD' | df_viz$protein == 'gbd1', 'group'] <- 'GBD'
df_viz[df_viz$protein == 'GRB2' | df_viz$protein == 'gab2', 'group'] <- 'GRB2'

# Replace NA group labels with '(+)'
df_viz$group <- na.replace(df_viz$group, replace = '(+)')

# Set group factor levels
df_viz$group <- 
  factor(df_viz$group, 
         levels = c('(+)', 'GBD', 'GRB2', 'PDZ', 'SH3'))

# Create Figure S2: visualize AUC for replicates by protein/group
FigS2 <- 
  ggplot(df_viz[df_viz$condition =='MTX' , ])+
  facet_grid(cols = vars(group), scales = 'free_x', space = 'free_x')+
  geom_jitter(aes(x = protein, y = auc, color = strain), 
              alpha = 0.8, position = position_jitterdodge(jitter.width = 0.2,  dodge.width = 0.4), size = 1.5)+
  stat_compare_means(aes(x =protein, y = auc, group = strain), 
                     hide.ns = TRUE, method = 't.test', label = 'p.signif', label.y = 210)+
  scale_color_manual(values = c('#51C3CC', '#FF8E32'), 
                     labels = c('F[3]', 'F[1,2]'))+
  ylab('PCA signal (corrected AUC)')+
  xlab('Protein')+
  theme_classic2()+
  gt+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'bottom')+
  guides(color = guide_legend(title = ''))

# Keep median information for further analysis
med_growth <- df_viz %>%
  dplyr::select(D12, D3, med_auc, med_dgr, orientation, group, condition, strain, protein) %>%
  unique()

# Filter for MTX condition only
med_growth <- med_growth[med_growth$condition == 'MTX', ]

# Get negative control (empty_empty plasmid) values
neg_growth <- med_growth[med_growth$D12 =='empty' & med_growth$D3 =='empty' & med_growth$condition =='MTX', ]

D12_neg <- unlist(neg_growth[neg_growth$strain=='BY4741-mVenus-D3', "med_auc" ])
D3_neg <- unlist(neg_growth[neg_growth$strain=='BY4741-mVenus-D12', "med_auc" ])

# Compute AUC ratio vs negative control for each strain
med_growth[med_growth$strain == 'BY4741-mVenus-D3', 'ratio_auc'] <- med_growth[med_growth$strain == 'BY4741-mVenus-D3', "med_auc"]/D12_neg
med_growth[med_growth$strain == 'BY4741-mVenus-D12', 'ratio_auc'] <- med_growth[med_growth$strain == 'BY4741-mVenus-D12', "med_auc"]/D3_neg

D12_neg <- unlist(neg_growth[neg_growth$strain=='BY4741-mVenus-D3', "med_dgr" ])
D3_neg <- unlist(neg_growth[neg_growth$strain=='BY4741-mVenus-D12', "med_dgr" ])

# Compute dgr ratio vs negative control for each strain
med_growth[med_growth$strain == 'BY4741-mVenus-D3', 'ratio_dgr'] <- med_growth[med_growth$strain == 'BY4741-mVenus-D3', "med_dgr"]/D12_neg
med_growth[med_growth$strain == 'BY4741-mVenus-D12', 'ratio_dgr'] <- med_growth[med_growth$strain == 'BY4741-mVenus-D12', "med_dgr"]/D3_neg

# Set DHFR fragment factor
med_growth$DHFR <- 
  factor(med_growth$strain, 
         levels = c('BY4741-mVenus-D3', 'BY4741-mVenus-D12'), 
         labels = c('F[1,2]', 'F[3]'))

# Assign type (PBD/peptide/empty) for visualization
med_growth[med_growth$D12 %in% PBD | med_growth$D3 %in% PBD, 'type'] <- 'PBD'
med_growth[med_growth$D12 %in% peptide | med_growth$D3 %in% peptide, 'type'] <- 'peptide'
med_growth$type <- na.replace(med_growth$type, '(+)')

# Visualize dgr ratios (main figure panel)
p2 <- 
  ggplot(med_growth)+
  facet_grid(cols = vars(type), scale = 'free_x', space = 'free')+
  geom_tile(aes(protein, DHFR, fill = ratio_dgr))+
  scale_fill_viridis_c(option = 'E', direction = 1)+
  theme_classic2()+
  xlab('Protein')+
  gt+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), 
        legend.position = 'right', 
        legend.title = element_text(size = 8))+
  guides(fill = guide_colorbar(title = 'Availability \nscore',
                               theme = theme(
                                 legend.key.width  = unit(1, "lines"),
                                 legend.key.height = unit(8, "lines")
                               )))

# Create figure 2 panel A (image) and combine with panel B (tile plot)
p1 <- 
  ggdraw()+
  draw_image('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig2A.png')

Fig2 <- 
  plot_grid(p1, p2, rel_widths = c(1, 2), labels = c('A', 'B'), 
            label_fontface = 'plain')

# Save figures
ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig2.svg', Fig2,
       height = 2.5, width = 6.8)

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/final_S1/FigS3.svg', 
       FigS2, height = 3, width = 7)
