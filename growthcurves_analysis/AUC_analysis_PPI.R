# Analysis of the PCA DHFR assay to quantify interaction strength

# Import required packages
library(ggplot2)
library(ggpubr)
library(magrittr)
library(stringr)
library(tidyverse)
library(gtools)
library(readxl)
library(cowplot)

# Set theme for graphs
gt <- theme_classic2()+
  theme(axis.title = element_text(size = 10, color = 'black'), 
        axis.text = element_text(size = 8, color = 'black'), 
        axis.line = element_line(linewidth = 0.4), 
        axis.ticks = element_line(linewidth = 0.8),
        legend.title = element_text(size = 10, color = 'black'), 
        legend.text = element_text(size = 8, color = 'black'))

# Import growth curve data
df_growth <- read.csv('~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/growthcurves_analysis/gc_data/condensed_PPI.csv')[, -1]

# Remove wells not inoculated (controls)
df_growth <- df_growth[df_growth$strain != 'control', ]

# Compute mean AUC and dgr for each strain (grouped by D12 and D3)
df_growth%<>%
  group_by(D12, D3)%>%
  mutate(med_auc = mean(auc), 
         med_dgr = mean(dgr))

# Prepare data for visualization
df_viz <- df_growth

# Define protein binding domains (PBD) and peptides
PBD <- c('PDZ',  'SH3', 'GBD','GRB2')
peptide <- c('pdz.1', 'pdz.2', 'pdz.3', 'sh3.1', 'sh3.2', 'sh3.3', 'sh3.4', 'gbd1', 'gab2')

# Classify orientation based on PBD and peptide presence
df_viz$orientation <- !(df_viz$D12 %in% peptide) & !(df_viz$D3 %in% PBD)

# Assign PBD and peptide columns for visualization
df_viz[ df_viz$D12 %in% PBD, 'PBD']<- df_viz[ df_viz$D12 %in% PBD, 'D12' ]
df_viz[ df_viz$D3 %in% PBD, 'PBD']<- df_viz[ df_viz$D3 %in% PBD, 'D3' ]

df_viz[ df_viz$D12 %in% peptide, 'peptide']<- df_viz[ df_viz$D12 %in% peptide, 'D12' ]
df_viz[ df_viz$D3 %in% peptide, 'peptide']<- df_viz[ df_viz$D3 %in% peptide, 'D3' ]

# Replace NA values with 'empty'
df_viz$PBD <- na.replace(df_viz$PBD, 'empty')
df_viz$peptide <- na.replace(df_viz$peptide, 'empty')

# Set orientation factor levels and labels
df_viz$orientation <- 
  factor(df_viz$orientation, 
         levels = c(TRUE, FALSE), 
         labels = c('PBD-F[1,2]_peptide-F[3]', 
                    'peptide-F[1,2]_PBD-F[3]'))

# Assign group based on PBD and peptide identity
df_viz[df_viz$PBD == 'SH3' | df_viz$peptide %in% c('sh3.1', 'sh3.2', 'sh3.3', 'sh3.4'), 'group'] <- 'SH3'
df_viz[df_viz$PBD == 'PDZ' | df_viz$peptide %in% c('pdz.1', 'pdz.2', 'pdz.3'), 'group'] <- 'PDZ'
df_viz[df_viz$PBD == 'GBD' | df_viz$peptide == 'gbd1', 'group'] <- 'GBD'
df_viz[df_viz$PBD == 'GRB2' | df_viz$peptide == 'gab2', 'group'] <- 'GRB2'

df_viz$group <- na.replace(df_viz$group, replace = '(-)')

df_viz$group <- 
  factor(df_viz$group, 
         levels = c('(-)', 'GBD', 'GRB2', 'PDZ', 'SH3'))

# Create combination label for plotting
df_viz$comb <- paste0(df_viz$PBD, '_', df_viz$peptide)

# Panel A: Visualize AUC by combination and orientation
FigS1A <- 
  ggplot(df_viz[order(df_viz$PBD), ], 
         aes(x=comb, y = auc))+
  facet_grid(cols = vars(group), scales = 'free_x', space = 'free_x')+
  geom_jitter(aes(color = orientation), 
              alpha = 0.7, position = position_jitterdodge(), size = 1.2)+
  stat_compare_means(aes(x = comb, y = auc, group = orientation), size = 2.5, label.y = 105,
                     method = 't.test', label = 'p.signif',hide.ns = TRUE)+
  scale_color_manual(values = c('#F5191CFF', '#36A5AAFF'), labels = c('F[1,2]-PBD_F[3]-peptide', 'F[1,2]-peptide_F[3]-PBD'))+
  ylab('PCA signal (corrected AUC)')+
  ylim(0, 110)+
  xlab('PBD_peptide combinations')+
  theme_classic2()+
    gt+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'bottom', axis.title.x = element_text( vjust = -0.5),
        axis.title.y = element_text( hjust = 1.1))+
  guides(color = guide_legend(title = ''))

# Panel: Scatter plot of dgr vs AUC, highlighting empty_empty controls in red
ggplot()+
  geom_point(data = df_viz[df_viz$comb !='empty_empty', ], 
            aes(x = dgr, y = auc), 
              alpha = 0.7, size = 2.5)+
  geom_point(data = df_viz[df_viz$comb =='empty_empty', ], 
             aes(x = dgr, y = auc), 
             alpha = 0.7, size = 2.5, color = 'red')+
  ylab('Corrected AUC')+
  xlab('Derivative growth rate')+
  theme_classic2()+
  gt+
  theme(
        legend.position = 'bottom')+
  guides(color = guide_legend(title = ''))

# Keep median information for each combination
df_viz%>%
  select(D12, D3, med_auc, med_dgr, orientation, PBD, peptide, group)%>%
  unique()->med_growth

# Compute AUC and dgr ratio vs empty_empty plasmid (negative control)
neg_growth <- med_growth[med_growth$D12 =='empty' & med_growth$D3 =='empty', ]

med_growth$ratio_auc <- med_growth$med_auc/neg_growth$med_auc
med_growth$ratio_dgr <- med_growth$med_dgr/neg_growth$med_dgr

# Assign the empty_empty observation in both orientations for completeness
x <- med_growth[med_growth$D12 == 'empty' & med_growth$D3 == 'empty', ]
x[, 'orientation'] <- 'peptide-F[1,2]_PBD-F[3]'

med_growth <- bind_rows(med_growth, x)

# Panel B: Visualize AUC ratios as heatmap for main orientation
p2 <- 
  ggplot(med_growth[med_growth$orientation == 'PBD-F[1,2]_peptide-F[3]', ])+
  geom_tile(aes(PBD, peptide, fill = ratio_auc))+
  scale_fill_viridis_c(option = 'F', direction = 1, end = 0.95)+
  theme_classic2()+
    gt+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
        legend.position = 'top', 
        legend.title=element_text(size = 8))+
  guides(fill = guide_colorbar(title = 'PPI score', 
                               theme = theme(
                                 legend.key.width  = unit(5, "lines"),
                                 legend.key.height = unit(1, "lines")
                               )))

# Panel B (supplementary): Heatmap for second orientation
s1 <- 
  ggplot(med_growth[med_growth$orientation == 'peptide-F[1,2]_PBD-F[3]', ])+
  geom_tile(aes(PBD, peptide, fill = ratio_auc))+
  scale_fill_viridis_c(option = 'F', direction = 1, end = 0.95)+
  theme_classic2()+
  gt+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        legend.position = 'bottom', legend.title=element_text(size = 7))+
  guides(fill = guide_colorbar(title = 'PPI score', 
                               theme = theme(legend.margin = margin(l = -40),
                                 legend.key.width  = unit(5, "lines"),
                                 legend.key.height = unit(1, "lines")
                               )))

## Compare PPI score vs affinity from reference data
affinity <- read_excel('~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/growthcurves_analysis/gc_data/ref_affinity.xlsx', 
                       col_names = c('Domain', 'motif', 'motif_seq', 'affinity'))
affinity$affinity <- as.numeric(affinity$affinity)

# Merge median growth data with affinity data
PPI_vs_affinity <- 
  merge(med_growth, 
        affinity, 
        by.x = c('PBD', 'peptide'), 
        by.y = c('Domain', 'motif'), 
        all.x = T)

# Remove GBD due to low abundance/high spurious binding in vivo
sub_viz <- PPI_vs_affinity[PPI_vs_affinity$orientation == 'PBD-F[1,2]_peptide-F[3]' & !is.na(PPI_vs_affinity$affinity)
                           , ]

# Panel C: Scatter plot of PPI score vs affinity

sub_viz[sub_viz$D12 %in% c('GBD'), 'color']<-'red'

sub_viz[!(sub_viz$D12 %in% c('GBD')), 'color']<-'black'

p3 <- 
  ggplot(sub_viz)+
  geom_point(aes(y = ratio_auc, x = affinity, shape = PBD, color = color), size = 3)+
  scale_shape_manual(values = c(15,16,17,18))+
  scale_color_manual(values = c('black', 'darkorange3'))+
  stat_cor(data = sub_viz[sub_viz$PBD != 'GBD', ], 
           (aes(y = ratio_auc, x = affinity)), method = 'spearman', label.y = 10, cor.coef.name = 'rho',  label.sep = "\n")+
  stat_cor((aes(y = ratio_auc, x = affinity)), method = 'spearman', color = 'darkorange3', label.y = 12, cor.coef.name = 'rho',  label.sep = "\n")+
  ylab(' PPI score ')+
  xlab(expression(paste('Affinity ',K[d], ' (', mu, 'M)')))+
  scale_x_continuous(transform = 'log10', breaks = c(0.1, 1, 10, 100, 1000),
                     labels = c(0.1, 1, 10, 100, 1000), limits = c(0.1, 1000))+
  theme_classic2()+
  gt+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
        legend.background = element_rect(fill = 'transparent', color = 'black'),
        legend.title = element_blank(), 
        legend.text = element_text(size = 6), 
        legend.spacing = unit(2, 'mm'), 
        legend.key.size = unit(4, 'mm'), 
        legend.margin = margin(t= 0.5, r = 0.7, b = 0.5, l =0.5, unit = 'mm'), 
        legend.key.spacing.x =unit(0, 'cm'))+
  guides(shape = guide_legend(nrow = 1, position = 'top'), 
         color = 'none')

# Panel A: Import image for figure
p1 <- 
  ggdraw()+
  draw_image('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig1A.png')

# Combine panels for main figure
Fig1 <- 
  plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1,1,1), align = 'h', axis = 'b', 
            labels = c('A', 'B', 'C'), label_fontface = 'plain')

library(svglite)

# Save main figure as SVG
ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig1.svg', Fig1,
      width = 7, height = 2.7, units = 'in')

# Combine panels for supplementary figure
FigS1 <- 
  plot_grid(FigS1A, s1, nrow = 1, rel_widths = c(1.5, 0.8), 
            labels = c('A', 'B'), label_fontface = 'plain', align = 'h',axis = 'tb')

# Save supplementary figure as PNG
ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/final_S1/FigS2.svg', FigS1,
       width = 8, height =3.5)
