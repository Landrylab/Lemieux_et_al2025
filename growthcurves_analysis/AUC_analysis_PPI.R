# Analysis of the PCA DHFR assay to quantify interaction  strength
#import package
library(ggplot2)
library(ggpubr)
library(magrittr)
library(stringr)
library(tidyverse)
library(gtools)
library(readxl)
library(cowplot)

# import data
df_growth <- read.csv('~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/data/condensed_PPI.csv')[, -1]

# remove wells not inoculated
df_growth <- df_growth[df_growth$strain != 'control', ]

# Compute median AUC and dgr for each strain

df_growth%<>%
  group_by(D12, D3)%>%
  mutate(med_auc = median(auc), 
         med_dgr = median(dgr))

# Visualize replicates
# Add columns to classify each observation and help visualization
df_viz <- df_growth

PBD <- c('PDZ',  'SH3', 'GBD','GRB2')
peptide <- c('pdz1', 'pdz2', 'pdz3', 'sh31', 'sh32', 'sh33', 'sh34', 'gbd1', 'gab2')

df_viz$orientation <- !(df_viz$D12 %in% peptide) & !(df_viz$D3 %in% PBD)

df_viz[ df_viz$D12 %in% PBD, 'PBD']<- df_viz[ df_viz$D12 %in% PBD, 'D12' ]
df_viz[ df_viz$D3 %in% PBD, 'PBD']<- df_viz[ df_viz$D3 %in% PBD, 'D3' ]

df_viz[ df_viz$D12 %in% peptide, 'peptide']<- df_viz[ df_viz$D12 %in% peptide, 'D12' ]
df_viz[ df_viz$D3 %in% peptide, 'peptide']<- df_viz[ df_viz$D3 %in% peptide, 'D3' ]

df_viz$PBD <- na.replace(df_viz$PBD, 'empty')
df_viz$peptide <- na.replace(df_viz$peptide, 'empty')

df_viz$orientation <- 
  factor(df_viz$orientation, 
         levels = c(TRUE, FALSE), 
         labels = c('PBD-F[1,2]_peptide-F[3]', 
                    'peptide-F[1,2]_PBD-F[3]'))

df_viz[df_viz$PBD == 'SH3' | df_viz$peptide %in% c('sh31', 'sh32', 'sh33', 'sh34'), 'group'] <- 'SH3'
df_viz[df_viz$PBD == 'PDZ' | df_viz$peptide %in% c('pdz1', 'pdz2', 'pdz3'), 'group'] <- 'PDZ'
df_viz[df_viz$PBD == 'GBD' | df_viz$peptide == 'gbd1', 'group'] <- 'GBD'
df_viz[df_viz$PBD == 'GRB2' | df_viz$peptide == 'gab2', 'group'] <- 'GRB2'

df_viz$group <- na.replace(df_viz$group, replace = '(-)')

df_viz$group <- 
  factor(df_viz$group, 
         levels = c('(-)', 'GBD', 'GRB2', 'PDZ', 'SH3'))


# change peptide labelling
sub_id <-  gsub('sh3', '', df_viz$peptide)
sh3 <- nchar(sub_id) == 1
sub_id[sh3]<- paste0('sh3.', sub_id[sh3])
df_viz$peptide <- sub_id

sub_id <-  gsub('pdz', '', df_viz$peptide)
pdz <- nchar(sub_id) == 1
sub_id[pdz]<- paste0('pdz.', sub_id[pdz])


df_viz$peptide <- sub_id

df_viz$comb <- paste0(df_viz$PBD, '_', df_viz$peptide)

FigS1A <- 
  ggplot(df_viz[order(df_viz$PBD), ], 
         aes(x=comb, y = auc))+
  facet_grid(cols = vars(group), scales = 'free_x', space = 'free_x')+
  geom_jitter(aes(color = orientation), 
              alpha = 0.7, position = position_jitterdodge(), size = 2.5)+
  stat_compare_means(aes(x = comb, y = auc, group = orientation),
                     method = 't.test', label = 'p.signif', hide.ns = TRUE)+
  scale_color_manual(values = c('#F5191CFF', '#36A5AAFF'))+
  ylab('Corrected area under\n the curve')+
  xlab('PBD_peptide combinations')+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14),
        axis.title.x = element_text(vjust = -1),
        strip.text = element_text(color = 'black', size = 12),
        legend.position = 'bottom')+
  guides(color = guide_legend(title = ''))

# Keep median information
med_growth <- unique(df_viz[, c(3:6, 9:14)])

# Compute AUC and dgr ratio vs empty_empty plasmid
neg_growth <- med_growth[med_growth$D12 =='empty' & med_growth$D3 =='empty', ]

med_growth$ratio_auc <- med_growth$med_auc/neg_growth$med_auc
med_growth$ratio_dgr <- med_growth$med_dgr/neg_growth$med_dgr

# assign the empty_empty observation in both orientation
x <- med_growth[med_growth$D12 == 'empty' & med_growth$D3 == 'empty', ]
x[, 'orientation'] <- 'peptide-F[1,2]_PBD-F[3]'

med_growth <- bind_rows(med_growth, x)

# Vizualise auc ratios
# keep only one orientation for main figure panel, show the second orientation in supplementary
p2 <- 
  ggplot(med_growth[med_growth$orientation == 'PBD-F[1,2]_peptide-F[3]', ])+
  geom_tile(aes(PBD, peptide, fill = ratio_auc))+
  scale_fill_viridis_c(option = 'F', direction = 1, breaks = c(0,1,2,3,4,5))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 12),
        legend.position = 'top')+
  guides(fill = guide_colorbar(title = 'PPI score', 
                               theme = theme(
                                 legend.key.width  = unit(10, "lines"),
                                 legend.key.height = unit(1, "lines")
                               )))


s1 <- 
  ggplot(med_growth[med_growth$orientation == 'peptide-F[1,2]_PBD-F[3]', ])+
  #facet_grid(cols =vars(orientation), scales = 'free', space = 'free', drop = T)+
  geom_tile(aes(PBD, peptide, fill = ratio_auc))+
  scale_fill_viridis_c(option = 'F', direction = 1, breaks = c(0,1,2,3,4,5))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 12),
        legend.position = 'bottom')+
  guides(fill = guide_colorbar(title = 'PPI score', 
                               theme = theme(
                                 legend.key.width  = unit(10, "lines"),
                                 legend.key.height = unit(1, "lines")
                               )))

## Compare PPI score vs affinity
affinity <- read_excel('~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/data/ref_affinity.xlsx', col_names = c('Domain', 'motif', 'motif_seq', 'affinity'))
affinity$affinity <- as.numeric(affinity$affinity)

PPI_vs_affinity <- 
  merge(med_growth, 
        affinity, 
        by.x = c('PBD', 'peptide'), 
        by.y = c('Domain', 'motif'), 
        all.x = T)

# remove GBD because of low abundance/high spurious binding in vivo showed by the abundance assay
sub_viz <- PPI_vs_affinity[PPI_vs_affinity$orientation == 'PBD-F[1,2]_peptide-F[3]' & !is.na(PPI_vs_affinity$affinity)
                           , ]

p3 <- 
  ggplot(sub_viz)+
  geom_point(aes(x = ratio_auc, y = affinity, shape = PBD, color= log10(affinity)), size = 3)+
  #stat_cor(aes(x = ratio_auc, y = affinity), label.x = 4, method = 'spearman')+
  scale_color_gradient(high = 'grey80', low = 'black')+
  scale_shape_manual(values = c(15,16,17,18))+
  xlab(' PPI score ')+
  ylab(expression(paste('Affinity .',K[d], ' (', mu, 'M)')))+
  scale_y_continuous(transform = 'log10', breaks = c(0.1, 1, 10, 100, 1000), labels = c(0.1, 1, 10, 100, 1000))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 8),
        axis.title = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 12),
        legend.background = element_rect(fill = 'transparent', color = 'black'),
        legend.title = element_blank(),
        legend.position.inside= c(0.75,0.87))+
  guides(shape = guide_legend(nrow = 2, position = 'inside'), 
         color = 'none')

# Create Fig1  
p1 <- 
  ggdraw()+
  draw_image('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig1A.png')

Fig1 <- 
  plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(0.8,1,1), 
            labels = c('A', 'B', 'C'), label_fontface = 'plain')

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig1.png', Fig1,
       width = 10, height =3.3)

# Create FigS1
FigS1 <- 
  plot_grid(FigS1A, s1, nrow = 1, rel_widths = c(1.5, 0.75), 
            labels = c('A', 'B'), label_fontface = 'plain', align = 'h',axis = 'tb')

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/FigS1.png', FigS1,
       width = 11, height =5)
