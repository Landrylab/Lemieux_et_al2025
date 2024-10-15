#import package
library(ggplot2)
library(ggpubr)
library(magrittr)
library(stringr)
library(tidyverse)
library(gtools)
library(readxl)
library(cowplot)
library(svglite)

# import data
df_growth <- read.csv('~/PL_projects/PL_papers/Scaffold_Letters/Data/condensed_Abundance.csv')[, -1]

# remove wells not inoculated
df_growth <- df_growth[df_growth$strain != 'control', ]

# Compute median AUC and dgr for each strain

df_growth%<>%
  group_by(D12, D3, condition)%>%
  mutate(med_auc = median(auc), 
         med_dgr = median(dgr))

# Visualize replicates
# Add columns to classify each observation and help visualization
df_viz <- df_growth

PBD <- c('PDZ',  'SH3', 'GBD','GRB2')
peptide <- c('pdz1', 'pdz2', 'pdz3', 'sh31', 'sh32', 'sh33', 'sh34', 'gbd1', 'gab2')

df_viz$orientation <- !(df_viz$D12 %in% peptide) & !(df_viz$D3 %in% PBD)

df_viz[ df_viz$D12 %in% PBD, 'protein']<- df_viz[ df_viz$D12 %in% PBD, 'D12' ]
df_viz[ df_viz$D3 %in% PBD, 'protein']<- df_viz[ df_viz$D3 %in% PBD, 'D3' ]

df_viz[ df_viz$D12 %in% peptide, 'protein']<- df_viz[ df_viz$D12 %in% peptide, 'D12' ]
df_viz[ df_viz$D3 %in% peptide, 'protein']<- df_viz[ df_viz$D3 %in% peptide, 'D3' ]

df_viz$protein <- na.replace(df_viz$protein, 'empty')

df_viz$orientation <- 
factor(df_viz$orientation, 
       levels = c(TRUE, FALSE), 
       labels = c('PBD-F[1,2] & peptide-F[3]', 
                  'peptide-F[1,2] & PBD-F[3]'))

df_viz[df_viz$protein == 'SH3' | df_viz$protein %in% c('sh31', 'sh32', 'sh33', 'sh34'), 'group'] <- 'SH3'
df_viz[df_viz$protein == 'PDZ' | df_viz$protein %in% c('pdz1', 'pdz2', 'pdz3'), 'group'] <- 'PDZ'
df_viz[df_viz$protein == 'GBD' | df_viz$protein == 'gbd1', 'group'] <- 'GBD'
df_viz[df_viz$protein == 'GRB2' | df_viz$protein == 'gab2', 'group'] <- 'GRB2'

df_viz$group <- na.replace(df_viz$group, replace = '(+)')

df_viz$group <- 
factor(df_viz$group, 
      levels = c('(+)', 'GBD', 'GRB2', 'PDZ', 'SH3'))

FigS2 <- 
ggplot(df_viz[df_viz$condition =='MTX' , ])+
  facet_grid(cols = vars(group), scales = 'free_x', space = 'free_x')+
 
  geom_jitter(aes(x = protein, y = auc, color = strain), 
              alpha = 0.8, width = 0.1, size = 2.5)+
  stat_compare_means(aes(x =protein, y = auc, group = strain), 
                     hide.ns = TRUE, method = 't.test', label = 'p.signif', label.y = 210)+
  scale_color_manual(values = c('#51C3CC', '#FF8E32'), 
                     labels = c('Protein-D[3]', 'Protein-D[1,2]'))+
  ylab('Area under the curve')+
  xlab('Protein')+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 12),
        legend.position = 'bottom')+
  guides(color = guide_legend(title = ''))

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/FigS2.png', FigS2,width = 8, height = 4)

# Keep median information
med_growth <- unique(df_viz[, c(3:6, 9:13)])

# Compute AUC and dgr ratio vs empty_empty plasmid and remove DMSO condition from analysis
med_growth <- med_growth[med_growth$condition == 'MTX', ]

neg_growth <- med_growth[med_growth$D12 =='empty' & med_growth$D3 =='empty' & med_growth$condition =='MTX', ]
  
D12_neg <- unlist(neg_growth[neg_growth$strain=='BY4741-mVenus-D3', "med_auc" ])
D3_neg <- unlist(neg_growth[neg_growth$strain=='BY4741-mVenus-D12', "med_auc" ])

med_growth[med_growth$strain == 'BY4741-mVenus-D3', 'ratio_auc'] <- med_growth[med_growth$strain == 'BY4741-mVenus-D3', "med_auc"]/D12_neg
med_growth[med_growth$strain == 'BY4741-mVenus-D12', 'ratio_auc'] <- med_growth[med_growth$strain == 'BY4741-mVenus-D12', "med_auc"]/D3_neg

D12_neg <- unlist(neg_growth[neg_growth$strain=='BY4741-mVenus-D3', "med_dgr" ])
D3_neg <- unlist(neg_growth[neg_growth$strain=='BY4741-mVenus-D12', "med_dgr" ])


med_growth[med_growth$strain == 'BY4741-mVenus-D3', 'ratio_dgr'] <- med_growth[med_growth$strain == 'BY4741-mVenus-D3', "med_dgr"]/D12_neg
med_growth[med_growth$strain == 'BY4741-mVenus-D12', 'ratio_dgr'] <- med_growth[med_growth$strain == 'BY4741-mVenus-D12', "med_dgr"]/D3_neg

med_growth$DHFR <- 
  factor(med_growth$strain, 
         levels = c('BY4741-mVenus-D3', 'BY4741-mVenus-D12'), 
         labels = c('F[1,2]', 'F[3]'))

med_growth[med_growth$D12 %in% PBD | med_growth$D3 %in% PBD, 'type'] <- 'PBD'
med_growth[med_growth$D12 %in% peptide | med_growth$D3 %in% peptide, 'type'] <- 'peptide'
med_growth$type <- na.replace(med_growth$type, '(+)')

# Vizualize auc ratios
# keep only one orientation for main figure panel, show the second orientation in supplementary
p2 <- 
ggplot(med_growth)+
  facet_grid(cols = vars(type), scale = 'free_x', space = 'free')+
  geom_tile(aes(protein, DHFR, fill = ratio_dgr))+
  scale_fill_viridis_c(option = 'E', direction = 1)+
  theme_classic2()+
  xlab('Protein')+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 12),
        axis.text.y = element_text(angle = 90, hjust = 0.5), 
        axis.title.y = element_blank(),
        legend.position = 'right')+
  guides(fill = guide_colorbar(title = 'Availability \nscore',
                               
                               theme = theme(
    legend.key.width  = unit(1.5, "lines"),
    legend.key.height = unit(10, "lines")
  )))

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig2XA.png', 
       width = 4, height = 5.5)  

p1 <- 
ggdraw()+
  draw_image('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Abundance PCA(4).png')

#Figure 2
Fig2 <- 
plot_grid(p1, p2, rel_widths = c(1, 2.2), labels = c('A', 'B'), 
          label_fontface = 'plain')

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig2.png', Fig2,
       height = 3, width = 10)
  