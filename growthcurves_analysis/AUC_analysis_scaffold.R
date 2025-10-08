# Script to analyse and plot the result of the Scaffold PCA
#import package
library(ggplot2)
library(ggpubr)
library(magrittr)
library(stringr)
library(viridisLite)
library(readxl)
library(cowplot)
library(svglite)
library(tidyverse)

#import data
df_growth <- read.csv('~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/growthcurves_analysis/gc_data/condensed_scaffold.csv')

# set the theme for graphics
gt <- theme(axis.text.x = element_text(), 
            axis.text = element_text(color = 'black', size = 8), 
            legend.text = element_text(color = 'black', size = 8),
            axis.title = element_text(color = 'black', size = 10),
            axis.title.x = element_text(vjust = -1),
            strip.text = element_text(color = 'black', size = 8))

# remove empty wells
df_growth <- df_growth[df_growth$strain !='empty', ]

df_growth$pep_com %<>%
  str_c(df_growth$D12, '_', df_growth$D3)

# compute the average auc and dgr
df_growth%<>%  
  group_by(pep_com, strain, estradiol) %>%
  mutate(med_dgr = mean(dgr), med_auc = mean(auc))


df_growth$estradiol <- factor(df_growth$estradiol, 
                              levels = c(0,5,10,20), 
                              labels = c(0,5,10,20))

# associate the scaffold structure to the strain
df_growth$scaffold <- 
  factor(df_growth$strain, 
         levels = c("BY4741", "PL0001", "PL0012", "PL0013", "PL0014", "PL0015", "PL0016"), 
         labels = c('no GEM TF & no scaffold', 'no scaffold', 'PDZ-SH3', 'PDZ-linker-SH3', 'SH3-PDZ', 'SH3-linker-PDZ', 'SH3-2xlinker-PDZ'))


# All auc values to SuppFigure

ggplot(df_growth)+
  facet_wrap(vars(pep_com),nrow = 3, ncol = 7)+
  geom_jitter(aes(y = auc, x = scaffold, color = estradiol), width = 0.1, size =1.3,alpha = 0.7)+
  scale_color_viridis_d(direction = -1)+
  ylab('PCA signal (corrected AUC)')+
  xlab('Scaffold system')+
  theme_classic2()+
  gt+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = 'bottom', 
        legend.title = element_text(size = 8))+
  guides(color = guide_legend(title = '[β-estradiol] (nM)'))

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/final_S1/FigS6.pdf', width = 7, height = 6)

# t-test to compare the PDZ-SH3 and SH3-PDZ scaffold with the sh3.1_pdz.1 peptide combination
x <- subset(df_growth, 
            scaffold %in% c('PDZ-SH3', 'SH3-PDZ') & 
              pep_com %in% c('sh3.1_pdz.1') &
              estradiol == 20)

compare_means(data = x, auc ~scaffold, method = 't.test', e)

# perform the ratio of auc and dgr vs the empty plasmid

df_ref <- unique(df_growth[df_growth$pep_com == 'empty_empty', c(2,3,4,10,13:16)])

colnames(df_ref) <- c("strain", "D12","D3", "estradiol", 
                      "pep_com","med_dgr_ref","med_auc_ref" , 'scaffold')

df_growth <-
  merge(df_growth, df_ref[,-c(2,3,5,8)], 
        by = c('strain', 'estradiol'))

df_growth%<>%
  mutate(r_auc = med_auc/med_auc_ref, 
         r_dgr = med_dgr/med_dgr_ref, 
         d_auc = med_auc - med_auc_ref, 
         d_dgr = med_dgr -med_dgr_ref)

res_growth <- unique(df_growth[, c(1,2,4,5,13:22)])

# Main figure for AUC ratio of plasmids vs empty in scaffold strain at 20nM estradiol

sub_pep <- c('empty_empty', 'pdz.1_sh3.1', 'pdz.1_sh3.3', 'pdz.3_sh3.1', 'pdz.3_sh3.3','sh3.1_pdz.1', 'sh3.1_pdz.3', 'sh3.3_pdz.1', 'sh3.3_pdz.3')

p2 <- 
  ggplot(res_growth[res_growth$pep_com %in% sub_pep & res_growth$estradiol ==20, ])+
  geom_tile(aes(y = pep_com, x = scaffold, fill = r_auc))+
  theme_classic2()+
  ylab('peptides F[1,2]_F[3]')+
  xlab('Scaffold system')+
  scale_fill_viridis_c(option = 'G', direction = -1)+
  gt+
  theme(axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1), 
        legend.position = 'top',
        legend.spacing = unit(0, 'mm'), 
        #legend.key.size = unit(4, 'mm'), 
        legend.margin = margin(t= 0, r = 0, b = 0, l =0, unit = 'mm'), 
        legend.key.spacing.x =unit(0, 'cm'),
        legend.title = element_text(size = 8), 
        axis.title.x = element_text(vjust = 3), 
        legend.justification.top = "left")+
  guides(fill = guide_colorbar(title = 'Scaffold \nscore',
                               theme = theme(
                                 legend.key.width  = unit(5, "lines"),
                                 legend.key.height = unit(0.85, "lines"), 
                                 
                               )))

#ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig3A.png', height = 4, width = 6)

# Comparison with affinity PBD-peptide
affinity <- read_excel('~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/growthcurves_analysis/gc_data/ref_affinity.xlsx')
affinity$affinity <- as.numeric(affinity$affinity)

# Classify Kd in category of affinity
affinity%>%
  filter(affinity<10)%>%
  mutate(type_affinity = 'high')->high_aff

affinity%>%
  filter(10<=affinity & affinity <100)%>%
  mutate(type_affinity = 'medium')->med_aff

affinity%>%
  filter(100<=affinity )%>%
  mutate(type_affinity = 'low')->low_aff

affinity <- bind_rows(high_aff,med_aff, low_aff)


res_affinity <- 
  merge(df_growth, 
        affinity[, c(2,4,5)], 
        by.x = 'D12', 
        by.y = 'peptide')

res_affinity <- 
  merge(res_affinity, 
        affinity[, c(2,4,5)], 
        by.x = 'D3', 
        by.y = 'peptide')

colnames(res_affinity)[c(23,25)]<- c('affinity_D12', 'affinity_D3')

res_affinity$sum_affinity <- res_affinity$affinity_D12 +res_affinity$affinity_D3



# plot the scaffold signal vs affinity categories
# orientation sh3.x-D12 and pdz.x-D3
# for PDZ-SH3 scaffold

sub_affinity <- res_affinity[grepl('sh3', res_affinity$D12) & !(res_affinity$strain %in% c('BY4741', 'PL0001')) & grepl('pdz', res_affinity$D3), ]

sub_affinity$comb_affinity <- str_c(sub_affinity$type_affinity.x, sub_affinity$type_affinity.y, sep = '_')

sub_affinity$comb_affinity <- 
  factor(sub_affinity$comb_affinity, 
         levels = c('high_high', 'medium_high', 'high_low', 'low_high', 'medium_low', 'low_low'))

p3 <- 
  ggplot(sub_affinity[sub_affinity$scaffold == 'PDZ-SH3', ])+
  facet_grid(cols = vars(scaffold))+
  geom_jitter(aes(x = comb_affinity, y = auc, color = estradiol), size = 1.2, alpha = 0.8, position = position_jitterdodge())+
  scale_color_viridis_d(option = 'D', direction = -1)+
  xlab('peptide affinity F[1,2]_F[3]')+
  ylab('PCA signal (corrected AUC)')+
  theme_classic2()+
  guides(color = guide_legend(title = '[β-estradiol] (nM)', 
                              nrow = 2, 
                              position = 'inside'))+
  gt+
  theme(axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1), 
        axis.title.x = element_text(hjust =  1.5),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.spacing = unit(0, 'mm'), 
        legend.key.size = unit(4, 'mm'), 
        legend.margin = margin(t= 0, r = 0, b = 0, l =0, unit = 'mm'), 
        legend.key.spacing.x =unit(0, 'cm'),
        legend.background = element_rect(fill = 'transparent', color = 'transparent'),
        legend.position.inside = c(0.7, 0.85))

# plot the scaffold signal vs affinity categories
# for all combinations

sub_affinity_1 <- res_affinity[grepl('sh3', res_affinity$D12) & !(res_affinity$strain %in% c('BY4741', 'PL0001')) & grepl('pdz', res_affinity$D3), ]
sub_affinity_2 <- res_affinity[grepl('pdz', res_affinity$D12) & !(res_affinity$strain %in% c('BY4741', 'PL0001')) & grepl('sh3', res_affinity$D3), ]

sub_affinity_1$orientation <- 'F[1,2]-sh3_F[3]-pdz'
sub_affinity_2$orientation <- 'F[1,2]-pdz_F[3]-sh3'

sub_affinity <- bind_rows(sub_affinity_1, sub_affinity_2)

sub_affinity$comb_affinity <- str_c(sub_affinity$type_affinity.x, sub_affinity$type_affinity.y, sep = '_')

sub_affinity$comb_affinity <- 
  factor(sub_affinity$comb_affinity, 
         levels = c('high_high', 'high_medium','medium_high', 'high_low', 'low_high', 'medium_low', 'low_medium', 'low_low'))


Supp5 <- 
  ggplot(sub_affinity)+
  facet_grid(cols = vars(scaffold), rows = vars(orientation))+
  geom_jitter(aes(x = comb_affinity, y = auc, color = estradiol), size =1.2, alpha = 0.8, position = position_jitterdodge())+
  scale_color_viridis_d(option = 'D', direction = -1)+
  xlab('peptide affinity F[1,2]_F[3]')+
  ylab('PCA signal(corrected AUC)')+
  theme_classic2()+
  guides(color = guide_legend(title = '[β-estradiol] (nM)', 
                              nrow = 1, 
                              position = 'bottom'))+
  gt+
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust =1), 
        axis.title.x = element_text(vjust = -1),
        legend.title = element_text(size = 8),
        legend.background = element_rect(fill = 'transparent', color = 'transparent'),
        legend.position.inside  = c(0.93, 0.85), 
        strip.text.y = element_text(size = 8),
        legend.text = element_text(size = 6), 
        legend.spacing = unit(0, 'mm'), 
        legend.key.size = unit(4, 'mm'), 
        legend.margin = margin(t= 0, r = 0, b = 0, l =0, unit = 'mm'), 
        legend.key.spacing.x =unit(0, 'cm'))

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/final_S1/FigS5.png', Supp5, width = 7, height = 4)


# Create figure 3

p1 <- 
  ggdraw()+
  draw_image('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig3A.png')

Fig3BC <-
  plot_grid( p3, p2, align = 'h', axis = 'b', rel_widths = c(1,1.2), labels = c( 'B', 'C'), label_fontface = 'plain')

Fig3 <- 
  plot_grid(p1, Fig3BC, nrow = 1, ncol=2,  rel_widths = c(1, 2.4),
            labels = c('A', '', ''), label_fontface = 'plain')

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig3.svg', Fig3, width = 7, height = 3.5)
