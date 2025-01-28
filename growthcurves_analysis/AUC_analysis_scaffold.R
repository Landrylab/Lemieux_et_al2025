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
df_growth <- read.csv('~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/data/condensed_scaffold.csv')

sub_id <-  gsub('sh3', '', df_growth$D12)
sh3 <- nchar(sub_id) == 1
sub_id[sh3]<- paste0('sh3.', sub_id[sh3])
df_growth$D12 <- sub_id

sub_id <-  gsub('sh3', '', df_growth$D3)
sh3 <- nchar(sub_id) == 1
sub_id[sh3]<- paste0('sh3.', sub_id[sh3])
df_growth$D3 <- sub_id

sub_id <-  gsub('pdz', '', df_growth$D12)
pdz <- nchar(sub_id) == 1
sub_id[pdz]<- paste0('pdz.', sub_id[pdz])
df_growth$D12 <- sub_id

sub_id <-  gsub('pdz', '', df_growth$D3)
pdz <- nchar(sub_id) == 1
sub_id[pdz]<- paste0('pdz.', sub_id[pdz])
df_growth$D3 <- sub_id


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
  geom_jitter(aes(y = auc, x = scaffold, color = estradiol), width = 0.1, size =2,alpha = 0.7)+
  scale_color_viridis_d(direction = -1)+
  ylab('Corrected area under \nthe curve')+
  xlab('Scaffold system')+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 12),
        legend.position = 'bottom')+
  guides(color = guide_legend(title = '[β-Estradiol] (nM)'))

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Supp_Fig3XA.png', width = 10, height = 8)


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
  theme(axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14),
        axis.title.x = element_text(vjust = 7),
        strip.text = element_text(color = 'black', size = 12),
        legend.position = 'top')+
  guides(fill = guide_colorbar(title = 'Scaffold \nscore',
                               theme = theme(
                                 legend.key.width  = unit(10, "lines"),
                                 legend.key.height = unit(1.5, "lines")
                               )))

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig3A.png', height = 4, width = 6)

# Comparison with affinity PBD-peptide
affinity <- read_excel('~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/data/ref_affinity.xlsx', col_names = c('PBD', 'peptide', 'peptide_seq', 'affinity'))
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
  geom_jitter(aes(x = comb_affinity, y = auc, color = estradiol), width =0.1)+
  scale_color_viridis_d(option = 'D', direction = -1)+
  xlab('peptide affinity F[1,2]_F[3]')+
  ylab('Corrected area under\n the curve')+
  theme_classic2()+
  guides(color = guide_legend(title = '[β-Estradiol] (nM)', 
                              nrow = 2))+
  theme(axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14),
        axis.title.x = element_text(vjust =  -1.4),
        strip.text = element_text(color = 'black', size = 12),
        legend.background = element_rect(fill = 'transparent', color = 'transparent'),
        legend.position = c(0.75, 0.8))
  
# plot the scaffold signal vs affinity categories
# for all combinations

sub_affinity_1 <- res_affinity[grepl('sh3', res_affinity$D12) & !(res_affinity$strain %in% c('BY4741', 'PL0001')) & grepl('pdz', res_affinity$D3), ]
sub_affinity_2 <- res_affinity[grepl('pdz', res_affinity$D12) & !(res_affinity$strain %in% c('BY4741', 'PL0001')) & grepl('sh3', res_affinity$D3), ]

sub_affinity_1$orientation <- 'sh3-F[1,2]_pdz-F[3]'
sub_affinity_2$orientation <- 'pdz-F[1,2]_sh3-F[3]'

sub_affinity <- bind_rows(sub_affinity_1, sub_affinity_2)

sub_affinity$comb_affinity <- str_c(sub_affinity$type_affinity.x, sub_affinity$type_affinity.y, sep = '_')

sub_affinity$comb_affinity <- 
  factor(sub_affinity$comb_affinity, 
         levels = c('high_high', 'high_medium','medium_high', 'high_low', 'low_high', 'medium_low', 'low_medium', 'low_low'))


Supp4 <- 
ggplot(sub_affinity)+
  facet_grid(cols = vars(scaffold), rows = vars(orientation))+
  geom_jitter(aes(x = comb_affinity, y = auc, color = estradiol), width =0.1)+
  scale_color_viridis_d(option = 'D', direction = -1)+
  xlab('peptide affinity F[1,2]_F[3]')+
  ylab('Corrected area under \nthe curve')+
  theme_classic2()+
  guides(color = guide_legend(title = '[β-Estradiol] (nM)', 
                              nrow = 2))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14),
        axis.title.x = element_text(vjust =  -1),
        strip.text = element_text(color = 'black', size = 12),
        legend.background = element_rect(fill = 'transparent', color = 'transparent'),
        legend.position = c(0.93, 0.85))

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/FigS4.png', Supp4, width = 12, height = 6)


# Create figure 3

p1 <- 
ggdraw()+
  draw_image('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig3A.png')
FigBC <- plot_grid(p2, p3, align = 'h', axis = 'b', rel_widths = c(1.2,1), labels = c('B', 'C'), label_fontface = 'plain')


Fig3 <- 
plot_grid(p1, FigBC, nrow = 1, ncol=2,  rel_widths = c(1, 2.4),
          labels = c('A', '', ''), label_fontface = 'plain')

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig3.png', Fig3, width = 12, height = 5)
