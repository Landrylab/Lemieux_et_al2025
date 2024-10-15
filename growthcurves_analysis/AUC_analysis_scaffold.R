#import package
library(ggplot2)
library(ggpubr)
library(magrittr)
library(stringr)
library(viridisLite)
library(readxl)
library(cowplot)
library(svglite)

#import data
df_growth <- read.csv('~/PL_projects/PL_papers/Scaffold_Letters/Data/condensed_SCAFFOLD.csv')

# remove empty wells
df_growth <- df_growth[df_growth$strain !='empty', ]

df_growth$pep_com %<>%
  str_c(df_growth$D12, '.', df_growth$D3)

# compute the median auc and dgr
df_growth%<>%  
  group_by(pep_com, strain, estradiol) %>%
  mutate(med_dgr = median(dgr), med_auc = median(auc))


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
  facet_wrap(vars(pep_com),nrow = 3)+
  geom_jitter(aes(y = auc, x = scaffold, color = estradiol), width = 0.1, size =2,alpha = 0.7)+
  scale_color_viridis_d(direction = -1)+
  ylab('Area under the \nthe curve')+
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

df_ref <- unique(df_growth[df_growth$pep_com == 'empty.empty', c(1,2,4,6,11:14)])

colnames(df_ref) <- c("strain", "estradiol", "D12","D3", 
                      "pep_com","med_dgr_ref","med_auc_ref" , 'scaffold')

df_growth <-
  merge(df_growth, df_ref[,-c(3,4,5,8)], 
                   by = c('strain', 'estradiol'))

df_growth%<>%
  mutate(r_auc = med_auc/med_auc_ref, 
         r_dgr = med_dgr/med_dgr_ref, 
         d_auc = med_auc - med_auc_ref, 
         d_dgr = med_dgr -med_dgr_ref)

res_growth <- unique(df_growth[, c(1,2,4,6,11:20)])

# Main figure for AUC ratio of plasmids vs empty in scaffold strain at 20nM estradiol

sub_pep <- c('empty.empty', 'pdz1.sh31', 'pdz1.sh33', 'pdz3.sh31', 'pdz3.sh33','sh31.pdz1', 'sh31.pdz3', 'sh33.pdz1', 'sh33.pdz3')

#p2 <- 
ggplot(res_growth[res_growth$pep_com %in% sub_pep & res_growth$estradiol ==20, ])+
  geom_tile(aes(y = pep_com, x = scaffold, fill = r_auc))+
  theme_classic2()+
  ylab('peptides F[1,2].F[3]')+
  xlab('Scaffold system')+
  scale_fill_viridis_c(option = 'G', direction = -1)+
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14),
        axis.title.x = element_text(vjust = 7),
        strip.text = element_text(color = 'black', size = 12),
        legend.position = 'right')+
  guides(fill = guide_colorbar(title = 'Scaffold \nscore',
                               theme = theme(
                                 legend.key.width  = unit(1.5, "lines"),
                                 legend.key.height = unit(10, "lines")
                               )))

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig3XA.png', height = 4, width = 6)

# Comparison with affinity PBD-peptide
affinity <- read_excel('~/PL_projects/docking/ref_affinity.xlsx', col_names = c('PBD', 'peptide', 'peptide_seq', 'affinity'))
affinity$affinity <- as.numeric(affinity$affinity)

res_affinity <- 
merge(res_growth, 
      affinity[, c(2,4)], 
      by.x = 'D12', 
      by.y = 'peptide')

res_affinity <- 
merge(res_affinity, 
      affinity[, c(2,4)], 
      by.x = 'D3', 
      by.y = 'peptide')

colnames(res_affinity)[15:16] <- c('affinity_D12', 'affinity_D3')

res_affinity$sum_affinity <- res_affinity$affinity_D12 +res_affinity$affinity_D3

# plot the scaffold score vs the sum of the affinity values
# orientation sh3.x-D12 and pdz.x-D3

sub_affinity <- res_affinity[grepl('sh3', res_affinity$D12) & !(res_affinity$strain %in% c('BY4741', 'PL0001')) & grepl('pdz', res_affinity$D3), ]

#p3 <- 
ggplot(sub_affinity, 
       aes(x = sum_affinity, y=r_auc, color = estradiol))+
  facet_grid(cols = vars(scaffold))+
 scale_x_log10()+
  geom_point(aes(x = sum_affinity, y=r_auc, color = estradiol), alpha = 0.7)+
  geom_smooth(aes(x = sum_affinity, y=r_auc, color = estradiol), alpha = 0.7,
              method = "nls", formula = y ~ a * exp(-S * x), 
              method.args = list(start = list(a = 3, S = 1)), se = FALSE)+
  scale_color_viridis_d(option = 'D', direction = -1)+
  xlab(expression(paste('Scaffold K'[d], ' (', mu, 'M)')))+
  ylab('Scaffold score')+
  theme_classic2()+
  guides(color = guide_legend(title = '[β-Estradiol] (nM)', 
                              nrow = 2))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1), 
        axis.text = element_text(color = 'black', size = 12), 
        legend.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 12),
        legend.background = element_rect(fill = 'transparent', color = 'transparent'),
        legend.position = c(0.92, 0.7))

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig3XB.png', height = 4, width = 6.5)

p1 <- 
ggdraw()+
  draw_image('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Scaffold PCA (Fig1A_Letters).png')
t <- 
plot_grid(p1,p2, rel_widths = c(1, 1.4), labels = c('A', 'B'), 
          label_fontface = 'plain')
Fig3 <- 
plot_grid(t, p3, rel_heights = c(1, 0.6), nrow = 2, ncol=1,
          labels = c('', 'C'), label_fontface = 'plain')

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig3.png', Fig3,
       width = 10, height = 7)
