# analysis of gc for PRM and their designed peptides

# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(stringr)
library(viridisLite)
library(paletteer)
library(grDevices)

# Set ggplot theme for all plots
gt <- theme_classic2()+
  theme(axis.title = element_text(size = 10, color = 'black'), 
        axis.text = element_text(size = 8, color = 'black'), 
        axis.line = element_line(linewidth = 0.4), 
        axis.ticks = element_line(linewidth = 0.8),
        legend.title = element_text(size = 10, color = 'black'), 
        legend.text = element_text(size = 8, color = 'black'))

# Set working directory
setwd('~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/')

# Load data
#all_data <- read_csv('pGD_variants_plates_1_9.csv')
condensed <- read_csv('growthcurves_analysis/gc_data/condensed_new_PBD-peptides.csv')

# Mark special case for D12 == '299' except plate 9
condensed[condensed$D12 == '299' & condensed$plate !=9,  'D12'] <- '299*'


# Extract peptide family from D3 column
pep <- unlist(condensed[!(condensed$D3 %in% c('empty', 'sh31')), 'D3'], '-')

pep_family <- as.numeric(matrix(unlist(strsplit(pep, '-')), ncol =2, byrow = T)[,1])

# Assign peptide family to condensed data
condensed[!(condensed$D3 %in% c('empty', 'sh31')), 'pep_family']<- as.character(pep_family)
condensed[condensed$D12 == '299*', 'pep_family'] <- '299*'

# Select true PPI pairs (matching family or controls)
condensed%>%
  filter(D12 == pep_family | D3 %in% c('empty', 'sh31') | D12 == 'empty') -> true_PPI

# Select false PPI pairs (non-matching family, BY4741 strain)
condensed%>%
  filter(D12 != pep_family & strain == 'BY4741')-> false_PPI

# label the method of design of each peptide
true_PPI[grepl('-1', true_PPI$D3), 'method']<-'PWM_top'
true_PPI[grepl('-2', true_PPI$D3), 'method']<-'HAL'
true_PPI[grepl('-3', true_PPI$D3), 'method']<-'NGS'
true_PPI[grepl('-4', true_PPI$D3), 'method']<-'NGS'
true_PPI[grepl('-5', true_PPI$D3), 'method']<-'NGS'
true_PPI[grepl('-6', true_PPI$D3), 'method']<-'NGS'
true_PPI[grepl('-7', true_PPI$D3), 'method']<-'NGS'

# Fill missing method as 'control'
true_PPI$method <- replace_na(true_PPI$method, 'control')

# Handle negative controls
sub_neg <- true_PPI[true_PPI$D12=='empty', ]
true_PPI[true_PPI$D12=='empty', 'method'] <- gsub('empty', 'neg.control', sub_neg$D12)

# Set factor levels for method
true_PPI$method <- 
factor(true_PPI$method, 
       levels = c('control', 'neg.control', 'HAL', 'NGS', 'PWM_top'), 
       labels = c(c('control', 'neg.control', 'HAL', 'NGS', 'PWM_top')))

# Remove NA and unwanted entries for clean analysis
true_PPI%>%
  filter(!is.na(pep_family))%>%
  filter(method != 'control')%>%
  filter(pep_family !='299*')%>%
  filter(!(D3 == '363-2' & plate ==7))-> true_PPI_clean

# Plot small scale PPI data for BY4741 strain
PPI_small_scale <- 
ggplot(true_PPI_clean[true_PPI_clean$strain == 'BY4741',])+
  facet_wrap(vars(pep_family), scales = 'free_x')+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 74), fill = 'grey90', alpha = 0.3)+
  geom_jitter(aes(x = D3, y = auc, color = method), width = 0.2, alpha=0.7)+
  theme_classic2()+
  scale_color_manual(values = c('grey50', paletteer_d("ggthemes::excel_Depth")[c(3,4,1)]))+
  ylab('PPI corrected auc')+
  xlab('F[3]-peptide')+
  gt+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.title.position = 'left', legend.position = 'bottom')+
  guides(color = guide_legend(title= 'peptide \nmethod'))

#saveRDS(PPI_small_scale, '~/PL_projects/PL_papers/Scaffold_Letters/Figures/FigSX.rds')
#ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/final_S1/FigS10.pdf', height = 8, width = 7)


# Re-assign negative control method for D3=='empty'
sub_neg <- true_PPI[true_PPI$D3=='empty', ]
true_PPI[true_PPI$D3=='empty', 'method'] <- gsub('empty', 'neg.control', sub_neg$D3)

# Set method factor levels and labels
true_PPI$method <- 
  factor(true_PPI$method,
         levels = c('control', 'neg.control', 'HAL', 'NGS', 'PWM_top'), 
         labels = c('control', 'peptide only', 'HAL', 'NGS', 'PWM_top'))

# Annotate PBD family types
ID <- c('25', '61','131', '152','214', '246', '250', '299', '299*', '304', '363', '366', '385', 'empty')
type <- c('WW', 'WW', 'SH3', 'SH3', 'SH3', 'SH3', 'SH3', 'SH3', 'SH3', 'SH3', 'PDZ', 'PDZ', 'PDZ', 'control')
ref_PBD <- tibble(ID, type)

# Merge PBD family info into true_PPI
true_PPI <- 
merge(true_PPI, ref_PBD, 
      by.x = 'D12', 
      by.y = 'ID', 
      all.x = TRUE)

## Plot abundance for PBD (PL19 strain)
Abun_d <- 
ggplot(true_PPI[true_PPI$strain %in% c('PL19') & true_PPI$D12 != '299*', ])+
  geom_vline(aes(xintercept = 37),  linetype='dotted', color = 'black', linewidth =0.5)+
  geom_point(aes(y = D12, x = auc, color = type), size=1,  alpha = 0.7)+
  scale_color_manual(values = c('grey50', '#156064', '#00c49a', '#e09f3e'))+
  theme_classic2()+
  xlim(0,100)+
  ylab('F[1,2]-PBD')+
  xlab('Avail. corrected auc')+
  gt+
  theme(axis.text = element_text(color='black'), 
        legend.position = 'top', 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0))+
  guides(color = guide_legend(title = 'PBD family'))

#saveRDS(Abun_d, '~/PL_projects/PL_papers/Scaffold_Letters/Figures/FigSY.rds')

#ref$name <- gsub(ref$name, pattern = 'PRM_0', replacement ='')

# Set method.ref for double empty controls
true_PPI[true_PPI$D12== 'empty' & true_PPI$D3 == 'empty', 'method.pep']<- 'control'

true_PPI[grepl('-1', true_PPI$D3), 'method.pep']<-'PWM_top'
true_PPI[grepl('-2', true_PPI$D3), 'method.pep']<-'HAL'
true_PPI[grepl('-3', true_PPI$D3), 'method.pep']<-'NGS'
true_PPI[grepl('-4', true_PPI$D3), 'method.pep']<-'NGS'
true_PPI[grepl('-5', true_PPI$D3), 'method.pep']<-'NGS'
true_PPI[grepl('-6', true_PPI$D3), 'method.pep']<-'NGS'
true_PPI[grepl('-7', true_PPI$D3), 'method.pep']<-'NGS'


# Set new_ID for missing values
#true_PPI[is.na(true_PPI$new_ID), 'new_ID'] <- 'empty'

## Plot abundance for peptides (PL17 strain)
Abun_p <- 
  ggplot(true_PPI[true_PPI$strain %in% c('PL17'), ])+
   geom_vline(aes(xintercept = 48),  linetype='dotted', color = 'black', linewidth =0.5)+
  geom_point(aes(y = D3, x = auc, color = method.pep), size= 1, alpha = 0.7)+
  scale_color_manual(values = c('grey50', paletteer_d("ggthemes::excel_Depth")[c(3,4,1)]))+
  theme_classic2()+
  xlim(0,100)+
  ylab('F[3]-peptide')+
  xlab('Avail. corrected auc')+
  gt+
  theme(
        legend.position = 'bottom', 
         legend.margin = margin(t = 0, r = 0, b = 0, l = 0))+
    guides(color = guide_legend(title = 'method', nrow = 2))
#saveRDS(Abun_p, '~/PL_projects/PL_papers/Scaffold_Letters/Figures/FigSZ.rds')

# Compute abundance score for PBD and peptide
data_abundance <- true_PPI[true_PPI$strain %in% c('PL17', 'PL19'), ]

## Reference values for normalization
PBD_ref <- mean(data_abundance[data_abundance$D12 == 'empty' & data_abundance$strain == 'PL19', 'auc'])
pep_ref <- mean(data_abundance[data_abundance$D3 == 'empty' & data_abundance$strain == 'PL17', 'auc'])

## Compute mean auc for each sequence
data_abundance%<>%
  group_by(strain, D12, D3)%>%
  mutate(mean_auc = mean(auc))

summary_abundance <- unique(data_abundance[, c('D3', 'D12', 'strain', 'mean_auc', 'method', 'pep_family')])

## Compute normalized scores for PBD and peptide
summary_abundance%>%
  filter(strain == 'PL19')%>%
  mutate(score = mean_auc/PBD_ref)-> PBD_score

summary_abundance%>%
  filter(strain == 'PL17')%>%
  mutate(score = mean_auc/pep_ref)-> pep_score

## Combine scores for matching PBD and peptide families
merge(PBD_score, pep_score, 
      by.x = 'D12',
      by.y = 'pep_family', 
      suffixes = c('.PBD', '.pep'))-> combined_abundance

## Sum scores for combined abundance
combined_abundance%<>%
  group_by(D12, D3.pep)%>%
  mutate(PBD_pep_score = sum(score.PBD, score.pep))

# Compute PPI score for BY4741 strain
data_PPI <- true_PPI_clean[true_PPI_clean$strain =='BY4741', ]

## Reference value for normalization
PPI_ref <- mean(true_PPI[true_PPI$D12 == 'empty' & true_PPI$D3 =='empty' & true_PPI$strain == 'BY4741', 'auc'])

## Compute mean auc and normalized PPI score
data_PPI%<>%
  group_by(D3, D12)%>%
  mutate(mean_auc = mean(auc))%>%
  select(D3, D12, method, mean_auc)

data_PPI%<>%
  group_by(D3,D12)%>%
  mutate(PPI_score = mean_auc/PPI_ref)
  
# Merge abundance and PPI scores for comparison
comparison <- 
merge(combined_abundance[, c('D12', 'score.PBD', 'D3.pep', 'score.pep', 'method.pep', 'PBD_pep_score')], 
      data_PPI[, c('D12', 'D3', 'PPI_score', 'method')], 
      by.x = c('D12', 'D3.pep'), 
      by.y = c('D12', 'D3'))

# Plot correlation between PPI score and combined abundance score
ggplot(comparison)+
  ggplot2::geom_point(aes(x = PPI_score, y =PBD_pep_score, color= method))+
  ggpubr::stat_cor(aes(x = PPI_score, y =PBD_pep_score), method = 'spearman')+
  gt

# Plot false PPI pairs for specificity check
ggplot(false_PPI[!is.na(false_PPI$pep_family), ])+
  facet_wrap(vars(D12), scales = 'free_x')+
  geom_jitter(aes(x = D3, y = auc), width = 0.2)+
  theme_classic2()+
  scale_color_manual(values = c('grey50', paletteer_d("ggthemes::excel_Depth")[c(1,3,4)]))+
  ylab('corrected auc')+
  ylim(0, 100)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.title = element_blank())+
  gt

# Create figure 1C: heatmaps for PBD family, method, PPI score, abundance
comparison%>%
  select(D12, D3.pep, score.PBD, score.pep, PPI_score, PBD_pep_score, method)%>%
  unique()->PPI_abund

PPI_abund$m <- 'method'

PPI_abund <- 
merge(PPI_abund, 
      ref_PBD, 
      by.x = 'D12', 
      by.y = 'ID')

library(ggpubr)
library("ggsci")
library("ggplot2")
library("gridExtra")

# Set factor levels for method and type
PPI_abund$method <- 
factor(PPI_abund$method, 
       labels = c('PWM_top', 'HAL', 'NGS'), 
       levels =  c('PWM_top', 'HAL', 'NGS'))

PPI_abund$type <- factor(PPI_abund$type, 
                         labels = c('PDZ', 'SH3', 'WW'), 
                         levels = c('PDZ', 'SH3', 'WW'))

# Prepare axis labels for heatmaps
lab_axis <- PPI_abund[, c('D12', 'type', 'D3.pep')]
lab_axis <- lab_axis[order(lab_axis$D3.pep), ]
lab_axis <- lab_axis[order(lab_axis$type), ]
lab_axis <- lab_axis[order(as.numeric(lab_axis$D12)), ]
lab_axis$axis <- 1:nrow(lab_axis)

PPI_abund <- left_join(PPI_abund, lab_axis)

# Heatmap for PBD family
pbd_fam <- 
  ggplot(PPI_abund)+
  geom_tile(aes(x = m, y = axis, fill =type))+
  scale_fill_manual(values = c( paletteer_d("ggthemes::Classic_Gray_5")[c(1,2,5)]))+
  theme_classic2()+
  labs(x = 'PBD \nfamily', y='PBD-peptide pairs')+
  scale_y_continuous(labels = PPI_abund$D3.pep, breaks = PPI_abund$axis, transform = 'reverse')+
  gt+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.4, color = 'black'),
        axis.line.x = element_line(linewidth = 0.4, color = 'black'),
        legend.position = 'bottom',
        legend.title = element_blank(),
        panel.border = element_rect(color = "transparent", fill = NA, linewidth = 0.4),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

l_fam <- ggpubr::get_legend(pbd_fam) 

# Heatmap for peptide method
PPI_abund[PPI_abund$method == 'PWM_top', 'peptide \ndesign']<-'max. predicted \npeptide'

PPI_abund[PPI_abund$method != 'PWM_top', 'peptide \ndesign']<-'other'

method_pep <- 
ggplot(PPI_abund)+
  geom_tile(aes(x = m, y = axis, fill =`peptide \ndesign`))+
  scale_fill_manual(values = c( paletteer_d("ggthemes::excel_Depth")[c(1,3)]))+
  scale_y_continuous(labels = PPI_abund$new_ID, breaks = PPI_abund$axis, transform = 'reverse')+
 theme_classic2()+
  labs(x = 'peptide\ndesign', y = 'PBD-peptide pairs')+
  gt+
  theme(axis.title.y = element_blank(), 
        axis.text= element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = 'bottom', 
        legend.title = element_blank(), 
        axis.line.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.4, color = 'black'),
        panel.border = element_rect(color = "transparent", fill = NA, size = 0.4),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

l_met <- ggpubr::get_legend(method_pep) 

# Heatmap for PPI score
ppi <- 
ggplot(PPI_abund)+
  geom_tile(aes(x = m, y = axis, fill =PPI_score))+
  geom_point(data = PPI_abund[PPI_abund$PPI_score>12, ], 
        aes(x = m, y = axis), shape= 4, color = 'black')+
  scale_fill_viridis_c( option = 'F', direction = 1, end = 0.95)+
  scale_y_continuous(labels = PPI_abund$new_ID, breaks = PPI_abund$axis, transform = 'reverse')+
  xlab('PPI score')+
  gt+
  theme(axis.title.y = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.4, color = 'black'),
        panel.border = element_rect(color = "transparent", fill = NA, size = 0.4),
        legend.position = 'bottom', legend.title = element_blank(), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

l_ppi <- ggpubr::get_legend(ppi) 

# Heatmap for PBD abundance score
abun_pbd <- 
ggplot(PPI_abund)+
  geom_tile(aes(x = m, y = axis, fill =score.PBD))+
  scale_fill_viridis_c( option = 'E', direction = 1, end = 0.95,  limits = c(0, 2.25), breaks = c(0, 0.5,  1, 1.5, 2))+
  scale_y_continuous(labels = PPI_abund$new_ID, breaks = PPI_abund$axis, transform = 'reverse')+
  xlab('Avail. score \nPBD')+
  gt+
  theme(axis.title.y = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line.y = element_blank(),
         axis.line.x = element_line(linewidth = 0.4, color = 'black'),
        panel.border = element_rect(color = "transparent", fill = NA, size = 0.4),
        legend.position = 'bottom', legend.title = element_blank(), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

l_abun <- ggpubr::get_legend(abun_pbd) 

# Heatmap for peptide abundance score
abun_pep <- 
ggplot(PPI_abund)+
  geom_tile(aes(x = m, y = axis, fill =score.pep))+
  scale_fill_viridis_c( option = 'E', direction = 1, end = 0.95,  limits = c(0, 2.25), breaks = c(0, 0.5,  1, 1.5, 2))+
  scale_y_continuous(labels = PPI_abund$new_ID, breaks = PPI_abund$axis, transform = 'reverse')+
  xlab('Avail. score \npeptide')+
  gt+
  theme(axis.title.y = element_blank(), 
      axis.text = element_blank(), 
      axis.ticks = element_blank(), 
      axis.line.y = element_blank(),
       axis.line.x = element_line(linewidth = 0.4, color = 'black'),
        panel.border = element_rect(color = "transparent", fill = NA, size = 0.4),
      legend.position = 'none',    plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

# Combine all heatmaps into one figure
library(cowplot)
top <- 
plot_grid(
          pbd_fam+theme(legend.position = 'none'),
          method_pep+theme(legend.position = 'none'), 
          ppi+theme(legend.position = 'none'), 
          abun_pbd+theme(legend.position = 'none'), 
          abun_pep, ncol = 5, align = 'h', axis= 'tb', rel_widths = c(0.8,0.4,0.6,0.6,0.6))
le <-
plot_grid( l_fam,   l_ppi,l_met,l_abun, cols = 2, rows = 2, align = 'hv', axis = 'tbrl')

plot_grid(top, 
          le, nrow = 2, rel_heights = c(1,0.2), align = 'v', axis = 'lr')

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig4C.svg', width = 4.5, height = 7)

# Combine abundance plots for supplementary figure
right_S<-
plot_grid(Abun_d+theme(legend.title = element_blank()),
         Abun_p+theme(legend.position = 'none'), ncol = 1, rel_heights = c(1,2.5), align = 'v',axis = 'lr', 
labels = c('B', 'C'), label_fontface = 'plain')

S<-
plot_grid(PPI_small_scale, right_S, ncol = 2, rel_widths = c(1.8,1), 
          labels = c('A', ''), label_fontface = 'plain')

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/final_S1/FigS10.pdf', S, width = 10, height = 8)

