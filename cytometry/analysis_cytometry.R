library(tidyverse)
library(ggpubr)
library(magrittr)

setwd(dir = '~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/cytometry/')

# Validation of induction with beta-estradiol
GFP_data <- read_csv("GFP_induction_data.csv" )

GFP_data$position <- as.factor(GFP_data$position)
GFP_data$rep <- as.factor(GFP_data$rep)
GFP_data$estradiol <- as.numeric(GFP_data$estradiol)

# Visualization of distribution of fluoresence per replicate
ggplot(GFP_data)+
  facet_grid(rows = vars(estradiol))+
  geom_density(aes(x = log10(GRN.B.HLin), color = rep))

GFP_data$GRN_normFSC <- 
  GFP_data$GRN.B.HLin/GFP_data$FSC.HLin

# Compute median of GFP signal (n=5000) per replicate
GFP_data%>%
  dplyr::group_by(position, estradiol)%>%
  dplyr::mutate(med_GFP = median(GRN.B.HLin),
         med_normGFP = median(GRN_normFSC), 
          na.rm = TRUE)%>%
  select(strain, position, estradiol, rep, med_GFP, med_normGFP)%>%
  unique()->simple_LP_GFP


SuppFig4<-
ggplot(simple_LP_GFP)+
  geom_smooth(aes(x=estradiol, y=(med_GFP)), color = 'green')+
  geom_point(aes(x=estradiol, y=(med_GFP)), alpha = 0.5, shape = 1)+
  ggpubr::theme_classic2()+
  scale_x_continuous(transform = 'log2',
                     breaks = c(0.5, 1, 2, 5, 10, 20, 30, 50))+
 # scale_y_continuous(transform = 'log2')+ 
  labs(x = expression(paste('[', beta, ' -estradiol] (nM)')), y = 'med. GFP signal (a.u.) ')+
  theme(axis.title = element_text(size = 10, color = 'black'), 
        axis.text = element_text(size = 8, color = 'black'), 
        axis.line = element_line(linewidth = 0), 
        axis.ticks = element_line(linewidth = 0.8),
        legend.title = element_text(size = 10, color = 'black'), 
        legend.text = element_text(size = 8, color = 'black'), 
        panel.border = element_rect(color = "black", 
                                    fill = NA,
                                    size = 0.8))

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/final_S1/FigS4.png', SuppFig4, width = 6, height = 4)

# Validation of effect of peptides on protein abundance(GFP)
GFP_pep<-
  read_csv('GFP_peptide_abundance.csv')

# Compute median of GFP signal (n=5000) per replicate
GFP_pep%>%
  dplyr::group_by(strain, position)%>%
  dplyr::mutate(med_GFP = median(GRN.B.HLin), na.rm = TRUE)%>%
  select(strain, position, rep, med_GFP)%>%
  unique()->simple_degron_GFP

ggplot(simple_degron_GFP)+
  geom_jitter(aes(x=strain, y=med_GFP), width = 0.1)+
  theme_classic()


# compare with growth curves measurements of binding availability
# for the peptides extracted from the PRM DB
gc_data <- 
  read_csv('../growthcurves_analysis/gc_data/condensed_new_PBD-peptides.csv')

gc_data%>%
  subset(D12 =='empty' & D3 == 'empty' & strain == 'PL17' )%>%
  dplyr::mutate(ref_auc = mean(auc))%>%
  select(ref_auc)%>%
  unlist()%>%
  unique-> ref_auc_1

# for the peptides previously used in synthetic scaffold
gc_data_scaffold <- read_csv('../growthcurves_analysis/gc_data/condensed_availability.csv')

gc_data_scaffold%>%
  subset(D12 =='empty' & D3 == 'empty' & strain == 'BY4741-mVenus-D12'& condition == 'MTX')%>%
  dplyr::mutate(ref_auc = mean(auc))%>%
  select(ref_auc)%>%
  unlist()%>%
  unique-> ref_auc_2

sub_Scaf <- 
  gc_data_scaffold[gc_data_scaffold$D3 %in% c('sh3.1', 'pdz.2') & gc_data_scaffold$condition == 'MTX', ]

# compute Avail. score for both growth curves experiments.
gc_data%<>%
  subset(strain == 'PL17')%>%
  dplyr::mutate(rel_auc= auc/ref_auc_1)

sub_Scaf%<>%
  dplyr::mutate(rel_auc= auc/ref_auc_2)


# merge both experiments

avail_sub <- 
  gc_data[gc_data$D3 %in% simple_degron_GFP$strain, ]

avail_sub <- 
  bind_rows(avail_sub, sub_Scaf)

# function to compute standard deviation per aa_sequence
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# compute summary of growth curve replicates
avail_sum <- 
  data_summary(avail_sub, 'rel_auc', 'D3')

# compute summary of cytometry replicates
gfp_sum <- 
  data_summary(simple_degron_GFP, 'med_GFP', 'strain')

# merge summaries of growth curves and cytometry measurements
comp <- 
  left_join(gfp_sum, 
            avail_sum, 
            join_by(strain ==  D3), 
            unmatched = 'drop')

# auc vs norm_sel_coef
Fig4D <- 
  ggplot(data = comp,
         aes(x = rel_auc,y =med_GFP))+ 
  geom_errorbarh(aes(xmin = rel_auc-sd.y,xmax = rel_auc+sd.y), color = 'darkgrey') + 
  geom_errorbar(aes(ymin = med_GFP-sd.x ,ymax = med_GFP+sd.x), color = 'darkgrey')+
  geom_point(size =0.8) + 
  stat_cor(method = 'pearson', cor.coef.name = 'rho', size = 4, label.sep = '\n', 
           label.y = 7, label.x = -1.5)+
  ylab('GFP signal (a.u.)')+
  xlab('Avail. score')+
  scale_x_continuous(transform = 'log2', 
                     breaks = scales::trans_breaks("log2", \(x) 2^x))+
  annotation_logticks(sides = 'b', 
                      short = unit(1,"mm"),
                      mid = unit(2,"mm"),
                      long = unit(3,"mm"))+
  theme_classic2()+
    theme(axis.title = element_text(size = 10, color = 'black'), 
          axis.text = element_text(size = 8, color = 'black'), 
          axis.line = element_line(linewidth = 0), 
          axis.ticks = element_line(linewidth = 0.8),
          legend.title = element_text(size = 10, color = 'black'), 
          legend.text = element_text(size = 8, color = 'black'), 
          panel.border = element_rect(color = "black", 
                                      fill = NA,
                                      size = 0.8))

ggsave('~/PL_projects/PL_papers/Scaffold_Letters/Figures/Fig4D.png', Fig4D, height = 2, width = 3)
