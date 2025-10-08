# generate peptides library from the selected PBD candidate
library(tidyverse)

## Variables to modify ##
# directory of the repository with the folder reference_data
path = '~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/PBD_selection/'
output_dir = './'

# selected manually from the Primary_ID of the PWM AAD divergence heatmap
manual_candidate <- 
  c('PRM_0250',
    'PRM_0304',
    'PRM_0131',
    'PRM_0152',
    'PRM_0246',
    'PRM_0299',	
    'PRM_0214', 
    'PRM_0363', 
    'PRM_0385', 
    'PRM_0366', 
    'PRM_0061',
    'PRM_0025')

################################

###### Functions ##########

# Function to identify the stronger binder for one PWM by selecting
# the aa with the highest score for each position
find_top_binder <- 
  function(pwm){
    max_pep <- colnames(pwm)[max.col(pwm, ties.method = 'first')]
    max_pep <-str_c(max_pep, collapse = '')
    return(max_pep)
  }

# read master prm_db
master <- read_csv(paste0(path,'reference_data/PRM_Master.csv'))

# keep information about manually seleted PBD
info_candidate <- master[master$Primary_ID %in% manual_candidate, ]

# Get the peptides validated by the Helisa (HAL) or in high abundance in the 
# phage display assay (NGS) for the selected PBD
abun_pep <- info_candidate[, 13:15]
NGS_supp <- strsplit(abun_pep$PDB_ID ,',')

NGS_supp <- lapply(NGS_supp, gsub, pattern = '#N/A', replacement = '')
NGS_supp <- lapply(NGS_supp, gsub, pattern = '^[a-z0-9]+', replacement = '')

NGS_pep <- tibble(PRM_ID = character(), peptide = character())

for(i in 1:length(NGS_supp)){
  seq <- NGS_supp[[i]][NGS_supp[[i]]!='']
  id <-info_candidate[i, 'Primary_ID']
  NGS_pep <- bind_rows(NGS_pep, bind_cols(PRM_ID = unlist(id), peptide = seq))
  
}

NGS_pep$method <- 'NGS'
top_NGS <- info_candidate[, c('NGS', 'Primary_ID')]
top_NGS$method <- 'NGS'
colnames(top_NGS)[1:2] <- c('peptide', 'PRM_ID')

NGS_pep <- 
rbind(top_NGS, NGS_pep)

HAL_pep <- tibble(peptide = info_candidate$HAL, PRM_ID = info_candidate$Primary_ID, method = 'HAL')
HAL_pep <- HAL_pep[HAL_pep$peptide !='#N/A', ]

# join all peptides for NGS ans HAL method for the candidate PBD
exp_pep <- rbind(NGS_pep, HAL_pep)


# Use the PWM to compute top binder for each PBD
PWM <- 
  readRDS(paste0(path, 'reference_data/PWM_updated.rds'))

candidate_PWM <- PWM[names(PWM) %in% info_candidate$Primary_ID]

candidate_top_binder <- 
  lapply(candidate_PWM, find_top_binder)

# add the top binder to the experimental peptides
id <- names(candidate_top_binder)
seq <- unlist(candidate_top_binder)
method <- 'PWM_top'

candidate_top_binder <- bind_cols(PRM_ID = id, peptide = seq, method = method)

exp_pep <- 
  bind_rows(exp_pep, candidate_top_binder)

# Add domain group to the peptide
exp_pep <- 
merge(exp_pep, 
      info_candidate[, c(1,4)], 
      by.x = 'PRM_ID', 
      by.y='Primary_ID', 
      order = F, 
      all.x = T, 
      all.y = T)


info_exp <- 
merge(info_candidate, 
      exp_pep, 
      by.x = c('Primary_ID', 'Domain_Group'), 
      by.y = c('PRM_ID', 'Domain_Group'))

info_exp <- 
  info_exp[, c(1,2,4,9,12,16,17)]

colnames(info_exp)[6] <- 'peptide_sequence'

# Save the library info
write_csv(info_exp, 
          paste0(output_dir, 'peptide_selection.csv'))

