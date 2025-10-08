# Compute distances between PWM with averaged AAD over the compared PWM position
# to select the most divergent PWM for PBD selection
library(tidyverse)

# The AAD for a single PWM position equals the
# sum of the absolute frequency deviations for all amino acids, divided by 20; AAD
# ranges from 0.0 (perfect agreement) to 0.1 (maximal divergence). 

################################
## Variables to modify ##
# directory of the repository
path = '~/PL_projects/PL_papers/Scaffold_Letters/Code/Scaffold_PCA/PBD_selection/'
output_dir= './'

# SH3, PDZ or WW
domain_group = 'SH3'
################################

setwd(path)
# import reference position weight matrices from PRM-DB (http://www.prm-db.org/download.php)
# stored in list of data frames 
PWM <- 
  read_rds(paste0(path,'reference_data/PWM_updated.rds'))

# import the reference database from PRM-DB 
master <- read_csv(paste0(path, 'reference_data/PRM_Master.csv'))

# make sure to remove PBDs with less than 30 peptides forming their position 
# weight matrices
master_filt <- master[master$Number_of_Logo_Peptides>30, ]

# select the matrix for the desired domain_group
pwm <- 
  PWM[names(PWM) %in% unlist(master_filt[master_filt$Domain_Group == domain_group, 1])]


# function to compute the frequency of each aa at each position from the pwm
pwm_to_freq <- 
function(pwm){
  tot_weigth <- rowSums(pwm)
  pro <- (apply(pwm, 2, function(x){x/tot_weigth}))
  pro <- as_tibble(pro)
  return(pro)
}


# Set the output object to compare minimal average AAD between 2 matrixes
div_pwm_matrix <- matrix(nrow = length(pwm), ncol = length(pwm))
colnames(div_pwm_matrix) <- names(pwm)
rownames(div_pwm_matrix) <- names(pwm)


# Compare all matrices against all matrices of the same domain group
for (i in 1:length(pwm)) {
  
  for(j in 1:length(pwm)){

    # transform scores to frequency for both PWMs
    pro1 <- pwm_to_freq(pwm[[i]])
    pro2 <- pwm_to_freq(pwm[[j]])
   
    # Select the longer PWM as reference
    if(nrow(pro1) > nrow(pro2)){
      ref <- pro1
      short <- pro2
    } else if (nrow(pro1) <= nrow(pro2)){
      ref <- pro2
      short <- pro1
    } 
    
    # Set the variable used to save all AAD values obtained from 
    # different iteration
    min_diff <- vector(mode = 'numeric', length = 1)
    # Set values used in the loop 
    r=1
    z=0
    # Set the number of all possible comparison between the 2 PWMs
    c <- sum(nrow(short)+nrow(ref)-1)
    
    # Loop through all possible comparison
    for(f in 1:c){
      
      # select the positions of the short PWM to use
      if (f<nrow(short)){
        y <- tail(1:nrow(short), f)
      } else if (f>=nrow(short) & f < nrow(ref)){
        y <-1:nrow(short)
      } else if (f>=nrow(ref)){
        y <- head(1:nrow(short), nrow(short)-z)
        z = z+1
      }
      sub_short <- short[y, ]
      
      # select the positions of the reference PWM to use
      if(f<nrow(short)){
        k <- f
        sub_ref <- ref[1:k, ]
      } else if(f==nrow(short)){
        k <- nrow(short)
        sub_ref <- ref[1:k, ]
      } else if (f > nrow(short) & f <= nrow(ref)){
        k <- tail(1:nrow(short)+r, nrow(ref)-r)
        r=r+1
        sub_ref <- ref[k, ]
      } else if(f > nrow(ref)){
        k <- tail(1:nrow(ref), nrow(sub_short))
        sub_ref <- ref[k, ]
      }
      
      # Compute the AAD between the 2 sub-PWM
      min_diff[f]<- mean(rowSums(abs(sub_ref-sub_short)/20))
    
      }
    
    # remove comparisons with less than half the PWMs overlapping
    b <- ceiling(nrow(short)*0.5)
    
    min_diff <- head(min_diff, -b)
    min_diff <- tail(min_diff, -b)
    
    # Keep the minimum AAD value for this combination of PWMs.
    div_pwm_matrix[i,j] <- min(min_diff, na.rm = T)

  }
}

# write output file for a domain group  
write_rds(div_pwm_matrix, paste0(output_dir, 'AAD_', domain_group,'.rds'))  

  
  
