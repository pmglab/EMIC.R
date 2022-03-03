
#Function: Compute the effective median value.

#Input parameters: 
#1) a sorted estimated causal effect sizes of the snps ("es" in below function);

#2) a sorted LD matrix (n*n) of the snps in a gene according to estimated causal effect size(see 1) above)("ld" in below function);

estimate_Effectively_independent_number <- function(es,ld) {
  
  ## The correlation of estimates at two loci
  Cor_es <- 0.4428*(ld^2) + 0.5665*abs(ld)
  diag(Cor_es)<-1
  
  ## Effectively independent number 
  eff_ind_num <- c()
  eff_ind_num[1] <- 1
  len <- dim(ld)[0]
  i <- 1
  
  for (i in 2:len){
  lamb <- eigen(Cor_es[1:i,1:i])$values
  num <- 0
  
  for (j in 1:i){
    if (lamb[j,j]>=1){
      num <- num + 1
    }else{
      num <- num + lamb[j,j]
    }
  }
  
  i <- i+1
  }
  
  eff_ind_num[i] <- num
  
  eff_median <- eff_ind_num/2

  index_i <- max(which(eff_ind_num<eff_median))
  index_j <- min(which(eff_ind_num>eff_median))
  
  M <- es[index_i] + (es[index_i]-es[index_j])*(eff_median - eff_ind_num[index_i])/(eff_ind_num[index_j] - eff_ind_num[index_i])
  
  return(M)
}
