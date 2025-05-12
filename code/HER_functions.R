# load packages
library(parallel)
library(dplyr)

# z-score normalization function
znorm <- function(dist){((dist - mean(dist))/sd(dist))}

# SE correcton factor
se_correct <- function(se,h2,n){
  factor <- sqrt(n) * (0.020225*h2 + 0.004225) + (-0.2352*h2 +0.9467)
  se_corrected <- se * factor
  return(se_corrected)
}

# Get pairwise products of trait values, as well as IBS relatedness from genotype probability-derived kinship matrix, for HE regression
pairwise_product2 <- function(t1,t2,k,type){
  #start <- Sys.time()
  kin <- k[[1]] # just need the kinship matrix for one autosome
  
  # make list of combinations to check
  ids1 <- names(t1)
  ids2 <- names(t2)
  
  kin_sub <- kin[rownames(kin) %in% ids1, colnames(kin) %in% ids2] # subset kinship matrix to only include individuals present in a sample of the data (this is written specifically for SE simulations).  When calculating cross-products normally, the two traits are already subset to only include phenotype information for individuals that are present in the kinship matrix. Thus, for normal h2 or GC calculations this line should not affect anything. 
  kin_sub <- kin_sub[match(names(t1),rownames(kin_sub)),]
  kin_sub <- kin_sub[,match(names(t2),colnames(kin_sub))]
  
  combs <- expand.grid(ids1,ids2)
  temp <- outer(t1,t2,FUN="*")
  if(type == "h"){
    ind <- which(upper.tri(temp,diag=F) , arr.ind = TRUE )
    
    out_df <- data.frame( col = dimnames(temp)[[2]][ind[,2]] ,
                          row = dimnames(temp)[[1]][ind[,1]] ,
                          prod = temp[ ind ] )
    
    ind2 <- which(upper.tri(kin_sub,diag=F) , arr.ind = TRUE )
    
    
    kin_sub2 <- data.frame( col = dimnames(kin_sub)[[2]][ind[,2]] ,
                            row = dimnames(kin_sub)[[1]][ind[,1]] ,
                            relatedness = kin_sub[ ind ] )
    
    out_df <- cbind(out_df,kin_sub2$relatedness)
    
  }else{
    out_df <-cbind(combs,as.vector(temp),as.vector(kin_sub))
  }
  
  colnames(out_df) <- c("id1","id2","prod","relatedness")
  out_df$id1 <- as.character(out_df$id1)
  out_df$id2 <- as.character(out_df$id2)
  out_df <- out_df[out_df$id1 != out_df$id2,] # exclude rows containing self-self comparisons
  
  rownames(out_df) <- 1:nrow(out_df)
  
  #end <- Sys.time()
  #end-start
  return(out_df)
}


# modified HE regression using the pairwise products of trait values...
# output is either narrow-sense heritability (h2) for a single trait... 
# or genetic correlation for a pair of traits

# pheno = dataframe of z-normalized residuals after fitting models that account for covariates in each trait of interest. 
#The first column of pheno should be mouse names/IDs, and each subsequent column should contain normalized residual trait values.
HEreg <- function(pheno,k,h2_df,type=c("h","gc"),traits,cores,correct=T){
  names <- colnames(pheno)[2:ncol(pheno)]
  if(is.null(traits)){
    if(type == "h"){
      pairs <- t(sapply(1:length(names),function(x){rep(names[x],times=2)}))
    }else{
      #herit_df <- read.table(h2_df,sep='\t',header=T)
      herit_df <- h2_df
      pairs <- t(combn(names,2))
    }
  }else{
    pairs <- traits
    herit_df <- h2_df
  }
  
  cor_df <- mclapply(1:nrow(pairs),function(x){
    #print(x)
    if(x %% 1000 == 0){print(x)}
    start <- Sys.time()
    # retreive names of traits used for HE regression
    trait1 <- pairs[x,1]
    trait2 <- pairs[x,2]
    
    # access trait values from phenotype dataframe
    T1 <- pheno[is.na(pheno[,trait1]) == F,trait1]
    names(T1) <- pheno[is.na(pheno[,trait1]) == F,1]
    T2 <- pheno[is.na(pheno[,trait2]) == F,trait2]
    names(T2) <- pheno[is.na(pheno[,trait2]) == F,1]
    
    n_1 <- length(T1)
    n_2 <- length(T2)
    
    # there are no non-na values, return na for the output
    if(length(T1) == 0 | length(T2) == 0){
      out <- c(trait1,trait2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,n_1,n_2)
      return(out)
    }
    
    # subset to only include individuals in the kinship matrix
    T1 <- T1[which(names(T1) %in% colnames(k[[1]]))]
    T2 <- T2[which(names(T2) %in% colnames(k[[1]]))]
    
    # pull shared values to calculate phenotypic correlation
    shared <- intersect(names(T1),names(T2))
    if(length(shared) >= 3){
      T1_sub <- T1[`shared`]
      T2_sub <- T2[`shared`]
      pcor <- cor(T1_sub,T2_sub)
      pcor_p <- cor.test(T1_sub,T2_sub)$p.value
    }else{
      pcor <- NA
      pcor_p <- NA
    }
    
    # get trait product and relatedness for each pair of animals
    #data <- pairwise_product(t1=T1,t2=T2,k=k)
    data <- pairwise_product2(t1=T1,t2=T2,k=k,type=type)
    
    fit <- lm(data$prod ~ data$relatedness) 
    gcor_p <- summary(fit)$coefficients[2,4]
    b <- fit$coefficients[2]
    
    ob <- b/(2*var(T1)) # because all traits are z-normalized in advance, the variance of both traits is going to be 1, so it doesn't matter which trait is used. 
    
    if(type == "h"){
      h2 <- ob
      h2_trait1 <- h2
      h2_trait2 <- h2
    }else{
      h2 <- NA
      h2_trait1 <- herit_df[herit_df$trait1 == trait1,"h2"]
      h2_trait2 <- herit_df[herit_df$trait1 == trait2,"h2"]
      se_h2_trait1 <- herit_df[herit_df$trait1 == trait1,"se_h2"]
      se_h2_trait2 <- herit_df[herit_df$trait1 == trait2,"se_h2"]
    }
    
    gcor <- ob/(sqrt(h2_trait1)*sqrt(h2_trait2))
    
    if(type == "gc"){
      #ss_reduction <-  1/(1+(2*ob*mean(data$relatedness*2)))
      
      b_error <- sqrt(sum(((data$prod - mean(data$prod))^2)) / (nrow(data) - 2)) / sqrt(sum(((data$relatedness - mean(data$relatedness))^2)))
      
      ob_error <- b_error/(2*var(T1))
      
      if(correct==T){
        ob_error <- se_correct(n=sqrt(n_1*n_2),h2=ob,se=ob_error)
      }else{
        ob_error <- ob_error
      }
      
      
      # source: http://www.met.rdg.ac.uk/~swrhgnrj/combining_errors.pdf
      #se_gc <- sqrt( (ob_error/ob)^2 + (.5*sqrt(se_h2_trait1) / sqrt(h2_trait1))^2 + (.5*sqrt(se_h2_trait2) / sqrt(h2_trait2))^2 ) #old, wrong version
      se_gc <- sqrt( (ob_error/ob)^2 + ((.5*sqrt(h2_trait1)*se_h2_trait1) / sqrt(h2_trait1))^2 + ((.5*sqrt(se_h2_trait2)*se_h2_trait2) / sqrt(h2_trait2))^2 )
      
      se_h2 <- NA
    }else{
      #b_error <- summary(fit)$coefficients[2,2]
      
      #b_error <- sqrt(sum(fit$residuals^2) / (sum( (data$relatedness - mean(data$relatedness))^2 ) * ((mean(c(length(T1),length(T2)))/2) - 2)))
      
      #b_error <- sqrt(sum(((data$prod - mean(data$prod))^2)) / ((mean(c(length(T1),length(T2)))/2) - 2)) / sqrt(sum(((data$relatedness - mean(data$relatedness))^2)))
      
      #ss_reduction <-  1/(1+(2*h2*mean(data$relatedness*2)))
      #ss_reduction <- 1
      
      b_error <- sqrt(sum(((data$prod - mean(data$prod))^2)) / (nrow(data) - 2)) / sqrt(sum(((data$relatedness - mean(data$relatedness))^2)))
      
      ob_error <- b_error/(2*var(T1))
      
      se_h2 <- ob_error
      
      if(correct==T){
        se_h2 <- se_correct(n=n_1,h2=h2,se=se_h2)
      }else{
        se_h2 <- se_h2
      }

      # source: http://www.met.rdg.ac.uk/~swrhgnrj/combining_errors.pdf
      # https://www.met.reading.ac.uk/~swrhgnrj/
      # https://scholar.google.com/citations?hl=en&user=-vUyGWsAAAAJ&view_op=list_works&sortby=pubdate
      se_gc <- sqrt( (ob_error/ob)^2 + ((.5*sqrt(h2)*se_h2) / sqrt(h2))^2 + ((.5*sqrt(se_h2)*se_h2) / sqrt(h2))^2)
      
      #se_gc <- NA
    }
    out <- c(trait1,trait2,gcor,gcor_p,pcor,pcor_p,h2,se_gc,se_h2,ob,h2_trait1,h2_trait2,n_1,n_2)
    
    end <- Sys.time()
    #print(end-start)
    return(out)
  },mc.cores=cores)
  
  cor_df <- do.call('rbind',cor_df)
  
  colnames(cor_df) <- c("trait1","trait2","gcor","gcor_p","pcor","pcor_p","h2","se_gc","se_h2","ob","h2_trait1","h2_trait2","n_1","n_2")
  
  cor_df <- data.frame(cor_df)
  for(i in 3:ncol(cor_df)){
    cor_df[,i] <- as.numeric(cor_df[,i])
  }
  for(i in 3:(ncol(cor_df)-2)){
    cor_df[,i] <- round(cor_df[,i],3)
  }
  return(cor_df)
}