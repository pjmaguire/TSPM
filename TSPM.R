#auther=Peter_Maguire
#contact=pmaguire@stanford.edu

library(edgeR)

#Turns off scientific notation in R
#options(scipen = 999)

TSPM <- function(counts, x1, x0, method = c("TMM", "lib_size"), offset, alpha.wh=0.05){
  #Description Not Updated Yet (Sorry)
  #Based off R code for "A Two-Stage Poisson Model for Testing RNA-Seq Data" by Paul L. Auer and R.W. Doerge
  
  ## Input:
  #counts: 	a matrix of RNA-Seq gene counts (genes are rows, samples are columns)
  #x1: 		a vector of treatment group factors (under the alternative hypothesis)
  #x0: 		a vector of treatment group factors (under the null hypothesis)
  #lib.size: 	a vector of RNA-Seq library sizes. This could simply be obtained
  #          	by specifying lib.size <- apply(counts,2,sum). It may also be any other
  #          	appropriate scaling factor.
  #alpha.wh:	the significance threshold to use for deciding whether a gene is overdispersed.
  #               Defaults to 0.05.
  
  ## Output:
  #log.fold.change:		a vector containing the estimated log fold changes for each gene
  #pvalues: 			a vector containing the raw p-values testing differential expression for each gene.
  #index.over.disp: 	a vector of integer values containing the indices of the over-dispersed genes.
  #index.not.over.disp:	a vector of integer values containing the indices of the non-over-dispersed genes.
  #padj:			a vector containing the p-values after adjusting for multiple testing using the 
  #				method of Benjamini-Hochberg
  
  
  ######## The main loop that fits the GLMs to each gene ########################
  
  # counts = total_counts[, HCM_CvP_columns]
  # x1=gl(n=2, k = length(HCM_CvP_columns)/2, labels = c("C", "T"))
  # x0=rep(1, times = length(HCM_CvP_columns))
  # method="TMM"
  # offset=TMM_norm[HCM_CvP_columns]
  # alpha.wh=0.05
  
  ### Initializing model parameters ####
  n <- dim(counts)[1]
  per.gene.disp <- NULL
  LRT <- NULL
  score.test <- NULL
  LFC <- NULL
  
  #Allows for library size to be used as the offset in the Poisson regression
  if(method == "lib_size"){
    offset = log(offset)
  }
  
  ###### Fitting the GLMs for each gene #################
  for(i in 1:n){
    ### Fit full and reduced models ###
    model.1 <- glm(as.numeric(counts[i,]) ~ x1, offset=offset, family=poisson)
    model.0 <- glm(as.numeric(counts[i,]) ~ x0, offset=offset, family=poisson)
    
    ### Obtain diagonals of Hat matrix from the full model fit ###
    hats <- hatvalues(model.1)
    
    ### Obtain Pearson overdispersion estimate ####
    per.gene.disp[i] <- sum(residuals(model.1, type="pearson")^2)/model.1$df.residual
    
    ### Obtain Likelihood ratio statistic ####
    LRT[i] <- deviance(model.0)-deviance(model.1)
    
    ### Obtain score test statistic ####
    score.test[i] <- 1/(2*length(counts[i,])) * sum(residuals(model.1, type="pearson")^2 - ((counts[i,] - hats*model.1$fitted.values)/model.1$fitted.values))^2
    
    ### Obtain the estimated log fold change ###
    LFC[i] <- -model.1$coef[2]
  }
  
  ## Initialize parameters for Working-Hotelling bands around the score TSs ###
  qchi <- qchisq(df=1, (1:n-0.5)/n)
  MSE <- 2
  UL <- NULL
  
  #### Obtain the upper boundary of the WH bands #######################################
  xbar <- mean(qchi)
  bottom <- sum((qchi-xbar)^2)
  top <- (qchi-xbar)^2
  s <- sqrt(MSE*(1/n) + (top/bottom))
  W <- sqrt(2*qf(df1=1, df2=n-1, p=1-(alpha.wh/n)))
  UL <- pmax(qchi + W*s,1)
  
  ###### Obtain the indices of the over-dispersed and not-over-dispersed genes, respectively ##########
  
  cutoff <- min(which(sort(score.test)-UL > 0))
  temp <- cutoff-1 + seq(cutoff:length(score.test))
  over.disp <- which(score.test %in% sort(score.test)[temp])
  not.over.disp <- setdiff(1:length(score.test), over.disp)
  
  ###### Compute p-values ####################################
  p.f <- pf(LRT[over.disp]/per.gene.disp[over.disp], df1=1, df2=model.1$df.residual, lower.tail=FALSE)
  p.chi <- pchisq(LRT[not.over.disp], df=1, lower.tail=FALSE)
  p <- NULL
  p[over.disp] <- p.f
  p[not.over.disp] <- p.chi
  
  ##### Adjust the p-values using the B-H method ####################
  p.bh.f <- p.adjust(p.f, method="BH")
  p.bh.chi <- p.adjust(p.chi, method="BH")
  final.p.bh.tagwise <- NULL
  final.p.bh.tagwise[over.disp] <- p.bh.f
  final.p.bh.tagwise[not.over.disp] <- p.bh.chi
  
  ### Output ###
  #Output results into a list of lists
  #list(log.fold.change=LFC, pvalues=p, index.over.disp=over.disp, index.not.over.disp=not.over.disp,
  #     padj=final.p.bh.tagwise)
  
  #Converts separate overdispersed indexes vectors into a single True/False vector
  true_false = vector()
  for(i in 1:length(LFC)){
    if(i %in% over.disp){
      true_false[i] = TRUE
    }
    else{
      true_false[i] = FALSE
    }
  }
  
  #Outputs the results in a dataframe
  data.frame(gene_id = row.names(counts), log_fold_change = LFC, p_values = p, over_dispersed = true_false,
             q_values = final.p.bh.tagwise)
}