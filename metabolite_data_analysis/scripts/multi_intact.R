#functions for performing Multi-INTACT

library(aod)
library(INTACT)
library(SQUAREM)

#pi1 function for chisq statistics

.pi1_fun_multi <- function(chisq_vec,df,lambda = 0.5){

  p_vec <- pchisq(chisq_vec,df = df,lower.tail = FALSE)

  p_vec <- p_vec[which(p_vec != 1)]

  pi0 <- length(which(p_vec > lambda))/(length(p_vec)*(1-lambda))

  pi0_max <-  0.99

  pi1 <- 1- min(pi0_max,pi0)

  return(pi1)
}


#Bayes factor functions (from Johnson 2005)

.f_bf <- function(f_vec, df1, df2){

        log10_bf <- (df1 + df2)/2 * log10((df2/df1 + f_vec)/(df2/df1 + 1)) - df1/2 * log10(f_vec)

        return(log10_bf)
}



.chisq_bf <- function(chisq,df){
	
	if(chisq >= df){

	log10_bf <- (df/2) * log10(df/chisq) + log10(exp((chisq - df)/2))
	
	}
	if(chisq < df){
	
	log10_bf <- 0
	
	}
      return(log10_bf)	      
}



#function to compute gene-level X2 statistic from individual-level data

multi_intact_chisq <- function(gene_name,gene_models_file,gene_sbams_file){

	model_wts_df <- read.table(gene_models_file,header=T,sep = '\t')
	model_wts <- as.matrix(model_wts_df[,3:ncol(model_wts_df)])	
	sbam <- read.table(gene_sbams_file,header=F)
        geno <- t(sbam[2:nrow(sbam),4:ncol(sbam)])
        geno[is.na(geno)] <- 0
	geno <- scale(geno)

	#compute predicted molecular gene product levels
	
	pred_gene_prod <- geno %*% model_wts

	pheno = as.numeric(sbam[1,4:ncol(sbam)])
	pheno_cent = pheno - mean(pheno)
	
	#chisq test

	fit <- lm(pheno ~ pred_gene_prod + 0)
        rst <- wald.test(Sigma = vcov(fit), b = coef(fit), Terms = 1:ncol(pred_gene_prod))
        chisq <- as.numeric(rst$result$chi2[1])

	out = c(gene_name, chisq)
	
	return(out)
}


#Function to compute log10(BF) from Servin & Stephens 2007
#Takes:
#Pdat: n x 2 matrix where first column is imputed expression; second is imputed protein 
#y: GWAS trait vector
#sigmas: a 2-vector with the first element prior E-->Y effect sd, followed by prior P-->Y effect sd
#Returns:
#BF_EP, BF_E, and BF_P

log10BF <- function(Pdat,y,sigmas){
	
	sigma_e = sigmas[1]
	sigma_p = sigmas[2]
	n=nrow(Pdat)
	
	X_EP = cbind(rep(1,n),Pdat)
	X_E = cbind(rep(1,n),Pdat[,1])
	X_P = cbind(rep(1,n),Pdat[,2])
	
	invnu_EP = diag(c(0,1/sigma_e^2,1/sigma_p^2))
	invnu_E = diag(c(0,1/sigma_e^2))
	invnu_P = diag(c(0,1/sigma_p^2))
	
	invOmega_EP = invnu_EP + t(X_EP) %*% X_EP
	invOmega_E = invnu_E + t(X_E) %*% X_E
	invOmega_P = invnu_P + t(X_P) %*% X_P
	
	B_EP = solve(invOmega_EP, t(X_EP) %*% cbind(y))
	B_E = solve(invOmega_E, t(X_E) %*% cbind(y))
	B_P = solve(invOmega_P, t(X_P) %*% cbind(y))
	invOmega0 = n
	
	BF_EP = -0.5*log10(det(invOmega_EP)) + 0.5*log10(invOmega0) - log10(sigma_e) - log10(sigma_p) -(n/2) * (log10( t(y- X_EP %*% B_EP) %*% y)- log10(t(y) %*% y - n*mean(y)^2))
	BF_E = -0.5*log10(det(invOmega_E)) + 0.5*log10(invOmega0) - log10(sigma_e)-(n/2) * (log10( t(y- X_E %*% B_E) %*% y)- log10(t(y) %*% y - n*mean(y)^2))
	BF_P = -0.5*log10(det(invOmega_P)) + 0.5*log10(invOmega0) - log10(sigma_p)-(n/2) * (log10( t(y- X_P %*% B_P) %*% y)- log10(t(y) %*% y - n*mean(y)^2))
	
	return(c(BF_EP,BF_E,BF_P))
}



#wrapper function to compute Stage 2 Bayes factors from individual-level data
#the first version will work with 2 gene product prediction models

multi_intact_stage2 <- function(gene_name,gene_models_file,gene_sbams_file,
				sigma_e_vec = c(0.4,0.8,1.6),
				sigma_p_vec = c(0.4,0.8,1.6)){
	
	model_wts_df <- read.table(gene_models_file,header=T,sep = '\t')
        model_wts <- as.matrix(model_wts_df[,3:ncol(model_wts_df)])
        sbam <- read.table(gene_sbams_file,header=F)
        geno <- t(sbam[2:nrow(sbam),4:ncol(sbam)])
        geno[is.na(geno)] <- 0
        geno <- scale(geno)

        #compute predicted molecular gene product levels

        pred_gene_prod <- geno %*% model_wts

        pheno <- as.numeric(sbam[1,4:ncol(sbam)])
	pheno_cent <- pheno - mean(pheno)

	sigma_grid <- cbind(sigma_e_vec,sigma_p_vec)

	#Apply BF computation over each sigma pair
	log10BF_grid = apply(X = sigma_grid,FUN = log10BF,MARGIN = 1, 
			     Pdat = pred_gene_prod,
			     y = pheno_cent)

	#Average over grid
	log10BFs <- apply(log10BF_grid,FUN=function(bf_vec){
				  x_star <- max(bf_vec)
			       	  logsumexp <- x_star + log10(sum(10^(bf_vec-x_star)))
			      	  out <- logsumexp - log10(length(bf_vec))
			      	  return(out)
			    	},
			    	MARGIN = 1)
	
        return(log10BFs)
}







#function that outputs Stage 1 posterior probabilities. Order of genes in chisq_vec must match the rows of coloc_rst. 

multi_intact_scan <- function(coloc_rst,chisq_vec = NULL,
			 chisq_vec_dof = n_gene_prod,	
			 prior_fun = linear,
			 t = 0.05, D = NULL,
			 xwas_priors = .pi1_fun_multi(chisq_vec,df=chisq_vec_dof),
			 xwas_BFs = NULL,
			 bf_type = "wakefield",
			 K = c(1,2,4,8,16),
			 glcp_aggreg = "max"){
	

	if (glcp_aggreg == "max"){

#		coloc_rst$glcp_aggregated <- apply(coloc_rst[,2:ncol(coloc_rst)],MARGIN = 1,FUN = max)
		coloc_rst$glcp_aggregated <- do.call(pmax, coloc_rst[,2:ncol(coloc_rst)])

	}

	if (glcp_aggreg == "one_minus_prod_one_minus"){
	
		coloc_rst$glcp_aggregated <- apply(coloc_rst[,2:ncol(coloc_rst)],
					       MARGIN = 1,
					       FUN = function(x){
						       out <- 1-prod(1-x)
					       return(out)})
		
	}

	if (!(glcp_aggreg %in% c("max","one_minus_prod_one_minus"))){

		stop("glcp_aggreg must be either max or one_minus_prod_one_minus")
		
	}


	 #Compute prior from GLCPs

        if (length(t) == 0 & length(D) == 0){
		
		prior <- prior_fun(GLCP = coloc_rst$glcp_aggregated,u = xwas_priors)
	
	}

	if (length(t) == 0 & length(D) > 0){

		prior <- prior_fun(GLCP = coloc_rst$glcp_aggregated,D = D,u = xwas_priors)
	
	}

	if (length(t) > 0 & length(D) == 0){

		prior <- prior_fun(GLCP = coloc_rst$glcp_aggregated,t = t,u = xwas_priors)

	}
	
	if (length(t) > 0 & length(D) > 0){

		prior <- prior_fun(GLCP = coloc_rst$glcp_aggregated,D = D,t = t,u = xwas_priors)
	
	}

	if (length(chisq_vec != 0)){

		coloc_rst$chisq <- chisq_vec

		coloc_rst$chisq_dof <- chisq_vec_dof

		coloc_rst$xwas_pval <- pchisq(chisq_vec,df = chisq_vec_dof,lower.tail = F)
		
		if (bf_type == "johnson"){

			coloc_rst$scan_log10bf <- vapply(coloc_rst$chisq, 
							   FUN = .chisq_bf, 
							   FUN.VALUE = numeric(1), df = chisq_vec_dof)

		}
		if (bf_type == "wakefield"){

			z_vec <- qnorm(coloc_rst$xwas_pval/2, lower.tail = F)

			bf_grid_log10 <- vapply(X=matrix(K),
						FUN.VALUE=numeric(length(z_vec)),
						FUN = function(K,z_vec){
							out <- 0.5*log10(1/(1+K)) + 
								(0.5*(z_vec^2) * (K/(1+K)))*log10(exp(1))
							return(out)
						},
				   z_vec=z_vec)

			coloc_rst$scan_log10bf <- apply(bf_grid_log10,FUN=function(bf_vec){
						  x_star <- max(bf_vec)
					    	  logsumexp <- x_star + log10(sum(10^(bf_vec-x_star)))
					    	  out <- logsumexp - log10(length(bf_vec))
					    	  return(out)
			       	   },
			       	   MARGIN = 1)

		}
		if (!(bf_type %in% c("johnson","wakefield"))){
		
			stop("bf_type must be either wakefield or johnson")
		
		}

		#compute log10(posterior odds)

		log10odds <- coloc_rst$scan_log10bf + log10((prior/(1-prior)))

		#compute posterior probability
    
		post <- 1/ (1 + 10^(-log10odds))

		coloc_rst$posterior <- post

		out <- coloc_rst[order(coloc_rst$posterior,decreasing=T),]

		return(out)


	}else{ #If XWAS Bayes factors are supplied, compute posteriors directly.

    	#Compute prior from GLCPs
    
	    	#compute log10(posterior odds)
    
		log10odds <- log10(xwas_BFs) + log10((prior/(1-prior)))

	    	#compute posterior probability

	    	post <- 1/ (1 + 10^(-log10odds))

		coloc_rst$posterior <- post

		out <- coloc_rst[order(coloc_rst$posterior,decreasing=T),]

		return(coloc_rst)
	
	}
		
}




#main function that outputs Stage 1 posterior probabilities and, if specified, Stage 2 BFs

multi_intact <- function(n_gene_prod,coloc_rst,chisq_vec = NULL,              
				stage2 = FALSE,
				stage2_rst = NULL,
				prior_fun = linear,
				t = 0.05, D = NULL,
			       	xwas_priors = .pi1_fun_multi(chisq_vec,df=n_gene_prod),
				xwas_BFs = NULL,
				alpha = 0.05,
				bf_type = "wakefield",
				glcp_aggreg = "max"){
	
	stage1_output = multi_intact_stage1(n_gene_prod = n_gene_prod,
					    coloc_rst = coloc_rst,
					    chisq_vec = chisq_vec,
					    prior_fun = prior_fun,
			       		    t = t,
			       		    D = D,
			       		    xwas_priors = xwas_priors,xwas_BFs = xwas_BFs,
					    bf_type = bf_type,
					    glcp_aggreg = glcp_aggreg)

	
	if(stage2 == FALSE){

		return(stage1_output)
	}

	if(stage2 == TRUE){
		
		if(n_gene_prod != 2){
		
			stop("Multi-INTACT Stage 2 currently only supports 2 predicted gene product types")
		}
		if(n_gene_prod == 2 & is.null(stage2_rst)){
		
			stop("stage2_rst must be supplied if stage2 == TRUE")
		
		}else{

			stage1_output$fdr_sig <- fdr_rst(posterior = stage1_output$posterior, 
						  alpha = alpha)$sig
		
			out <- merge(stage1_output, stage2_rst,by="gene")

			#out$BF_E[out$fdr_sig==FALSE] <- NA
			#out$BF_P[out$fdr_sig==FALSE] <- NA
			#out$BF_EP[out$fdr_sig==FALSE] <- NA

			return(out[order(out$posterior,decreasing=T),])
		}
		
	}else{
	
		stop("stage2 must be set to TRUE or FALSE")
	}

}


fdr_rst2 <- function (posterior, alpha = 0.05)
{
    gene_num <- seq(1, length(posterior))
    lfdr_sort = sort(1 - posterior)
    FDR = cumsum(lfdr_sort)/(1:length(lfdr_sort))
    thresh = 1 - lfdr_sort[max(which(FDR <= alpha))]
    rej_gene = as.numeric(gene_num[which(posterior >= thresh)])
    out_tmp <- rep(FALSE, length(posterior))
    out_tmp[rej_gene] <- TRUE
    out <- data.frame(posterior = posterior, sig = out_tmp)
    return(out)
}



###################################
##EM with coloc, estimating pi0 with qvalue
###################################

bf.weighted_sum_coloc<-function(w,bf,fp_coloc,i){
        K <- length(w)
        fp_coloc <- fp_coloc[i]
        bf.sum <- 0
        bf.gene <- bf[((i-1)*K+1):(i*K)]
        bf.m <- max(bf.gene)
        bf.null <- bf.gene[1]
        bf.alt <- bf.gene[2:length(bf.gene)]
        bf.sum <- sum(c((w[1]*fp_coloc + 1 - fp_coloc)*exp(bf.null-bf.m),
                        w[2:K]*fp_coloc*exp(bf.alt-bf.m)))
        return(bf.m+log(bf.sum))
}

#bf.weighted_sum_coloc(w = pi_start, bf = bf_df$BF, fp_coloc = df$fp_coloc,i = 1)

bf.loglik_coloc<-function(w,bf,fp_coloc){
  K<-length(w)
  n<-length(bf)/K
  loglik<-0
  sumloglik<-0
  for (i in 1:n){
          bfs <- bf[((i-1)*K+1):(i*K)]
          med <- c(exp(log(w[1]*fp_coloc[i]+ 1 - fp_coloc[i]) + bfs[1] - bf.weighted_sum_coloc(w,bf,fp_coloc,i)),
                   exp(log(w[2:K]*fp_coloc[i]) + bfs[2:K] - bf.weighted_sum_coloc(w,bf,fp_coloc,i)))
          loglik <- sum(c(log(w[1]*fp_coloc[i]+ 1 - fp_coloc[i])*med[1], log(w[2:K])*med[2:K]))
          #loglik <- sum(c(log(w[1]*fp_coloc[i]+ 1 - fp_coloc[i])*med + bf[1]*med, log(w[2:K])*med + bf[2:K]*med))
    sumloglik<-sumloglik+loglik
  }
  return(sumloglik)
}

#bf.loglik_coloc(w = pi_start, bf = bf_df$BF, fp_coloc = df$fp_coloc)


bf.em_coloc_pi0<-function(w,bf,fp_coloc){
  K<-length(w)
  wnew<-rep(NA,K)
  n<-length(bf)/K
  #E-step
  ei<-matrix(NA,n,K)
  for (i in 1:n){
          bfs <- bf[((i-1)*K+1):(i*K)]
          ei[i,] <- c(exp(log(w[1]*fp_coloc[i]+ 1 - fp_coloc[i]) + bfs[1] - bf.weighted_sum_coloc(w,bf,fp_coloc,i)),
          exp(log(w[2:K]*fp_coloc[i]) + bfs[2:K] - bf.weighted_sum_coloc(w,bf,fp_coloc,i)))
  }
  #M-step
  wnew[1] <- w[1]
  diff <- 1 - wnew[1]
  tmp <- colSums(ei[,2:K],na.rm = T)
  scaling <- sum(tmp)/diff
  #if(scaling != 0){
  wnew[2:K] <- tmp/scaling
  #}
  #if(scaling == 0){
  #wnew[2:K] <- tmp
  #}
  return(wnew)
}

#bf.em_coloc_pi0(w = pi_start, bf = bf_df$BF, fp_coloc = df$fp_coloc)


class_posteriors<-function(w,bf,fp_coloc){
  K<-length(w)
  wnew<-rep(NA,K)
  n<-length(bf)/K
  #E-step
  ei<-matrix(NA,n,K)
  for (i in 1:n){
          bfs <- bf[((i-1)*K+1):(i*K)]
          ei[i,] <- c(exp(log(w[1]*fp_coloc[i]+ 1 - fp_coloc[i]) + bfs[1] - bf.weighted_sum_coloc(w,bf,fp_coloc,i)),
          exp(log(w[2:K]*fp_coloc[i]) + bfs[2:K] - bf.weighted_sum_coloc(w,bf,fp_coloc,i)))
  }
  return(ei)
}

wakefield_bf_z_ln <- function(z_vec, K = c(1,2,4,8,16)){

     bf_grid_log <- vapply(X=matrix(K),FUN.VALUE=numeric(length(z_vec)),
                           FUN = function(K,z_vec){
      out <- 0.5*log(1/(1+K)) + (0.5*(z_vec^2) * (K/(1+K)))
      return(out)
    },
    z_vec=z_vec)


    #Perform Bayesian model averaging using log sum exp trick

    bf_log <- apply(bf_grid_log,FUN=function(bf_vec){
      x_star <- max(bf_vec)
      logsumexp <- x_star + log(sum(exp(bf_vec-x_star)))
      out <- logsumexp - log(length(bf_vec))
      return(out)
    },
    MARGIN = 1)

    return(bf_log)
}

