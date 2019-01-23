library(knitr)


pairwise_association_test = function(mx2,cdataNames,edataNames,dataName='out',dir_out = 'output'
  , methods = rep('lm',length(edataNames)), n_decimal = 3) {

  n_cdata = length(cdataNames)
  n_edata = length(edataNames)
  mx_pval = matrix(NA, n_cdata, n_edata)
  rownames(mx_pval) = cdataNames
  colnames(mx_pval) = edataNames
  
  mx_beta = mx_pval
  mx_qval = mx_pval
  
  
  pdf(file=sprintf('%s/scatterplot_%s.pdf',dir_out,dataName),paper='special',width=6,height=8)
  layout(matrix(1:6,3,2,byrow=TRUE))
  
  
  for(j in 1:length(edataNames)){
  
    for(i in 1:length(cdataNames)){
      # i = 1; j =1;
      edataName = edataNames[j]
      cdataName = cdataNames[i]
      
      ### need to fix this so that the cvrt name will be displayed in the resultsing tables
      model_design = sprintf('%s ~ %s', edataName, cdataName)

      if(methods[j] == 'lm'){
      	#lm1 = lm(model_design, na.action = na.omit,data=mx2)
      	lm1 = lm(mx2[,edataName] ~ mx2[,cdataName], na.action = na.omit)

  	  }
      sum1 = summary(lm1)
      plot(mx2[,cdataName],mx2[,edataName],xlab=cdataName,ylab=edataName)
      if(class(mx2[,cdataName]) == 'numeric'){
	      if(var(mx2[,cdataName], na.rm = TRUE) !=0){
	        abline(lm1)
	      }
	  }

      ### change the row.names

      
      ### if p-value valid such (when var(X) == 0, p-value not exist)
      if(dim(sum1$coefficients)[1] == 2){
      	rownames(sum1$coefficients) = c('Intercept',cdataName)
        mx_pval[i,j] = sum1$coefficients[2,4]
        mx_beta[i,j] = sum1$coefficients[2,1]

      } else if(class(mx2[,cdataName]) == 'factor'){
      	levelNames = gsub('mx2\\[, cdataName\\]','',rownames(sum1$coefficients))[-1]
      	rownames(sum1$coefficients) = c('Intercept',paste(cdataName,'.',levelNames,sep=''))
        f <- sum1$fstatistic 
        mx_pval[i,j] = pf(f[1], f[2], f[3], lower=FALSE)
        mx_beta[i,j] = NA

      }

      ### add R2, confidence interval and df
      cfs = confint(lm1)
      #cfs = round(cfs, 3)
      df1 = sum1$df[2]
      var_X = c()

      idx_keep = which(!is.na(mx2[,edataName]) & !is.na(mx2[,cdataName]))
      if(is.numeric(mx2[,cdataName]) ) {
      	var_X = var(mx2[idx_keep,cdataName], na.rm = TRUE)
      	var_X = c(0,var_X)
      } else{
      	xx = mx2[idx_keep,cdataName]
      	#subdf$letters <- factor(subdf$letters)
      	xx <- factor(xx)
      	X1 = model.matrix(~xx)
      	var_X = apply(X1, 2, var, na.rm = TRUE)
      }

      coefs = sum1$coefficients[,1]

      if(length(coefs) == 1) {
      	coefs = c(coefs, 0)
      	sum1$coefficients = rbind(sum1$coefficients,NA)
      }

      R2 = coefs^2 * var_X/(sum(coefs^2 * var_X) + sum1$sigma^2)
      #R2 = round(R2,3)

      #cfs[,1] = format(cfs[,1], digits= 3)
      #cfs[,2] = format(cfs[,2], digits= 3)
    
      #sum1$coefficients[,1] = format(sum1$coefficients[,1], digits = 3)
      #R2 = format(R2, digits = 3)

      #if(flag_round == TRUE){
      #  mx_out = cbind(sum1$coefficients, round(cfs,n_decimal), round(R2, n_decimal))
      #} else{
        mx_out = cbind(sum1$coefficients, cfs, R2)
        for(jj in 1:dim(mx_out)[2]) {
          mx_out[,jj] = format(as.numeric(mx_out[,jj]), digits = n_decimal)
        }

      #}
      #mx_beta[i,j] = sum1$coefficients[2,1]
      #rownames(sum1$coefficients) = c('Intercept',cdataName)
      print(kable(mx_out,caption=paste(dataName,': ', edataName,' vs ',cdataName, ', df=',df1,sep='')))
      
    }
    mx_qval[,j] = p.adjust(mx_pval[,j],method='BH')
  
  }


  fileout_pval = sprintf('%s/pval_%s.csv',dir_out,dataName)
  write.csv(mx_pval,file=fileout_pval, row.names=TRUE)
  
  fileout_beta = sprintf('%s/beta_%s.csv',dir_out,dataName)
  write.csv(mx_beta,file=fileout_beta, row.names=TRUE)
  
  fileout_qval = sprintf('%s/qval_%s.csv',dir_out,dataName)
  write.csv(mx_qval,file=fileout_qval, row.names=TRUE)
  
  
  dev.off()


}



library(survival)
pairwise_association_test_surv = function(mx2,cdataNames,edataNames,dataName='out',dir_out = 'output', methods = rep('lm',length(edataNames))) {
  
  n_cdata = length(cdataNames)
  n_edata = length(edataNames)
  mx_pval = matrix(NA, n_cdata, n_edata)
  rownames(mx_pval) = cdataNames
  colnames(mx_pval) = edataNames
  
  mx_beta = mx_pval
  mx_qval = mx_pval
  
  
  pdf(file=sprintf('%s/scatterplot_%s.pdf',dir_out,dataName),paper='special',width=6,height=8)
  layout(matrix(1:6,3,2,byrow=TRUE))
  
  
  for(j in 1:length(edataNames)){
    
    for(i in 1:length(cdataNames)){
      # i = 1; j =1;
      edataName = edataNames[j]
      cdataName = cdataNames[i]
      
      ### need to fix this so that the cvrt name will be displayed in the resultsing tables
      model_design = sprintf('%s ~ %s', edataName, cdataName)
      
      if(methods[j] == 'lm'){
        #lm1 = lm(model_design, na.action = na.omit,data=mx2)
        lm1 = lm(mx2[,edataName] ~ mx2[,cdataName], na.action = na.omit)
        
      } else if(methods[j] == 'survival'){
        status = as.numeric(mx2[,edataName] != 12)
        lm1 = coxph(Surv(mx2[,edataName], status) ~ mx2[,cdataName])
        sum1 = summary(lm1)
      }
      sum1 = summary(lm1)
      plot(mx2[,cdataName],mx2[,edataName],xlab=cdataName,ylab=edataName)
      if(var(mx2[,cdataName], na.rm = TRUE) !=0){
        abline(lm1)
      }
      
      ### change the row.names
      
      
      ### if p-value valid such (when var(X) == 0, p-value not exist)
      if(methods[j] == 'survival'){
        rownames(sum1$coefficients) = c(cdataName)
        mx_pval[i,j] = sum1$coefficients[5]
        mx_beta[i,j] = sum1$coefficients[1]
        
      } else if(dim(sum1$coefficients)[1] == 2){
        rownames(sum1$coefficients) = c('Intercept',cdataName)
        mx_pval[i,j] = sum1$coefficients[2,4]
        mx_beta[i,j] = sum1$coefficients[2,1]
        
      } else if(class(mx2[,cdataName]) == 'factor'){
        levelNames = gsub('mx2\\[, cdataName\\]','',rownames(sum1$coefficients))[-1]
        rownames(sum1$coefficients) = c('Intercept',paste(cdataName,'.',levelNames,sep=''))
        f <- sum1$fstatistic 
        mx_pval[i,j] = pf(f[1], f[2], f[3], lower=FALSE)
        mx_beta[i,j] = NA
        
      }
      
      #mx_beta[i,j] = sum1$coefficients[2,1]
      #rownames(sum1$coefficients) = c('Intercept',cdataName)
      print(kable(sum1$coefficients,caption=paste(edataName,' vs ',cdataName,sep='')))
      
    }
    mx_qval[,j] = p.adjust(mx_pval[,j],method='BH')
    
  }
  fileout_pval = sprintf('%s/pval_%s.csv',dir_out,dataName)
  write.csv(mx_pval,file=fileout_pval, row.names=TRUE)
  
  fileout_beta = sprintf('%s/beta_%s.csv',dir_out,dataName)
  write.csv(mx_beta,file=fileout_beta, row.names=TRUE)
  
  fileout_qval = sprintf('%s/qval_%s.csv',dir_out,dataName)
  write.csv(mx_qval,file=fileout_qval, row.names=TRUE)
  
  
  dev.off()
  
  
}


### currently for some reason, it is not working and keep asking about j not exist
### weird !!!
pairwise_association_test_lmm_surv = function(mx2,cdataNames,edataNames, randomName, cvrtNames=''
								, dataName='out',dir_out = 'output', methods = rep('lmm',length(edataNames))) {

  n_cdata = length(cdataNames)
  n_edata = length(edataNames)
  mx_pval = matrix(NA, n_cdata, n_edata)
  rownames(mx_pval) = cdataNames
  colnames(mx_pval) = edataNames
  
  mx_beta = mx_pval
  mx_qval = mx_pval
  
  #for(j in 1:length(edataNames)) {

  for(j in 1:length(edataNames)) {
    for(i in 1:length(cdataNames)) {
      # i = 1; j = 1;
      cvrt_all = c(cdataNames[i], cvrtNames)
      cvrt_model = paste(cvrt_all, collapse=' + ')
      model_design = sprintf('%s ~ %s + (1|%s)', edataNames[j], cvrt_model, randomName)
        
      #lmm_fit <- lmer(mx2[,edataNames[j]] ~ mx2[,cdataNames[i]]+episode+(1|SUBID), data = mx2)
      if(methods[j] == 'lmm'){
      	lmm_fit <- lmer(model_design, data = mx2)
      	coefs <- data.frame(coef(summary(lmm_fit)))

	      cvrtNames_dummy = c()

	      for(cvrtName in cvrtNames){
	      	if(class(mx2[,cvrtName]) == 'factor') {
	      		cvrtNames_dummy = c(cvrtNames_dummy, paste(cvrtName,'.',levels(mx2[,cvrtName])[-1],sep=''))
	      	} else{
	      		cvrtNames_dummy = c(cvrtNames_dummy, cvrtName)
	      	}

	      }
	      rowNames_coeff = c('Intercept',cdataNames[i],cvrtNames_dummy)
	      rownames(coefs) = rowNames_coeff
	    
	      # use normal distribution to approximate p-value
	      coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
	        
	      mx_pval[i,j] = coefs[2,4]
	      mx_beta[i,j] = coefs[2,1]

      } else if(methods[j] == 'survival'){
        status = as.numeric(mx2[,edataNames[j]] != 12)

        lm1 = coxme(Surv(mx2[,edataNames[j]], status) ~ mx2[,cdataNames[i]] + episode + (1|SUBID), data=mx2)
        #sum1 = summary(lm1)
        mx_pval[i,j] = anova(lm1)[2,4]
        mx_beta[i,j] = fixef(lm1)[1]
        coefs = anova(lm1)
        rownames(coefs) = c("",cdataNames[i],'episode')
      }
      
      print(kable(coefs,caption = paste(edataNames[j], " VS ", cdataNames[i])))
      
    }
    mx_qval[,j] = p.adjust(mx_pval[,j],method='BH')
  }
  
  fileout_pval = sprintf('%s/pval_%s.csv',dir_out,dataName)
  write.csv(mx_pval,file=fileout_pval, row.names=TRUE)
  
  fileout_beta = sprintf('%s/beta_%s.csv',dir_out,dataName)
  write.csv(mx_beta,file=fileout_beta, row.names=TRUE)
  
  fileout_qval = sprintf('%s/qval_%s.csv',dir_out,dataName)
  write.csv(mx_qval,file=fileout_qval, row.names=TRUE)

}




### linear mixed effect model and add coxme for repeated measure
pairwise_association_test_lmm = function(mx2,cdataNames,edataNames, randomName, cvrtNames=''
                                         , dataName='out',dir_out = 'output', methods = rep('lmm',length(edataNames))) {
  
  n_cdata = length(cdataNames)
  n_edata = length(edataNames)
  mx_pval = matrix(NA, n_cdata, n_edata)
  rownames(mx_pval) = cdataNames
  colnames(mx_pval) = edataNames
  
  mx_beta = mx_pval
  mx_qval = mx_pval
  
  for(j in 1:length(edataNames)){
    for(i in 1:length(cdataNames)){
      # i = 1; j = 1;
      cvrt_all = c(cdataNames[i], cvrtNames)
      cvrt_model = paste(cvrt_all, collapse=' + ')
      model_design = sprintf('%s ~ %s + (1|%s)', edataNames[j], cvrt_model, randomName)
      
      #lmm_fit <- lmer(mx2[,edataNames[j]] ~ mx2[,cdataNames[i]]+episode+(1|SUBID), data = mx2)
      lmm_fit <- lmer(model_design, data = mx2)
      
      coefs <- data.frame(coef(summary(lmm_fit)))
      
      cvrtNames_dummy = c()
      
      for(cvrtName in cvrtNames){
        if(class(mx2[,cvrtName]) == 'factor') {
          cvrtNames_dummy = c(cvrtNames_dummy, paste(cvrtName,'.',levels(mx2[,cvrtName])[-1],sep=''))
        } else{
          cvrtNames_dummy = c(cvrtNames_dummy, cvrtName)
        }
        
      }
      rowNames_coeff = c('Intercept',cdataNames[i],cvrtNames_dummy)
      rownames(coefs) = rowNames_coeff
      
      # use normal distribution to approximate p-value
      coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
      
      mx_pval[i,j] = coefs[2,4]
      mx_beta[i,j] = coefs[2,1]
      

      

      ### add R2, confidence interval and df
      sum1 = summary(lmm_fit)
      Pval = coefs$p.z
      sum1$coefficients = cbind(sum1$coefficients, Pval)
      cfs = confint(lmm_fit, method="Wald")[-c(1,2),]
      cfs = round(cfs, 3)
      df1 = sum1$devcomp$dims[['N']] - sum1$devcomp$dims[['q']] - length(rowNames_coeff)
      var_X = c()


      idx_keep = which(!is.na(mx2[,edataNames[j]]) & !is.na(mx2[,cdataNames[i]]))
      #cvrtNames_X = c(cdataNames[i],cvrtNames_dummy)

      xx = mx2[idx_keep,]
      fm1 = as.formula(paste('~',cvrt_model))
      X1 = model.matrix(fm1, data = xx)
      var_X = apply(X1, 2, var, na.rm = TRUE)

      # if(is.numeric(mx2[,cvrtNames_X]) ) {
      # 	var_X = var(mx2[idx_keep,cvrtNames_X], na.rm = TRUE)
      # 	var_X = c(0,var_X)
      # } else{
      # 	xx = mx2[idx_keep,cdataName]
      # 	#subdf$letters <- factor(subdf$letters)
      # 	xx <- factor(xx)
      # 	X1 = model.matrix(~xx)
      # 	var_X = apply(X1, 2, var, na.rm = TRUE)
      # }

      R2 = coefs[,1]^2 * var_X/(sum(coefs[,1]^2 * var_X) + sum1$sigma^2)
      R2 = round(R2,3)
      mx_out = cbind(sum1$coefficients, cfs, R2)
      
      print(kable(mx_out,caption = paste(dataName,': ',edataNames[j], " VS ", cdataNames[i], ', df=',df1,sep='')))
      
    }
    mx_qval[,j] = p.adjust(mx_pval[,j],method='BH')
  }
  
  fileout_pval = sprintf('%s/pval_%s.csv',dir_out,dataName)
  write.csv(mx_pval,file=fileout_pval, row.names=TRUE)
  
  fileout_beta = sprintf('%s/beta_%s.csv',dir_out,dataName)
  write.csv(mx_beta,file=fileout_beta, row.names=TRUE)
  
  fileout_qval = sprintf('%s/qval_%s.csv',dir_out,dataName)
  write.csv(mx_qval,file=fileout_qval, row.names=TRUE)
  
}





cal_estimated_n_testing = function(mx1){
  ## the input matrix is p by p matrix of correlation
  mat = mx1
  
  
  n_marker = dim(mx1)[1]
  ##mat is pairwise haplotypic Pearson's correlation coefficients(linkage disequilibrium (LD) measure) matrix
  M<- n<- n_marker # the number of marker
  
  #mat<-matrix(NA,nrow=n,ncol=n)
  # mat[upper.tri(mat)]<- c(.957,.365,.367,.743,.343,.362,.027,.029,.133,.139,.12,.39,.408,.758,.721,.714,.049,.557,.586,.466,.425)
  # diag(mat)<-rep(1,n)
  # mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]
  # colnames(mat)<-rownames(mat)<-c("rs7085104","rs3740400","rs3740393","rs3740390","rs11191439","rs10748835","rs1046778")
  
  ##estimate Meff from Cheverud-Nyholt method
  #sq.mat<-sqrt(mat)
  #eig.val<-eigen(sq.mat)
  
  eig.val <- eigen(mat)
  
  var.eig<-var(eig.val$values) ##variance of eigenvalues of the pair correlation matrix
  Meff1<-1+(M-1)*(1-(var.eig)/M)  ##effective number of test using formula formula (1) in the paper
  
  ##estimate Meff from alternative formula for Cheverud-Nyholt method
  Meff2<-1+(sum(1-mat))/M   ## effective number of test using formula (7) in the paper
  
  ##estimate Neff from Valentina
  alpha<-0.05
  alpha<-0.001
  
  r<-kappa<-rep(NA,n)
  for(j in 2:n){
    r[j]<-max(abs(mat[1:(j-1),j]))
    ##kappa[j], number between 0 and 1 measuring degree of statistical independece of the test
    ##          for the jth marker from the tests for the preceding markers
    kappa[j]<-sqrt(1-r[j]^(-1.31*log10(alpha)))
  }
  Keff<-1+sum(kappa,na.rm=T)  ##estimated effective number of tests
  PI<-1-(1-alpha)^Keff        ##the exact overall type I error probability for the given individual significance level,0.05
  Neff=log(1-PI)/log(1-alpha) ##effective number of test usign formulae (2) through (6)
  
  #print(Neff);print(Meff1);print(Meff2)
  mx_out = matrix(c(Neff, Meff1, Meff2),nrow=1)
  rownames(mx_out) = c('Estimated Number of Testing')
  colnames(mx_out) = c('Neff','Meff1','Meff2')
  return(mx_out)
  
}



