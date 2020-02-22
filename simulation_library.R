library(pscl)
library(lmtest)

# snv simulation
SNV.simulation <- function(sample.size, zero.rate, gene.mu, is.snv=TRUE){
  if(is.snv){
    # Gene
    zero.size <- ceiling(sample.size * zero.rate)
    nonzero.size <- sample.size - zero.size
    gene.zero <- rep(0, zero.size)
    gene.nonzero <- rpois(nonzero.size, gene.mu)
    
    #SNV
    snv.zero <- rep(0,zero.size)
    error.rate <- 1/gene.nonzero
    error <- (runif(nonzero.size) < error.rate)
    snv.nonzero <- as.integer(!error)
    snv <- c(snv.zero, snv.nonzero)
  }
  else{
    snv <- rep(0,sample.size)
  }   
  return(sample(snv))
}

# gene simulation
gene.simulation <- function(sample.size,zero.rate, gene.mu){
  zero.size <- ceiling(sample.size * zero.rate)
  nonzero.size <- sample.size - zero.size
  gene.zero <- rep(0, zero.size)
  gene.nonzero <- rpois(nonzero.size, gene.mu)
  gene <- c(gene.zero,gene.nonzero)
  return(sample(gene))
}

# snv-gene-relationship simulation
related.pair.simulation <- function(sample.issnv.size, sample.notsnv.size, zero.snv.rate, zero.gene.rate, 
                                    snv.gene.mu,gene.with.snv.mu,gene.without.snv.mu){
  issnv <- SNV.simulation(sample.issnv.size,zero.snv.rate,snv.gene.mu,is.snv=TRUE)
  notsnv <- SNV.simulation(sample.notsnv.size,zero.snv.rate,snv.gene.mu,is.snv=FALSE)
  snv <- c(issnv,notsnv)
  
  gene.with.snv <- gene.simulation(sample.issnv.size,zero.gene.rate, gene.with.snv.mu)
  gene.without.snv <- gene.simulation(sample.notsnv.size, zero.gene.rate, gene.without.snv.mu)
  gene <- c(gene.with.snv,gene.without.snv)
  label <- c(rep(1,sample.issnv.size), rep(0,sample.notsnv.size))
  
  simulation.data <- data.frame(gene,snv, label)
  return(simulation.data)
}

# pvalue calculation with zinb model
calculate.pvalue <- function(df){
  # here df is a dataframe with 2 cols: gene, snv
  # if all genes are non-zero,back to poisson regression
  ##print(sum(df$gene==0))
  if(sum(df$gene==0)>0){
    zinb.model <- suppressWarnings(try(zeroinfl(formula = gene ~ snv, data = df, dist = "negbin"),silent = TRUE))
    m0 <- suppressWarnings(try(zeroinfl(formula = gene ~ 1|snv , data = df, dist = "negbin"),silent = TRUE))
    if (class(zinb.model)=='try-error' | class(m0)=='try-error'){
      pvalue <- NA
    }else{
      pvalue <- waldtest(m0,zinb.model)[2,'Pr(>Chisq)']
    }
    return(pvalue)
  }else{
    ##print('yes')
    poisson.model <- try(glm(formula = gene ~ snv, family="poisson", data=df),silent = TRUE)
    m0 <- try(glm(formula = gene ~ 1, family="poisson", data=df),silent = TRUE)
    message("Since all genes are non-zero, back to Poisson regression.")
    if ('try-error' %in% class(poisson.model)| 'try-error' %in% class(m0)){
      pvalue <- NA
    }else{
      pvalue <- waldtest(m0,poisson.model)[2,'Pr(>F)']
    }
    return(pvalue)
  }
}

# pvalue calculation with log-linear model
calculate.lm.pvalue <- function(df){
    lm.model <- suppressWarnings(try(lm(formula = log(gene+0.1) ~ snv, data=df),silent = TRUE))
    if ('try-error' %in% class(lm.model)){
        pvalue.lm <- NA
    }else{
        pvalue.lm <- anova(lm.model)['snv','Pr(>F)']
    }
    return(pvalue.lm)
}

# for compare different test methods for zinb model
calculate.pvalue.different.test <- function(df){
  # here df is a dataframe with 2 cols: gene, snv
  zinb.model <- try(zeroinfl(formula = gene ~ snv, data = df, dist = "negbin"),silent = TRUE)
  m0 <- try(zeroinfl(formula = gene ~ 1|snv , data = df, dist = "negbin"),silent = TRUE)
  if (class(zinb.model)=='try-error' | class(m0)=='try-error'){
    pvalue <- NA
    pvalue.chi.test <- NA
  }else{
    pvalue <- waldtest(m0,zinb.model)[2,'Pr(>Chisq)']
    pvalue.chi.test <- pchisq(2 * (logLik(zinb.model) - logLik(m0)), df = 1, lower.tail=FALSE)
  }
  result <- c(pvalue, pvalue.chi.test)
  names(result) <- c('pvalue', 'pvalue.chi.test')
  return(result)
}


calculate.poisson.pvalue <- function(df){
  # here df is a dataframe with 2 cols: gene, snv
  poisson.model <- try(glm(formula = gene ~ snv, family="poisson", data=df),silent = TRUE)
  m0 <- try(glm(formula = gene ~ 1, family="poisson", data=df),silent = TRUE)
  if ('try-error' %in% class(poisson.model)| 'try-error' %in% class(m0)){
    pvalue <- NA
  }else{
    pvalue <- waldtest(m0,poisson.model)[2,'Pr(>F)']
  }
  return(pvalue)
}

# for comparing different test methods for poisson regression
calculate.poisson.pvalue.different.test <- function(df){
  # here df is a dataframe with 2 cols: gene, snv
  poisson.model <- try(glm(formula = gene ~ snv, family="poisson", data=df),silent = TRUE)
  m0 <- try(glm(formula = gene ~ 1, family="poisson", data=df),silent = TRUE)
  if ('try-error' %in% class(poisson.model)| 'try-error' %in% class(m0)){
    pvalue <- NA
    pvalue.chi.test <- NA
  }else{
    pvalue <- waldtest(m0,poisson.model)[2,'Pr(>F)']
    #pvalue.chi.test <- pchisq(poisson.model$deviance, df=poisson.model$df.residual, lower.tail=FALSE)
    pvalue.chi.test <- pchisq(2 * (logLik(poisson.model) - logLik(m0)), df = 1, lower.tail=FALSE)
  }
  result <- c(pvalue, pvalue.chi.test)
  names(result) <- c('pvalue', 'pvalue.chi.test')
  return(result)
}

# for comparing log-linear and poisson regression
calculate.poisson.loglinear.zinb.pvalue <- function(df){
  # here df is a dataframe with 2 cols: gene, snv
  poisson.model <- try(glm(formula = gene ~ snv, family="poisson", data=df),silent = TRUE)
  m0 <- try(glm(formula = gene ~ 1, family="poisson", data=df),silent = TRUE)
  
  lm.model <- try(lm(formula = log(gene+0.01) ~ snv, data=df),silent = TRUE)
  
  zinb.model <- try(zeroinfl(formula = gene ~ snv, data = df, dist = "negbin"),silent = TRUE)
  m0.zinb <- try(zeroinfl(formula = gene ~ 1|snv , data = df, dist = "negbin"),silent = TRUE)
  if (class(zinb.model)=='try-error' | class(m0.zinb)=='try-error'){
    pvalue.zinb <- NA
  }else{
    pvalue.zinb <- waldtest(m0.zinb,zinb.model)[2,'Pr(>Chisq)']
  }
  if ('try-error' %in% class(poisson.model)| 'try-error' %in% class(m0)){
    pvalue.poi <- NA
  }else{
    pvalue.poi <- waldtest(m0,poisson.model)[2,'Pr(>F)']
  }
  
  if ('try-error' %in% class(lm.model)){
    pvalue.lm <- NA
  }else{
    pvalue.lm <- anova(lm.model)['snv','Pr(>F)']
  }
  result <- c(pvalue.poi, pvalue.lm,pvalue.zinb)
  names(result) <- c('pvalue.poi', 'pvalue.lm','pvalue.zinb')
  return(result)
}
if(FALSE){
# for comparing log-linear and zinb regression
calculate.loglinear.zinb.pvalue <- function(df){
  # here df is a dataframe with 2 cols: gene, snv
  lm.model <- try(lm(formula = log(gene+0.1) ~ snv, data=df),silent = TRUE)
  zinb.model <- try(zeroinfl(formula = gene ~ snv, data = df, dist = "negbin"),silent = TRUE)
  m0.zinb <- try(zeroinfl(formula = gene ~ 1|snv , data = df, dist = "negbin"),silent = TRUE)
  if (class(zinb.model)=='try-error' | class(m0.zinb)=='try-error'){
    pvalue.zinb <- NA
  }else{
    pvalue.zinb <- waldtest(m0.zinb,zinb.model)[2,'Pr(>Chisq)']
  }
  
  if ('try-error' %in% class(lm.model)){
    pvalue.lm <- NA
  }else{
    pvalue.lm <- anova(lm.model)['snv','Pr(>F)']
  }
  result <- c(pvalue.lm,pvalue.zinb)
  names(result) <- c('pvalue.lm','pvalue.zinb')
  return(result)
}
}

# for comparing log-linear and zinb regression
calculate.loglinear.zinb.pvalue <- function(df){
    # here df is a dataframe with 2 cols: gene, snv
    pvalue.lm <- calculate.lm.pvalue(df)
    pvalue.zinb <- calculate.pvalue(df)
    result <- c(pvalue.lm,pvalue.zinb)
    names(result) <- c('pvalue.lm','pvalue.zinb')
    return(result)
}
