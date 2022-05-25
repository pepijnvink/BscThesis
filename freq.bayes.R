# R function to simulate data for frequentist and Bayesian t-tests
freq.bayes <- function(sample.size = 10, iterations = 10, mean.low = 100, sd = 10, difference = 0){
  col_names <- c("difference", "n", "mean1", "mean2", "sd", "y1", "y2", "diff_obs", "p", "ci.low","ci.high", "HDI.low", "HDI.high", "BayesFactor", "BFa", "p.decision", "BF1.decision", "BF3.decision", "BF10.decision", "BF30.decision", "BF.3.cor", "BF.10.cor", "BF.30.cor", "BF.3.ninc", "BF.10.ninc", "BF.30.ninc", "BF.evidence", "post.null", "post.alt", "HDI.decision", "CI.est.correct", "HDI.est.correct", "ci.width", "hdi.width")
  
  save <- matrix(data = NA, nrow = iterations, ncol = 34) # create matrix for data to be added to
  for(i in 1:iterations){
    diff <- difference
    n <- sample.size
    
    mean2 <- mean.low # end up with a positive (or null) mean difference
    mean1 <- mean2 + diff
    
    y1 <- rnorm(n = n, mean = mean1, sd = sd) # simulate first sample
    y2 <- rnorm(n = n, mean = mean2, sd = sd) # simulate second sample
    
    result <- t.test(y1, y2, var.equal = TRUE) # frequentist t-test
    p <- result$p.value # extract p-value
    ci <- result$conf.int # extract CI
    
    mean_y1 <- mean(y1) # sample mean of group 1
    mean_y2 <- mean(y2) # sample mean of group 2
    
    ci.low <- ci[1] # split up CI into its lower and upper bounds
    ci.high <- ci[2]
    
    mean_diff <- mean(y1) - mean(y2) # sample difference
    
    bfa <- ttestBF(y1, y2) # Bayesian t-test
    
    bf0 <- 1/bfa # compute BF0a
    
    post <- bayestestR::hdi(bfa) # extract HDI of the posterior distribution
    
    HDI.low <- post$CI_low # split up HDI into lower and upper bound
    HDI.high <- post$CI_high
    
    BFa <- extractBF(bfa)$bf # extract Bayes factor a0
    
    BF0 <- extractBF(bf0)$bf # extract Bayes factor 0a
    
    p_decision <- case_when(p < 0.05 ~ "alt", # label p value decision
                           TRUE ~ "null")
    
    BF.3 <- case_when(BF0 >= 3 ~ "null", # label Bayes factor decisions
                      BF0 < 3 & BFa < 3 ~ "U",
                      BFa >= 3 ~ "alt")
    
    BF.10 <- case_when(BF0 >= 10 ~ "null",
                       BF0 < 10 & BFa < 10 ~ "U",
                       BFa >= 10 ~ "alt")
    
    BF.30 <- case_when(BF0 >= 30 ~ "null",
                       BF0 < 30 & BFa < 30 ~ "U",
                       BFa >= 30 ~ "alt")
    
    BF.1 <- case_when(BF0 > 1 ~ "null",
                      TRUE ~ "alt")
    
    BF.3.dich.cor <- case_when(BFa > 3 & diff != 0 ~ "correct", # label if BF decision is correct
                               BF0 > 3 & diff == 0 ~ "correct",
                               TRUE ~ "incorrect")
    BF.10.dich.cor <- case_when(BFa > 10 & diff != 0 ~ "correct",
                                BF0 > 10 & diff == 0 ~ "correct",
                                TRUE ~ "incorrect")
    
    BF.30.dich.cor <- case_when(BFa > 30 & diff != 0 ~ "correct",
                                BF0 > 30 & diff == 0 ~ "correct",
                                TRUE ~ "incorrect")
    
    BF.3.dich.nincor <- case_when(BF0 < 3 & diff != 0 ~ "correct",
                                  BFa < 3 & diff == 0 ~ "correct",
                                  TRUE ~ "incorrect")
    
    BF.10.dich.nincor <- case_when(BF0 < 10 & diff != 0 ~ "correct",
                                   BFa < 10 & diff == 0 ~ "correct",
                                   TRUE ~ "incorrect")
    
    BF.30.dich.nincor <- case_when(BF0 < 30 & diff != 0 ~ "correct",
                                   BFa < 30 & diff == 0 ~ "correct",
                                   TRUE ~ "incorrect")
    
    BF.evidence <- case_when(BF0 < 1/100 ~ -5, # quantify degree of evidence using BF based on Jeffrey's (1961)
                             BF0 >= 1/100 & BF0 < 1/30 ~ -4,
                             BF0 >= 1/30 & BF0 < 1/10 ~ -3,
                             BF0 >= 1/10 & BF0 < 1/3 ~ -2,
                             BF0 >= 1/3 & BF0 < 1 ~ -1,
                             BF0 == 1 ~ 0,
                             BF0 > 1 & BF0 <= 3 ~ 1,
                             BF0 > 3 & BF0 <= 10 ~ 2,
                             BF0 > 10 & BF0 <= 30 ~ 3,
                             BF0 > 30 & BF0 <= 100 ~ 4,
                             BF0 > 100 ~ 5)
    post.alt <- 1/(BF0+1) # compute PMP's if prior odds are 1
    
    post.null <- 1/(BFa+1)
    
    HDI_decision <- case_when(HDI.low <= 0 & HDI.high >= 0 ~ "null", # label decision using HDI
                              TRUE ~ "alt")
    
    HDI_correct_est <- case_when(HDI.low <= diff & HDI.high >= diff ~ "correct", # label if the HDI contains the population difference
                                 TRUE ~ "incorrect")
    
    CI_correct_est <- case_when(ci.low <= diff & ci.high >= diff ~ "correct", # label if the CI contains the population difference
                                TRUE ~ "incorrect")
    
    HDI.width <- HDI.high - HDI.low # compute width of the HDI
    
    CI.width <- ci.high - ci.low # compute width of the CI
    
    vec <- c(diff, n, mean1, mean2, sd, mean_y1, mean_y2, mean_diff, p, ci.low, ci.high, HDI.low, HDI.high, BF0, BFa, p_decision, BF.1, BF.3, BF.10, BF.30, BF.3.dich.cor, BF.10.dich.cor, BF.30.dich.cor, BF.3.dich.nincor, BF.10.dich.nincor, BF.30.dich.nincor, BF.evidence,  post.null, post.alt, HDI_decision, CI_correct_est, HDI_correct_est, CI.width, HDI.width)
    
    save[i,] <- vec # add variables to the dataframe
  }
  colnames(save) <- col_names # name the columns
  save # makes sure it is printed
}