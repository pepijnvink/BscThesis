interval_inclusion <- function(D){
  it <- nrow(D)
  save <- data.frame(ncol = 2, nrow = it, stringsAsFactors = TRUE)
  for(i in 1:it){
    if(i == it){
      diff <- D[1, "y1"] - D[1, "y2"]
    } else {
      diff <- D[i+1, "y1"] - D[i+1, "y2"]
    }
    ci_low <- D[i, "ci.low"]
    ci_high <- D[i, "ci.high"]
    hdi_low <- D[i, "HDI.low"]
    hdi_high <- D[i, "HDI.high"]
    
    save[i, 1] <- case_when(ci_low <= diff & ci_high >= diff ~ "included",
                            TRUE ~ "not")
    save[i, 2] <- case_when(hdi_low <= diff & hdi_high >= diff ~ "included",
                            TRUE ~ "not")
  }
  colnames(save) <- c("CI.Included", "HDI.Included")
  save$CI.Included <- as_factor(save$CI.Included)
  save$HDI.Included <- as_factor(save$HDI.Included)
  save
}