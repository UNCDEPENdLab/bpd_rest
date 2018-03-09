meanreplace <- function(df, ID, varprefix="SPECC", subscaleitems) {
  subjrow <- df$SPECC_ID == ID
  # browser()
  dfitems <- subset(df, select=paste0(varprefix, subscaleitems))
  dfitems[dfitems == 99] <- NA
  subjitems <- dfitems[subjrow,]
  whichmiss <- names(subjitems)[which(is.na(subjitems))]
  obsmean <- mean(unlist(subjitems), na.rm=TRUE)
  df[subjrow, whichmiss] <- obsmean
  df
}
