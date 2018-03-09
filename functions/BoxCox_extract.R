BoxCox_extract <- function(data, lambdas, plot = TRUE){
  Box <- boxcox(data ~ 1, lambda = lambdas)
  Cox = data.frame(Box$x, Box$y)            # Create a data frame with the results
  
  Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] # Order the new data frame by decreasing y
  Cox2[1,]    #display lambda with highest log-likelihood
  lambda = Cox2[1, "Box.x"] #extract lambda
  
  transformed <- (data ^lambda -1)/lambda
  
  if(plot) {
    before <- qplot(data)
    print(before)
    
    x <- qplot(transformed)
    print(x)
  }
  return(transformed)
}