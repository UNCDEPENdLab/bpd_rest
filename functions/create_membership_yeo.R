getwd()
setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/data/")
schaef.master <- read.csv("schaefer422_masterlookup.csv")
head(schaef.master)

membership.yeo <- array(NA, nnodes)
  
for(i in 1:nnodes){
if(schaef.master$X7_networks[i] == "Vis") {membership.yeo[i] <- 1}
  else{
    if(schaef.master$X7_networks[i] == "SomMot") {membership.yeo[i] <- 2}
    else{
      if(schaef.master$X7_networks[i] == "DorsAttn") {membership.yeo[i] <- 3}
      else{
        if(schaef.master$X7_networks[i] == "SalVentAttn") {membership.yeo[i] <- 4}
        else{
          if(schaef.master$X7_networks[i] == "Limbic") {membership.yeo[i] <-5}
          else{
            if(schaef.master$X7_networks[i] == "FPN") {membership.yeo[i] <- 6}
            else{
              membership.yeo[i] <- 7
            }
          }
        }
      }
    }
  }
}

write.csv(membership.yeo, file = "membership.yeo.csv", row.names = FALSE, col.names = FALSE)
