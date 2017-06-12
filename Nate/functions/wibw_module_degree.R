#adapted function from igraph
wibw_module_degree <- function (g, community_attr="community") {
  stopifnot(is_igraph(g))
  if (is.null(get.vertex.attribute(g, community_attr))) { 
    warning("Cannot find community vertex attribute: ", community_attr) 
    return(NULL)
  }
  memb <- get.vertex.attribute(g, community_attr)
  modules <- unique(memb)
  N <- length(modules)
  A <- as_adj(g, sparse = FALSE, names = FALSE)
  z_within <- Ki <- rep(0, nrow(A))
  z_between <- Bi <- rep(0, nrow(A))
  Ksi <- sigKsi <- rep(0, nrow(A))
  Bsi <- sigBsi <- rep(0, nrow(A))
  for (S in modules) {
    x_wi <- A[memb == S, ] %*% (memb == S) #from S to S
    x_bw <- A[memb == S, ] %*% (memb != S) #from S to other nodes
    Ki[memb == S] <- x_wi; Ksi[memb==S] <- mean(x_wi); sigKsi[memb==S] <- sd(x_wi)
    Bi[memb == S] <- x_bw; Bsi[memb==S] <- mean(x_bw); sigBsi[memb==S] <- sd(x_bw)
  }
  z_within <- (Ki - Ksi)/sigKsi
  z_within <- ifelse(!is.finite(z_within), 0, z_within)
  z_between <- (Bi - Bsi)/sigBsi
  z_between <- ifelse(!is.finite(z_between), 0, z_between)
  return(data.frame(module=memb, Ki=Ki, Ksi=Ksi, sigKsi=sigKsi, z_within=z_within, Bi=Bi, Bsi=Bsi, sigBsi=sigBsi, z_between=z_between))
}