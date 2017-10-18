gateway_coeff_NH <- function (g, memb, centr = c("btwn.cent", "degree", "strength")) 
{
  stopifnot(is_igraph(g))
  if ("degree" %in% vertex_attr_names(g)) {
    Ki <- V(g)$degree
  }
  else {
    Ki <- degree(g)
  }
  centr <- match.arg(centr)
  if (centr == "btwn.cent") {
    if ("btwn.cent" %in% vertex_attr_names(g)) {
      cent <- V(g)$btwn.cent
    }
    else {
      cent <- betweenness(g, weights=1/E(g)$weight)
    }
  }
  else if (centr == "degree") {
    cent <- Ki
  }
  else if (centr == "strength") {
    if ("strength" %in% vertex_attr_names(g)) {
      cent <- V(g)$strength
    }
    else {
      cent <- strength(g)
    }
  }
  N <- max(memb)
  Cn <- max(vapply(seq_len(N), function(x) sum(cent[which(memb == 
                                                            x)]), numeric(1)))
  A <- as_adj(g, sparse = FALSE, names = FALSE)
  Kis <- vapply(seq_len(N), function(x) colSums(A * (memb == 
                                                       x)), numeric(length(Ki)))
  M <- max(which(table(memb) > 1))
  Kjs <- matrix(0, nrow = N, ncol = N)
  Kjs[1:M, 1:M] <- vapply(seq_len(M), function(x) colSums(Kis[which(memb == 
                                                                      x), 1:M]), numeric(M))
  barKis <- Cis <- matrix(0, nrow = length(Ki), ncol = N)
  for (i in which(Ki > 0)) {
    barKis[i, ] <- Kis[i, ]/Kjs[, memb[i]]
    Vi <- which(A[, i] == 1)
    Cis[i, ] <- vapply(seq_len(N), function(x) sum(cent[Vi[memb[Vi] == 
                                                             x]]), numeric(1))
  }
  barCis <- Cis/Cn
  gis <- 1 - barKis * barCis
  G <- 1 - ((1/Ki^2) * rowSums(Kis^2 * gis^2, na.rm = TRUE))
  return(G)
}