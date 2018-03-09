setwd("~/Desktop/")
ts_file <- read.table("001ra_07dec2013_schaefer422_ts_file.txt")

library(corpcor)
x <- cor.shrink(ts_file)
hist(x)

y <- cor(ts_file)
hist(y)
id <- diag(422)

lambda <- 0.0228

r_star <- lambda*id + (1-lambda)*y

hist(r_star)

#r_star == x

#all.equal(r_star, x)

diffmat <- r_star - x
hist(diffmat[lower.tri(diffmat)])

pearson_m_shrink <- r_star - y
hist(pearson_m_shrink[lower.tri(pearson_m_shrink)])






r_star <- atanh(r_star)
y <- atanh(y)

sep=""
setwd("~/Desktop/")

tempcor_file <- read.table("001ra_07dec2013_schaefer422_cor.shrink.txt",sep=sep, header=FALSE)
pearson_file <- read.table("001ra_07dec2013_schaefer422_pearson.txt", sep = sep, header = FALSE)

# diffmat_1 <- r_star - tempcor_file
# hist(diffmat_1[lower.tri(diffmat_1)])
# 
# 
# diffmat_2 <- pearson_file - y
# hist(diffmat_2[lower.tri(diffmat_2)])
diag(pearson_file) <-1
hist(as.matrix(pearson_file))

diag(tempcor_file) <- 1

range(tempcor_file)
hist(as.matrix(tempcor_file))

subj1 <- allmats[1,,]

diff <- pearson_file - tempcor_file
hist(diff[lower.tri(diff)])

pdf("corcomparison.pdf", width=10, height=8)
for (i in 1:dim(allmats)) {
  pearson <- allmats_[i,,]
  shrink <- allmats[i,,]
  diffmat <- pearson-shrink
  allcorvec <- data.frame(p=pearson[lower.tri(pearson)], s=shrink[lower.tri(shrink)], d=diffmat[lower.tri(diffmat)]) %>%
    gather()
  ggplot(allcorvec, aes(x=value)) + geom_histogram() + facet_wrap(~key, scales="free")
}
dev.off()

divmat <- shrink/pearson

