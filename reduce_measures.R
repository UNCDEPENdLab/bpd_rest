library(ggplot2)
library(reshape2)
library(lattice)
library(reshape2)
library(psych)

nodemeasures <- list.files("/Volumes/Serena/SPECC/Neil/bpd_rest/neil/stats_output_v2/", pattern='^[01]',full.names=TRUE)
#nodemeasures = nodemeasures[1:2]

allsubjs <- c()
for (f in nodemeasures) {
  d <- read.csv(f, header=TRUE)
  d_wide <- dcast(d,roi ~ edge_def+parameter+statistic, value.var='value')
  d_wide$id <- basename(f)
  allsubjs <- rbind(allsubjs, d_wide)
}

bysubj <- split(allsubjs, f=allsubjs$id)
byroi <- split(allsubjs, f=allsubjs$roi)

graph_measures <- grep('HARD_0.*',names(allsubjs),value=TRUE)
graph_measures <- c(graph_measures,grep('SOFT.*(pagerank|eigenvector|local).*',names(allsubjs),value=TRUE))
#graph_measures <- names(allsubjs[!names(allsubjs) %in% c("id", "roi")])

pr <- prcomp(allsubjs[,graph_measures], scale.=TRUE)
sum(pr$sdev)
cumvariance <- (cumsum((pr$sdev)^2) / sum(pr$sdev^2))
print(cumvariance)
pr$rotation #varimax or promax rotation if you need to interpret this

####

#promax for oblique, varimax for orthogonal rotation
f1 <- principal(allsubjs[,graph_measures], nfactors=3, rotate="varimax")
f1

#with thresholded loadings
print(f1$loadings, cutoff=0.25)

pattype = ifelse(allsubjs$id > "0" & allsubjs$id < "1",0,1) # 0 = patient, 1 = control (consistent with file name, not with common understanding

#add the PC scores onto the original dataset

dataset <- cbind(allsubjs, f1$scores,pattype)

components = colnames(f1$scores)
colmeans_control = colMeans(dataset[dataset$pattype==1,components])
colmeans_pat = colMeans(dataset[dataset$pattype==0,components])
diff.observed = colmeans_pat-colmeans_control

len_p = length(dataset[dataset$pattype==0,components[1]])
len_c = length(dataset[dataset$pattype==1,components[1]])

num_permutations = 1000
diff.random = matrix(,num_permutations,length(components))
for (i in 1 : num_permutations){
	for (j in 1 : length(components))
	{
		shuffled = sample(dataset[,components[j]],length(dataset[,components[j]]))
		#a.random = shuffled[1:len_c]
		#b.random = shuffled[(len_c)+1:length(dataset[,components[j]])]
		a.random = head(shuffled,len_c)
		b.random = tail(shuffled,len_p)
		diff.random[i,j] = mean(a.random) - mean(b.random)
	}
}

pvalue = colSums(abs(diff.random) >= abs(diff.observed)) / num_permutations
print (pvalue)
