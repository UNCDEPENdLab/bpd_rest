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

permutation.test <- function(data,components,num_permutations = 1000) {
	# assumes separator variable is pattype
	colmeans_control = colMeans(data[data$pattype==1,components])
	colmeans_pat = colMeans(data[data$pattype==0,components])
	diff.observed = colmeans_pat-colmeans_control

	len_p = length(data[data$pattype==0,components[1]])
	len_c = length(data[data$pattype==1,components[1]])

	num_permutations = 1000
	diff.random = matrix(,num_permutations,length(components))
	for (i in 1 : num_permutations){
		for (j in 1 : length(components))
		{
			shuffled = sample(data[,components[j]],length(data[,components[j]]))
			#a.random = shuffled[1:len_c]
			#b.random = shuffled[(len_c+1):length(data[,components[j]])]
			a.random = head(shuffled,len_c)
			b.random = tail(shuffled,len_p)
			diff.random[i,j] = mean(a.random) - mean(b.random)
		}
	}
	pvalue = colSums(abs(diff.random) >= abs(diff.observed)) / num_permutations
	return(pvalue)
}

all_rois = unique(dataset$roi)
values = matrix(,length(all_rois),length(components))
colnames(values) = components
for (i in 1:length(all_rois)){
	data = dataset[dataset$roi == all_rois[i],]
	values[i,]=(permutation.test(data,components))
}
values = cbind(all_rois,values)
values[values[,'PC1']<0.05 | values[,'PC2'] < 0.05 | values[,'PC3'] < 0.05,]

graph.hist <- function(one,two,breaks=10) {
	p1 = hist(one,breaks=breaks,plot=FALSE)
	p2 = hist(two,breaks=breaks,plot=FALSE)
	plot(p1,col=rgb(0,0,1,1/4))
	plot(p2,col=rgb(1,0,0,1/4),add=T)
}

graph.hist(dataset[dataset$pattype==1,'PC1'],dataset[dataset$pattype==0,'PC2'],breaks=100)
graph.hist(dataset[dataset$pattype==1 & dataset$roi == 1,'PC1'],dataset[dataset$pattype==0 & dataset$roi == 1,'PC2'],breaks=10)

