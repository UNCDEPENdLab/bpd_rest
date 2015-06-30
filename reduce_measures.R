library(ggplot2)
library(reshape2)
library(lattice)
library(reshape2)
library(psych)

nodemeasures <- list.files("/Volumes/Serena/SPECC/Neil/bpd_rest/neil/stats_output_v3/", pattern='^[01]',full.names=TRUE)
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

#pr <- prcomp(allsubjs[,graph_measures], scale.=TRUE)
#sum(pr$sdev)
#cumvariance <- (cumsum((pr$sdev)^2) / sum(pr$sdev^2))
#print(cumvariance)
#pr$rotation #varimax or promax rotation if you need to interpret this

####
# Not needed, as it turns out, as psych->principal appears to scale by default
allsubjs.scaled = cbind(scale(allsubjs[,graph_measures]),allsubjs[,c('id','roi')])

#promax for oblique, varimax for orthogonal rotation
f1 <- principal(allsubjs.scaled[,graph_measures], nfactors=3, rotate="varimax",scores=TRUE)
#f1 <- principal(allsubjs.scaled[,graph_measures], nfactors=3, rotate="promax")
f1

#with thresholded loadings
print(f1$loadings, cutoff=0.25)

pattype = ifelse(allsubjs.scaled$id > "0" & allsubjs.scaled$id < "1",0,1) # 0 = patient, 1 = control (consistent with file name, not with common understanding

#add the PC scores onto the original dataset

dataset <- cbind(allsubjs.scaled, f1$scores,pattype)

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
values[values[,components[1]]<0.05 | values[,components[2]] < 0.05 | values[,components[3]] < 0.05,]

graph.hist <- function(one,two,breaks=10) {
	p1 = hist(one,breaks=breaks,plot=FALSE)
	p2 = hist(two,breaks=breaks,plot=FALSE)
	xlim <- range(p1$breaks,p2$breaks)
	ylim <- range(0,p1$density,p2$density)
	plot(p1,col=rgb(0,0,1,1/4),xlim = xlim, ylim = ylim,freq = FALSE)
	plot(p2,col=rgb(1,0,0,1/4),add=TRUE,xlim=xlim, ylim=ylim,freq=FALSE)
}
graph.density <- function(one,two,labelOne = "One",labelTwo="Two",main = "Density") {
	p1 = density(one)
	p2 = density(two)
	xlim <- range(p1$x,p2$x)
	ylim <- range(0,p1$y,p2$y)
	oneCol <- rgb(0,0,1,0.2)
	twoCol <- rgb(1,0,0,0.2)
	plot(p1,xlim=xlim,ylim = ylim,panel.first = grid(),main=main)
	polygon(p1,density = -1, col = oneCol)
	polygon(p2,density = -1, col = twoCol)
	legend('topleft',c(labelOne,labelTwo),fill = c(oneCol,twoCol),bty = 'n', border = NA)
}



library(data.table)
outlierReplace = function(dataframe, cols, rows, newValue = NA) {
	if (any(rows)) {
		set(dataframe, rows, cols, newValue)
	}
}

outlierReplace(dataset,'PC3',which(dataset$PC3 > 6),6) # INVESTIGATE FURTHER
outlierReplace(dataset,'PC2',which(dataset$PC2 > 6),6) # Due to local_efficiency outliers in the power-law edge definition group

pat = dataset[dataset$pattype == 0,]
control = dataset[dataset$pattype == 1,]

results = matrix(,length(all_rois),length(components))
for (i in 1 : length(all_rois)){
	for (j in 1: length(components)){
		t1 = control[control$roi == all_rois[i],components[j]]
		t2 = pat[pat$roi == all_rois[i],components[j]]
		t = t.test(t1,t2)
		results[i,j] = t$p.value
	}
}
colnames(results) = components
results = cbind(all_rois,results)
alpha = 0.05 / length(results)
interesting_rows = results[results[,'PC1']<alpha | results[,'PC2'] < alpha | results[,'PC3'] < alpha,]

if ( FALSE){
graph.hist(dataset[dataset$pattype==1,'PC1'],dataset[dataset$pattype==0,'PC1'],breaks=100)
graph.hist(dataset[dataset$pattype==1 & dataset$roi == 1,'PC2'],dataset[dataset$pattype==0 & dataset$roi == 1,'PC2'],breaks=10)
graph.hist(dataset[dataset$pattype==1 & dataset$roi == 253,'PC1'],dataset[dataset$pattype==0 & dataset$roi == 253,'PC1'],breaks=20)
graph.hist(dataset[dataset$pattype==0 & dataset$roi == 263,'PC2'],dataset[dataset$pattype==1 & dataset$roi == 263,'PC2'],breaks=20)
graph.density(dataset[dataset$pattype==0 & dataset$roi == 263,'PC2'],dataset[dataset$pattype==1 & dataset$roi == 263,'PC2'])
t.test(pat[pat$roi == 263,'SOFT_10_local_efficiency'],control[control$roi==263,'SOFT_10_local_efficiency'])
var.test(pat[pat$roi == 263,'PC2'],control[control$roi==263,'PC2'])
graph.density(dataset$SOFT_14_local_efficiency,dataset[dataset$PC2>=6,'SOFT_14_local_efficiency']) # outliers on PC2
graph.density(dataset$HARD_0.95_betweenness_binary,dataset[dataset$PC3>=6,'HARD_0.95_betweenness_binary']) # outliers on PC3

for (i in 1:length(interesting_rows[,'all_rois'])){
	par(mfrow=c(3,1))
	roi = interesting_rows[i,'all_rois']
	for (j in 1:length(components)){
		graph.density(control[control$roi==roi,components[j]],pat[pat$roi==roi,components[j]],labelOne="control",labelTwo="patient",main=sprintf("%d: %s, p=%f",roi,components[j],interesting_rows[i,components[j]]))
	}
	readline()
}
dev.off()

# test actual graph measures for each of the interesting rows to see which ones are contributing to statistical significance
for (i in 1:length(interesting_rows[,'all_rois'])){
	roi = interesting_rows[i,'all_rois']
	for (j in 1: length(graph_measures)){
		t1 = control[control$roi == roi,graph_measures[j]]
		t2 = pat[pat$roi == roi,graph_measures[j]]
		t = t.test(t1,t2)
		if (t$p.value < alpha){
			print(sprintf("p-val: %f roi: %d, graph measure %s",t$p.value,roi,graph_measures[j]))
		}
	}
}

outliers = which(dataset$PC2>= 6 | dataset$PC3>=6)
subset_measures = grep('HARD_0.95|SOFT_12',graph_measures,value=TRUE)
for (i in 1:length(outliers)){
	par(mfrow=c(4,3))
	row = dataset[outliers[i],]
	roi = row$roi
	for (j in 1:length(subset_measures)){
		#graph.density(control[control$roi==roi,subset_measures[j]],pat[pat$roi==roi,subset_measures[j]],labelOne="control",labelTwo="patient",main=sprintf("%d: %s",roi,subset_measures[j]))
		graph.density(dataset[,subset_measures[j]],dataset[dataset$roi==roi,subset_measures[j]],labelOne="dataset",labelTwo="roi",main=sprintf("%d: %s",roi,subset_measures[j]))
	}
	for (j in 1:length(components)){
		#graph.density(control[control$roi==roi,components[j]],pat[pat$roi==roi,components[j]],labelOne="control",labelTwo="patient",main=sprintf("%d: %s, p=%f",roi,components[j],row[components[j]]))
		graph.density(dataset[,components[j]],dataset[dataset$roi==roi,components[j]],labelOne="dataset",labelTwo="roi",main=sprintf("%d: %s, p=%f",roi,components[j],row[components[j]]))
	}
	readline()
}
dev.off()

}
