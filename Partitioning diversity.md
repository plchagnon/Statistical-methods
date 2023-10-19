# Partitioning diversity following Laliberté et al. 2020 Ecol lett

As outlined by Laliberté et al. 2020 for the case of spectral diversity, it can be useful/insightful to partition components of diversity ($\gamma, \beta \text{ and } \alpha$) when rows in the dataset can be ascribed to discrete groups. This way, one can look at local row contribution to $\gamma$ diversity, feature contribution to $\gamma$ diversity, group contribution to $\beta$-diversity, feature contribution to $\beta$-diversity, and feature contribution to local $\alpha$-diversity in a given group. This may highlight, for example, some features that are very variable among groups (high FCD $_{\beta}$), but quite stable within groups (i.e., low average FCD $_\alpha$ for all groups). 

Below, by taking the ``mite`` dataset in the R library ``vegan``, we start by defining heuristically 4 groups of samples (rows) more similar to each other than to any other sample in the metacommunity, and partition diversity into $\alpha, \beta$ and $\gamma$ components, along with local and feature contributions, following [Laliberté et al. 2020](https://doi.org/10.1111/ele.13429).

```r

rm(list=ls())
### Import a dataset
library(vegan)
data(mite)
y=mite[-67,]

### Say we cluster this metacommunity into 4 groups
groups=kmeans(y,centers=4)$cluster

### Partitioning function:

part=function(y,groups){

	#gamma
	sdis=function(x){(x-mean(x))^2}
	sij=apply(y,2,sdis)
	SSgamma=sum(sij)
	SDgamma=SSgamma/(nrow(y)-1)
	SSgamma_j=colSums(sij)
	FCDgamma=SSgamma_j/SSgamma
	LCDgamma=rowSums(sij)/SSgamma

	#beta
	ypred=matrix(nr=max(unique(groups)),nc=ncol(y))
	for(i in 1:max(unique(groups))){
	ypred[i,]=colMeans(y[groups==i,])}
	gr_dis=function(x){(x-colMeans(y))^2}
	skj=t(apply(ypred,1,gr_dis))
	gr_size=numeric(length(unique(groups)))
	for(i in 1:length(gr_size)){gr_size[i]=nrow(y[groups==i,])}
	SSbeta_k=rowSums(apply(skj,2,function(x){x*gr_size}))
	SSbeta=sum(SSbeta_k)
	SDbeta=SSbeta/(nrow(y)-1)
	LCDbeta_k=SSbeta_k/SSbeta
	SSbeta_j=colSums(apply(skj,2,function(x){x*gr_size}))
	FCDbeta_j=SSbeta_j/SSbeta

	#alpha
	sub=SSalpha_jk=list(length(unique(groups)))
		for(i in 1:length(unique(groups))){
		sub[[i]]=y[groups==i,]
		SSalpha_jk[[i]]=colSums(t(apply(sub[[i]],1,function(x){(x-ypred[i,])^2})))
		}

	SSalpha_k=as.numeric(lapply(SSalpha_jk,sum))
	SDalpha_k=SSalpha_k/as.numeric(lapply(sub,nrow))
	SSalpha=sum(SSalpha_k)
	FCDalpha_jk=do.call(rbind,lapply(SSalpha_jk,function(x){x/sum(x)}))
	
	res=list(SSgamma,SSalpha,SSbeta,SDgamma,SDalpha_k,SDbeta,LCDgamma,LCDbeta_k,FCDgamma,FCDalpha_jk,FCDbeta_j)
	names(res)=c("##### SS_gamma #####","##### SS_alpha #####","##### SS_beta #####","##### SD_gamma #####","##### SD_alpha_per_group #####","##### SD_beta #####","##### LCD_gamma #####","##### LCD_beta_k(group contrib to regional beta-div) #####","##### FCD_gamma #####","##### FCD to alpha diversity within groups #####","##### FCD to beta-diversity in the region #####")
	return=res}


RES=part(y,groups)

```

This code is available [here](./Partitioning%20diversity.R) as an R script.