########### HERE CODE USING BRAY-CURTIS INDEX OF DISSIMILARITY INSTEAD OF WHITTAKER AS IN POISOT 2012

rm(list=ls())
rich_p=20
rich_f=20
reg_p=paste0("plant",1:rich_p)
reg_f=paste0("fungus",1:rich_f)
n=4
min_sp_per_comp=7
nw=list(n)
for(i in 1:n){
	rows=sample(min_sp_per_comp:rich_f,1)
	cols=sample(min_sp_per_comp:rich_p,1)
	nw[[i]]=matrix(rbinom(rows*cols,10,.1),nc=cols)
	rownames(nw[[i]])=sample(reg_p,nrow(nw[[i]]),replace=F)
	colnames(nw[[i]])=sample(reg_f,ncol(nw[[i]]),replace=F)}


#### the "nw" list will be the input for the function below:


#### Function to change a network matrix to a list format
nw_to_tab=function(x){as.data.frame(cbind(rep(rownames(x),ncol(x)),rep(colnames(x),each=nrow(x)),c(x)))}

#### Add missing spp to pair of networks

complete=function(focal,other){

	### add absent spp to nw1
	m1=focal
	j_1=colnames(other)[!colnames(other)%in%colnames(focal)]
	m1=cbind(focal,matrix(rep(0,length(j_1)*nrow(focal)),nc=length(j_1)))
	colnames(m1)=c(colnames(focal),j_1)
	i_1=rownames(other)[!rownames(other)%in%rownames(focal)]
	m1=rbind(m1,matrix(rep(0,length(i_1)*ncol(m1)),nc=ncol(m1)))
	rownames(m1)=c(rownames(focal),i_1)

return=m1}



beta=function(nw1,nw2,binary){

### TOTAL DISSIMILARITY

	if(binary==TRUE){
		nw1[nw1>0]<-1
		nw2[nw2>0]<-1}

	c1=complete(nw1,nw2)
	c2=complete(nw2,nw1)

	t1=nw_to_tab(c1)
	t2=nw_to_tab(c2)
	
	t1=t1[order(t1[,1],t1[,2]),]
	t2=t2[order(t2[,1],t2[,2]),]
	
	x=apply(cbind(t1[,3],t2[,3]),2,as.numeric)
	B_tot=sum(abs(apply(x,1,diff)))/sum(apply(x,1,sum))

### REWIRING DISSIMILARITY

	s1=nw1[rownames(nw1)%in%rownames(nw2),colnames(nw1)%in%colnames(nw2)]
	s2=nw2[rownames(nw2)%in%rownames(nw1),colnames(nw2)%in%colnames(nw1)]
	s1=s1[rownames(s2),colnames(s2)]
	
	x=apply(cbind(c(s1),c(s2)),2,as.numeric)
	B_rw=sum(abs(apply(x,1,diff)))/sum(apply(x,1,sum))

	B_turn=B_tot-B_rw

	res=list(B_tot,B_rw,B_turn)
	names(res)=c("B_tot","B_rw","B_turn")
	return=res}




### Now function to compute network beta-diversity at regional scale (all pairwise comparisons)

beta_reg=function(nw,binary){
	
	B_tot=matrix(nr=length(nw),nc=length(nw))
	rownames(B_tot)=colnames(B_tot)=paste0("nw",1:length(nw))	
	B_rw=B_turn=B_tot

	for(i in 1:length(nw)){for(j in 1:length(nw)){
	
		if(j>=i){next}
		
	x=beta(nw[[i]],nw[[j]],binary)
	B_tot[i,j]=x[[1]]
	B_rw[i,j]=x[[2]]
	B_turn[i,j]=x[[3]]}}
	
	res=list(B_tot,B_rw,B_turn)
	names(res)=c("B_tot","B_rw","B_turn")
	return=res}


### Now function to compare a list of individual networks to their counterpart in the whole regional combined network (termed B'_os in Poisot et al. 2012):

counter=function(nw,binary){

reg_p=unique(unlist(lapply(nw,rownames)))
reg_f=unique(unlist(lapply(nw,colnames)))

regional=matrix(nr=length(reg_p),nc=length(reg_f))
rownames(regional)=reg_p[order(reg_p)]
colnames(regional)=reg_f[order(reg_f)]		

comp=list(length(nw))
for(i in 1:length(nw)){
	comp[[i]]=complete(nw[[i]],regional)
	comp[[i]]=comp[[i]][rownames(regional),colnames(regional)]}

reg=Reduce('+',comp)

res=numeric(length(nw))
for(i in 1:length(nw)){

	nw1=nw[[i]]
	nw2=reg[rownames(nw1),colnames(nw1)]

if(binary==TRUE){
	nw1[nw1>0]<-1
	nw2[nw2>0]<-1}

	t1=nw_to_tab(nw1)
	t2=nw_to_tab(nw2)
	
	t1=t1[order(t1[,1],t1[,2]),]
	t2=t2[order(t2[,1],t2[,2]),]
	
	x=apply(cbind(t1[,3],t2[,3]),2,as.numeric)
	res[i]=sum(abs(apply(x,1,diff)))/sum(apply(x,1,sum))}
	return=res}






########## NOW TROJELSGAARD 2015 Proc Roy Soc B Method:

beta_troj=function(nw1,nw2,binary){

### TOTAL DISSIMILARITY

	if(binary==TRUE){
		nw1[nw1>0]<-1
		nw2[nw2>0]<-1}

	c1=complete(nw1,nw2)
	c2=complete(nw2,nw1)
	
	unshared_p=c(rownames(nw1)[!rownames(nw1)%in%rownames(nw2)],rownames(nw2)[!rownames(nw2)%in%rownames(nw1)])
	unshared_f=c(colnames(nw1)[!colnames(nw1)%in%colnames(nw2)],colnames(nw2)[!colnames(nw2)%in%colnames(nw1)])

	t1=nw_to_tab(c1)
	t2=nw_to_tab(c2)
	
	t1=t1[order(t1[,1],t1[,2]),]
	t2=t2[order(t2[,1],t2[,2]),]
	
	x=cbind(t1[,1:2],apply(cbind(t1[,3],t2[,3]),2,as.numeric))
	I_tot=sum(abs(apply(x[,3:4],1,diff)))
	I_rw=sum(abs(apply(x[!x[,1]%in%unshared_p&!x[,2]%in%unshared_f,3:4],1,diff)))
	I_p=sum(abs(apply(x[x[,1]%in%unshared_p&!x[,2]%in%unshared_f,3:4],1,diff)))
	I_f=sum(abs(apply(x[!x[,1]%in%unshared_p&x[,2]%in%unshared_f,3:4],1,diff)))
	I_pf=sum(abs(apply(x[x[,1]%in%unshared_p&x[,2]%in%unshared_f,3:4],1,diff)))

	
	res=list(c(I_tot,I_tot/I_tot),c(I_rw,I_rw/I_tot),c(I_p,I_p/I_tot),c(I_f,I_f/I_tot),c(I_pf,I_pf/I_tot))
	names(res)=c("I_tot","I_rw","I_plant","I_fungal","I_plant+fungal")
	return=res}







		

