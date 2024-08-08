##### Basic function to calculate beta_sor between two networks:
sor=function(nw1,nw2){
	x=cbind(c(nw1),c(nw2))
	return=1-((2*sum(apply(x,1,min)))/sum(apply(x,1,sum)))
	}


##### Basic function to compare with and without intrasp variation
beta_int=function(nw1,nw2,sp1,sp2){

# keep only intersect of spp between two nws
nw1=nw1[sp1%in%intersect(sp1,sp2),]
nw2=nw2[sp2%in%intersect(sp1,sp2),]
sp1=sp1[sp1%in%intersect(sp1,sp2)]
sp2=sp2[sp2%in%intersect(sp1,sp2)]

# remove any empty row/column
nw1=nw1[rowSums(nw1)>0,colSums(nw1)>0]
nw2=nw2[rowSums(nw2)>0,colSums(nw2)>0]

# Calculate raw sorensen with intrasp variation
sor_intra=sor(nw1,nw2)

# Build nw without intrasp variation
spp=levels(as.factor(intersect(sp1,sp2)))
nw1_no=matrix(nr=length(spp),nc=ncol(nw1))
nw2_no=matrix(nr=length(spp),nc=ncol(nw2))
for(i in 1:length(spp)){
	nw1_no[i,]=colMeans(nw1[sp1==spp[i],])
	nw2_no[i,]=colMeans(nw2[sp2==spp[i],])
	}
sor_no=sor(nw1_no,nw2_no)
	
res=list(sor_intra,sor_no)
names(res)=c("Beta_total","Beta_no_intrasp_variation")
return=res}
	


nw1=matrix(round(rlnorm(400,.5,2)),nc=20)
nw2=matrix(round(rlnorm(400,.5,2)),nc=20)
sp1=rep(c("A","B","C","D"),each=5)
sp2=rep(c("B","C","D","E"),each=5)
print(beta_int(nw1,nw2,sp1,sp2))


