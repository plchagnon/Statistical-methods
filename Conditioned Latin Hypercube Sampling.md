# Conditioned Latin Hypercube Sampling in R
Introduced by [Minasny & McBratney (2006)](https://doi.org/10.1016/j.cageo.2005.12.009), the conditioned latin hypercube sampling (cLHS) is a useful strategy to take advantage of available data a priori (e.g., spectral data, satellite imagery, veris soil surveys, etc.) to guide a sampling scheme aiming at best capturing the inherent variability in a given landscape. This technique is an extension of the latin hypercube sampling, itself being an extension of the latin square. At its very simplest expression, the latin square simply describes a two-dimensional scheme where a given attribute (e.g., inclusion in a sampling design) is **only represented once** in each stratum for dimension 1 and dimension 2.   

Say you have, for example, an experiment with 9 plots distributed as a 3 x 3 grid in the field, and three different treatments to distribute, then a latin square would be one where treatment A is only found once in each column AND each row of the grid. In principle, this readily extends to situations where dimensions are >2. We can design sampling schemes that sample each stratum (equivalent to rows and columns of your latin square) only once for all dimensions. In a situation where you wished to guide a sampling design based on soil spectral surveys, each wavelength for which you have an absorbance/reflectance data value for all samples is a dimension in your hypercube. And you would wish to include in your sampling design soils that represent distinct strata for this wavelength (i.e., that have higher or lower reflectance values for this wavelength). However, while this can work in theory, it is very probable that some combinations of strata for different dimensions are poorly represented in your samples. For example, it may well be that there are simply no soils that have high reflectance for $\lambda=450$ nm and low reflectance for $\lambda=382$ nm. This is where the cLHS sampling becomes interesting: it finds an optimal trade-offs between (1) best matching the latin hypercube standards and (2) dealing with samples that have heterogeneous frequency distributions for the different dimensions you have data for. You end up defining a sampling ensemble $s$ that results from an optimization procedure (using simulated annealing) that best matches theoretical expectations from the latin hypercube.


This process works in few easy steps:


1. You create an initial partition of your data points into $s$ samples that will be analyzed and the rest of the dataset (i.e., the reserve)
2. You calculate a cost ( $C$ ) function representing how much your partition deviates from expectations of the latin hypercube
3. You conduct a simulated annealing algorithm to :
    - switch sample(s) in and out of your ensemble $s$
    - re-calculate $C$
    - Accept or reject the switch above

After such a process, you end-up with an optimized partition of samples with a $s$ ensemble that best fits the latin hypercube. Let's do a test run with some imaginary samples and a few dimensions for which we have data to characterize them.

  
<br>

```R
rm(list=ls())
### Assume a dataset "dat" with categorical and non categorical variables, for N samples:
N=80
dat=data.frame(rnorm(N,200,5),rnorm(N,20,3),rnorm(N,10,1),rnorm(N,400,30),sample(c("A","B"),N,replace=T),sample(c("C","D"),N,replace=T))
categ=c(0,0,0,0,1,1)      ##### vector specifying what variables are categorical
names(dat)=paste("var",1:ncol(dat),sep=" ")
```

<br>

Then, imagine you have sufficient funds to sample and characterize thoroughly only 12 samples $n=12$ out of the 80 samples overall. This defines also the number of strata into which you will separate your full data vectors for all continuous variables included in it. Here, you have 4 of them in the example (4 data vectors) with 80 values in each of them. So you will want to first randomly select 12 samples out of your 80 samples:
  
<br> 

```R
n=12
s=dat[sample(1:nrow(dat),n,replace=F),]
```

<br>

Then you will want to distinguish continuous from categorical variables, and separate continuous variables' vectors into $n$ equiprobable strata. This will be important for the calculation of your cost function $C$:
  
<br>

```R
dat_con=as.matrix(dat[,categ==0])		### continuous variables in the full dataset
dat_cat=as.matrix(dat[,categ==1])		### categorical variables in the full dataset
s_con=as.matrix(s[,categ==0])			### continuous variables in the ensemble "s"	
s_cat=as.matrix(s[,categ==1])			### categorical variables in the ensemble "s"
f=function(x){quantile(x,probs=seq(0,1,length.out=n+1))}
quan=apply(as.matrix(dat_con),MARGIN=c(2),FUN=f)
```

<br>

Next step is to define a R function to compute our cost function $C$. Minasny and McBratney (2006) described this cost function using three additive components:

- $C_1$ is related to the representation of each data vector *stratum*, or quantile, in the tentative sampling scheme. In a perfect latin hypercube, each of the $n$ strata should be represented exactly once. As this may not be strictly possible with any given subset of the full dataset, the objective si to calculate the deviation between the "ideal" scenario and the actual strata representation in the realized sampling scheme.
- $C_2$ is somewhat analogous to $C_1$ for categorical variables. In an ideal scenario, the relative representation of each class/category for each variable should be the same in the full dataset vs. the expected $n$ samples to be part of the sampling scheme.
- $C_3$ evaluates the correlation patterns among the continuous variables, to determine, as for $C_2$, if the expected $n$ selected samples mirror patterns found in the full dataset. 

These three components can be weighted as desired by user, if one of the components was deemed more significant in a given case scenario. Additional cost components could also be added as needed. For example, constraints related to the inclusion of specific samples in the $n$ subset (e.g., cost related to ease of sampling, travel distances, etc.). In the simplest expression of the cost function formulated by Minasny and McBratney (2006), the cost function is calculated as:

$C=C_1+C_2+C_3$

<br>

```R
Obj=function(dat,n,categ,w1,w2,w3,quan,s){

## "dat" is the initial dataset, with all samples as rows and all ancillary variables as columns
## "n" is the number of datapoints to be included in the sampling design
## "categ" is a vector specifying is each variable (column in "dat") is categorical (1) or not (0)
## "w1","w2","w3" are weights to be put on Objective functions 1-2-3 for the overall objective function
## "quan" = for all continuous variables in "dat", the quantile distributions in "n" strata
## "s" is the portion of dataset put into the sampling bin

dat_con=as.matrix(dat[,categ==0])		### continuous variables in the full dataset
dat_cat=as.matrix(dat[,categ==1])		### categorical variables in the full dataset
s_con=as.matrix(s[,categ==0])			### continuous variables in the ensemble "s"	
s_cat=as.matrix(s[,categ==1])			### categorical variables in the ensemble "s"

#### For each variable in "dat", count the number of samples in "s" found in each of the "n" strata:
counts=matrix(nr=n^2,nc=ncol(s_con)+2)
colnames(counts)=c("sample","stratum",colnames(s_con))
counts[,1]=rep(1:n,each=n)
counts[,2]=rep(1:n,n)
for(i in 1:nrow(counts)){
	for(j in 3:ncol(counts)){
		counts[i,j]=as.numeric(s_con[counts[i,1],j-2]>=quan[counts[i,2],j-2]&s_con[counts[i,1],j-2]<=quan[counts[i,2]+1,j-2])}}

#### Summarize results above into a matrix "d" summarizing counts per stratum per variable
d=matrix(nr=n,nc=ncol(s_con))
for(i in 1:nrow(d)){
	d[i,]=colSums(counts[counts[,2]==i,3:(ncol(dat)-sum(categ))])}
colnames(d)=paste("var",1:(ncol(dat)-sum(categ)),sep=" ")
rownames(d)=paste("stratum",1:n,sep=" ")

#### Likewise, summarize results from "counts" into a matrix summarizing counts for each sample for each stratum (will serve to determine which sample from "s" to swap out of "s" to try improving the subsample "s")
e=matrix(nr=n,nc=n)
for(i in 1:nrow(e)){
	temp=counts[counts[,1]==i,3:ncol(dat)]
	for(j in 1:ncol(e)){e[i,j]=rowSums(temp)[j]}}
colnames(e)=paste("stratum",1:n,sep=" ")
rownames(e)=paste("sample",1:n,sep=" ")

#### Calculate cost function #1
O1=sum(abs(d-1))

#### Determine what stratum deviates most from the perfect latin hypercube, and what sample is most represented in this stratum
discrep_stratum=sample(1:n,1,prob=as.numeric(rowSums(d)==max(rowSums(d))))
rem=rownames(s)[sample(1:n,1,prob=as.numeric(e[,discrep_stratum]==max(e[,discrep_stratum])))]

#### Calculate discrepancies in relative representation of each category for all categorical variables, between "s" and "dat", to calculate cost function #2
if(sum(categ)==0){Discrep=0}else{
Discrep=numeric(ncol(s_cat))
for(i in 1:ncol(dat_cat)){
	nlev=length(levels(as.factor(dat_cat[,i])))
	temp1=temp2=numeric(nlev)
	for(j in 1:nlev){
		temp1[j]=length(s_cat[,i][s_cat[,i]==levels(as.factor(dat_cat[,i]))[j]])
		temp2[j]=length(dat_cat[,i][dat_cat[,i]==levels(as.factor(dat_cat[,i]))[j]])}}
Discrep=temp1/sum(temp1)-temp2/sum(temp2)}
O2=sum(abs(Discrep))

#### Calculate cost function 3 with continuous variables
O3=sum(abs(cor(s_con)-cor(dat_con)))
O=w1*O1+w2*O2+w3*O3

#### Return the cost function itself, and the sample that should be removed from "s" in the next iteration of the optimization algorithm (to be used later on...)
return=as.numeric(c(O,rem))}

```

<br>  

Now that we have our cost function in hand, we can integrate it into an optimization algorithm that uses simulated annealing to solve such NP-hard problem. Simulated annealing is a generic optimization strategy that uses a cooling factor to alter the conservatism level in the optimality search criterion/criteria as iterations proceed. In the early phase of the algorithm, the cooling factor is still "hot", and allows for non-optimal movements to be made (here, swaps of samples between the selected $n$ samples and the full dataset that could increase the overall cost function $C$). As the algorithm progresses, such movements will become increasingly improbable, up to a point where they will receive a probability of 0 ("cold"/late phase of the algorithm). We must thus first determine a shape for the Temperature $\times$ time/iterations curve (typically an exponential decay). In the example below, we will use an initial Temperature of 100 and a decay factor of 0.995, such that:

$T_t=T_{in}(0.995)^t$


```R
#### Set initial paremeters (Tin, cooling constant/factor, nb of runs to be made, and the number of runs to consider to determine if the cost function is stagnating or not (based on the coefficient of variation in the last "min_stationary" runs...
Tin=100
const=.995
nrun=50000
min_stationary=10000

####Calculate the cost function with the initial data partition "s"
Oin=Obj(dat=dat,n=n,categ=categ,w1=1,w2=1,w3=1,quan=quan,s=s)
  
#### Then create a vector to store results during the algorithm, and start the loop:
Oval=numeric(nrun)
for(i in 1:nrun){

  T=Tin*const^i
  
  ####Then create a swap in the ensemble "s", removing the most "problematic" sample from it and replacing this sample with a random draw from the rest of the full dataset. Calculate the cost function from this alternative data partition called "s2".
	s2=s[!rownames(s)%in%Oin[2],]
	s2=rbind(s2,dat[sample(rownames(dat[!rownames(dat)%in%rownames(s2),]),1),])
	Ofin=Obj(dat=dat,n=n,categ=categ,w1=1,w2=1,w3=1,quan=quan,s=s2)
	
  #### Determine the impact of the swap on the cost function. Swaps that improve results (i.e., reduce cost function) will be accepted with probability = 1. Swaps that increase cost function could be accepted, based on the actual difference between cost functions for the two alternative "s" and "s2" considered, and on the cooling factor. 
  delta=Ofin[1]-Oin[1]
	metro=exp(-delta/T)
	rand=runif(1,0,1)	
	if(rand<metro){
		s=s2
		Oval[i]=Ofin[1]
    Oin=Ofin}
	if(rand>metro){
		s=s
		Oval[i]=Oin[1]}	

  ####If the cost function has been stationary (cv < 0.001, user-defined) in the last 10 000 steps (user-defined), stop the algorithm. Alternatively, keep going with the loop...
	if(i>min_stationary){cv=sd(Oval[(i-min_stationary):i])/mean(Oval[(i-min_stationary):i])}
	if(exists("cv")==FALSE){cv=1}
	if(cv<.001)break
	}

#### We can plot the results of the algorithm
plot(Oval[1:nrun]~c(1:nrun),type="l")

```

<br>

The last line can be useful to visually confirm that (1) the cost function was indeed stationary enough during the last phases of the algorithm, confirming that the number of time steps was sufficient, and (2) the algorithm was "liberal" enough during the early time steps, to ensure that various possibilities were explored even if they could involve non-optimal swaps. This is the key to avoid getting trapped in local optima, and hoping to find the global minimum for this cost function. It is also wise to re-run the algorithm more than once, as this process is heuristic: optimal solution could change from run to run. Large variations in the cost function $C$ for different runs would be a symptom of poor algorithm parameters. 
