# Kennard-Stone algorithm

The Kennard-Stone algorithm was developed in the late 1960's. The objective of this algorithm is to draw $n$ samples from a set $N$ such that samples in $n$ best represent the variability originally present in $N$. 

This can be especially useful, for example, if you want to distribute among experimental treaments some units (e.g., plots) for which you have *a priori* data. For each treatment level, you will certainly want to allocate units in such a way that most of this *a priori* known variability among plots is optimally represented in plots of all treatment levels. You would not want, for example, to test the effect of of a nitrogen fertilizer, but knowing that clay content is spatially variable in your field, and allocate all high-N treatments to a region of the field where clay content is high... You would want both high-N and low-N plots to be found in lower and higher clay content conditions. This can be extended to situations where you know more than clay content, i.e., when you have data for not only 1, but $p$ variables. 

Another frequent use of the Kennard-Stone algorithm is to determine which samples in statistical models will go to (1) calibration of the model vs (2) validation of the model. For example, if you want to calibrate a model predicting soil OM content based on vis-NIR spectra, you will want both you calibration and validation sets of samples to best represent the overall "spectral variability" found in your whole spectral library. Thus, in an ideal world, you would **not** allocate samples to calibration *vs.* validation only through random draws.

<br>
<br>

## How does the algorithm work?
---------------
The algorithm is very simple and works in an iterative manner. Let's present a case scenario where you would want to select $n$ samples from a pool of $N$ to include into a calibration set ($k$) to build a statistical model.For all of your $N$ samples, you have data on $p$ descriptors.

<br>

1. Select the two most distant samples in your set $N$, and put one into $k$. This samples becomes $k_1$. 
2. At eash step $i$, you find, among the remaining samples in $N\notin k$, the one that is furthest away from $k_{i-1}$, the last element you added to $k$ at step $(i-1)$.
3. You go on until you reach an ensemble $k$ of size $n$.


<br>
<br>

## Implementation in R
---------------
A very handy function already exists in the R package ``prospectr`` (function ``kenStone``). Alternatively, you can go on this way:

```R
rm(list=ls())

p=data.frame(rnorm(100,10,2),rnorm(100,4,1),rnorm(100,200,34))
sc.p=scale(p)

d=as.matrix(dist(sc.p))

n=80

k=numeric(n)
k[1]=sample(1:nrow(d),1)


for(i in 2:n){
	d=d[,!colnames(d)%in%k]
	x=d[k[i-1],]
	x=x[x>0]
	k[i]=names(x)[x==min(x)][1]}

	
```

A handy feature in the ``kenStone`` function is that you can choose to calculate distance after dimensionality reduction (if $p$ is very large), and calculate distances based only on principal components capturing at least a certain proportion (user-defined) of the total variance.