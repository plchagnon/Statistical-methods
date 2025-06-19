# PHATE in Python

PHATE is a dimensionality reduction method introduced by [Moon et al. *2019*](https://www.nature.com/articles/s41587-019-0336-3). It compares favorably to other commonly used method such as the classical PCA or *t*-SNE (stochastic neighbour embedding). The former could be thought to place too much emphasis on global fitting (i.e., each principal component finding the linear dimension most aligned with the data globally) leaving behind some details related to "close neighbourhood". On the other hand *t*-SNE could be thought to focus overly on local fitting, better suited to elucidate the presence of well-separated **clusters**, but lacking the ability to depict **global**, or overall, transitions and structure in the data. 

Below, I present an example usage of PHATE (for **P**otential of **H**eat diffusion for **A**ffinity-based **T**ransition **E**mbedding) within Python. To allow usage from within Python, I simply ran ``pip install --user phate``, but had to first downgrade ``numpy`` to an older version:

```{Python}
pip install numpy==1.24.4
pip install --user phate
```

Once installed, one can simply open a Powershell, run Python, and go:

```{Python}
import phate
```

<br><br>

## Case example with an imported table
Let's say we would like to conduct a dimensionality reduction on a **microbial metacommunity**, that would be stored in a tabular format in the object ``comm.txt``. 

<br>

**Import the table**
<div style="margin-top: -0.4cm;"></div>

```Python
### From within Python
import pandas as pd
data=pd.read_csv("comm.txt",sep="\t",index_col=0)
data
```

<br>

The last line allows to confirm via printing the overview of the object and summary, that you have properly imported the correct table, with proper number of rows and columns. Then, we can further proceed with the PHATE embedding per se:
 <br> 

**Define PHATE operator**
<div style="margin-top: -0.4cm;"></div>

```Python
### Import phate
import phate

### Define the phate "operator" with wanted parameters
phate_op=phate.PHATE(
        n_components = 2,
        knn = 4,
        mds_solver = 'smacof',
        knn_dist = 'euclidean',
        mds_dist = 'euclidean',
        n_jobs = -2,
        verbose = 1)
```
<br>

All the parameters to be selected are listed and explained at [this page](https://phate.readthedocs.io/en/stable/). Here, we have selected two PHATE dimensions, $k=4$ nearest neighbours to build kernels, a MDS solver that is slightly slower but more optimal than the other ``sdg`` alternative, Euclidean distance as the metric to compute pairwise distance for both knn kernels and MDS computation, and jobs to be distributed on all but one cores on the computer (``n_jobs=-2``). Verbose not being equal to 0 means that all output will be printed on the console.

Then we can fit this embedding on our data object:

**Compute PHATE**
<div style="margin-top: -0.4cm;"></div>

```Python
fit = phate_op.fit_transform(data)

### Get it as a data frame
fit_out=pd.DataFrame(fit,index=data.index,columns=["Dim1","Dim2"])

### Export it
fit_out.to_csv("phate_out.txt",sep="\t")
```

<br><br>

### To come

* Smoothing: Diffusion matrix (rows $\times$ rows) can be multiplied with original matrix to modify it according to expected geometry and diffusion probabilities across rows. This is a method thought to **reduce noise** in the detection of columns! Genes, in Moon et al. 2019. Could it help explore plant-microbial networks?? Where associated microbes are a noisy thing to detect?? Maybe **modules** in such networks would make more biological sense if investigated on "smoothed" networks?