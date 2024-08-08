# Interindividual variation in network $\beta$-diversity

Consider a network of interactions among $n$ plant individuals and $f$ fungal species, as an adjacency matrix $N_{n,f}$. These $n$ plant individuals each have a species affiliation $S$, with $p$ species in total.

If we have two similar networks $N_1$ and $N_2$, comprising the same $i$ species, and we want to compute network $\beta$-diversity based on species level means, we can create, for each respective network $N$, a new $p \times f$ network where each row represents the average interactions of plant species $i$ with the fungal taxon $j$:

$$N'_{ij}=\frac{\sum_{S\in S_i}{N_{Sj}}}{n_i}$$

where 

$S_i =$ the set of rows in N that belong to species $i$
<br>$n_i =$ the number of plant individuals belonging to species $i$ 

<br>These two networks averaged by plant species, $N'_1$ and $N'_2$ can then be easily compared for $\beta$-diversity computation using Sorensen's index:

$$D=1- \frac{2\:\:\sum_{i=1}^p\sum_{j=1}^k|N'_{1,ij}\cap N'_{2,ij}|\:}{\ \ \sum_{i=1}^p\sum_{j=1}^k|N'_{1,ij}|+|N'_{2,ij}|\ \ }$$



