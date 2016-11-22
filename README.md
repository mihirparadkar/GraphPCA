# GraphPCA
Implements graph-regularized PCA as described in Jiang et al. (Graph-Laplacian PCA: Closed-form Solution and Robustness)


http://www.cv-foundation.org/openaccess/content_cvpr_2013/papers/Jiang_Graph-Laplacian_PCA_Closed-Form_2013_CVPR_paper.pdf

Use the function
```julia
GLPCA(X, β, L, k)
```
to construct an object of type GLPCA, which also solves for the low-rank components.

Use
```julia
reconstruct(g::GLPCA)
```
to get the reconstruction
