using GraphPCA
using Base.Test

# X matrix
X1 = 10*randn(20,4)
X2 = X1 + randn(20,4)
X = [X1 X2]

#Construct the Laplacian matrix
L = eye(8)
for i in 1:4
  j = i + 4
  L[i,j] = -1
  L[j,i] = -1
end

λ = eigmax(Symmetric(X'*X))
η = eigmax(Symmetric(L))
eet = ones(8,8)
n = 8

#Pure PCA
XtX = Symmetric(I - X'*X/λ)
Dpca,Vpca = eig(XtX, 1:4)
Qpca = Vpca
Upca = X*Qpca
Xpca = Upca*Qpca'

#Pure Laplacian Embedding
Ln = Symmetric(L/η)
Dn, Vn = eig(Ln, 1:4)
Qlap = Vn
Ulap = X*Qlap
Xlap = Ulap*Qlap'

#Laplacian Embedding with Offset
# Lnp = Symmetric(L/η + eet/n)
# Dnp,Vnp = eig(Lnp, 1:4)
# Qlapp = Vnp
# Ulapp = X*Qlapp
# Xlapp = Ulapp*Qlapp'

GLPCA0 = GLPCA(X, 0, L, 4)
GLPCA1 = GLPCA(X, 1, L, 4)
@test GLPCA1.U*GLPCA1.Q' ≈ Xlap
@test GLPCA0.U*GLPCA0.Q' ≈ Xpca
