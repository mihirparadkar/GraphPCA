module GraphPCA

include("utils.jl")

# package code goes here
export GLPCA, reconstruct, reconstruction_error

"""
Retrieves the matrix for finding Y in graph PCA
"""
function getG(X::AbstractMatrix, β::Number, L::AbstractMatrix)
  if β < 0 || β > 1
    throw(ArgumentError("β must be between 0 and 1"))
  end
  n = size(L,2)
  XtX = Symmetric(X'*X)
  L = Symmetric(L)
  #eet = Symmetric(ones(n,n))
  @assert n == size(XtX,2)
  #Largest eigenvalue of XtX to normalize
  λ = eigmax(XtX)
  η = eigmax(L)
  #G = (1 - β)*(I - XtX/λ) + β*(L/η + eet/n)
  G = XtX
  mul!(G, 1/λ)
  sub!(I, G)
  mul!(G, 1 - β)

  #Don't want to copy L
  H = deepcopy(L)
  mul!(H, 1/η)
  #mul!(eet, 1/n)
  #add!(H, eet)
  mul!(H, β)
  add!(G,H)
  G
end

type GLPCA
  U::Matrix{Float64}
  Q::Matrix{Float64}
  X::Matrix{Float64}
  L::Matrix{Float64}
  β::Number
  k::Int
end

"""
Fits a gLPCA model
"""
function GLPCA(X::AbstractMatrix, β::Number, L::AbstractMatrix, k::Int)
  g = getG(X, β, L)
  F = eigfact(g, 1:k)
  Q = F[:vectors]
  U = X*Q
  GLPCA(U, Q, X, L, β, k)
end

"""
Reconstructs the matrix from low-rank components U and Q
"""
function reconstruct(g::GLPCA)
  g.U*g.Q'
end

function reconstruction_error(g::GLPCA)
  g.X - reconstruct(g)
end

end # module
