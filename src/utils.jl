function mul!(S::Symmetric, a::Number)
  n = size(S,2)
  for j in 1:n
    for i in 1:j
      S.data[i,j] *= a
    end
  end
  S
end

function sub!(U::UniformScaling, S::Symmetric)
  n = size(S,2)
  mul!(S, -1)
  for i in 1:n
    S.data[i,i] += U.λ
  end
  S
end

function add!(S::Symmetric, T::Symmetric)
  n = size(S,2)
  if size(S) != size(T)
    throw(ArgumentError("Matrices must be same size"))
  end
  for j in 1:n
    for i in 1:j
      S.data[i,j] += T.data[i,j]
    end
  end
  S
end
