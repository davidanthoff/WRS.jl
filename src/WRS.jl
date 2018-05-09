__precompile__()
module WRS

using Distributions

export hd

# This implements (3.16) from (2017) p. 71-72
function hd(X, q=0.5, issorted=false)
    n = length(X)

    a = (n+1)*q
    b = (n+1)*(1-q)

    β = Beta(a, b)

    sorted_X = issorted ? X : sort(X)
    
    θ = 0.
    for i=1:n
        w = cdf(β, i/n) - cdf(β, (i-1)/n)
        θ += w*sorted_X[i]
    end

    return θ
end

end # module
