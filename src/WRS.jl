__precompile__()
module WRS

using Distributions
using Distributed
using DataFrames

export hd, pb2gen, to_be_named

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

function pb2gen(x, y, est, q=0.5; alpha=0.05, nboot=2000)
    bootstrapped_diff_est = pmap(1:nboot) do i
        sampled_values_x = sample(x, length(x), replace=true)
        sampled_values_y = sample(y, length(y), replace=true)

        sort!(sampled_values_x)
        sort!(sampled_values_y)

        return est(sampled_values_x, q, true) - est(sampled_values_y, q, true)
    end
    sort!(bootstrapped_diff_est)

    low = round(Int, (alpha/2)*nboot)+1
    up = nboot-low

    est_x = est(x, q)
    est_y = est(y, q)
    est_diff = est_x - est_y
    ci = (bootstrapped_diff_est[low], bootstrapped_diff_est[up])

    A = count(i->i<0,bootstrapped_diff_est)
    C = count(i->i==0,bootstrapped_diff_est)    
    p_hat_star = A/nboot+0.5*C/nboot
    pvalue = 2*(min(p_hat_star,1-p_hat_star))

    se = var(bootstrapped_diff_est)

    return est_x, est_y, est_diff, ci, pvalue, se
end

function to_be_named(x, y, est, quantiles=[0.05, 0.25, 0.5, 0.75, 0.95]; alpha=0.05, nboot=2000)
    results = DataFrame(q=Float64[], est_x=Float64[], est_y=Float64[], est_diff=Float64[], ci_low=Float64[], ci_up=Float64[], pvalue=Float64[], se=Float64[], signif=Bool[])
    for q in quantiles
        est_x, est_y, est_diff, ci, pvalue, se = pb2gen(x, y, est, q; alpha=alpha, nboot=nboot)
        # push!(results, q, est_x, est_y, est_diff, ci[1], ci[2], pvalue, se, !(ci[1]<0 && ci[2]>0))
        push!(results, (q, est_x, est_y, est_diff, ci[1], ci[2], pvalue, se, pvalue < alpha))
    end
    return results
end

end # module
