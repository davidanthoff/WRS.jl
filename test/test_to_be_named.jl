# example of to_be_named

using Distributions

include("../src/WRS.jl")
using .WRS

n = 10000
d = Normal()

# generate data
sample_a = rand(d, n)
sample_b = rand(d, n)

# run comparison
results = to_be_named(sample_a, sample_b, hd)
