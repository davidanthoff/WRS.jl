using WRS
using Distributions
using Test

n = 10000
d = Normal()

# generate data
sample_a = rand(d, n)
sample_b = rand(d, n)

# run comparison
results = pb2gen(sample_a, sample_b)

@test length(results) == 5
