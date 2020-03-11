## Getting started
Julia

## Usage

#### Julia: persisterModel construction
```julia
mutable struct persisterParams
    alpha::Float64
    beta::Float64
    tau::Float64
    delta::Array{Int64,2}
    labels::Vector{String}
end

alpha = 0.1
beta = 0.1
tau = 2
delta = [
    -1 1 1;
    1 -1 0
]
labels = ["Growing", "Non-growing"]

mod = persisterParams(
    alpha,
    beta,
    tau,
    delta,
    labels
)

X, T, tsteps = run_gillespie(mod)
```
#### CLI: not implemented yet
```bash
$ julia gillespie.jl
```
