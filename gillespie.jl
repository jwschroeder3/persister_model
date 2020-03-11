using Plots

"""Gillespie algorithm.

    Args:
        x (ndarray(int)): Initial counts for growing, non-growing cells
        alpha (float): switch constant from growing to non-growing
        beta (float): switch constant from non-growing to growing
        tau (float): division constant for growing cells
        delta (ndarray): 2D diffusion matrix for system
        N (int): number of iterations for Gillespie

    returns:
        X (ndarray(int)): Trace of component counts for each iteration.
        T (ndarray(float)): Time adjusted trace of time during simulation.
        tsteps (ndarray(float)): Time weight trace; duration of time spent in each state.
"""
function gillespie(x, alpha, beta, tau, delta, N)
    
    # Initialisation
    t = 0
    T = zeros(N)
    tsteps = zeros(N)
    X = zeros((size(delta)[1], N))
    t_arr = rand(N)
    choice_arr = rand(N)

    # Simulation
    for i in 1:N
        # Determine rates
        rates = [
            x[1] * alpha,
            x[2] * beta,
            x[2] * tau
        ]
        summed = sum(rates)

        # Determine WHEN state change occurs
        t_step = (-1) / summed * log(t_arr[i])
        t = t + t_step
        T[i] = t
        tsteps[i] = t_step

        # Determine WHICH reaction occurs with relative propabilities
        reac = sum(isless.(cumsum(rates/summed), choice_arr[i])) + 1 # add one to get to 1-indexed julia
        x = x + delta[:, reac]
        X[:, i] = x
    end
    return (X, T, tsteps)
end


function plot_gillespie(model, T, X)
    plot_array = Any[]
    plot_number = size(X)[1]
    for i in 1:plot_number
        push!(
            plot_array,
            plot(
                T,
                X[i,:],
                label=model.labels[i]
            )
        )
    end
    plot(plot_array..., layout=(length(plot_array), 1))
end

# """Computes theoretical means, variances, and covariances for system.  
#     returns: tuple of means, variances, and covar
# """
# function compute_theoretical(model)
#     # Theoretical states
#     x1 = model.alpha * model.tau1
#     x2 = model.lambd * model.tau2 * x1
#     x3 = model.lambd * model.tau3 * x1

#     # Theoretical variance and covariances
#     v1 = 1 / x1
#     v2 = 1 / x2 + model.tau1 / (x1 * (model.tau1 + model.tau2))
#     v3 = 1 / x3 + model.tau1 / (x1 * (model.tau1 + model.tau3))

#     cov12 = model.tau1 / (x1 * (model.tau1 + model.tau2))
#     cov13 = model.tau1 / (x1 * (model.tau1 + model.tau3))
#     cov23 = (
#         model.tau1
#         * (2 * model.tau2 * model.tau3 + model.tau1 * (model.tau2 + model.tau3))
#     ) / (
#         x1
#         * (model.tau1 + model.tau2)
#         * (model.tau1 + model.tau3)
#         * (model.tau2 + model.tau3)
#     )

#     return ([x1, x2, x3], [v1, v2, v3], [cov12, cov13, cov23])
# end

"""Runs Gillespie algorithm using model parameters.
    Args:
    ------------------

"""
function run_gillespie(model; N=10000, start_near_ss=false)
    
    # Optional flag for starting simulation from predicted steady states when testing. 
    # false by default.
    # if start_near_ss
        # x = compute_theoretical(model)[1]
    # else
    x = [
        1000000,
        10000    
    ]
    # end
    X, T, tsteps = gillespie(
        ceil.(Int,x),
        model.alpha,
        model.beta,
        model.tau,
        model.delta,
        N,
    )
    return(X, T, tsteps)
end

# """Calculates difference in fluxes for each component.
#     returns: List containing the flux balances for each component.
# """
# function get_fluxes(model,X,tsteps)
#     R1p = sum(
#         fill(model.alpha, length(X[1])[1]) * tsteps
#     ) / sum(tsteps)
#     R1m = sum((X[1] / model.tau1) * tsteps) / sum(tsteps)
#     R2p = sum((model.lambd * X[1]) * tsteps) / sum(tsteps)
#     R2m = sum((X[2] / model.tau2) * tsteps) / sum(tsteps)
#     R3p = R2p
#     R3m = sum((X[3] / model.tau3) * tsteps) / sum(tsteps)
#     return [R1p - R1m, R2p - R2m, R3p - R3m]
# end


"""
Composit type for storing attributes of the dynamics of mRNA for two genes regulated by a transcription factor.
    
    Attributes:
        alpha (float): switch constant from growing to non-growing
        beta (float): switch constant from non-growing to growing
        tau (float): doubling constant for growing population
        delta (ndarray): 2D diffusion matrix for system
        labels (list(str)): labels for each component in system
"""
mutable struct persisterParams
    alpha::Float64
    beta::Float64
    tau::Float64
    delta::Array{Int64,2}
    labels::Vector{String}
end

alpha = 0.1
beta = 0.2
tau = 200
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

X, T, tsteps = run_gillespie(mod, N=1000000)

plot_gillespie(mod, T, X)