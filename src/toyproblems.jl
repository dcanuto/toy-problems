module toyproblems

using MAT
using Distributions

# solver parameters
include("solverparams.jl")

# type definitions
include("buildall.jl")

# data assimilation
include("builderrors.jl")
include("paramwalk.jl")

export SolverParams
export ModelParams
export buildall
export builderrors
export paramwalk!

end
