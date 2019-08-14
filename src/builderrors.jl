type Errors
    a::Float64 # parameter distribution smoothing
    h::Float64

    odev::Vector{Float64} # σ for observation, param distributions
    pdev::Vector{Float64}

    lb::Vector{Float64} # lower/upper bounds for truncated normal param distrs.
    ub::Vector{Float64}

    function Errors(nparams=9)
        this = new();
        δ = 0.99; # ∃ (0, 1, typically within 0.95-0.99), lower values = less parameter dispersion
        this.a = (3*δ-1)/(2*δ);
        this.h = sqrt.(1-this.a^2);
        this.odev = [0.1,0.1,0.1,0.1];
        this.pdev = [0.1,0.1,0.1]; # no params
        this.lb = -Inf*ones(nparams);
        this.ub = Inf*ones(nparams);
        this.lb[1] = 1e-3; # R1
        this.lb[2] = 1e-3; # R2
        this.lb[3] = 1e-3; # R3
        return this
    end
end
