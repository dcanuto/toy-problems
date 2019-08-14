type CVSystem
    Pb::Vector{Float64}
    Q1::Vector{Float64}
    Q2::Vector{Float64}
    Q3::Vector{Float64}
    t::Vector{Float64}
    sparams::SolverParams
    mparams::ModelParams

    function CVSystem(nvar=0,old=Dict("a"=>0),restart="no")
        this = new()
        if restart == "no"
            this.mparams = ModelParams();
        elseif restart == "yes"
            mparams = old["system"]["mparams"];
            this.mparams = ModelParams(mparams,restart);
        end
        this.Pb = [];
        this.Q1 = [];
        this.Q2 = [];
        this.Q3 = [];
        this.t = [];
        this.sparams = SolverParams(nvar);
        return this
    end
end
