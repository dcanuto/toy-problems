type ModelParams
    R1::Float64
    R2::Float64
    R3::Float64

    function ModelParams(old=Dict("a"=>0),restart="no")
        this = new()
        if restart == "no"
            this.R1 = 1.;
            this.R2 = 1.;
            this.R3 = 2.;
        elseif restart == "yes"
            this.R1 = old["R1"];
            this.R2 = old["R2"];
            this.R3 = old["R3"];
        end
        return this
    end
end
