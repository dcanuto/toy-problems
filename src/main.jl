using toyproblems

function main()
    # option flags
    rstflag = "no"

    # initialization
    tf = 30; # total number of cardiac cycles
    ensemblesize = 50;
    if rstflag == "no"
        y0 = [[rand(Distributions.Normal(1,0.1));0;0;0] for i=1:ensemblesize]; # random
        nvar = [4 for i=1:ensemblesize];
        systems = pmap((a1)->toyproblems.CVSystem(a1),nvar);
    elseif rstflag == "yes"
        fnames = ["system_$i.mat" for i=1:ensemblesize];
        vars = [MAT.matread(fnames[i]) for i=1:ensemblesize];
        nvar = [4 for i=1:ensemblesize];
        rst = [rstflag for i=1:ensemblesize];
        systems = pmap((a1,a2,a3)->toyproblems.CVSystem(a1,a2,a3),nvar,vars,rst);
        y0 = [[vars[1]["system"]["Pb"][end];vars[1]["system"]["Q1"][end];
            vars[1]["system"]["Q2"][end];vars[1]["system"]["Q3"][end]] for i=1:ensemblesize];
    end
    t0 = [0 for i=1:ensemblesize];
    sparams = [systems[i].sparams for i=1:ensemblesize];
    mparams = [systems[i].mparams for i=1:ensemblesize];
    to = Vector{Float64}[];
    yo = Vector{Vector{Float64}}[];
    append!(yo,[[zeros(1)]]);
    append!(to,[zeros(1)]);

    # data assimilation setup
    # allocators for ensemble augmented state, measurements
    nparams = 3;
    X = [zeros(nvar[1]) for i in (1:length(systems))];
    θ = [zeros(nparams) for i in (1:length(systems))];
    Y = [zeros(4) for i in (1:length(systems))];

    # parameter scalings
    θs = zeros(nparams);
    for i = 1:length(systems)
        θs[1] += mparams[i].R1;
        θs[2] += mparams[i].R2;
        θs[3] += mparams[i].R3;
    end
    θs /= ensemblesize;

    # allocators for state, forecast measurement mean
    xhat = zeros(nvar[1]);
    θhat = zeros(nparams);
    yhat = zeros(4);

    # allocators for measurement replicates
    yi = [zeros(4) for i in 1:length(systems)];

    # output variables for ensemble avg. state, meas., params. w/ times
    tout = Float64[];
    xout = Vector{Float64}[];
    xoutv = Vector{Float64}[];
    yout = Vector{Float64}[];
    youtv = Vector{Float64}[];
    innov = Vector{Float64}[];
    ioutv = Vector{Float64}[];
    θout = Vector{Float64}[];
    θoutv = Vector{Float64}[];
    Pθout = Vector{Float64}[];
    Pθoutv = Vector{Float64}[];
    Pxout = Vector{Float64}[];
    Pxoutv = Vector{Float64}[];j
    lbx = Vector{Float64}[];
    ubx = Vector{Float64}[];
    lbxv = Vector{Float64}[];
    ubxv = Vector{Float64}[];

    # distributions for state, measurement, parameters
    error = toyproblems.Errors(nparams);

    # RTPS allocators
    p = 0.5; # relaxation amount, ∃ [0.5, 0.95]
    c = zeros(length(systems));
    σb = zeros(length(xhat));
    σa = zeros(length(xhat));
    σtb = zeros(length(θhat));
    σta = zeros(length(θhat));

    # "true" system
    R1true = 1.;
    R2true = 1.;
    R3true = 2.;

    # solver loop
    tic()
    while to[1][end] <= tf
        println("Assimilating data at t = $(to[1][end])")

        # non-dimensional measurement
        P1true = sin(2*pi*to[1][end]);
        Pbtrue = P1true/(R1true*(1/R1true+1/R2true+1/R3true));
        Q1true = (P1true - Pbtrue)/R1true;
        Q2true = Pbtrue/R2true;
        Q3true = Pbtrue/R3true;
        y = [P1true,Pbtrue,Q2true,Q3true];

        # generate measurement replicates
        for i = 1:length(systems)
            yi[i][1] = y[1] + rand(Distributions.Normal(0,error.odev[1]));
            yi[i][2] = y[2] + rand(Distributions.Normal(0,error.odev[2]));
            yi[i][3] = y[3] + rand(Distributions.Normal(0,error.odev[3]));
            yi[i][4] = y[4] + rand(Distributions.Normal(0,error.odev[4]));
        end

        # forecast parameters and means
        for i = 1:length(systems)
            θ[i][1] = mparams[i].R1;
            θ[i][2] = mparams[i].R2;
            θ[i][3] = mparams[i].R3;
        end

        # forecast parameters (dimensional)
        # println("Forecast parameters:")
        θ = [WK3.paramwalk!(error,θs,θ[i]) for i in 1:length(systems)];

        # non-dimensionalize parameters
        for i = 1:length(systems)
            θ[i] = θ[i]./θs;
        end

        # RTPS prior standard deviation
        for i in 1:length(θ[1])
            for j in 1:length(systems)
                c[j] = θ[j][i];
            end
            σtb[i] = std(c;corrected=true);
        end

        # parameters back into model
        for i = 1:ensemblesize
            mparams[i].R1 = θ[i][1]*θs[1];
            mparams[i].R2 = θ[i][2]*θs[2];
            mparams[i].R3 = θ[i][3]*θs[3];
        end

        # forecast state w/ forecast parameters (single time step)
        for i = 1:ensemblesize
            P1 = sin(2*pi*to[i][end]);
            Q1true = (P1true - Pbtrue)/R1true;
            Q2true = Pbtrue/R2true;
            Q3true = Pbtrue/R3true;
            append!(systems[i].Pb,P1/(mparams[i].R1*(1/mparams[i].R1+1/mparams[i].R2+1/mparams[i].R3)
            append!(systems[i].Q1,(P1-systems[i].Pb[end])/mparams[i].R1)
            append!(systems[i].Q2,systems[i].Pb[end]/mparams[i].R2)
            append!(systems[i].Q3,systems[i].Pb[end]/mparams[i].R3)
        end

        # vector of forecast state vectors
        X = [[systems[i].Pb[end],systems[i].Q1[end],systems[i].Q2[end],systems[i].Q3[end]] for i in 1:ensemblesize];

        # vector of forecast measurements
        Y = [[systems[i].Pb[end],systems[i].Q1[end],systems[i].Q2[end],systems[i].Q3[end]] for i in 1:ensemblesize];

        # forecast mean state, parameters, measurement
        xhat = mean(X);
        yhat = mean(Y);
        θhat = mean(θ);
        # println("Normalized ̂x, first forecast: $xhat")
        # println("Normalized ̂y, first forecast: $yhat")
        # println("Normalized ̂θ, first forecast: $θhat")

        # forecast params./meas. cross covariance, meas. covariance
        Pty = zeros(length(θhat),length(yhat))
        Pyy = zeros(length(yhat),length(yhat))
        for i = 1:ensemblesize
            Pty += *((θ[i] .- θhat),(Y[i] .- yhat)');
            Pyy += *((Y[i] .- yhat),(Y[i] .- yhat)');
        end
        Pty ./= (ensemblesize);
        Pyy ./= (ensemblesize);

        # add meas. noise to meas. covariance (allows invertibility)
        Pyy += diagm([(error.odev[1])^2,(error.odev[2])^2,(error.odev[3])^2,(error.odev[4])^2],0);

        # parameter Kalman gain
        K = Pty*inv(Pyy);
        # println("Parameter Kalman gain:")
        # display(K)

        # parameter analysis step
        for i = 1:ensemblesize
            # println("ith measurement replicate: $(yi[i])")
            # println("ith forecast measurement: $(Y[i])")
            θ[i][:] += K*(yi[i] .- Y[i]);
        end

        # println("Normalized θ after analysis: $(mean(θ))")

        # RTPS parameter covariance inflation
        for i in 1:length(θ[1])
            for j in 1:length(systems)
                c[j] = θ[j][i];
            end
            σta[i] = std(c;corrected=true);
        end
        θhat = mean(θ);
        for i = 1:ensemblesize
            θ[i] .= θ[i] .+ p.*((σtb.-σta)./σta).*(θ[i].-θhat);
        end

        # println("Normalized θ after RTPS: $(mean(θ))")

        # analysis parameters back into ensemble members
        for i = 1:ensemblesize
            mparams[i].R1 = θ[i][1]*θs[1];
            mparams[i].R2 = θ[i][2]*θs[2];
            mparams[i].R3 = θ[i][3]*θs[3];
            if mparams[i].R1 <= error.lb[1]
                println("Warning: analysis R1 below lower bound for member $i.
                    Setting to lower bound of $(error.lb[1]).")
                mparams[i].R1 = error.lb[1];
                θ[i][1] = error.lb[1]/θs[1];
            end
            if mparams[i].R2 <= error.lb[2]
                println("Warning: analysis R2 below lower bound for member $i.
                    Setting to lower bound of $(error.lb[2]).")
                mparams[i].R2 = error.lb[2];
                θ[i][2] = error.lb[2]/θs[2];
            end
            if mparams[i].R3 <= error.lb[3]
                println("Warning: analysis R3 below lower bound for member $i.
                    Setting to lower bound of $(error.lb[3]).")
                mparams[i].R3 = error.lb[3];
                θ[i][3] = error.lb[3]/θs[3];
            end
            systems[i].mparams = mparams[i];
        end

        # recalculate and output parameter averages (non-dimensional)
        θhat = mean(θ);
        append!(θout,[θhat])
        # println("Normalized analysis parameter averages: $θhat")

        # analysis parameter variance (for post-processing)
        Ptt = zeros(length(θhat))
        for i = 1:ensemblesize
            Ptt += (θ[i] .- θhat).^2;
        end
        Ptt ./= ensemblesize;
        append!(Pθout,[Ptt])

        # corrected forecast w/ analysis parameters
        for i = 1:ensemblesize
            P1 = sin(2*pi*to[i][end]);
            Q1true = (P1true - Pbtrue)/R1true;
            Q2true = Pbtrue/R2true;
            Q3true = Pbtrue/R3true;
            append!(systems[i].Pb,P1/(mparams[i].R1*(1/mparams[i].R1+1/mparams[i].R2+1/mparams[i].R3)
            append!(systems[i].Q1,(P1-systems[i].Pb[end])/mparams[i].R1)
            append!(systems[i].Q2,systems[i].Pb[end]/mparams[i].R2)
            append!(systems[i].Q3,systems[i].Pb[end]/mparams[i].R3)
        end

        # vector of analysis state vectors
        X = [[systems[i].Pb[end],systems[i].Q1[end],systems[i].Q2[end],systems[i].Q3[end]] for i in 1:ensemblesize];

        # vector of analysis measurements
        Y = [[systems[i].Pb[end],systems[i].Q1[end],systems[i].Q2[end],systems[i].Q3[end]] for i in 1:ensemblesize];

        # for i in 1:ensemblesize
        #     println("Forecast values, member $i: $(X[i][1:end])")
        # end

        # analysis mean state, measurement
        xhat = mean(X);
        yhat = mean(Y);
        # println("Normalized ̂x, second forecast: $xhat")
        # println("Normalized ̂y, second forecast: $yhat")

        # analysis state/measurement, innovation outputs
        append!(yout,[yhat])
        append!(innov,[yhat-y])
        append!(xout,[xhat])

        # state variance (for post-processing)
        Pxx = zeros(length(xhat))
        for i = 1:ensemblesize
            Pxx += (X[i] .- xhat).^2;
        end
        Pxx ./= ensemblesize;
        append!(Pxout,[Pxx])

        # state 2-sigma quantiles (for post-processing)
        for i = 1:length(xhat)
            xq = Float64[];
            for j = 1:ensemblesize
                push!(xq,X[j][i])
            end
            q = quantile(xq,[0.025,0.975]);
            if i == 1
                append!(lbx,[ones(1)*q[1]])
                append!(ubx,[ones(1)*q[2]])
            else
                push!(lbx[end],q[1])
                push!(ubx[end],q[2])
            end
        end

        # output measurement time
        push!(tout,to[1][end])
        for i = 1:ensemblesize
            append!(systems[i].t,to[i][end])
        end

        # set up new measurement time
        if to[1][end] < tf
            tfnew = to[1][end] + systems[1].sparams.dtsav;
            for i = 1:ensemblesize
                append!(to[i][end],tfnew)
            end
        end
    end
    toc()

    # reshape ensemble averages into vectors of individual time series
    for i in 1:length(xout[1]) # i indexes state variables
        xo = Float64[];
        Pxo = Float64[];
        lb = Float64[];
        ub = Float64[];
        for j in 1:length(xout) # j indexes time steps
            push!(xo,xout[j][i])
            push!(Pxo,Pxout[j][i])
            push!(lb,lbx[j][i])
            push!(ub,ubx[j][i])
        end
        append!(xoutv,[xo])
        append!(Pxoutv,[Pxo])
        append!(lbxv,[lb])
        append!(ubxv,[ub])
    end

    for i in 1:length(θout[1]) # i indexes state variables
        θo = Float64[];
        Pθo = Float64[];
        for j in 1:length(θout) # j indexes time steps
            push!(θo,θout[j][i])
            push!(Pθo,Pθout[j][i])
        end
        append!(θoutv,[θo])
        append!(Pθoutv,[Pθo])
    end

    for i in 1:length(yout[1])
        yo = Float64[];
        io = Float64[];
        for j in 1:length(yout)
            push!(yo,yout[j][i])
            push!(io,innov[j][i])
        end
        append!(youtv,[yo])
        append!(ioutv,[io])
    end

    # save for post-processing
    vnames = ["t" "x" "Px" "theta" "Pt" "Pv" "lb" "ub"];
    for i in 1:length(vnames)
        file = MAT.matopen("$(vnames[i]).mat","w");
        if vnames[i] == "t"
            write(file,"t",tout)
        elseif vnames[i] == "x"
            write(file,"x",xoutv)
        elseif vnames[i] == "Px"
            write(file,"Px",Pxoutv)
        elseif vnames[i] == "theta"
            write(file,"theta",θoutv)
        elseif vnames[i] == "Pt"
            write(file,"Pt",Pθoutv)
        elseif vnames[i] == "Pv"
            write(file,"Pv",Pvout)
        elseif vnames[i] == "lb"
            write(file,"lb",lbxv)
        elseif vnames[i] == "ub"
            write(file,"ub",ubxv)
        end
        close(file)
    end

    # output
    return systems,tout,xoutv,youtv,ioutv,θoutv,Pθoutv,Pxoutv,lbxv,ubxv
end
