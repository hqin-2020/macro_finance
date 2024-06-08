ENV["JULIA_PKG_OFFLINE"] = "true"
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Optim
using Roots
using NPZ
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--Delta"
            help = "time step"
            arg_type = Float64
            default = 1.  
        "--gamma"
            help = "risk aversion"
            arg_type = Float64
            default = 8.0
        "--rho"
            help = "inverse IES"
            arg_type = Float64
            default = 1.0   
        "--alpha"
            help = "productivity"
            arg_type = Float64
            default = 0.1844
        "--kappa"
            help = "capital elasticity of substitution"
            arg_type = Float64
            default = 0.0
        "--zeta"
            help = "capital weight"
            arg_type = Float64
            default = 0.5  
        "--beta1"
            help = "capital 1 persistence"
            arg_type = Float64
            default = 0.01
        "--beta2"
            help = "capital 2 persistence"
            arg_type = Float64
            default = 0.01
        "--action_name"
            help = "action name"
            arg_type = String
            default = "publish"
    end
    return parse_args(s)
end

@show parsed_args = parse_commandline()
Delta                = parsed_args["Delta"]
gamma                = parsed_args["gamma"]
rho                  = parsed_args["rho"]
alpha                = parsed_args["alpha"]
kappa                = parsed_args["kappa"]
zeta                 = parsed_args["zeta"]
beta1                = parsed_args["beta1"]
beta2                = parsed_args["beta2"]
action_name          = parsed_args["action_name"]
include("utils_twocapitals_two_dimensions.jl")

outputdir = "./output/"*action_name*"/Delta_"*string(Delta)*"/beta1_"*string(beta1)*"_beta2_"*string(beta2)*"/kappa_"*string(kappa)*"_zeta_"*string(zeta)*"/gamma_"*string(gamma)*"_rho_"*string(rho)*"_alpha_"*string(alpha)*"/"
isdir(outputdir) || mkpath(outputdir)

## Calibration
delta = 0.01
betaz = 0.056
eta1 = -0.04
eta2 = -0.04
phi1 = 8.0
phi2 = 8.0
smean = 6.30e-06;
sigma_k1 = [sqrt(2)*0.92, .0, .4] * sqrt(smean) * sqrt(12)
sigma_k2 = [0, sqrt(2)*0.92, .4] * sqrt(smean) * sqrt(12)
sigma_z =  [0, 0, 5.7] * sqrt(smean) * sqrt(12)

## Construct state space grid
II, JJ = trunc(Int,1001), trunc(Int,201);
if kappa == 0.0
    rmax = 6.0;
else 
    rmax = 1.0;
end
rmin = -rmax;
zmax = 1.0;
zmin = -zmax;

maxit = 200000;     # maximum number of iterations in the HJB loop
crit  = 10e-6;      # criterion HJB loop

# Initialize model objects 
baseline1 = Baseline(zeta, kappa, betaz, sigma_z, beta1, sigma_k1, delta);
baseline2 = Baseline(zeta, kappa, betaz, sigma_z, beta2, sigma_k2, delta);
technology1 = Technology(alpha, phi1, eta1);
technology2 = Technology(alpha, phi2, eta2);
model = TwoCapitalEconomy(baseline1, baseline2, technology1, technology2);

## Initialize grid and FDM
grid = Grid_rz(rmin, rmax, II, zmin, zmax, JJ);
params = FinDiffMethod(maxit, crit, Delta);

## Initialize value function and make a guess for consumption-capital ratio
if rho == 1.0
    preloadV0 = -3*ones(grid.I, grid.J)
    preloadcons = 0.03*ones(grid.I, grid.J)
else
    preload_rho = 1.0
    if kappa == 0.0
        if beta1 == 0.0
            preload_Delta = 0.1       
        elseif beta1 == 0.04
            preload_Delta = 1.0
        end
    else
        preload_Delta = 1.0
    end
    preload_alpha = 0.1844
    preloaddir = "./output/"*action_name*"/Delta_"*string(preload_Delta)*"/beta1_"*string(beta1)*"_beta2_"*string(beta2)*"/kappa_"*string(kappa)*"_zeta_"*string(zeta)*"/gamma_"*string(gamma)*"_rho_"*string(preload_rho)*"_alpha_"*string(preload_alpha)*"/"
    preload = npzread(preloaddir*"res.npz")
    println("preload location : "*preloaddir)
    preloadV0 = preload["V"]
    preloadcons = preload["cons"]
end

## Solve the model
times = @elapsed begin
A, V, val, d1_F, d2_F, d1_B, d2_B, h1_F, h2_F, hz_F, h1_B, h2_B, hz_B,
        mu_1_F, mu_1_B, mu_r_F, mu_r_B, mu_z_F, mu_z_B, V0, Vr, Vr_F, Vr_B, Vz_B, Vz_F, cF, cB, Vz, rr, zz, dr, dz =
        value_function_twocapitals(gamma, rho, model, grid, params, preloadV0, preloadcons, beta1);
end
println("Convegence time (minutes): ", times/60)
g = stationary_distribution(A, grid)

## Control variables
mu_1 = (mu_1_F + mu_1_B)/2.;
mu_r = (mu_r_F + mu_r_B)/2.;
mu_z = (mu_z_F + mu_z_B)/2.;
h1 = (h1_F + h1_B)/2.;
h2 = (h2_F + h2_B)/2.;
hz = (hz_F + hz_B)/2.;
d1 = (d1_F + d1_B)/2;
d2 = (d2_F + d2_B)/2;

r = range(rmin, stop=rmax, length=II);    # capital ratio vector
rr = r * ones(1, JJ);
IJ = II*JJ;
k1a = zeros(II,JJ)
k2a = zeros(II,JJ)
for i=1:IJ
    p = rr[i];
    k1a[i] = (1-zeta + zeta*exp.(p*(1-kappa))).^(1/(kappa-1));
    k2a[i] = ((1-zeta)*exp.(p*((kappa-1))) + zeta).^(1/(kappa-1));
end
d1k = d1.*k1a
d2k = d2.*k2a
c = alpha*ones(II,JJ) - d1k - d2k

## Save results
results = Dict(
# Parameters    
"delta" => delta, "betaz"=> betaz,  "beta1" => beta1, "beta2" => beta2, "eta1" => eta1, "eta2" => eta2, 
"sigma_k1" => sigma_k1, "sigma_k2" => sigma_k2, "sigma_z" =>  sigma_z, 
"alpha" => alpha, "kappa" => kappa,"zeta" => zeta, "phi1" => phi1, "phi2" => phi2,  
"gamma" => gamma, "rho" => rho,
# Grid
"I" => II, "J" => JJ, "rmax" => rmax, "rmin" => rmin, "zmax" => zmax, "zmin" => zmin, "rr" => rr, "zz" => zz,
# Iteration
"maxit" => maxit, "crit" => crit, "Delta" => Delta,
# Results
"cons" => c, "h1" => h1, "h2" => h2, "hz" => hz,"d1" => d1, "d2" => d2,"k1a" => k1a, "k2a"=> k2a,
"V0" => V0, "V" => V, "Vr" => Vr, "Vz" => Vz, "val" => val,"dr" => dr, "dz" => dz,
"mu_1" => mu_1, "mu_r" => mu_r, "mu_z" => mu_z,
"g" => g)
npzwrite(outputdir*"res.npz", results)
