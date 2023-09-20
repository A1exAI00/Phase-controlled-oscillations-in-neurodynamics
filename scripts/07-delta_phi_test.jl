using StaticArrays
using CairoMakie

include("../src/neurodynamics_integration.jl")
include("../src/misc_tools.jl")
include("../src/time_series_tools.jl")


########################################################################
########################################################################

DEFAULT_PLOT_RES = (1000, 800)
DEFAULT_SAVING_PATH = "generated"
PLOT_FILENAME = "delta_phi_test"
DEFAULT_PX_PER_UNIT_PNG = 2

########################################################################
########################################################################


a, ϵ, I = 0.01, 0.02, 0.01
reference_param = SA[a, ϵ, I, 0, 0, 0]
param = SA[a, ϵ, I, 0, 0, 0]

U₀ = SA[0.14996701966490192, 0.021055242172979355]

t₀, t₁ = 0, 1000
t_SPAN = [t₀, t₁]


########################################################################
########################################################################


reference_sol = neurodynamics_integrate(reference_param, U₀, t_SPAN)
reference_u_sol = reference_sol[1,:]
#reference_v_sol = reference_sol[2,:]
reference_t_sol = reference_sol.t

T = mesure_T(reference_u_sol, reference_t_sol)
println("T=$(T)")


########################################################################
########################################################################


N_test = 50
Δt = T/N_test

φ = zeros(N_test)
t_start = [Δt*i for i = 1:N_test]

for i = 1:N_test
    U₀_curr = reference_sol(T-t_start[i])
    curr_sol = neurodynamics_integrate(param, U₀_curr, t_SPAN)
    curr_u_sol = curr_sol[1,:]
    curr_t_sol = curr_sol.t
    φ[i] = φ_time_series(curr_u_sol, curr_t_sol, reference_u_sol, reference_t_sol, true)
end


########################################################################
########################################################################


fig = Figure(resolution=DEFAULT_PLOT_RES)

ax = Axis(fig[1, 1], 
    title="Δφ with no spike : Aₛₚ=$(0), T=$(T)", 
    xlabel="φ_start", ylabel="φ_end")

lines!(ax, 2π*t_start/T, φ)

savingpath = joinpath(DEFAULT_SAVING_PATH, PLOT_FILENAME*"$(time_ns()).png")
save(savingpath, fig, px_per_unit=DEFAULT_PX_PER_UNIT_PNG)
