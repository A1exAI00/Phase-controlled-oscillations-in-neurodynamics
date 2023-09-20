using StaticArrays
using CairoMakie

include("../src/neurodynamics_integration.jl")
include("../src/misc_tools.jl")
include("../src/time_series_tools.jl")


########################################################################
########################################################################
# General program settings

DEFAULT_PLOT_RES = (1000, 800)
DEFAULT_SAVING_PATH = "../generated"
PLOT_FILENAME = "phase_sync"
DEFAULT_PX_PER_UNIT_PNG = 2


########################################################################
########################################################################
# System and spike parameters for reference oscillation

a, ϵ, I = 0.01, 0.02, 0.01
reference_param = SA[a, ϵ, I, 0, 0, 0]


########################################################################
########################################################################
# Initial values and time span for simulation

u₀, v₀ = 0.14996701966490192, 0.021055242172979355
U₀ = SA[u₀, v₀]

t₀, t₁ = 0, 1200
t_SPAN = [t₀, t₁]




########################################################################
########################################################################
# reference oscillation integration

reference_sol = neurodynamics_integrate(reference_param, U₀, t_SPAN)
reference_u_sol = reference_sol[1,:]
reference_v_sol = reference_sol[2,:]
reference_t_sol = reference_sol.t

T = mesure_T(reference_u_sol, reference_t_sol)
println("T=$(T)")


########################################################################
########################################################################
# Spike parameters for test oscillations

# numper of test oscillations
N_oscillations = 100
Δt = T/N_oscillations

spike_duration_period = 0.4
periods_before_spike = 5
Aₛₚ, tₛₚ, τₛₚ = -1.1, periods_before_spike*T, spike_duration_period*T # 2, periods_before_spike*T, spike_duration_period*T
param = SA[a, ϵ, I, Aₛₚ, tₛₚ, τₛₚ]
f_spike(t) = if (t<tₛₚ) || (t>tₛₚ+τₛₚ) 0 else Aₛₚ end


########################################################################
########################################################################


fig = Figure(resolution=DEFAULT_PLOT_RES)

ax1 = Axis(fig[1, 1], 
    title="spike : Aₛₚ=$(Aₛₚ), tₛₚ=$(tₛₚ), τₛₚ=$(spike_duration_period)*T, T=$(T)", 
    xlabel="t", ylabel="u")
ax2 = Axis(fig[2, 1], 
    title="reference oscillation : a=$(a), ϵ=$(ϵ) I=$(I)", 
    xlabel="t", ylabel="u")

lines!(ax1, reference_t_sol, f_spike.(reference_t_sol))

#lines!(ax2, reference_t_sol, reference_u_sol, linewidth=4, color=:red)
for i = 1:N_oscillations
    U₀_curr = reference_sol(Δt*i)
    curr_sol = neurodynamics_integrate(param, U₀_curr, t_SPAN)
    curr_u_sol = curr_sol[1,:]
    curr_t_sol = curr_sol.t
    lines!(ax2, curr_t_sol, curr_u_sol)
end

savingpath = joinpath(DEFAULT_SAVING_PATH, PLOT_FILENAME*"$(time_ns()).png")
save(savingpath, fig, px_per_unit=DEFAULT_PX_PER_UNIT_PNG)
