using StaticArrays
using CairoMakie

include("../src/neurodynamics_integration.jl")
include("../src/misc_tools.jl")


########################################################################
########################################################################
# General program variables

DEFAULT_PLOT_RES = (1000, 800)
DEFAULT_SAVING_PATH = "../generated"
PLOT_FILENAME = "phase_reset"
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

t₀, t₁ = 0, 1000
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
# Spike parameters for test oscillation

spike_duration_period = 4
periods_before_spike = 5
Aₛₚ, tₛₚ, τₛₚ = 1, periods_before_spike*T, spike_duration_period*T
param = SA[a, ϵ, I, Aₛₚ, tₛₚ, τₛₚ]
f_spike(t) = if (t<tₛₚ) || (t>tₛₚ+τₛₚ) 0 else Aₛₚ end


sol = neurodynamics_integrate(param, U₀, t_SPAN)
u_sol = sol[1,:]
v_sol = sol[2,:]
t_sol = sol.t

debug_display(sol, false)


########################################################################
########################################################################


fig = Figure(resolution=DEFAULT_PLOT_RES)

ax1 = Axis(fig[1, 1], 
    title="spike : Aₛₚ=$(Aₛₚ), tₛₚ=$(tₛₚ), τₛₚ=$(spike_duration_period)*T, T=$(T)", 
    xlabel="t", ylabel="u")
ax2 = Axis(fig[2, 1], 
    title="reference oscillation : a=$(a), ϵ=$(ϵ) I=$(I)", 
    xlabel="t", ylabel="u")

lines!(ax1, t_sol, f_spike.(t_sol))

lines!(ax2, reference_t_sol, reference_u_sol, linewidth=4, color=:red)
lines!(ax2, t_sol, u_sol)

savingpath = joinpath(DEFAULT_SAVING_PATH, PLOT_FILENAME*"$(time_ns()).png")
save(savingpath, fig, px_per_unit=DEFAULT_PX_PER_UNIT_PNG)
