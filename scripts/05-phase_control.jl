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
PLOT_FILENAME = "phase_control"
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

t₀, t₁ = 0, 2000
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

spike_duration_period = 0.4
periods_before_spike = 5

tₛₚ, τₛₚ = periods_before_spike*T, spike_duration_period*T
Aₛₚ_start, Aₛₚ_end, Aₛₚ_N = -3, 3, 2000
Aₛₚ_span = range(Aₛₚ_start, Aₛₚ_end, Aₛₚ_N)

reference_t_maxes = calc_maxes(reference_u_sol, reference_t_sol)
reference_t_first_max = first_greater_then(reference_t_maxes, tₛₚ)



########################################################################
########################################################################


φ_array = zeros(Aₛₚ_N)
for (i, Aₛₚ) in enumerate(Aₛₚ_span)
    global φ_array
    param = SA[a, ϵ, I, Aₛₚ, tₛₚ, τₛₚ]
    
    println("Aₛₚ=$(Aₛₚ)")
    
    sol = neurodynamics_integrate(param, U₀, t_SPAN)
    u_sol = sol[1,:]
    t_sol = sol.t
    
    t_maxes = calc_maxes(u_sol, t_sol)
    t_first_max = first_greater_then(t_maxes, (tₛₚ+τₛₚ)+10)
    
    φ = mod(2*π*(t_first_max-reference_t_first_max)/T, 2*π)
    #φ = (t_first_max-reference_t_first_max)/T
    
    φ_array[i] = φ
end


########################################################################
########################################################################


fig = Figure(resolution=DEFAULT_PLOT_RES)

ax = Axis(fig[1, 1], 
    title="Control of phase : a=$(a), ϵ=$(ϵ) I=$(I); T=$(T), τₛₚ=$(τₛₚ)", 
    xlabel="Aₛₚ", ylabel="φ*")

lines!(ax, Aₛₚ_span, φ_array)

savingpath = joinpath(DEFAULT_SAVING_PATH, PLOT_FILENAME*"$(time_ns()).png")
save(savingpath, fig, px_per_unit=DEFAULT_PX_PER_UNIT_PNG)
