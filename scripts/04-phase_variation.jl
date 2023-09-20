using StaticArrays
using CairoMakie
using Statistics: std

include("../src/neurodynamics_integration.jl")
include("../src/misc_tools.jl")
include("../src/time_series_tools.jl")


########################################################################
########################################################################
# General program settings

DEFAULT_PLOT_RES = (1000, 800)
DEFAULT_SAVING_PATH = "../generated"
PLOT_FILENAME = "phase_deviation"
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

t₀, t₁ = 0, 1500
t_SPAN = [t₀, t₁]

N_oscillations = 100
Δt = estimated_T/N_oscillations


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
Aₛₚ_start, Aₛₚ_end, Aₛₚ_N = -3, 3, 500 # -3, 3, 500
Aₛₚ_span = range(Aₛₚ_start, Aₛₚ_end, Aₛₚ_N)

#f_spike(t) = if (t<tₛₚ) || (t>tₛₚ+τₛₚ) 0 else Aₛₚ end

reference_t_maxes = calc_maxes(reference_u_sol, reference_t_sol)
reference_t_first_max = first_greater_then(reference_t_maxes, tₛₚ)


########################################################################
########################################################################


φ_deviation = zeros(Aₛₚ_N)
for (i, Aₛₚ) in enumerate(Aₛₚ_span)
    global φ_deviation
    param = SA[a, ϵ, I, Aₛₚ, tₛₚ, τₛₚ]
    
    progress("Aₛₚ=$(Aₛₚ)", true)
    
    φ_arr = zeros(N_oscillations)
    for j = 1:N_oscillations
        global U₀
        U₀ = reference_sol(j*Δt)
        sol = neurodynamics_integrate(param, U₀, t_SPAN)
        u_sol = sol[1,:]
        t_sol = sol.t
        
        t_maxes = calc_maxes(u_sol, t_sol)
        t_first_max = first_greater_then(t_maxes, (tₛₚ+τₛₚ))

        # TODO: use the time_series_tools.jl function for φ

        #φ = mod(2*π*(t_first_max-reference_t_first_max)/T, 2*π)
        φ = (t_first_max-reference_t_first_max)/T
        φ_arr[j] = φ
    end
    φ_deviation[i] = std(φ_arr)
end


########################################################################
########################################################################


fig = Figure(resolution=DEFAULT_PLOT_RES)

ax = Axis(fig[1, 1], 
    title="Standart deviation of phase : a=$(a), ϵ=$(ϵ) I=$(I); T=$(T), τₛₚ=$(spike_duration_period)*T", 
    xlabel="Aₛₚ", ylabel="δφ")

lines!(ax, Aₛₚ_span, φ_deviation)

savingpath = joinpath(DEFAULT_SAVING_PATH, PLOT_FILENAME*"$(time_ns()).png")
save(savingpath, fig, px_per_unit=DEFAULT_PX_PER_UNIT_PNG)
