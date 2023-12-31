#=
Фазовая перестановка под действием спайка.
После прекращения спайка система приходит к предельному циклу, 
но будет фазовое смещение от опорного сигнала.

Скрипт генерирует график осцилляций 2 нейронов: один - опорный, 
второй - на который воздействует спайк.
=#

#########################################################################################

using StaticArrays
using CairoMakie

include("../src/neurodynamics_integration.jl")
include("../src/misc_tools.jl")
include("../src/time_series_tools.jl")

#########################################################################################

DEFAULT_PLOT_RES = (1000, 800)
DEFAULT_SAVING_PATH = "generated"
PLOT_FILENAME = "02_phase_reset"
DEFAULT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Параметры системы
a, ϵ, I = 0.01, 0.02, 0.01

# Начальные условия опорного сигнала
u₀, v₀ = 0.14996701966490192, 0.021055242172979355

# Время интегрирования
t₀, t₁ = 0, 200

#########################################################################################

reference_param = SA[a, ϵ, I, 0, 0, 0]
U₀ = SA[u₀, v₀]
t_SPAN = [t₀, t₁]

#########################################################################################

# Интегрирование опорного сигнала
reference_sol = neurodynamics_integrate(reference_param, U₀, t_SPAN)
reference_u_sol = reference_sol[1,:]
reference_v_sol = reference_sol[2,:]
reference_t_sol = reference_sol.t

T = mesure_T(reference_u_sol, reference_t_sol)
println("T=$(T)")

#########################################################################################

# Параметры спайка
spike_duration_period = 0.01
periods_before_spike = 1
Aₛₚ, tₛₚ, τₛₚ = 1, periods_before_spike*T, spike_duration_period*T
f_spike(t) = if (t<tₛₚ) || (t>tₛₚ+τₛₚ) 0 else Aₛₚ end

param = SA[a, ϵ, I, Aₛₚ, tₛₚ, τₛₚ]

sol = neurodynamics_integrate(param, U₀, t_SPAN)
u_sol = sol[1,:]
v_sol = sol[2,:]
t_sol = sol.t

debug_display(sol, false)

#########################################################################################

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
