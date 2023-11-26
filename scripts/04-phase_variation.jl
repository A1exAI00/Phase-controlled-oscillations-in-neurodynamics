#=
В зависимости от амплитуды спайка можно достигнуть разной степени синхронизации.

Скрипт генерирует график СКО фазы после фазовой переустановки в зависимости от амплитуды спайка.

TODO: Есть баг, что программа не понимает после какого времени считать систему "фазово переустановленной". 
Поэтому для некоторых значений амплитуды значения СКО не верные.
=#

#########################################################################################

using StaticArrays
using CairoMakie
using Statistics: std

include("../src/neurodynamics_integration.jl")
include("../src/misc_tools.jl")
include("../src/time_series_tools.jl")


#########################################################################################

DEFAULT_PLOT_RES = (1000, 800)
DEFAULT_SAVING_PATH = "generated"
PLOT_FILENAME = "04_phase_deviation"
DEFAULT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Параметры опорной системы
a, ϵ, I = 0.01, 0.02, 0.01

# Начальные условия опорного сигнала
u₀, v₀ = 0.14996701966490192, 0.021055242172979355

# Время интегрирования
t₀, t₁ = 0, 1500

# Количество нейронов, по которым будет бить спайк
N_oscillations = 100

# Параметры спайка
spike_duration_period = 0.4
periods_before_spike = 5
Aₛₚ_start, Aₛₚ_end, Aₛₚ_N = -3, 3, 500

#########################################################################################

reference_param = SA[a, ϵ, I, 0, 0, 0]
U₀ = SA[u₀, v₀]
t_SPAN = [t₀, t₁]
Aₛₚ_span = range(Aₛₚ_start, Aₛₚ_end, Aₛₚ_N)

#########################################################################################

# Интегрирование опорной системы
reference_sol = neurodynamics_integrate(reference_param, U₀, t_SPAN)
reference_u_sol = reference_sol[1,:]
reference_v_sol = reference_sol[2,:]
reference_t_sol = reference_sol.t

T = mesure_T(reference_u_sol, reference_t_sol)
println("T=$(T)")

#########################################################################################

Δt = T/N_oscillations
tₛₚ, τₛₚ = periods_before_spike*T, spike_duration_period*T

reference_t_maxes = times_of_max(reference_u_sol, reference_t_sol)
reference_t_first_max = first_greater_then(reference_t_maxes, tₛₚ)

#########################################################################################

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
        
        t_maxes = times_of_max(u_sol, t_sol)
        t_first_max = first_greater_then(t_maxes, (tₛₚ+τₛₚ))

        # TODO: use the time_series_tools.jl function for φ

        #φ = mod(2*π*(t_first_max-reference_t_first_max)/T, 2*π)
        φ = (t_first_max-reference_t_first_max)/T
        φ_arr[j] = φ
    end
    φ_deviation[i] = std(φ_arr)
end

#########################################################################################

fig = Figure(resolution=DEFAULT_PLOT_RES)

ax = Axis(fig[1, 1], 
    title="Standart deviation of phase : a=$(a), ϵ=$(ϵ) I=$(I); T=$(T), τₛₚ=$(spike_duration_period)*T", 
    xlabel="Aₛₚ", ylabel="δφ")

lines!(ax, Aₛₚ_span, φ_deviation)

savingpath = joinpath(DEFAULT_SAVING_PATH, PLOT_FILENAME*"$(time_ns()).png")
save(savingpath, fig, px_per_unit=DEFAULT_PX_PER_UNIT_PNG)
