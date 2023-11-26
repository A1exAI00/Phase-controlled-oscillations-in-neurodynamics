#=
TODO
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
PLOT_FILENAME = "08_bif_diagram"
DEFAULT_PX_PER_UNIT_PNG = 2

#########################################################################################

# Параметры опорной системы
a, ϵ, I = 0.01, 0.02, 0.01

# Начальные условия опорного сигнала
u₀, v₀ = 0.14996701966490192, 0.021055242172979355

# Время интегрирования
t₀, t₁ = 0, 204.49914703598859

# Параметры спайка
spike_duration_period = 0.01
periods_before_spike = 1
Aₛₚ_start, Aₛₚ_end, Aₛₚ_N = -1.0, 3.5, 150

N_starting_φ = 10

# Количество итераций в построении бифуркационной диаграммы
N_iterations = 1000

# Количество последних точек, которые будут отмечены на бифуркационной диаграмме
N_show_last_iterations = 100 # 100

#########################################################################################

reference_param = SA[a, ϵ, I, 0, 0, 0]
U₀ = SA[u₀, v₀]
t_SPAN = [t₀, t₁]
Aₛₚ_span = range(Aₛₚ_start, Aₛₚ_end, Aₛₚ_N)
starting_φ_span = range(0, 2π, N_starting_φ)

#########################################################################################

reference_sol = neurodynamics_integrate(reference_param, U₀, t_SPAN)
reference_u_sol = reference_sol[1,:]
reference_t_sol = reference_sol.t

T = mesure_T(reference_u_sol, reference_t_sol)
println("T=$(T)")

#########################################################################################

tₛₚ, τₛₚ = periods_before_spike*T, spike_duration_period*T

#########################################################################################


φ_end = Vector{Float64}()
A = Vector{Float64}()
for (i,Aₛₚ) in enumerate(Aₛₚ_span)
    global φ_end, A
    param = SA[a, ϵ, I, Aₛₚ, tₛₚ, τₛₚ]
    curr_φ_end = Vector{Float64}()
    for (j,φ) in enumerate(starting_φ_span) 
        U₀_curr = reference_sol(T-φ/(2π)*T)
        for k in 1:N_iterations
            curr_sol = neurodynamics_integrate(param, U₀_curr, t_SPAN)
            U₀_curr = curr_sol[:,end]  # curr_sol(t₁)
            if k > (N_iterations-N_show_last_iterations)
                φ = φ_time_series(curr_sol[1,:], curr_sol.t, reference_u_sol, reference_t_sol, true)
                push!(curr_φ_end, φ)
            end
        end
    end
    
    append!(φ_end, curr_φ_end)
    append!(A, Aₛₚ*ones(length(curr_φ_end)))

    progress("$i/$Aₛₚ_N "*get_easy_time(), true)
end

#########################################################################################

fig = Figure(resolution=DEFAULT_PLOT_RES)

ax = Axis(fig[1, 1], 
    title="bif diagram : tₛₚ=$(tₛₚ), τₛₚ=$(spike_duration_period)*T, T=$(T)", 
    xlabel="Aₛₚ", ylabel="φ_end")

scatter!(ax, A, φ_end, markersize=2, color=:black) # 0.5

savingpath = joinpath(DEFAULT_SAVING_PATH, PLOT_FILENAME*"$(time_ns()).png")
save(savingpath, fig, px_per_unit=DEFAULT_PX_PER_UNIT_PNG)