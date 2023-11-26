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
PLOT_FILENAME = "06_phase_space_mapping"
DEFAULT_PX_PER_UNIT_PNG = 2

#########################################################################################

function step_points(points, t_SPAN_frame)
    global param
    for j = 1:n_points
        U₀ = SVector{2}(points[j,:])
        sol = neurodynamics_integrate(param, U₀, t_SPAN_frame)
        points[j,:] = sol(t_SPAN_frame[2])
    end
    return points
end

#########################################################################################

# Параметры опорной системы
a, ϵ, I = 0.01, 0.02, 0.01

# Начальные условия опорного сигнала
u₀, v₀ = 0.14996701966490192, 0.021055242172979355

# Количество кадров в gif-файле
n_frames = 1000

# Количество изображающех точек в фазовом пространстве
n_points = 100

# Параметры спайка
spike_duration_period = 0.4
periods_before_spike = 3
Aₛₚ = 0.5

# Время интегрирования
t₀, t₁ = 0, 500

# Границы фазового пространства
x_min, x_max = -0.6, 1.2
y_min, y_max = -0.2, 0.35 

#########################################################################################

reference_param = SA[a, ϵ, I, 0, 0, 0]
U₀ = SA[u₀, v₀]
Δt_frame = t₁/n_frames
t_SPAN = [t₀, t₁]
t_SPAN_frame = [t₀, Δt_frame]

x_slow_movement = range(x_min, x_max, 100)
f_slow_movement(x) = x*(x-a)*(1-x)
y_slow_movement = f_slow_movement.(x_slow_movement)

#########################################################################################

# Интегрирование опорной системы
reference_sol = neurodynamics_integrate(reference_param, U₀, t_SPAN)
reference_u_sol = reference_sol[1,:]
reference_v_sol = reference_sol[2,:]
reference_t_sol = reference_sol.t

T = mesure_T(reference_u_sol, reference_t_sol)
println("T=$(T)")

#########################################################################################

Δt = T/n_points
tₛₚ, τₛₚ = periods_before_spike*T, spike_duration_period*T
param = (a, ϵ, I, Aₛₚ, tₛₚ, τₛₚ)

#########################################################################################

# Заполнить начальные условия для колебаний
points = zeros(n_points,2)
for i = 1:n_points
    points[i,:] = reference_sol(i*Δt)
end

#########################################################################################

fig = Figure(resolution=DEFAULT_PLOT_RES)

ax = Axis(fig[1, 1], 
    title="Mapping : a=$(a), ϵ=$(ϵ) I=$(I); T=$(T), τₛₚ=$(spike_duration_period)*T, Aₛₚ=$(Aₛₚ)", 
    xlabel="Aₛₚ", ylabel="φ*")

limits!(ax, x_min, x_max, y_min, y_max)


savingpath = joinpath(DEFAULT_SAVING_PATH, PLOT_FILENAME*"$(time_ns()).gif")
record(fig, savingpath, 1:n_frames; framerate=30) do i
    global t_SPAN_frame, points
    progress("frame #$(i)/$(n_frames)", true)
    
	empty!(ax)
    
    # integrate points
    points = step_points(points, t_SPAN_frame)
    
    # slow movement line
    lines!(ax, x_slow_movement, y_slow_movement, color=:red)
    
    # phase points
    x_points = points[:,1]
    y_points = points[:,2]
    scatter!(ax, 
        x_points, y_points, 
        color=:blue)
    
    t_SPAN_frame += Δt_frame*[1,1]
end