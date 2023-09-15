using StaticArrays
using CairoMakie

include("../src/neurodynamics_integration.jl")
using .neurodynamics_integration


########################################################################
########################################################################


DEFAULT_PLOT_RES = (1000, 800)
DEFAULT_SAVING_PATH = "../generated"
PLOT_FILENAME = "phase_spase"
DEFAULT_PX_PER_UNIT_PNG = 2


########################################################################
########################################################################


a, ϵ, I = 0.01, 0.02, 0.01
Aₛₚ, tₛₚ, τₛₚ = 0, 0, 0
param = SA[a, ϵ, I, Aₛₚ, tₛₚ, τₛₚ]

# reference oscillation initial val
u₀, v₀ = 0.14996701966490192, 0.021055242172979355
U₀ = SA[u₀, v₀]

# time span
t₀, t₁ = 0, 1000
t_SPAN = [t₀, t₁]


########################################################################
########################################################################

# streamplot
f(u,p,t) = Point2f(neurodynamics_system(u,p,t))
f₀(u,p) = f(u,p,0)
f₀ₚ(u) = f₀(u,param)


########################################################################
########################################################################


# limit cycle
sol = neurodynamics_integrate(param, U₀, t_SPAN)
u_sol = sol[1,:]
v_sol = sol[2,:]

########################################################################
########################################################################


fig = Figure(resolution=DEFAULT_PLOT_RES)
ax = Axis(fig[1, 1], 
    title="Phase space : a=$(a), ϵ=$(ϵ) I=$(I)", 
    xlabel="U", ylabel="V")
    
hlines!(ax, 0, color=:black)
vlines!(ax, 0, color=:black)
scatter!(ax, I, I*(I-a)*(1-I))
lines!(ax, u_sol, v_sol, color=:red, linewidth=3) # linewidth=1.5 default
streamplot!(ax, f₀ₚ, -0.5..1.1, -0.2..0.35, colormap=:magma)

savingpath = joinpath(DEFAULT_SAVING_PATH, PLOT_FILENAME*"$(time_ns()).png")
save(savingpath, fig, px_per_unit=DEFAULT_PX_PER_UNIT_PNG)
