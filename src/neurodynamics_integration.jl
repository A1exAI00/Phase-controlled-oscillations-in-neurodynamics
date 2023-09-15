module neurodynamics_integration
    using Statistics: mean
    using OrdinaryDiffEq, StaticArrays
    
########################################################################
########################################################################

    
    const RELTOL::Float64, ABSTOL::Float64 = 1e-5, 1e-5
    const MAXITERS::Int64 = Int(1e7)
    const ALG = Rodas5P
    const X_BOUNDS, Y_BOUNDS, Z_BOUNDS = (-3.0,3.0), (-3.0,3.0), (0.0,1e10)


########################################################################
########################################################################

    function neurodynamics_system(u,p,t)
        #z, γ, I₀ = p
        a, ϵ, I, Aₛₚ, tₛₚ, τₛₚ = p
        Iₜ(t) = if (t<tₛₚ) || (t>tₛₚ+τₛₚ) 0 else Aₛₚ end
        return SA[u[1]*(u[1]-a)*(1-u[1]) - u[2], 
                ϵ*(u[1]-I-Iₜ(t))]
    end
    
    function neurodynamics_system_inplace(du, u,p,t)
        #z, γ, I₀ = p
        a, ϵ, I, Aₛₚ, tₛₚ, τₛₚ = p
        Iₜ(t) = if (t<tₛₚ) || (t>tₛₚ+τₛₚ) 0 else Aₛₚ end
        du = SA[u[1]*(u[1]-a)*(1-u[1]) - u[2], 
                ϵ*(u[1]-I-Iₜ(t))]
    end


########################################################################
########################################################################

    function neurodynamics_integrate(param, U₀, t_span, check_success=true)
        prob = ODEProblem(neurodynamics_system, U₀, t_span, param)
        sol = solve(prob, ALG(), reltol=RELTOL, abstol=ABSTOL, maxiters=MAXITERS)
        
        if (check_success && sol.retcode!=:Success) 
            return NaN
        end
        
        return sol
    end
    
########################################################################
########################################################################
    
    export neurodynamics_integrate, neurodynamics_system
end
