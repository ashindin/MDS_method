include("../01_perturbation_model/perturbation_model.jl")

function sqrt_epsilon_fun(ω, z, t, z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)
    N_model = N_model_fun(z, t, z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)        
    ωₚ = 2pi*Ne_to_fpe(N_model)
    if ωₚ > ω 
        return NaN
    end
    return sqrt(1 - ωₚ^2/ω^2)
end

function sqrt_eps_integral(z_step, ω, t, z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)
    z = z1
    S = sqrt_epsilon_fun(ω, z, t, z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)/2    
    while true        
        z = z + z_step
        s = sqrt_epsilon_fun(ω, z, t, z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)
        if isnan(s)
            S = S - sqrt_epsilon_fun(ω, z-z_step, t, z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)/2
            break
        end
        S = S + s
    end
    return S*z_step    
end

# sqrt_eps_integral(1000., 45e5, 5., z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)

function phase_variations(ω, t1, t2, z_step, z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)
    c = 299792458.
    phase2 = 2ω/c * sqrt_eps_integral(z_step, ω, t2, z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)
    phase1 = 2ω/c * sqrt_eps_integral(z_step, ω, t1, z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)
    return phase2 - phase1
end

# f_axe = 4e6:1000.:5e6
# ω_axe = 2pi*f_axe
# t_axe = 0.0:0.1:5.
# z_step = 10.
# DPhi = zeros(length(t_axe),length(f_axe))

# for i=2:length(t_axe)
#     println(i)   
#     for j=1:length(f_axe)        
#         # println(i, " ",j)        
#         DPhi[i,j] = phase_variations(ω_axe[j], t_axe[i-1], t_axe[i], z_step, z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)
#     end
# end

# using Plots


# i=12; plot(f_axe/1e6, DPhi[i,:])
# title(string("t=",t_axe[i]))
# ylims((0,15))