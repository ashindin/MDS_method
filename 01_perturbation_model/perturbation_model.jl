function Ne_to_fpe(Ne)
    # convert electron concentration to electron plasma frequency
    return 8978.18235285877*sqrt(abs(Ne)) # модуль нужен, чтобы небыло предупреждений при построении графиков из-за взятия корня от -0.
end

function fpe_to_Ne(fpe)
# convert electron plasma frequency to electron concentration
    return 1.24057537429491*fpe^2*1e-8
end

function heaviside(t)
    if t==0
        return 0
    end
    return 0.5 * (sign(t) + 1)
 end

function N0_fun(z, Nmax, z1, z2, zmax)
    N0 = Nmax*(z-z1)*(z-z2)/(zmax-z1)/(zmax-z2)*(heaviside(z-z1)-heaviside(z-z2))
    return N0
end

function dN_max_fun(z, dN, z0, d_z)
    return -dN*exp(-(z-z0)^2/2/d_z^2)
end

function func_f(t,t_r,t_d,tau_r,tau_d)    
    return (1-exp(-(t-t_r)/tau_r))*(heaviside(t-t_r)-(heaviside(t-t_d)))+(1-exp(-(t_d-t_r)/tau_r))*exp(-(t-t_d)/tau_d)*heaviside(t-t_d)
end

function N_model_fun(z, t, z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)
    N0 = N0_fun(z, Nmax, z1, z2, zmax)
    dN_max = dN_max_fun(z, dN, z0, d_z)
    f = func_f(t,t_r,t_d,tau_r,tau_d)
    return N0 + dN_max*f
end

# z_axe = z1:z2
# t_axe = 0:2:14

# N_model = zeros(length(z_axe),length(t_axe))
# for i=1:length(z_axe)
#     for j=1:length(t_axe)
#         N_model[i,j] = N_model_fun(z_axe[i], t_axe[j], z1, z2, zmax, dN, z0, d_z, t_r, t_d, tau_r, tau_d)        
#     end
# end

# using Plots
# plot(N_model,z_axe/1000,lw=3)
# ylims!((z0/1000-10,z0/1000+10))
# xlims!((0.85*N_pump,1.15*N_pump))