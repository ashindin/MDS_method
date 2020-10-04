#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double Ne_to_fpe(double Ne)
{
    return 8978.18235285877*sqrt(abs(Ne));
}
double fpe_to_Ne(double fpe)
{
    return 1.24057537429491*pow(fpe,2)*1e-8;
}
double heaviside(double x, double a)
{
    if(x<=0) return a;
    return 1.;
}
double N0_fun(double z)
{
    double Nmax = 310143.8435737275, z1 = 150e3, z2 = 350e3, zmax = 250e3;
    double N0 = Nmax*(z-z1)*(z-z2)/(zmax-z1)/(zmax-z2)*(heaviside(z-z1, 0.)-heaviside(z-z2, 0.));
    return N0;
}
double DN(long int T, double z, double dt)
{
    double z0 = 210809.259786272;
    double dN = 26250.574920080297;
    double d_z = 3000.;
    double t_r = 0., t_d = 6., tau_r = 10., tau_d = 10.;
    double DN_max = -dN*(exp(-pow(z-z0,2))/2/pow(d_z,2));
    double f = (1-exp(-(T*dt-t_r)/tau_r))*(heaviside(T*dt-t_r, 0.)-(heaviside(T*dt-t_d, 0.)))+
    (1-exp(-(t_d-t_r)/tau_r))*exp(-(T*dt-t_d)/tau_d)*heaviside(T*dt-t_d, 0.);
    return DN_max*f;
}
double omega_fun(long int T, double z, double dt)
{
    return 2*M_PI*Ne_to_fpe(N0_fun(z)+DN(T,z,dt));
}

int main()
{
    double c = 299792458;
    double scale_ex = 4*M_PI*c*1e-7, scale_dx = 1/c;
    double dt = 8e-9, ddx = dt*(2*c);
    long int KE = 500000/ddx;

    double pulse_freq1 = 4.4*1e6, pulse_freq2 = 4.6*1e6;
    double pulse_start_time = 52.5e-6, pulse_end_time = 152.5e-6;
    
    long int pulse_start_time_ind = pulse_start_time/dt;
    long int pulse_end_time_ind = pulse_end_time/dt;
    
    // int z1 = 150000, z2 = 350000, zmax = 250000;
    // double fmax = 5e6, Nmax = fpe_to_Ne(fmax);

    double *h_axe = malloc(KE * sizeof(double));
    // double *n_axe = malloc(KE * sizeof(double));    
    // double *fp_axe = malloc(KE * sizeof(double));

    for(int i = 0; i < KE; i++)
    {
        h_axe[i] = 0. + i*ddx;
        // n_axe[i] = Nmax*(h_axe[i]-z1)*(h_axe[i]-z2)/(zmax-z1)/(zmax-z2)*
        // (heaviside(h_axe[i]-z1, 0.)-heaviside(h_axe[i]-z2, 0.));
        // fp_axe[i] = Ne_to_fpe(n_axe[i]);
    }

    double p_tstart = 2.5e-6, p_tend = 102.5e-6;
    int pT_end = 200e-6/dt;

    double *Pulse = malloc(pT_end * sizeof(double));
    
    for(int i = 0; i < pT_end; i++)
    {
        Pulse[i] = 0.33/2/scale_ex*(sin(2*M_PI*pulse_freq1*i*dt)+sin(2*M_PI*pulse_freq2*i*dt))*
        heaviside(i*dt-p_tstart, 0.)*(1-heaviside(i*dt-p_tend, 0.));
    }
    
    double pPeriod = 100e-3;
    int pPeriod_ind = pPeriod/dt;

    double f_pump = 4.6e6, N_pump = fpe_to_Ne(f_pump);

    int num_of_pulses = 601;
    long long *p_start_inds = malloc(num_of_pulses * sizeof(long long));
    for(int i = 0; i < num_of_pulses; i++)
    {
        p_start_inds[i] = 0 + i*12500000LL; // 12500000 - 100 ms period
    }

    double *dx = malloc(KE * sizeof(double));
    double *ex = malloc(KE * sizeof(double));
    double *hy = malloc(KE * sizeof(double));
    double *sx = malloc(KE * sizeof(double));
    double *sxm1 = malloc(KE * sizeof(double));
    double *sxm2 = malloc(KE * sizeof(double));
    double *omega = malloc(KE * sizeof(double));

    long long T_end = 125000000LL*60LL + 5LL*125000LL; // 60.005 s

    // printf("%i\n", KE);
    // printf("%i %i", pulse_start_time_ind, pulse_end_time_ind);
    // printf("%f\n", Nmax);
    // printf("%f\n", h_axe[KE-1]);
    // printf("%f\n", h_axe[1]);
    printf("%lld\n", p_start_inds[0]);
    printf("%lld\n", p_start_inds[1]);
    printf("%lld\n", p_start_inds[num_of_pulses-1]);
    return 0;
} 