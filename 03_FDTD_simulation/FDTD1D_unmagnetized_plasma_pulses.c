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
double DN(long long T, double z, double dt)
{
    double z0 = 210809.259786272;
    double dN = 26250.574920080297;
    double d_z = 3000.;
    double t_r = 0., t_d = 6., tau_r = 10., tau_d = 10.;
    double DN_max = -dN * (  exp( -pow(z-z0,2.)/2/pow(d_z,2.) )  );
    double f = (1-exp(-(T*dt-t_r)/tau_r))*(heaviside(T*dt-t_r, 0.)-(heaviside(T*dt-t_d, 0.)))+
    (1-exp(-(t_d-t_r)/tau_r))*exp(-(T*dt-t_d)/tau_d)*heaviside(T*dt-t_d, 0.);
    return DN_max*f;
}
double omega_fun(long long T, double z, double dt)
{
    return 2*M_PI*Ne_to_fpe(N0_fun(z)+DN(T,z,dt));
}

int main()
{
    double c = 299792458;
    double scale_ex = 4*M_PI*c*1e-7;
    double dt = 8e-9, ddx = dt*(2*c);
    long int KE = 500000/ddx;
    double pulse_freq1 = 4.4*1e6, pulse_freq2 = 4.6*1e6;
    
    double *dx = malloc(KE * sizeof(double));
    double *ex = malloc(KE * sizeof(double));
    double *hy = malloc((KE-1) * sizeof(double));
    double *sx = malloc(KE * sizeof(double));
    double *sxm1 = malloc(KE * sizeof(double));
    double *sxm2 = malloc(KE * sizeof(double));
    double *omega = malloc(KE * sizeof(double));

    double *h_axe = malloc(KE * sizeof(double));
    for(int i = 0; i < KE; i++) h_axe[i] = 0. + i*ddx;
    
    double p_tstart = 2.5e-6, p_tend = 102.5e-6;
    
    int pT_end = 200e-6/dt;
    double *Pulse = malloc(pT_end * sizeof(double));
    
    // FILE *dbg_write_ptr;
    // char dbg_filename[80];
    FILE *out_write_ptr;
    char out_filename[80];

    double nu = 1000.;

    for(int i = 0; i < pT_end; i++)
    {
        Pulse[i] = 0.33/2/scale_ex*(sin(2*M_PI*pulse_freq1*i*dt)+sin(2*M_PI*pulse_freq2*i*dt))*
        heaviside(i*dt-p_tstart, 0.)*(1-heaviside(i*dt-p_tend, 0.));
    }

    int num_of_pulses = 601;
    // int num_of_pulses = 3;
    
    // Pulse start indexes
    long long *p_start_inds = malloc(num_of_pulses * sizeof(long long));
    printf("%s\n", "Pulse start indexes:");
    for(int i = 0; i < num_of_pulses; i++)
    {
        p_start_inds[i] = 0 + i*12500000LL; // 12500000 - 100 ms period
        printf("%d - %lld\n", i, p_start_inds[i]);
    }
    printf("%s", "\n");
    
    // double *E_out = malloc(125000000 * sizeof(double)); // 1.0 s
    // double *E_out = malloc(25625000 * sizeof(double)); // 0.205 s
    double *E_out = malloc(5*125000 * sizeof(double)); // 0.005 s

    // for(int pulse_num = 0; pulse_num < num_of_pulses; pulse_num++)  
    for(int pulse_num = 1; pulse_num < num_of_pulses; pulse_num++)    
    {
        printf("Pulse number = %03d\n", pulse_num);
        // Arrays initialization
        for(int i = 0; i<KE; i++)
        {
            omega[i] = omega_fun(p_start_inds[pulse_num], h_axe[i], dt);
            dx[i]    = 0.;
            ex[i]    = 0.;        
            sx[i]    = 0.;
            sxm1[i]  = 0.;
            sxm2[i]  = 0.;
        }
        for(int i=0; i<KE-1; i++) hy[i] = 0.;
        double dx_old_1 = 0., dx_old_N_1 = 0.;
        
        sprintf(out_filename,"E_out_%012lld.bin",p_start_inds[pulse_num]);
        E_out[0] = ex[5];        

        for(long long T = p_start_inds[pulse_num]+1; T < p_start_inds[pulse_num] + 5*125000; T++)
        {
            dx_old_1 = dx[1];
            dx_old_N_1 = dx[KE-2];
        
            for (int i = 1; i < KE-1; i++) dx[i] = dx[i] - 0.5 * (hy[i] - hy[i-1]);

            if(T-p_start_inds[pulse_num] < pT_end) dx[5] = dx[5] + Pulse[T-p_start_inds[pulse_num]];

            dx[0]    = dx_old_1   + (-1./3.) * (dx[1]    - dx[0]   );  // bc0
            dx[KE-1] = dx_old_N_1 + (-1./3.) * (dx[KE-2] - dx[KE-1]);  // bcN

            for (int i = 0; i < KE; i++) ex[i] = dx[i]-sx[i];
                    
            E_out[T%(5*125000)] = ex[5]; // rewrite it!

            for (int i = 0; i < KE; i++) 
            {
                sx[i] = (1+exp(-nu*dt))*sxm1[i] - exp(-nu*dt) * sxm2[i] + pow(omega[i],2.)*dt/nu*(1-exp(-nu*dt))*ex[i];
                sxm2[i] = sxm1[i];
                sxm1[i] = sx[i];            
            }

            for (int i = 0; i < KE-1; i++) hy[i] = hy[i] - 0.5 * (ex[i+1] - ex[i]);
        }
        out_write_ptr = fopen(out_filename, "wb");
        fwrite(E_out,sizeof(E_out), 5*125000, out_write_ptr);
        fclose(out_write_ptr);

    }    

    return 0;
} 