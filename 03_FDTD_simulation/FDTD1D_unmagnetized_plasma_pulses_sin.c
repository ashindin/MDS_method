#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{  
    FILE *out_write_ptr1, *out_write_ptr2;
    char out_filename1[80], out_filename2[80];

    int num_of_pulses = 601;  
    
    double f0 = 4500000.;
    double dt = 8e-9;

    // Pulse start indexes
    long long *p_start_inds = malloc(num_of_pulses * sizeof(long long));
    printf("%s\n", "Pulse start indexes:");
    for(int i = 0; i < num_of_pulses; i++)
    {
        p_start_inds[i] = 0 + i*12500000LL; // 12500000 - 100 ms period
        printf("%d - %lld\n", i, p_start_inds[i]);
    }
    printf("%s", "\n");    
    
    double *Sin_out = malloc(5*125000 * sizeof(double)); // 0.005 s
    double *Cos_out = malloc(5*125000 * sizeof(double)); // 0.005 s
    
    for(int pulse_num = 0; pulse_num < num_of_pulses; pulse_num++)    
    {
        printf("Pulse number = %03d\n", pulse_num);
        
        sprintf(out_filename1,"sin_%012lld.bin",p_start_inds[pulse_num]);     
        sprintf(out_filename2,"cos_%012lld.bin",p_start_inds[pulse_num]);     

        for(long long T = p_start_inds[pulse_num]; T < p_start_inds[pulse_num] + 5*125000; T++)
        {                   
            Sin_out[T%(5*125000)] = sin(2*M_PI*f0*T*dt);
            Cos_out[T%(5*125000)] = cos(2*M_PI*f0*T*dt);
        }
        out_write_ptr1 = fopen(out_filename1, "wb");
        out_write_ptr2 = fopen(out_filename2, "wb");
        fwrite(Sin_out,sizeof(Sin_out), 5*125000, out_write_ptr1);
        fwrite(Cos_out,sizeof(Cos_out), 5*125000, out_write_ptr2);
        fclose(out_write_ptr1);
        fclose(out_write_ptr2);
    }    
    return 0;
} 