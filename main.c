#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#define CROW 4 // Number of cells in a row
#define M 4 // Number of particles in a cell
#define m 1.0
#define Kb 1.0
#define epsilon 1.0
#define PI 3.1415926535897932384626433
#define BINS 150 // number of bins for the histogram of g(r)
#define LOGARITHMIC_VOLUME_MOVES 0

#define N (CROW*CROW*CROW*M) // Total number of particles (constant)
double P = 0.0; // Pressure (constant) 
double T = 0.0; // Temperature (constant)

double L = 0.0; // Cube Side (variable)
#define V (L*L*L) // Volume (variable)

double delta_s = 0.0; // Jitter (position move)
double delta_V = 0.0; // Jitter (volume move)
double delta_l = 0.0; // Jitter (log volume move)


// (x/L, y/L, z/L)
typedef struct Particle {
    double x; double y; double z; 
} Particle;

typedef struct vec3 {
    double x; double y; double z;
} vec3;

typedef struct Measure {
    double temp; // temperature
    double comp; // compressibilità
    double ener; // energia
    double vol;  // volume
    double rho;  // volume
} Measure;

// Lennard Jones potential energy separated in the two components (V6, V12). Vtot = V6 + V12
typedef struct LennardJonesPotential {
    double V6;  // part that goes like 1/r6
    double V12; // part that goes like 1/r12
} LennardJonesPotential;

LennardJonesPotential add_potentials(LennardJonesPotential P1, LennardJonesPotential P2) {
    return (LennardJonesPotential) { P1.V6 + P2.V6, P1.V12 + P2.V12 };
}



//                      HISTOGRAM OF PARTICLE DISTANCES

double* histogram;

void add_to_distribution_histogram(Particle S[]) {
    double bin_width = 1.0 / BINS;
    for (int i = 1; i < N; i++) {
        for (int j = 1; j < N; j++) {
            if(i == j) continue;
            vec3 Rij = (vec3) {
                1.0 * (  (S[i].x - S[j].x) - rint(S[i].x - S[j].x)  ),
                1.0 * (  (S[i].y - S[j].y) - rint(S[i].y - S[j].y)  ),
                1.0 * (  (S[i].z - S[j].z) - rint(S[i].z - S[j].z)  ),
            };
            double r = sqrt(Rij.x*Rij.x + Rij.y*Rij.y + Rij.z*Rij.z);
            int bin_index = r / bin_width;
            if(r > sqrt(3)*1.0) { printf("Distanza fuori scala %f > %f\n", r, 1.0);  }
            if(bin_index > BINS - 1) { 
                printf("Indice fuori scala %d > %d\n", bin_index, BINS);
            } else {
                histogram[bin_index] += 1; 
            }
        }
    }
}

void save_distribution_histogram(FILE* file, unsigned int steps) {
    double bin_width = 1.0 / BINS;
    // Normalize and print to file
    for (int j = 1; j < BINS; j++) {
        double dV = 4*PI/3.0*((j+1)*(j+1)*(j+1) - j*j*j)*bin_width*bin_width*bin_width;
        histogram[j] /= dV * (m*N/V) * N * (steps/2);
        fprintf(file, "%10.5e %10.5e\n", (j + 0.0) / BINS, histogram[j]);
    }
}




//                      SCATTER PLOT OF PARTICLES' POSITIONS

void save_scatter_plot(Particle S[]) {
    FILE* scatter = fopen("data/scatter.dat", "w+");
    vec3 center = (vec3) {0.5, 0.5, 0.5};
    for(int i = 0; i < N; i++) {
        vec3 Ri = (vec3) {
            (S[i].x - center.x) - rint(S[i].x - 0.5),
            (S[i].y - center.y) - rint(S[i].y - 0.5),
            (S[i].z - center.z) - rint(S[i].z - 0.5),
        };
        fprintf(scatter, "%f %f %f\n", Ri.x, Ri.y, Ri.z);
    }
}


//                      LENNARD JONES

LennardJonesPotential get_interaction_lennard_jones(Particle S[], int i, int j) {
    // Calculate V_cut
    vec3 Rij = (vec3) {
        L * (  (S[i].x - S[j].x) - rint(S[i].x - S[j].x)  ),
        L * (  (S[i].y - S[j].y) - rint(S[i].y - S[j].y)  ),
        L * (  (S[i].z - S[j].z) - rint(S[i].z - S[j].z)  ),
    };
    double r2 = Rij.x*Rij.x + Rij.y*Rij.y + Rij.z*Rij.z;
    if(r2 < (L/2)*(L/2)) { // truncate forces after r = L/2
        // Calculate V(L/2)
        double lennard_jones_L_halves_6;
        double lennard_jones_L_halves_12;
        {
            double sigma2 = 1.0;
            double sigma_su_L_mezzi_2 = sigma2 / (L/2) / (L/2);
            double sigma_su_L_mezzi_6 = sigma_su_L_mezzi_2 * sigma_su_L_mezzi_2 * sigma_su_L_mezzi_2;
            lennard_jones_L_halves_6 = -4*epsilon*sigma_su_L_mezzi_6;
            lennard_jones_L_halves_12 = 4*epsilon*sigma_su_L_mezzi_6*sigma_su_L_mezzi_6;
        }
        double sr6 = pow(r2, -3.0);
        return (LennardJonesPotential) {
            -4*epsilon*sr6 - lennard_jones_L_halves_6,
            4*epsilon*sr6*sr6 - lennard_jones_L_halves_12
        };
    }

    return (LennardJonesPotential) {0, 0};
}

LennardJonesPotential get_potential_energy_lennard_jones(Particle S[]) {
    // Calculate total potential energy
    LennardJonesPotential potential_energy = (LennardJonesPotential) {0, 0};
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < i; j++) {
            potential_energy = add_potentials(potential_energy, get_interaction_lennard_jones(S, i, j));
        }
    }
    return potential_energy;
}

double get_forces_lennard_jones(Particle S[], vec3 F[]) {
    double sigma2 = 1.0;

    double W_forces = 0.0;
    memset(F, 0, N*sizeof(vec3)); // zero out forces
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < i; j++) {
            vec3 Rij = (vec3) {
                L * (  (S[i].x - S[j].x) - rint(S[i].x - S[j].x)  ),
                L * (  (S[i].y - S[j].y) - rint(S[i].y - S[j].y)  ),
                L * (  (S[i].z - S[j].z) - rint(S[i].z - S[j].z)  ),
            };
            double r2 = Rij.x*Rij.x + Rij.y*Rij.y + Rij.z*Rij.z;
            if(r2 < (L/2.0)*(L/2.0)) { // truncate forces after r = L/2
                double sr2 = sigma2 / r2;
                double sr6 = sr2 * sr2 * sr2;
                double factor = 24.0*epsilon*sr6*(2.0*sr6 - 1.0) / r2; 
                vec3 Fij = (vec3) { Rij.x * factor, Rij.y * factor, Rij.z * factor };
                F[i].x += Fij.x;
                F[i].y += Fij.y;
                F[i].z += Fij.z;
                F[j].x -= Fij.x;
                F[j].y -= Fij.y;
                F[j].z -= Fij.z;
                W_forces += Fij.x*Rij.x + Fij.y*Rij.y + Fij.z*Rij.z;
            }
        }
    }
    W_forces /= N;

    return W_forces;
}




//                      MEASURES

FILE* file_measures = NULL;
Measure measure_metropolis(Particle S[], vec3 F[], double time, double potential_energy, double W_forces) {
    // Compressibilità
    double compressibility = 1 + W_forces / 3.0 / Kb / T;

    // 1. Step     2. PotentialEnergyPerParticle     3. Compressibility
    if(file_measures) fprintf(file_measures, "%10.10e %10.10e %10.10e %10.10e\n", time, potential_energy / N, compressibility, V);

    return (Measure) { T, compressibility, potential_energy / N, V, m*N/V };
}

Measure averages(Measure* measures, int start, int end) {
    Measure avg = (Measure) {0, 0, 0};
    for(int i = start; i < end; i++) { 
        avg.temp += measures[i].temp; 
        avg.comp += measures[i].comp; 
        avg.ener += measures[i].ener; 
        avg.vol  += measures[i].vol; 
        avg.rho  += measures[i].rho; 
    }
    avg.temp /= end - start;
    avg.comp /= end - start;
    avg.ener /= end - start;
    avg.vol  /= end - start;
    avg.rho  /= end - start;
    return avg;
}

Measure variances(Measure* measures, Measure avg, int start, int end) {
    Measure var = (Measure) {0, 0, 0};
    for(int i = start; i < end; i++) { 
        var.temp += (measures[i].temp - avg.temp)*(measures[i].temp - avg.temp); 
        var.comp += (measures[i].comp - avg.comp)*(measures[i].comp - avg.comp); 
        var.ener += (measures[i].ener - avg.ener)*(measures[i].ener - avg.ener); 
        var.vol  += (measures[i].vol  - avg.vol )*(measures[i].vol  - avg.vol ); 
        var.rho  += (measures[i].rho  - avg.rho )*(measures[i].rho  - avg.rho ); 
    }
    var.temp /= end - start;
    var.comp /= end - start;
    var.ener /= end - start;
    var.vol  /= end - start;
    var.rho  /= end - start;
    return var;
}

void averages_and_plateau(Measure *measures, int steps, int k) {
    // Global Averages (on the second half of the data)
    Measure avg = averages(measures, steps/2, steps);
    Measure var = variances(measures, avg, steps/2, steps); // Meaningful only with velocity verlet

    printf("Temperatura = %f ± %f\n", avg.temp, sqrt(var.temp));
    printf("Compressibilità = %f ± %f\n", avg.comp, sqrt(var.comp));
    printf("Energia = %f ± %2.5e\n", avg.ener, sqrt(var.ener));
    printf("Volume = %f ± %2.5e\n", avg.vol, sqrt(var.vol));
    printf("Densità = %f ± %2.5e\n", avg.rho, sqrt(var.rho));

    // Plateau file
    char varAvgBFile_name[20]; 
    sprintf(varAvgBFile_name, "data/plateau_%02d.dat", k);
    FILE* varAvgBFile = fopen(varAvgBFile_name, "w+");

    // Variances for metropolis
    int start = steps / 2;
    int NN = steps - start; // Size of the data
    for(int B = 1; B < NN / 2; B++) {  // Loop over the number of points per block
        int NB = NN / B; // The number of blocks
        Measure* avgB = malloc(sizeof(Measure) * NB); // Averages in the blocks 
        for(int b = 0; b < NB; b++){
            avgB[b] = averages(measures, start + B*b, start + B*(b + 1));
        }
        Measure varAvgB = variances(avgB, averages(avgB, 0, NB), 0, NB); // Variances in the blocks
        if(varAvgBFile) fprintf(varAvgBFile, "%d %f %f\n", B, sqrt(varAvgB.ener / NB), sqrt(varAvgB.comp / NB));
        free(avgB);
    }

    fclose(varAvgBFile);
}




//                      INITIALIZATIONS

void initialize_constants(double simulation[4]) {
    P = simulation[0]; // Pressure
    T = simulation[1]; // Temperature
    double rho0 = simulation[2]; // Initial density
    delta_s = simulation[3]; // Random Jump Size (particle move)
    delta_V = simulation[4]; // Random Jump Size (volume move)
    delta_l = simulation[5]; // Random Jump Size (volume move)

    // TODO: initialize pressure P!?

    L = CROW * pow(M / rho0, 1.0/3.0); // rho = N/L3 = P3*M/L3
    
    printf("\nCROW = %d, N = %d, L0 = %f, T = %f, rho0 = %f, P = %f\n", CROW, N, L, T, rho0, P);
}

void initialize_positions(Particle S0[]) {
    double a = 1.0 / (double)CROW;
    
    vec3 CC[1]  = { (vec3){0, 0, 0} };
    vec3 BCC[2] = { (vec3){0, 0, 0}, (vec3){0.5, 0.5, 0.5} };
    vec3 FCC[4] = { (vec3){0, 0, 0}, (vec3){0.5, 0.5, 0}, (vec3){0.5, 0, 0.5}, (vec3){0, 0.5, 0.5} };
    vec3* CELL_TYPES[] = {NULL, CC, BCC, NULL, FCC};

    unsigned int index = 0;
    for(int nx = 0; nx <= CROW - 1; nx++) {
    for(int ny = 0; ny <= CROW - 1; ny++) {
    for(int nz = 0; nz <= CROW - 1; nz++) {
        for(int p = 0; p < M; p++) {
            S0[index].x = a * nx + a*CELL_TYPES[M][p].x; 
            S0[index].y = a * ny + a*CELL_TYPES[M][p].y; 
            S0[index].z = a * nz + a*CELL_TYPES[M][p].z;
            index += 1;            
        }
    }}}
}




//                      SYSTEM MOVES

// Get random number uniformly distributed between min and max
double sample_uniform(double min, double max) {
    double x = rand()/(RAND_MAX + 1.0);
    return min + x*(max - min);
}

LennardJonesPotential move_one_particle(Particle S0[], int k, vec3 delta) {
    LennardJonesPotential interaction_before = (LennardJonesPotential) {0, 0};
    LennardJonesPotential interaction_after = (LennardJonesPotential) {0, 0};
    for(int i = 0; i < N; i++) {
        if(i == k) continue;
        interaction_before = add_potentials(interaction_before, get_interaction_lennard_jones(S0, k, i));
    }
    S0[k].x += delta.x;
    S0[k].y += delta.y;
    S0[k].z += delta.z;
    for(int i = 0; i < N; i++) {
        if(i == k) continue;
        interaction_after = add_potentials(interaction_after, get_interaction_lennard_jones(S0, k, i));
    }

    // return the change in potential energy
    return (LennardJonesPotential) {
        interaction_after.V6  - interaction_before.V6,
        interaction_after.V12 - interaction_before.V12,
    };
}




//                      MONTE CARLO SIMULATION

Measure* metropolis(Particle S0[], unsigned int steps) {
    vec3* F = malloc(N * sizeof(vec3));

    // Calcolo le forze (serve per la compressibilità)
    double W_forces = get_forces_lennard_jones(S0, F);

    // Calcolo l'energia potenziale
    LennardJonesPotential V0 = get_potential_energy_lennard_jones(S0);

    // Inizializzo le misure e misuro
    Measure* measures = calloc(steps + 1, sizeof(Measure));
    measures[0] = measure_metropolis(S0, F, 0.0, V0.V6 + V0.V12, W_forces);

    int position_jumps = 0;
    int volume_jumps = 0;
    printf("Step: 0/%d", steps);
    for(int i = 1; i <= steps; i++) {
        if(i % 100 == 0) { 
            printf("\rStep: %d/%d, Acceptance: %.2f", i, steps, (position_jumps + volume_jumps) / (double)i); fflush(stdout); 
        }

        int k = floor(sample_uniform(0.0, N+1));
        if(k == N || k == N + 1) {
            if(LOGARITHMIC_VOLUME_MOVES) {
                // MOSSA NEL LOGARITMO DEL VOLUME (TODO)
            } else {
                // MOSSA DI VOLUME
                double epsilon_volume = sample_uniform(-0.5, 0.5);
                double nuovo_volume = fmax(0.0, V + epsilon_volume * delta_V);
                double differenza_volumi = nuovo_volume - V;
                double rapporto_volumi = nuovo_volume / V;
                double differenza_potenziali = 
                    V0.V6 *(powl(rapporto_volumi,  -6.0  / 3.0) - 1) + 
                    V0.V12*(powl(rapporto_volumi,  -12.0 / 3.0) - 1); 
                double jump_probability = fmin(1, 
                    powl(rapporto_volumi, N) * exp(-(differenza_potenziali + P*differenza_volumi)/Kb/T)
                );
                if(sample_uniform(0, 1) < jump_probability) {
                    // Accetto il salto
                    L = powl(nuovo_volume, 1.0/3.0);
                    V0 = (LennardJonesPotential) {
                        V0.V6  * powl(rapporto_volumi,  -6.0 / 3.0),
                        V0.V12 * powl(rapporto_volumi, -12.0 / 3.0),
                    };
                    W_forces = get_forces_lennard_jones(S0, F);
                    volume_jumps += 1;
                } else {
                    // Rifiuto il salto
                }
            }
        }else{
            // MOSSA DI POSIZIONE
            vec3 delta = (vec3) { 
                sample_uniform(-0.5, 0.5) * delta_s, 
                sample_uniform(-0.5, 0.5) * delta_s, 
                sample_uniform(-0.5, 0.5) * delta_s 
            };
            LennardJonesPotential change_in_potential_energy = move_one_particle(S0, k, delta);
            double diff = change_in_potential_energy.V6 + change_in_potential_energy.V12;
            double jump_probability = fmin(1, exp(-diff/Kb/T));
            if(sample_uniform(0, 1) < jump_probability) {
                // Accetto il salto
                V0 = add_potentials(V0, change_in_potential_energy);
                W_forces = get_forces_lennard_jones(S0, F);
                position_jumps += 1;
            } else {
                // Rifiuto il salto
                S0[k].x -= delta.x;
                S0[k].y -= delta.y;
                S0[k].z -= delta.z;
            }
        }

        // Misura delle osservabili
        measures[i] = measure_metropolis(S0, F, i, V0.V6 + V0.V12, W_forces);
        if(histogram && i >= steps/2) add_to_distribution_histogram(S0);
    }
    printf("\rAcceptance: %f, position_jumps/volume_jumps: %f\t\t\t\t\n", (volume_jumps + position_jumps) / ((double) steps), position_jumps /  (double)volume_jumps); 

    free(F);

    return measures;
}

int main() {
    const int SIM = 1;
    double simulations[SIM][6] = {
    //   P       T       rho0      delta_s  delta_V  delta_l
    //  {1.0,    1.1,    0.001,    0.20,      10,    1.00},
    //  {1.0,    1.1,    0.010,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    0.100,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    0.200,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    0.250,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    0.300,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    0.400,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    0.500,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    0.600,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    0.700,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    0.750,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    0.800,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    1.000,    0.01,    0.01,    0.01},
    //  {1.0,    1.1,    1.200,    0.02,      10,    1.00}, 
        {4.0,    1.1,    1.200,    0.02,      10,    1.00}, 
    };

    Particle* S0 = malloc(N * sizeof(Particle));
    histogram = malloc(BINS * sizeof(double));
    unsigned int steps = 1e5;
    for(int k = 0; k < SIM; k++){

        // Set file to which measurements are written (measure_metropolis)
        char file_measures_name[40]; 
        sprintf(file_measures_name, "data/measures_%02d.dat", k);
        file_measures = fopen(file_measures_name, "w+");

        // Reset histogram
        memset(histogram, 0, BINS * sizeof(double)); 

        // Initialize T, delta_s, delta_V, delta_l, L
        initialize_constants(simulations[k]); 
        // Place particles in the grid
        initialize_positions(S0);

        Measure* measures = metropolis(S0, steps);

        // Calculate averages and produce plateau graph for error
        averages_and_plateau(measures, steps, k);

        free(measures);
        fclose(file_measures);

    }


    // Finalize and export histogram
    FILE* histogram_file = fopen("data/histogram.dat", "w+");
    save_distribution_histogram(histogram_file, steps);
    fclose(histogram_file);

    // Scatter plot of particles' positions
    save_scatter_plot(S0);

    return 0;
}