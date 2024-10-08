#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#define CROW 5 // Number of cells in a row
#define M 4 // Number of particles in a cell
#define m 1.0
#define Kb 1.0
#define epsilon 1.0
#define PI 3.1415926535897932384626433
#define HISTOGRAM_SAMPLE_RATE 500
#define BINS 150 // number of bins for the histogram of g(r)
#define MAX_DENSITY_BIN 1.0
#define DENSITY_BINS 300 // number of bins for the histogram of rho during the simulation
#define LOGARITHMIC_VOLUME_MOVES 1

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
    double ener; // energia per particella (meno energia cinetica)
    double vol;  // volume
    double rho;  // densità
    double enth;  // entalpia totale (meno energia cinetica)
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

double* histogram = NULL;

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
    int samples = steps / HISTOGRAM_SAMPLE_RATE;
    // Normalize and print to file
    for (int j = 1; j < BINS; j++) {
        double dV = 4*PI/3.0*((j+1)*(j+1)*(j+1) - j*j*j)*bin_width*bin_width*bin_width;
        histogram[j] /= dV * (m*N/V) * N * (samples/2);
        fprintf(file, "%10.5e %10.5e\n", j / (double)BINS, histogram[j]);
    }
}




//                      HISTOGRAM OF SYSTEM DENSITIES

int* density_histogram = NULL;

void add_to_density_histogram(double density) {
    double bin_width = MAX_DENSITY_BIN / DENSITY_BINS;
    if(density < MAX_DENSITY_BIN & density > 0) {
        density_histogram[(int) floor(density / bin_width)] += 1;
    }else{
        printf("Density expected to be between 0.0 and %f", MAX_DENSITY_BIN);
        exit(1);
    }
}

void save_density_histogram(FILE* file) {
    double bin_width = MAX_DENSITY_BIN / DENSITY_BINS;
    for (int j = 0; j < DENSITY_BINS; j++) {
        fprintf(file, "%10.5e %d\n", j * bin_width, density_histogram[j]);
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
Measure measure_metropolis(Particle S[], double time, double potential_energy, double W_forces) {
    double rho = m * N / V;
    
    // Compressibilità
    double compressibility = 1 + W_forces / 3.0 / Kb / T;
    compressibility = P / rho / Kb / T;

    // 1. Step     2. PotentialEnergyPerParticle     3. Compressibility     4. Volume
    if(file_measures) fprintf(file_measures, "%10.10e %10.10e %10.10e %10.10e\n", time, potential_energy / N, compressibility, V);

    double enthalpy = (potential_energy + P * V) / N; // entalpia per particella

    return (Measure) { T, compressibility, potential_energy / N, V, rho, enthalpy };
}

Measure averages(Measure* measures, int start, int end) {
    Measure avg = (Measure) {0, 0, 0, 0, 0, 0};
    for(int i = start; i < end; i++) { 
        avg.temp += measures[i].temp; 
        avg.comp += measures[i].comp; 
        avg.ener += measures[i].ener; 
        avg.vol  += measures[i].vol; 
        avg.rho  += measures[i].rho; 
        avg.enth += measures[i].enth; 
    }
    avg.temp /= end - start;
    avg.comp /= end - start;
    avg.ener /= end - start;
    avg.vol  /= end - start;
    avg.rho  /= end - start;
    avg.enth /= end - start;
    return avg;
}

Measure variances(Measure* measures, Measure avg, int start, int end) {
    Measure var = (Measure) {0, 0, 0, 0, 0, 0};
    for(int i = start; i < end; i++) { 
        var.temp += (measures[i].temp - avg.temp)*(measures[i].temp - avg.temp); 
        var.comp += (measures[i].comp - avg.comp)*(measures[i].comp - avg.comp); 
        var.ener += (measures[i].ener - avg.ener)*(measures[i].ener - avg.ener); 
        var.vol  += (measures[i].vol  - avg.vol )*(measures[i].vol  - avg.vol ); 
        var.rho  += (measures[i].rho  - avg.rho )*(measures[i].rho  - avg.rho ); 
        var.enth += (measures[i].enth - avg.enth)*(measures[i].enth - avg.enth); 
    }
    var.temp /= end - start;
    var.comp /= end - start;
    var.ener /= end - start;
    var.vol  /= end - start;
    var.rho  /= end - start;
    var.enth /= end - start;
    return var;
}

void averages_and_plateau(Measure *measures, int steps, int k, Measure* avg, Measure* var, double* cp) {
    // Global Averages (on the second half of the data)
    *avg = averages(measures, steps/2, steps);
    *var = variances(measures, *avg, steps/2, steps); // Meaningful only with velocity verlet
    *cp = 1.5 + var->enth / Kb / T / T * N;

    printf("Temperatura = %f ± %f\n", avg->temp, sqrt(var->temp));
    printf("Compressibilità = %f ± %f\n", avg->comp, sqrt(var->comp));
    printf("Energia = %f ± %2.5e\n", avg->ener, sqrt(var->ener));
    printf("Volume = %f ± %2.5e\n", avg->vol, sqrt(var->vol));
    printf("Densità = %f ± %2.5e\n", avg->rho, sqrt(var->rho));
    printf("Calore specifico (cP/kB) = 5/2 + %f = %f\n", *cp - 2.5, *cp);

    // Plateau file
    char varAvgBFile_name[20]; 
    sprintf(varAvgBFile_name, "data/plateau_%02d.dat", k);
    FILE* varAvgBFile = fopen(varAvgBFile_name, "w+");

    // Variances for metropolis
    int start = steps / 2;
    int NN = steps - start; // Size of the data
    for(int B = 1; B < (80000 < NN/2 ? 80000 : NN/2); B++) {  // Loop over the number of points per block
        int NB = NN / B; // The number of blocks
        Measure* avgB = malloc(sizeof(Measure) * NB); // Averages in the blocks 
        for(int b = 0; b < NB; b++){
            avgB[b] = averages(measures, start + B*b, start + B*(b + 1));
        }
        Measure varAvgB = variances(avgB, averages(avgB, 0, NB), 0, NB); // Variances in the blocks
        if(varAvgBFile) fprintf(varAvgBFile, "%d %10.10e %10.10e %10.10e\n", B, sqrt(varAvgB.ener / NB), sqrt(varAvgB.comp / NB), sqrt(varAvgB.vol / NB));
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

    L = CROW * pow(M / rho0, 1.0/3.0); // rho = N/L3 = P3*M/L3
    
    if( LOGARITHMIC_VOLUME_MOVES) printf("CROW = %d, N = %d, L0 = %f, T = %f, rho0 = %f, P = %f [LOG]\n", CROW, N, L, T, rho0, P);
    if(!LOGARITHMIC_VOLUME_MOVES) printf("CROW = %d, N = %d, L0 = %f, T = %f, rho0 = %f, P = %f [NOLOG]\n", CROW, N, L, T, rho0, P);
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
    double W_forces = 0.0; // get_forces_lennard_jones(S0, F);

    // Calcolo l'energia potenziale
    LennardJonesPotential V0 = get_potential_energy_lennard_jones(S0);

    // Inizializzo le misure e misuro
    Measure* measures = calloc(steps + 1, sizeof(Measure));
    measures[0] = measure_metropolis(S0, 0.0, V0.V6 + V0.V12, W_forces);

    int position_jumps = 0;
    int volume_jumps = 0;
    int position_steps = 0;
    int volume_steps = 0;
    printf("Step: 0/%d", steps);
    double last_position_acceptance = 0.5;
    double last_volume_acceptance = 0.5;
    for(int i = 1; i <= steps; i++) {
        if(i % 5000 == 0) { 
            double acceptance = (volume_jumps + position_jumps) / ((double) i);
            double position_acceptance = position_jumps /  (double) position_steps;
            double volume_acceptance = volume_jumps / (double) volume_steps;
            printf("\rStep: %d/%d, Acceptance: %f, Position Acceptance: %f (%3.1e), Volume Acceptance: %f (%3.1e)\t\t\t\t", 
                i, steps,
                acceptance,
                position_acceptance, delta_s, 
                volume_acceptance, LOGARITHMIC_VOLUME_MOVES ? delta_l : delta_V
            ); 
            fflush(stdout); 
            if(i < steps/2) {
                delta_s *= 1.0 + fmax(-0.5, (position_acceptance - 0.5)) // forza di richiamo
                               + fmax(-0.5, fmin(0.2, (position_acceptance - last_position_acceptance) * 2)); // attrito
                delta_V *= 1.0 + fmax(-0.5, (volume_acceptance - 0.5)) // forza di richiamo
                               + fmax(-0.5, fmin(0.2, (volume_acceptance - last_volume_acceptance) * 2)); // attrito
                delta_l *= 1.0 + fmax(-0.5, (volume_acceptance - 0.5)) // forza di richiamo
                               + fmax(-0.5, fmin(0.2, (volume_acceptance - last_volume_acceptance) * 2)); // attrito
                delta_s = fmin(1.0, delta_s);
            }
            last_position_acceptance = position_acceptance;
            last_volume_acceptance = volume_acceptance;
        }

        int k = floor(sample_uniform(0.0, N + 1));
        if(k == N || k == N + 1) {
            volume_steps += 1;
            if(LOGARITHMIC_VOLUME_MOVES) {
                // MOSSA NEL LOGARITMO DEL VOLUME
                double epsilon_l = sample_uniform(-0.5, 0.5);
                double lx = log(V);
                double ly = lx + epsilon_l * delta_l;
                double nuovo_volume = exp(ly);
                double differenza_volumi = nuovo_volume - V;
                double rapporto_volumi = nuovo_volume / V;
                double differenza_potenziali = 
                    V0.V6 *(powl(rapporto_volumi,  -6.0  / 3.0) - 1) + 
                    V0.V12*(powl(rapporto_volumi,  -12.0 / 3.0) - 1); 
                double jump_probability = fmin(1, 
                    powl(rapporto_volumi, N + 1) * exp(-(differenza_potenziali + P*differenza_volumi)/Kb/T)
                );
                if(sample_uniform(0, 1) < jump_probability) {
                    // Accetto il salto
                    L = powl(nuovo_volume, 1.0/3.0);
                    V0 = (LennardJonesPotential) {
                        V0.V6  * powl(rapporto_volumi,  -6.0 / 3.0),
                        V0.V12 * powl(rapporto_volumi, -12.0 / 3.0),
                    };
                    // W_forces = get_forces_lennard_jones(S0, F);
                    volume_jumps += 1;
                } else {
                    // Rifiuto il salto
                }
            } else {
                // MOSSA DI VOLUME
                double epsilon_volume = sample_uniform(-0.5, 0.5);
                double nuovo_volume = fmax(1e-8, V + epsilon_volume * delta_V);
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
                    // W_forces = get_forces_lennard_jones(S0, F);
                    volume_jumps += 1;
                } else {
                    // Rifiuto il salto
                }
            }
        }else{
            // MOSSA DI POSIZIONE
            position_steps += 1;
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
                // W_forces = get_forces_lennard_jones(S0, F);
                position_jumps += 1;
            } else {
                // Rifiuto il salto
                S0[k].x -= delta.x;
                S0[k].y -= delta.y;
                S0[k].z -= delta.z;
            }
        }

        // Misura delle osservabili
        measures[i] = measure_metropolis(S0, i, V0.V6 + V0.V12, W_forces);
        if(histogram && i >= steps/2 && i % HISTOGRAM_SAMPLE_RATE == 0) add_to_distribution_histogram(S0);
        if(density_histogram && i >= steps/2) add_to_density_histogram(m*N/V);
    }
    printf("\rAcceptance: %f, Position Acceptance: %f, Volume Acceptance: %f\t\t\t\t\t\t\t\t \n", 
        (volume_jumps + position_jumps) / ((double) steps), 
        position_jumps /  (double) position_steps, 
        volume_jumps / (double) volume_steps
    ); 

    free(F);

    return measures;
}

int main() {
    const int SIM = 53;
    double simulations[SIM][6] = {

    //   P            T       rho0      delta_s  delta_V  delta_l
    //  {  0.001,     1.0,    0.001,    0.700,     10,    0.35},
    //  {  0.050,     1.0,    0.300,    0.065,     10,    0.08},
    //  {  0.500,     1.0,    0.700,    0.040,     10,    0.05},
    //  {  5.000,     1.0,    1.000,    0.02,      10,    0.03},
    //  {  50,        1.0,    1.300,    0.02,      10,    0.02},
    //  {  80,        1.0,    1.400,    0.02,      10,    0.02},

    //   P             T       rho0      delta_s  delta_V  delta_l
    //  {  0.001,      1.0,    0.001,    0.700,     10,    0.35},
    //  {  0.001,      1.0,    0.010,    0.700,     10,    0.35},
    //  {  0.001,      1.0,    0.100,    0.700,     10,    0.35},
    //  {  0.001,      1.0,    0.300,    0.700,     10,    0.35},
    //  {  0.001,      1.0,    0.600,    0.700,     10,    0.35},
    //  {  0.001,      1.0,    0.900,    0.700,     10,    0.35},
    //  {  0.001,      1.0,    1.200,    0.700,     10,    0.35},

    //     P           T       rho0      delta_s  delta_V  delta_l
      //  {  0.001,      2.0,    0.049,    0.700,     10,    0.35},
    //  {  0.002,      2.0,    0.099,    0.700,     10,    0.35},
    //  {  0.005,      2.0,    0.025,    0.500,     10,    0.30},
    //  {  0.010,      2.0,    0.050,    0.450,     10,    0.29},
    //  {  0.020,      2.0,    0.010,    0.450,     10,    0.28},
    //  {  0.040,      2.0,    0.020,    0.350,     10,    0.27},
    //  {  0.080,      2.0,    0.042,    0.350,     10,    0.26},
    //  {  0.090,      2.0,    0.047,    0.300,     10,    0.25},
    //  {  0.100,      2.0,    0.053,    0.300,     10,    0.24},
    //  {  0.110,      2.0,    0.060,    0.250,     10,    0.23},
    //  {  0.150,      2.0,    0.082,    0.250,     10,    0.22},
      //  {  0.200,      2.0,    0.113,    0.200,     10,    0.21},
    //  {  0.250,      2.0,    0.144,    0.200,     10,    0.20},
    //  {  0.300,      2.0,    0.175,    0.150,     10,    0.19},
      //  {  0.350,      2.0,    0.203,    0.150,     10,    0.18},
      //  {  0.400,      2.0,    0.232,    0.100,     10,    0.17},
    //  {  0.450,      2.0,    0.261,    0.100,     10,    0.16},
      //  {  0.500,      2.0,    0.290,    0.050,     10,    0.15},
    //  {  0.550,      2.0,    0.325,    0.050,     10,    0.14},
      //  {  0.600,      2.0,    0.330,    0.050,     10,    0.13},
    //  {  0.650,      2.0,    0.374,    0.050,     10,    0.12},
    //  {  0.700,      2.0,    0.375,    0.050,     10,    0.11},
      //  {  0.750,      2.0,    0.406,    0.050,     10,    0.08},
    //  {  0.800,      2.0,    0.396,    0.045,     10,    0.08},
    //  {  0.850,      2.0,    0.428,    0.045,     10,    0.08},
    //  {  0.900,      2.0,    0.429,    0.040,     10,    0.08},
    //  {  0.950,      2.0,    0.445,    0.040,     10,    0.08},
    //  {  1.000,      2.0,    0.453,    0.035,     10,    0.07},
    //  {  1.050,      2.0,    0.468,    0.035,     10,    0.07},
    //  {  1.100,      2.0,    0.485,    0.030,     10,    0.07},
    //  {  1.150,      2.0,    0.482,    0.045,     10,    0.07},
      //  {  1.200,      2.0,    0.512,    0.045,     10,    0.06},
    //  {  1.250,      2.0,    0.524,    0.050,     10,    0.06},
    //  {  1.300,      2.0,    0.513,    0.055,     10,    0.06},
      //  {  1.400,      2.0,    0.526,    0.060,     10,    0.06},
    //  {  1.500,      2.0,    0.548,    0.065,     10,    0.05},
      //  {  1.600,      2.0,    0.555,    0.065,     10,    0.05},
    //  {  1.750,      2.0,    0.578,    0.065,     10,    0.05},
    //  {  2.000,      2.0,    0.586,    0.065,     10,    0.05},
      //  {  2.250,      2.0,    0.631,    0.088,     10,    0.05},
    //  {  2.500,      2.0,    0.643,    0.020,     10,    0.05},
    //  {  2.750,      2.0,    0.679,    0.065,     10,    0.05},
    //  {  3.000,      2.0,    0.670,    0.100,     10,    0.05},
      //  {  3.500,      2.0,    0.713,    0.100,     10,    0.05},
    //  {  4.000,      2.0,    0.736,    0.110,     10,    0.05},
    //  {  4.500,      2.0,    0.740,    0.031,     10,    0.05},
    //  {  5.000,      2.0,    0.752,    0.032,     10,    0.05},
    //  {  6.000,      2.0,    0.788,    0.120,     10,    0.05},
      //  {  7.000,      2.0,    0.805,    0.010,     10,    0.05},
    //  {  8.000,      2.0,    0.817,    0.005,     10,    0.05},
    //  {  8.500,      2.0,    0.817,    0.005,     10,    0.05},
    //  {  9.000,      2.0,    0.844,    0.005,     10,    0.05},
    //  {  9.500,      2.0,    0.844,    0.005,     10,    0.05},
    //  { 10.000,      2.0,    0.855,    0.005,     10,    0.05},
    //  { 10.500,      2.0,    0.855,    0.005,     10,    0.05},
    //  { 11.000,      2.0,    0.855,    0.005,     10,    0.05},
    //  { 11.500,      2.0,    0.855,    0.005,     10,    0.05},
    //  { 11.750,      2.0,    0.855,    0.005,     10,    0.05},
    //  { 12.000,      2.0,    0.882,    0.005,     10,    0.05},
    //  { 12.250,      2.0,    0.882,    0.005,     10,    0.05},
    //  { 12.500,      2.0,    0.882,    0.005,     10,    0.10},
    //  { 12.750,      2.0,    0.882,    0.005,     10,    0.10},
      //  { 13.000,      2.0,    0.882,    0.005,     10,    0.10},
      // { 13.150,      2.0,    0.925,    0.005,     10,    0.10},
    //  { 13.500,      2.0,    1.050,    0.005,     10,    0.10},
    //  { 13.750,      2.0,    1.050,    0.005,     10,    0.10},
    //  { 14.000,      2.0,    1.050,    0.005,     10,    0.10},
      //  { 14.250,      2.0,    1.050,    0.005,     10,    0.10},
      //  { 14.500,      2.0,    1.050,    0.005,     10,    0.10},
    //  { 14.750,      2.0,    1.100,    0.005,     10,    0.10},
    //  { 15.000,      2.0,    0.950,    0.005,     10,    0.10},
    //  { 15.250,      2.0,    1.100,    0.005,     10,    0.10},
    //  { 15.500,      2.0,    1.100,    0.005,     10,    0.10},
      //  { 15.750,      2.0,    1.100,    0.005,     10,    0.10},
    //  { 16.000,      2.0,    1.100,    0.005,     10,    0.10},    
    //  { 17.000,      2.0,    1.100,    0.005,     10,    0.10},    
    //  { 18.000,      2.0,    1.150,    0.005,     10,    0.10},
    //  { 19.000,      2.0,    1.150,    0.005,     10,    0.10},
      //  { 20.000,      2.0,    1.150,    0.005,     10,    0.10},
    //  { 21.000,      2.0,    1.150,    0.005,     10,    0.10},
    //  { 22.000,      2.0,    1.150,    0.005,     10,    0.10},
    //  { 23.000,      2.0,    1.150,    0.005,     10,    0.10},
    //  { 24.000,      2.0,    1.200,    0.005,     10,    0.10},
    //  { 25.000,      2.0,    1.200,    0.005,     10,    0.10},
      //  { 26.000,      2.0,    1.200,    0.005,     10,    0.10},
    //  { 27.000,      2.0,    1.200,    0.005,     10,    0.10},
    //  { 28.000,      2.0,    1.200,    0.005,     10,    0.10},
      //  { 32.000,      2.0,    1.200,    0.005,     10,    0.10},
      //  { 40.000,      2.0,    1.200,    0.005,     10,    0.10},
    //  { 50.000,      2.0,    1.200,    0.005,     10,    0.10},
    //  { 60.000,      2.0,    1.300,    0.005,     10,    0.10},
      //  { 70.000,      2.0,    1.300,    0.005,     10,    0.10},
    //  { 80.000,      2.0,    1.300,    0.005,     10,    0.10},
      //  { 120.000,     2.0,    1.500,    0.005,     10,    0.10},
    //

    //  { 12.750,      2.0,    0.882,    0.020,     10,    0.04},

    //  { 10.000,      2.0,    0.9,    0.005,     10,    0.10},


    //  {  2.00/1.1,       1.0,    0.900,    0.02,      10,    0.03},

    //  {  0.0011*2.0/1.1,     2,  0.001,    0.800,     10,    2*0.30},
    //  {  0.0105059*2.0/1.1,  2,  0.010,    0.600,     10,    2*0.25},
    //  {  0.0650441*2.0/1.1,  2,  0.100,    0.300,     10,    2*0.12},
    //  {  0.0619119*2.0/1.1,  2,  0.200,    0.070,     10,    2*0.09},
    //  {  0.0353067*2.0/1.1,  2,  0.300,    0.065,     10,    2*0.08},
    //  {  0.0069784*2.0/1.1,  2,  0.400,    0.060,     10,    2*0.07},
    //  {  0.0001*2.0/1.1,     2,  0.500,    0.055,     10,    2*0.06},
    //  {  0.0618*2.0/1.1,     2,  0.600,    0.050,     10,    2*0.05},
    //  {  0.528*2.0/1.1,      2,  0.700,    0.040,     10,    2*0.05},
    //  {  1.76*2.0/1.1,       2,  0.800,    0.025,     10,    2*0.05},
    //  {  4.84*2.0/1.1,       2,  1.000,    0.02,      10,    2*0.03},
    //  { 21.6*2.0/1.1,        2,  1.200,    0.02,      10,    2*0.02},

    //  {  2.00*2.0/1.1,       2,  0.900,    0.02,      10,    2*0.03},

    //  {0.000994868,	1.0, 1.00e-03, 0.1, 10, 0.1 }, 
    //  {0.00946293,	1.0, 1.00e-02, 0.1, 10, 0.1 }, 
    //  {0.03705045,	1.0, 5.00e-02, 0.1, 10, 0.1 }, 
    //  {0.0440148,	    1.0, 1.00e-01, 0.1, 10, 0.1 }, 
    //  {0.03948945,	1.0, 1.50e-01, 0.1, 10, 0.1 }, 
    //  {0.0278508,	    1.0, 2.00e-01, 0.1, 10, 0.1 }, 
    //  {0.017992975,	1.0, 2.50e-01, 0.1, 10, 0.1 }, 
    //  {-0.01238076,	1.0, 3.00e-01, 0.1, 10, 0.1 }, 
    //  {-0.04029655,	1.0, 3.50e-01, 0.1, 10, 0.1 }, 
    //  {-0.0758,	    1.0, 4.00e-01, 0.1, 10, 0.1 }, 
    //  {-0.12302775,	1.0, 4.50e-01, 0.1, 10, 0.1 }, 
    //  {-0.1726835,	1.0, 5.00e-01, 0.1, 10, 0.1 }, 
    //  {-0.18889915,	1.0, 5.50e-01, 0.1, 10, 0.1 }, 

      //  { 1.60, 1.0, 0.9, 0.005, 10, 0.10},
      //  { 1.80, 1.0, 0.9, 0.005, 10, 0.10},
      //  { 1.90, 1.0, 0.9, 0.005, 10, 0.10},
      //  { 1.92, 1.0, 0.9, 0.005, 10, 0.10},
      //  { 1.93, 1.0, 0.9, 0.005, 10, 0.10},
      //  { 1.94, 1.0, 0.9, 0.005, 10, 0.10},
      //  { 1.95, 1.0, 0.9, 0.005, 10, 0.10},
      //  { 2.00, 1.0, 0.9, 0.005, 10, 0.10},
      //  { 2.05, 1.0, 0.9, 0.005, 10, 0.10},
      //  { 2.10, 1.0, 0.9, 0.005, 10, 0.10},


      //  { 13.5, 2.0, 0.97, 0.005, 10, 0.10},
      //  { 13.7, 2.0, 0.97, 0.005, 10, 0.10},
      //  { 13.9, 2.0, 0.97, 0.005, 10, 0.10},
      //  { 13.2, 2.0, 0.97, 0.005, 10, 0.10},
      //  { 13.5, 2.0, 0.97, 0.005, 10, 0.10},
      //  { 13.8, 2.0, 0.97, 0.005, 10, 0.10},
      //  { 13.9, 2.0, 0.97, 0.005, 10, 0.10},
      //  { 14.0, 2.0, 0.97, 0.005, 10, 0.10},
      //  { 14.1, 2.0, 0.97, 0.005, 10, 0.10},
      //  { 14.2, 2.0, 0.97, 0.005, 10, 0.10},
      //  { 14.3, 2.0, 0.97, 0.005, 10, 0.10},


     // DA QUI...
     // { 0.001, 1.0, 0.001, 0.005, 10, 0.10},
     // { 0.001, 1.0, 1.000, 0.005, 10, 0.10},
     // { 0.002, 1.0, 0.001, 0.005, 10, 0.10},
     // { 0.002, 1.0, 1.000, 0.005, 10, 0.10},

     // { 0.010, 1.0,  0.001, 0.005, 10, 0.10},
     // { 0.010, 1.0,  0.2  , 0.005, 10, 0.10},
     // { 0.010, 1.0,  1.0  , 0.005, 10, 0.10},
     // { 0.020, 1.0,  0.001, 0.005, 10, 0.10},
     // { 0.020, 1.0,  0.2  , 0.005, 10, 0.10},
     // { 0.020, 1.0,  1.0  , 0.005, 10, 0.10},
     // { 0.040, 1.0,  0.001, 0.005, 10, 0.10},
     // { 0.040, 1.0,  0.2  , 0.005, 10, 0.10},
     // { 0.040, 1.0,  1.0  , 0.005, 10, 0.10},
     // { 0.080, 1.0,  0.001, 0.005, 10, 0.10},
     // { 0.080, 1.0,  0.2  , 0.005, 10, 0.10},
     // { 0.080, 1.0,  1.0  , 0.005, 10, 0.10},
     // { 0.160, 1.0,  0.001, 0.005, 10, 0.10},
     // { 0.160, 1.0,  0.2  , 0.005, 10, 0.10},
     // { 0.160, 1.0,  1.0  , 0.005, 10, 0.10},

     // { 0.2, 1.0, 1.0, 0.005, 10, 0.10},
     // { 0.4, 1.0, 1.0, 0.005, 10, 0.10},
     // { 0.6, 1.0, 1.0, 0.005, 10, 0.10},
     // { 0.8, 1.0, 1.0, 0.005, 10, 0.10},
     // { 1.0, 1.0, 1.0, 0.005, 10, 0.10},
     // { 1.2, 1.0, 1.0, 0.005, 10, 0.10},
     // { 1.4, 1.0, 1.0, 0.005, 10, 0.10},
     // { 1.6, 1.0, 1.0, 0.005, 10, 0.10},
     // { 1.8, 1.0, 1.0, 0.005, 10, 0.10},
     // { 2.0, 1.0, 1.0, 0.005, 10, 0.10},
     // { 2.2, 1.0, 1.0, 0.005, 10, 0.10},
     // { 2.4, 1.0, 1.0, 0.005, 10, 0.10},

     // { 0.2, 1.0, 0.6, 0.005, 10, 0.10},
     // { 0.4, 1.0, 0.6, 0.005, 10, 0.10},
     // { 0.6, 1.0, 0.6, 0.005, 10, 0.10},
     // { 0.8, 1.0, 0.6, 0.005, 10, 0.10},
     // { 1.0, 1.0, 0.6, 0.005, 10, 0.10},
     // { 1.2, 1.0, 0.6, 0.005, 10, 0.10},
     // { 1.4, 1.0, 0.6, 0.005, 10, 0.10},
     // { 1.6, 1.0, 0.6, 0.005, 10, 0.10},
     // { 1.8, 1.0, 0.6, 0.005, 10, 0.10},
     // { 2.0, 1.0, 0.6, 0.005, 10, 0.10},
     // { 2.2, 1.0, 0.6, 0.005, 10, 0.10},
     // { 2.4, 1.0, 0.6, 0.005, 10, 0.10},

     // { 2.6, 1.0, 0.6, 0.005, 10, 0.10},
     // { 2.8, 1.0, 0.6, 0.005, 10, 0.10},
     // { 3.2, 1.0, 0.6, 0.005, 10, 0.10},
     // { 3.6, 1.0, 0.6, 0.005, 10, 0.10},
     // { 4.0, 1.0, 0.6, 0.005, 10, 0.10},
     // { 6.0, 1.0, 0.6, 0.005, 10, 0.10},
     // { 8.0, 1.0, 0.6, 0.005, 10, 0.10},
     // { 16.0, 1.0, 0.6, 0.005, 10, 0.10},
     // { 32.0, 1.0, 0.6, 0.005, 10, 0.10},
     // { 64.0, 1.0, 0.6, 0.005, 10, 0.10},

     // { 0.2, 1.0, 1.0, 0.005, 10, 0.10},
     // { 0.4, 1.0, 1.0, 0.005, 10, 0.10},
     // { 0.6, 1.0, 1.0, 0.005, 10, 0.10},
     // { 0.8, 1.0, 1.0, 0.005, 10, 0.10},
     // { 1.0, 1.0, 1.0, 0.005, 10, 0.10},
     // { 1.2, 1.0, 1.0, 0.005, 10, 0.10},
     // { 1.4, 1.0, 1.0, 0.005, 10, 0.10},
     // { 1.6, 1.0, 1.0, 0.005, 10, 0.10},
     // { 1.8, 1.0, 1.0, 0.005, 10, 0.10},
     // { 2.0, 1.0, 1.0, 0.005, 10, 0.10},
     // { 2.2, 1.0, 1.0, 0.005, 10, 0.10},
     // { 2.4, 1.0, 1.0, 0.005, 10, 0.10},

     // { 2.6, 1.0, 1.0, 0.005, 10, 0.10},
     // { 2.8, 1.0, 1.0, 0.005, 10, 0.10},
     // { 3.2, 1.0, 1.0, 0.005, 10, 0.10},
     // { 3.6, 1.0, 1.0, 0.005, 10, 0.10},
     // { 4.0, 1.0, 1.0, 0.005, 10, 0.10},
     // { 6.0, 1.0, 1.0, 0.005, 10, 0.10},
     // { 8.0, 1.0, 1.0, 0.005, 10, 0.10},
     // A QUI... simulare nolog


     // { 0.0450, 2.0,  0.001, 0.005, 10, 0.10}, // un picco solo
     // { 0.0010, 2.0,  0.001, 0.005, 10, 0.10}, // un picco solo
     // { 0.6, 2.0,  1.0  , 0.005, 10, 0.10}, // un picco solo

     // { 2.0, 1.0, 0.84, 0.005, 10, 0.10},
     // { 2.2, 1.0, 0.84, 0.005, 10, 0.10},
     // { 2.4, 1.0, 0.85, 0.005, 10, 0.10},
     // { 2.6, 1.0, 0.85, 0.005, 10, 0.10},
     // { 2.8, 1.0, 0.86, 0.005, 10, 0.10},
     // { 3.2, 1.0, 0.87, 0.005, 10, 0.10},
     // { 3.6, 1.0, 0.88, 0.005, 10, 0.10},
     // { 4.0, 1.0, 0.89, 0.005, 10, 0.10},
     
     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // breaks the loop

    };

    FILE* file_averages = fopen("data/averages.dat", "w+");
 
    Particle* S0 = malloc(N * sizeof(Particle));
    
    if(SIM == 1) histogram = malloc(BINS * sizeof(double));
    if(SIM == 1) density_histogram = malloc(DENSITY_BINS * sizeof(int));

    unsigned int steps = 1e7;
    for(int k = 0; k < SIM; k++){
        if(simulations[k][1] == 0.0) break;
        printf("\n(%d/%d) ", k + 1, SIM);

        // Set file to which measurements are written (measure_metropolis)
        // char file_measures_name[40]; 
        // sprintf(file_measures_name, "data/measures_%02d.dat", k);
        // file_measures = fopen(file_measures_name, "w+");
        file_measures = NULL;

        // Reset histogram
        if(histogram) memset(histogram, 0, BINS * sizeof(double)); 
        // if(density_histogram) memset(density_histogram, 0, DENSITY_BINS * sizeof(int)); 

        // Initialize T, delta_s, delta_V, delta_l, L
        initialize_constants(simulations[k]); 
        // Place particles in the grid
        initialize_positions(S0);

        Measure* measures = metropolis(S0, steps);

        // Calculate averages and produce plateau graph for error
        Measure avg; Measure var; double cp;
        averages_and_plateau(measures, steps, k, &avg, &var, &cp);
        fprintf(file_averages, "%d %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e\n", 
            k, avg.temp, avg.comp, avg.ener, avg.vol, avg.rho, P, cp
        );

        free(measures);
        fclose(file_measures);

    }

    fclose(file_averages);

    // Finalize and export distribution histogram
    if(histogram) {
        FILE* histogram_file = fopen("data/histogram.dat", "w+");
        save_distribution_histogram(histogram_file, steps);
        fclose(histogram_file);
    }

    // Export density histogram
    if(density_histogram) {
        FILE* density_histogram_file = fopen("data/density_histogram.dat", "w+");
        save_density_histogram(density_histogram_file);
        fclose(density_histogram_file);
    }

    // Scatter plot of particles' positions
    if(SIM == 1) save_scatter_plot(S0);

    return 0;
}