// Define global parameters of the program
bool print_all;
bool output_discrete_steps;
bool output_k; 

ofstream outputFile; // output file stream 
ifstream inputFile; 
string outputFilename;

// Parameters read in when program is started
int L_energy, L_mu, L_kx, L_time; 
double t_min, t_max;
double A0, A0x, A0y, delta; 
double tau, sigma, omega_pump;

double J, filling, V0, omega_D; 
double temperature_init;
double T1, T2;

// Parameters calculated from initialized parameters
int N; // total number of points in BZ
double mu; // chemical potential
double rho0; // DOS integrated over region [-omegaD, omegaD] around mu
double Delta0; // initial gap

double *spins; // table that stores the values of the spins

