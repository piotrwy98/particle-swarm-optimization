#define _USE_MATH_DEFINES

#include <iostream>
#include <omp.h>
#include <math.h>
#include <iomanip>

using namespace std;

double LOWER_BOUND = -32.768;
double UPPER_BOUND = 32.768;

struct Particle
{
    double* position;
    double* best_position;
    double best_value = DBL_MAX;
    double* velocity;
};

int get_cmd_option(char** begin, char** end, const string& option, int defaultValue)
{
    char** itr = find(begin, end, option);

    if (itr != end && ++itr != end)
    {
        int value = atoi(*itr);
        return value != 0 ? value : defaultValue;
    }

    return defaultValue;
}

double ackley(double positions[], int n)
{
    double firstPart = 0.0;
    double secondPart = 0.0;

    #if SERIAL
        for (int i = 0; i < n; i++)
        {
            firstPart += positions[i] * positions[i];
            secondPart += cos(2 * M_PI * positions[i]);
        }
    #endif

    #if PARALLEL
        #pragma omp parallel for reduction(+:firstLoopResult, secondLoopResult)
        for (int i = 0; i < n; i++)
        {
            firstPart += positions[i] * positions[i];
            secondPart += cos(2 * M_PI * positions[i]);
        }
    #endif

    return -20.0 * exp(-0.2 * sqrt(1.0 / n * firstPart)) - 
        exp(1.0 / n * secondPart) + 20.0 + DBL_EPSILON;
}

void pso(int d, int m, int c1, int c2, int v, int i, int s, int threads)
{
    double start = omp_get_wtime();
    double w = 0.5;
    int k = 1;
    Particle* particles = new Particle[m];
    double best_general_value = DBL_MAX;
    double* best_general_position = new double[d];

    #if SERIAL
        // pętla do inicjalizacji pozycji i prędkości
        for (int particle = 0; particle < m; particle++)
        {
            particles[particle].position = new double[d];
            particles[particle].velocity = new double[d];

            for (int dimension = 0; dimension < d; dimension++)
            {
                // initialize position x_id randomly within permissible range
                particles[particle].position[dimension] = (UPPER_BOUND - LOWER_BOUND) * ((double) rand() / (double) RAND_MAX) + LOWER_BOUND;

                // initialize velocity v_id randomly within permissible range
                particles[particle].velocity[dimension] = v * ((double) rand() / (double) RAND_MAX);
            }
        }

        do
        {
            // ---------------------------------------------------------------------------------------------
            
            // 1 wersja
            for (int particle = 0; particle < m; particle++)
            {
                double value = ackley(particles[particle].position, d);
                
                if (value < particles[particle].best_value)
                {
                    particles[particle].best_value = value;
                    particles[particle].best_position = particles[particle].position;

                    if (value < best_general_value)
                    {
                        best_general_value = value;
                        best_general_position = particles[particle].position;
                    }
                }             
            }

            for (int particle = 0; particle < m; particle++)
            {
                for (int dimension = 0; dimension < d; dimension++)
                {
                    double rand_1 = ((double) rand() / (double) RAND_MAX);
                    double rand_2 = ((double) rand() / (double) RAND_MAX);
                    
                    // update prędkości danej próbki per wymiar
                    // v_id(k+1) = w * v_id(k) + c1 * rand_1 * (p_best_id - x_id) + c2 * rand2 * (g_best_d - x_id)
                    
                    particles[particle].velocity[dimension] = w * particles[particle].velocity[dimension] + 
                        c1 * rand_1 * (particles[i].best_value - particles[particle].position[dimension]) + 
                        c2 * rand_2 * (best_general_value - particles[particle].position[dimension]);
                    
                    // update pozycji danej próbki per wymiar
                    // x_id(k+1) = x_id(k) + v_id(k+1)
                    particles[particle].position[dimension] += particles[particle].velocity[dimension];
                }
            }

            // ---------------------------------------------------------------------------------------------
            
            // 2 wersja
            /*for (int j = 0; j < m; j++)
            {
                for (int jj = 0; jj < d; jj++)
                {
                    // rand musi być z zakresu (0,1)

                    // update prędkości danej próbki per wymiar
                    // v_id(k+1) = w * v_id(k) + c1 * rand_1(p_best_id - x_id) + c2 * rand2(g_best_d - x_id)
                }

                // update pozycji danej próbki
                // x_i(k+1) = x_i(k) + v_i(k+1)

                // obliczenie wartości funkcji

                // jeżeli wartość jest mniejsza niż p_best 
                    // ustaw wartość jako p_best
                    // jeżeli wartość jest mniejsza niż g_best  
                        // ustaw wartość jako g_best 
            }*/

            // ---------------------------------------------------------------------------------------------

            // increase counter
            k++;
        }
        while (k <= i);
    #endif

    #if PARALLEL
        omp_set_num_threads(threads);

        #pragma omp parallel for
        for (int particle = 0; particle < m; particle++)
        {
            particles[particle].position = new double[d];
            particles[particle].velocity = new double[d];

            #pragma omp parallel for
            for (int dimension = 0; dimension < d; dimension++)
            {
                particles[particle].position[dimension] = (UPPER_BOUND - LOWER_BOUND) * ((double)rand() / (double)RAND_MAX) + LOWER_BOUND;
                particles[particle].velocity[dimension] = v * ((double)rand() / (double)RAND_MAX);
            }
        }

        // czy tutaj też omp? jest jakaś dyrektywa do WHILE czy należy zamienić na FOR?
        do
        {
            #pragma omp parallel for
            for (int particle = 0; particle < m; particle++)
            {
                double value = ackley(particles[particle].position, d);

                if (value < particles[particle].best_value)
                {
                    particles[particle].best_value = value;
                    particles[particle].best_position = particles[particle].position;

                    if (value < best_general_value)
                    {
                        best_general_value = value;
                        best_general_position = particles[particle].position;
                    }
                }
            }

            #pragma omp parallel for
            for (int particle = 0; particle < m; particle++)
            {
                #pragma omp parallel for
                for (int dimension = 0; dimension < d; dimension++)
                {
                    double rand_1 = ((double)rand() / (double)RAND_MAX);
                    double rand_2 = ((double)rand() / (double)RAND_MAX);

                    particles[particle].velocity[dimension] = w * particles[particle].velocity[dimension] + 
                        c1 * rand_1 * (particles[i].best_value - particles[particle].position[dimension]) +
                        c2 * rand_2 * (best_general_value - particles[particle].position[dimension]);
                    particles[particle].position[dimension] += particles[particle].velocity[dimension];
                }
            }

            k++;
        }
        while (k <= i);
    #endif

    double end = omp_get_wtime();
    cout << "# TIME:     " << fixed << setprecision(15) << end - start << endl;
    cout << "  VALUE:    " << best_general_value << endl;
    cout << "  POSITION: (";

    for (int dimension = 0; dimension < d; dimension++)
    {
        cout << best_general_position[dimension] << (dimension != d - 1 ? ", " : "");
    }

    cout << ")" << endl;
}

int main(int argc, char* argv[])
{
    srand(time(0));

    int d = get_cmd_option(argv, argv + argc, "-D", 2);
    int m = get_cmd_option(argv, argv + argc, "-m", 8);
    int c1 = get_cmd_option(argv, argv + argc, "-c", 2);
    int c2 = get_cmd_option(argv, argv + argc, "-C", 2);
    int v = get_cmd_option(argv, argv + argc, "-V", 60);
    int i = get_cmd_option(argv, argv + argc, "-i", 1);
    int s = get_cmd_option(argv, argv + argc, "-s", 1);
    int threads = get_cmd_option(argv, argv + argc, "-t", 6);

    pso(d, m, c1, c2, v, i, s, threads);
}
