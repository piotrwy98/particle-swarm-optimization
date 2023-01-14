#include <iostream>
#include <omp.h>

using namespace std;

int getCmdOption(char** begin, char** end, const string& option, int defaultValue)
{
    char** itr = find(begin, end, option);

    if (itr != end && ++itr != end)
    {
        int value = atoi(*itr);
        return value != 0 ? value : defaultValue;
    }

    return defaultValue;
}

void pso(int d, int m, int c1, int c2, int v, int i, int s)
{
    double start = omp_get_wtime();
    double w = 0.5;

    #if SERIAL
        cout << "pso_serial" << endl;

        // do przetrzymywania pozycji per próbka potrzebujemy listy list
        // długość zewnętrznej listy = ilość cząstek m
        // długość wewnętrznej listy = ilość wymiarów d
        // np. [ [x1,y1,z1], [x2,y2,z2], [x3,y3,z3]]

        // do przetrzymywania najlepszych pozycji per próbka potrzebujemy takiej samej tablicy jak ta wyżej

        // do przetrzymywania najlepszych wartości per próbka potrzebujemy tablicy (długość talicy = ilość cząstek m)

        // do przetrzymywania najlepszej wartości globalnej potrzebujemy jednej zmiennej

        // do przetrzymywania najlepszych współrzędnych potrzebujemy listy
        // długość listy = ilość wymiarów d

        // pętla do inicjalizacji pozycji i prędkości
        for (int j = 0; j < m; j++)
        {
            for (int jj = 0; jj < d; jj++)
            {
                // initialize position x_id randomly within permissible range
                // initialize velocity v_id randomly within permissible range

                // tutaj można też wybrać najlepszą pozycję globalną
            }
        }

        int k = 1;
        do
        {
            // ---------------------------------------------------------------------------------------------
            
            // 1 wersja
            for (int j = 0; j < m; j++)
            {
                // obliczenie wartości funkcji
                
                // jeżeli wartość jest mniejsza niż p_best 
                    // ustaw wartość jako p_best
                    // jeżeli wartość jest mniejsza niż g_best  
                        // ustaw wartość jako g_best                
            }

            for (int j = 0; j < m; j++)
            {
                for (int jj = 0; jj < d; jj++)
                {
                    // rand musi być z zakresu (0,1)
                    
                    // update prędkości danej próbki per wymiar
                    // v_id(k+1) = w * v_id(k) + c1 * rand_1(p_best_id - x_id) + c2 * rand2(g_best_d - x_id)
                    
                    // update pozycji danej próbki per wymiar
                    // x_id(k+1) = x_id(k) + v_id(k+1)
                }
            }

            // ---------------------------------------------------------------------------------------------
            
            // 2 wersja
            for (int j = 0; j < m; j++)
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
            }

            // ---------------------------------------------------------------------------------------------

            // increase counter
            k = k + 1;
        }
        while (k <= i);

    #endif

    #if PARALLEL
        cout << "pso_parallel" << endl;
        omp_set_num_threads(d);
    #endif

    double end = omp_get_wtime();
    cout << "# " << end - start << endl;
}

int main(int argc, char* argv[])
{
    int d = getCmdOption(argv, argv + argc, "-D", 2);
    int m = getCmdOption(argv, argv + argc, "-m", 8);
    int c1 = getCmdOption(argv, argv + argc, "-c", 2);
    int c2 = getCmdOption(argv, argv + argc, "-C", 2);
    int v = getCmdOption(argv, argv + argc, "-V", 60);
    int i = getCmdOption(argv, argv + argc, "-i", 1);
    int s = getCmdOption(argv, argv + argc, "-s", 1);

    pso(d, m, c1, c2, v, i, s);
}
