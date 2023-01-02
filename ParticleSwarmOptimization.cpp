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

    #if SERIAL
        cout << "pso_serial" << endl;
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
