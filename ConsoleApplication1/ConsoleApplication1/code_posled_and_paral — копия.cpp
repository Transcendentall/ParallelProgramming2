#include <iostream>
#include <chrono>
#include <concepts>
#include <omp.h>
//#include <cstring.h>

using namespace std;

const int N = 100000000;

template<class F> requires std::is_invocable_r_v <F, double>
    double integral_rr(double a, double b, F f)
    {
        double sum = 0, dx = (b - a) / N;
        unsigned P = omp_get_num_procs();
#pragma omp parallel
        {
            unsigned T = omp_get_num_threads();
            unsigned t = omp_get_thread_num();
            for (unsigned k = 0; t + k * T < N; ++k)
                sum += f(a + (t + k * T) * dx);
        }
        return sum * dx;


    }

double integrate_seq(double a, double b, double(*f)(double))
{
    double res = 0;
    double dx = (-a + b) / N;
    for (int i = 0; i < N; ++i)
    {
        res += f(a + i * dx);
    }
    return dx * res;
}

double integrate_par(double a, double b, double(*f)(double))
{
    double res = 0;
    double dx = (-a + b) / N;
#pragma omp parallel for reduction(+:res)
    for (int i = 0; i < N; ++i)
    {
        res += f(a + i * dx);
    }
    return dx * res;
}


double integrate_par_rr_old(double a, double b, double(*f)(double))
{
        double sum = 0, dx = (b - a) / N;
        unsigned P = omp_get_num_procs();
#pragma omp parallel
        {
            unsigned T = omp_get_num_threads();
            unsigned t = omp_get_thread_num();
            for (unsigned k = 0; t + k * T < N; ++k)
                sum += f(a + (t + k * T) * dx);
        }
        return sum * dx;
}

double integrate_par_rr_new(double a, double b, double(*f)(double))
{
    double sum = 0, dx = (b - a) / N;
    unsigned P = omp_get_num_procs();
    double* partial_sums = (double*)calloc(P, sizeof(double));
    sum = 0;
#pragma omp parallel
    {
        unsigned T = omp_get_num_threads();
        ++sum = 1;
        unsigned t = omp_get_thread_num();
        for (unsigned k = 0; t + k * T < N; ++k)
            sum += f(a + (t + k * T) * dx);
    }
    for (unsigned t = 0; t < P; ++t)
        sum += partial_sums[t];
    free(partial_sums);
    return sum * dx;
}


double integrate_par_rr_new_new(double a, double b, double(*f)(double))
{
    double sum = 0, dx = (b - a) / N;
    unsigned P = omp_get_num_procs();
    double* partial_sums = (double*)calloc(P, sizeof(double));
    sum = 0;
#pragma omp parallel
    {
        unsigned T = omp_get_num_threads();
        ++sum = 1;
        unsigned t = omp_get_thread_num();
        for (unsigned k = 0; t + k * T < N; ++k)
            sum += f(a + (t + k * T) * dx);
    }
    for (unsigned t = 0; t < P; ++t)
        sum += partial_sums[t];
    free(partial_sums);
    return sum * dx;
}


struct result_t
{
    double value, milliseconds;
};

// result_t run_experiment(double(*integrate)(double a, double b, double(*f)(double)));

int main()
{
    auto f = [](double x) { return x * x; };

    //   auto r_seq=run_experiment(integrate_seq, -1,1,f);
    //   auto r_par=run_experiment(integrate_par, -1,1,f);



    auto tm_seq = std::chrono::steady_clock::now();
    auto r_seq = integrate_seq(-1, 1, f);
    cout << "SEQ Time:    " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - tm_seq).count() << endl;

    auto tm_par = std::chrono::steady_clock::now();
    auto r_par = integrate_par(-1, 1, f);
    cout << "PAR Time:    " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - tm_par).count() << endl;

    auto tm_par_rr_old = std::chrono::steady_clock::now();
    auto r_par_rr_old = integrate_par_rr_old(-1, 1, f);
    cout << "PAR_RR_old Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - tm_par_rr_old).count() << endl;

    auto tm_par_rr_new = std::chrono::steady_clock::now();
    auto r_par_rr_new = integrate_par_rr_new(-1, 1, f);
    cout << "PAR_RR_new Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - tm_par_rr_new).count() << endl;
    /*
    auto tm_par_rr_new_new = std::chrono::steady_clock::now();
    auto r_par_rr_new_new = integrate_par_rr_new(-1, 1, f);
    cout << "PAR_RR_new_new Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - tm_par_rr_new_new).count() << endl;
    */
}
