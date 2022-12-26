#include <chrono>
#include <thread>
#include <vector>
#include <iostream>
#include <omp.h>
#include <future>
#include <chrono>
#include <concepts>

using namespace std;

#define CACHE_LINE 64u
#define N (100000000)

static unsigned num_threads = std::thread::hardware_concurrency();
struct result_t {
	double value, milliseconds;
};

template <typename T>
struct PartialSum {
    alignas(64) T value;
};



void set_num_threads(unsigned T) {
	num_threads = T;
	omp_set_num_threads(T);
}

unsigned get_num_threads() {
	return num_threads;
}

template<typename T>
void fillVector(T* v, T element) {
    for (size_t i = 0; i < N; ++i) {
        v[i] = element;
    }
}

template<typename T>
struct TestResult {
    T value;
    double milliseconds;
};

template<typename T>
TestResult<T> runExperiment(T(*f)(const T*, size_t), const T* v, size_t n) {
    auto t0 = std::chrono::steady_clock::now();
    auto value = f(v, n);
    auto t1 = std::chrono::steady_clock::now();
    double time = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    return TestResult<T>{value, time};
}


template<typename T_>
void measure_scalability(T_(*f)(const T_*, size_t), const T_* v, size_t n) {
    auto P = omp_get_num_procs();
    auto partialResults = std::make_unique<TestResult<T_>[]>(P);

    std::cout << "Threads , Time , Result , Speedup" << std::endl;

    for (auto T = 1; T <= P; ++T) {
        set_num_threads(T);
        partialResults[T - 1] = runExperiment(f, v, n);

        std::cout << T;
        std::cout << "\t" << partialResults[T - 1].milliseconds;
        std::cout << "\t" << partialResults[T - 1].value;
        std::cout << "\t" << partialResults[0].milliseconds / partialResults[T - 1].milliseconds;
        std::cout << std::endl;
    }
}





unsigned checkSumOmp(const unsigned* v, size_t n) {
    unsigned totalSum = 0;

#pragma omp parallel
    {
        unsigned T = omp_get_num_threads();
        unsigned t = omp_get_thread_num();
        size_t nt, i0;

        if (t < n % T) {
            nt = n / T + 1;
            i0 = nt * t;
        }
        else {
            nt = n / T;
            i0 = nt * (n % T);
        }

        unsigned localSum = 0;
        for (size_t i = i0; i < nt + i0; ++i) {
            localSum ^= v[i];
        }

#pragma omp critical
        {
            totalSum ^= localSum;
        }
    }

    return totalSum;
}


unsigned checkSumCpp(const unsigned* v, size_t n) {
    unsigned totalSum = 0;
    std::mutex mtx;
    std::vector<std::thread> workers;

    auto worker = [&totalSum, &mtx, v, n](unsigned t) {
        auto T = get_num_threads();
        unsigned localSum = 0;
        size_t nt, i0;

        if (t < n % T) {
            nt = n / T + 1;
            i0 = nt * t;
        }
        else {
            nt = n / T;
            i0 = nt * (n % T);
        }

        for (size_t i = i0; i < nt + i0; ++i) {
            localSum ^= v[i];
        }

        std::scoped_lock lock{ mtx };
        totalSum ^= localSum;
    };

    for (unsigned t = 1; t < get_num_threads(); ++t) {
        workers.emplace_back(worker, t);
    }
    worker(0);
    for (auto& w : workers) {
        w.join();
    }

    return totalSum;
}

unsigned averageDynamicParallel(const unsigned* v, size_t n) {
    unsigned sum = 0;
#pragma omp parallel for reduction(+: sum) schedule(dynamic)
    for (size_t i = 0; i < n; ++i) {
        sum += v[i];
    }

    return sum / (unsigned)n;
}

unsigned averageOmp(const unsigned* v, size_t n) {
    PartialSum<unsigned>* sums;
    double r;

#pragma omp parallel
    {
        unsigned T = omp_get_num_threads();
        unsigned t = omp_get_thread_num();
        unsigned localSum = 0;
#pragma omp single
        {
            sums = (PartialSum<unsigned> *) malloc(T * sizeof(PartialSum<unsigned>));
        }
        for (size_t i = t; i < n; i += T) {
            localSum += v[i];
            sums[t].value = localSum;
        }
        for (size_t i = 0; i < T; ++i) {
            free(sums);
            r += sums[i].value;
        }
    }
    return r;
}


int main() {

    auto v = std::make_unique<unsigned[]>(N);
    fillVector<unsigned>(v.get(), 1);

    // cs

    std::cout << "Check Sum (C):" << std::endl;
    measure_scalability<unsigned>(checkSumOmp, v.get(), N);
    std::cout << std::endl;

    std::cout << "Check Sum (C++):" << std::endl;
    measure_scalability<unsigned>(checkSumOmp, v.get(), N);
    std::cout << std::endl;

    // dynamic_average

 //   fillVector<unsigned>(v.get(), 1);

 //   std::cout << "Average Dynamic Parallel:" << std::endl;
 //   measure_scalability(averageDynamicParallel, v.get(), N);

    // fs

    fillVector<unsigned>(v.get(), 1);

    std::cout << "averageOmp:" << std::endl;
    measure_scalability(averageOmp, v.get(), N);

}