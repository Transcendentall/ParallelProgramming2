#include <chrono>
#include <thread>
#include <vector>
#include <iostream>
#include <omp.h>

#define CACHE_LINE 64u
#define N (100000000)

static unsigned num_threads = std::thread::hardware_concurrency();
struct result_t {
	double value, milliseconds;
};

union partial_sum_t {
	double value;
	alignas(double) char pd[64];
};

void set_num_threads(unsigned T) {
	num_threads = T;
	omp_set_num_threads(T);
}

unsigned get_num_threads() {
	return num_threads;
}

void fillVector(double* v, size_t n) {
	for (size_t i = 0; i < n; ++i) {
		v[i] = 1.0;
	}
}

result_t
run_experiment(double (*average)(const double*, size_t), const double* v, size_t n) {
	auto tm1 = std::chrono::steady_clock::now();
	double value = average(v, n);
	auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - tm1).count();
	result_t res{ value, (double)time };
	return res;
}

void measure_scalability(auto averageFunction) {
	auto P = omp_get_num_procs();
	auto partial_res = std::make_unique<result_t[]>(P);
	auto v = std::make_unique<double[]>(N);
	fillVector(v.get(), N);
	for (auto T = 1; T <= P; ++T) {
		set_num_threads(T);
		partial_res[T - 1] = run_experiment(averageFunction, v.get(), N);
		auto speedup = partial_res[0].milliseconds / partial_res[T - 1].milliseconds;
		std::cout << "Threads = " << T << std::endl;
		std::cout << "Time    = " << partial_res[T - 1].milliseconds << std::endl;
		std::cout << "Value   = " << partial_res[T - 1].value << std::endl;
		std::cout << "Speedup = " << speedup << std::endl << std::endl;
	}
}

double average_par_static(const double* v, size_t m) {
	unsigned T;
	partial_sum_t* sums;
	double result = 0;
#pragma omp parallel shared(T, sums)
	{
		unsigned t = omp_get_thread_num();
		double local_sum = 0;
#pragma omp single
		{
			T = (unsigned)omp_get_num_threads();
			sums = (partial_sum_t*)malloc(T * sizeof(partial_sum_t));
		}
#pragma omp parallel for schedule(static)
		for (long long i = t; i < m; i += T) {
			local_sum = v[i];
		}
		sums[t].value = local_sum;
	}

	for (size_t i = 0; i < T; ++i) {
		result += sums[i].value;
	}

	free(sums);
	return result;
}

double average_par_dynamic(const double* v, size_t m) {
	unsigned T;
	partial_sum_t* sums;
	double result = 0;
#pragma omp parallel shared(T, sums)
	{
		unsigned t = omp_get_thread_num();
		double local_sum = 0;
#pragma omp single
		{
			T = (unsigned)omp_get_num_threads();
			sums = (partial_sum_t*)malloc(T * sizeof(partial_sum_t));
		}
#pragma omp parallel for schedule(dynamic)
		for (long long i = t; i < m; i += T) {
			local_sum = v[i];
		}
		sums[t].value = local_sum;
	}

	for (size_t i = 0; i < T; ++i) {
		result += sums[i].value;
	}

	free(sums);
	return result;
}

int main() {
	std::cout << "AverageParStatic:" << std::endl;
	measure_scalability(average_par_static);
	std::cout << "AverageParDynamic:" << std::endl;
	measure_scalability(average_par_dynamic);
}



/*

Записи 12.12.2022

double randomize(word seed, int *V, size_t n, int a, int b)
{

	word s0 = seed;
	static const auto lut = get_lut(std::thread::hardware_concurrency());
	auto T = get_num_threads();
	auto partial_sums = std::make_unique<partial_sum_t[]>(T);
	std::vector<std::thread> thr;



	auto worker = [T,n,a,b] (unsigned t)
	{
		int local_sum = 0;
		size_t st = lut[t].a * s0 + lut[t].b;
		for (unsigned k=t; k<n; k+=t)
		{
			V[k] = a + int(st % (b-a+1));
			local_sum += V[k];
			st = lut[t].a * st + lut[t].b;
		}
		partial sums[t].value = local_sum;
	};




	for (size_t step = 1, n_step = step * 2; step < T
	}



	for (unsigned t=1; t < get_num_threads(); ++t)
	{
		thr.emplace_back(worker, t);
	}
	worker(0);
	for (auto &w; thr)
	{
		w.join();
	}
	return partial_sums[0].value / n;
}




*/

/*

Динамический параллелизм (параллелизм на задачах, task-based parallelism)

_dyn

// # pragma omp parallel schedule (static)
// # pragma omp parallel schedule (dynamic)

double* average_stat(const double *V, size_t n)
{
	double s=0;
	// # pragma omp parallel for schedule (static)
	for (size_t i=0; i<n; ++i)
		s+=V[i];
	return s/n;
}



*/