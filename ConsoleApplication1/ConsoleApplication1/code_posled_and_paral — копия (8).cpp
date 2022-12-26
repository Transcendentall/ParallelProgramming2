#include <chrono>
#include <thread>
#include <vector>
#include <iostream>
#include <omp.h>
#include <algorithm>
#include <memory>
#include <immintrin.h>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <barrier>
#include <stdint.h>

typedef int16_t WORD;

#define CACHE_LINE 64u

#include <cstdint>
#include <iostream>


void srand(int seed);
typedef uint64_t word;
word rand(word s);


struct lut_row
{
	WORD A, B;
};



auto get_lut(unsigned T)
{
	WORD A = UINT64_C(6364136223846793005);
	WORD B = 1;
	auto lut = std::make_unique<lut_row[]>(T + 1);
	lut[0].A = 1;
	lut[0].B = 0;
	for (unsigned i = 1; i < T; i++)
	{
		lut[i].A = lut[i - 1].A * A;
		lut[i].B = A * lut[i - 1].B + B;
		return lut;
	}
}



double randomize(WORD seed, int* V, size_t n, int a, int b)
{
	WORD s0 = seed;
	size_t T = omp_get_num_procs();
	// static const auto y = get_lut(omp_get_num_procs());
	static auto lut = get_lut(T);
# pragma omp parallel
	{
		size_t t = omp_get_thread_num();
		size_t st = lut[t].A * s0 + lut[t].B;
		for (unsigned k = t; k < n; k += T)
		{
			V[k] = a + int(st % (b - a + 1));
			st = lut[T].A * st + lut[T].B;
		}
	}
	auto get_lut(unsigned T);
}


int main() {
	size_t n = 100;
	int* V = new int[n];
	randomize(456, V, n, 0, 1);
	for (int i = 0; i < n; i++) {
		std::cout << V[i];
	}

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