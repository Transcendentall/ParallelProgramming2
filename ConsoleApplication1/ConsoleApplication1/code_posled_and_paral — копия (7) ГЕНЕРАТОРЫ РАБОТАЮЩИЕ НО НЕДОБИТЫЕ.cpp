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