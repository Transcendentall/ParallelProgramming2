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

#define CACHE_LINE 64u

/*
class my_barrier
{
	std::condition_variable cv; std::mutex mtx;
	const unsigned T; bool gen;
	unsigned wT = 0;
public:
	my_barrier(unsigned threads) : T(threads) { }
//	void arrive_and_wait();
};
*/

/*
class condition_variable
{
public:
	void wait(std::unique_lock<std::mutex> & lock);
	template <class Pr> requires std::is_invokable_r<bool,Pr>
		void wait(std::unique_lock<std::mutex>& lock, Pr pred)
		{
			while (pred())
				this->wait(lock);
		}
};
*/

// _ _ __ ___ _____ _  _  __ __ _ ____  ___ ___ ______  

static unsigned num_threads = std::thread::hardware_concurrency();
void set_num_threads(unsigned T) {
	num_threads = T;
	omp_set_num_threads(T);
}

unsigned get_num_threads() {
	return num_threads;
}


#include <mutex>
class my_barrier 
{
	std::condition_variable cv; std::mutex mtx;
	const unsigned T; bool gen;
	unsigned wT = 0;
public:
	void arrive_and_wait();
};

int N = 1000;

# if !defined(__GNUC__) || __GNUC__ >= 11
# include <barrier>

typedef std::barrier<> barrier;
# else
typedef my_barrier barrier;
# endif

void fillVector(double* v, size_t n) {
	for (size_t i = 0; i < n; ++i) {
		v[i] = 1.0;
	}
}

struct result_t {
	double value, milliseconds;
};

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

# include <atomic>
int main(int argc, char** argv)
{

	std::cout << "AveragePar1:" << std::endl;
	measure_scalability(average_par);



	auto P = omp_get_num_procs();
	auto partial_res = std::make_unique<result_t[]>(P);
	auto v = std::make_unique<double[]>(N);
	fillVector(v.get(), N);
	std::cout << "AverageCsCpp:" << std::endl;

}



/*
void Latch::arrive_and_wait() // чисто одноразовый барьер
{
	std::unique_lock l { mtx };
	if (++wT < T )
	{

	while (wT < T)
		cv.wait(l);
	}
	else
	{
		cv.notify_all();
	}
}
*/


void my_barrier::arrive_and_wait()
{
	std::unique_lock l{ mtx };
	if (++wT < T)
	{
		bool my_gen = gen;
		while (wT < T)
			cv.wait(l);
	}
	else
	{
		cv.notify_all();
		wT = 0;
		gen=!gen;
	}
}



union partial_sum_t {
	double value;
	alignas(double) char pd[64];
};


double average_par (const double *V, size_t n)
{
	std::vector<std::thread> workers;
	auto T = get_num_threads();
	auto partial_sums = std::make_unique<partial_sum_t[]>(T);
	std::vector<std::thread> thr;
	barrier bar(T);
	auto worker = [&partial_sums, T, &bar, V, n](unsigned t)
	{
		double local_sum = 0;
		for (size_t step = 1, next_step = step + step; step < T; step = next_step, next_step += next_step)
		{
			if (t % next_step == 0 && t + step < T)
			{
				partial_sums[t].value += partial_sums[t + step].value;
			}
			bar.arrive_and_wait();
		}
		// partial_sums[t].value = local_sum;
	}; // worker
	for (unsigned t = 0; t < get_num_threads(); ++t)
		workers.emplace_back(worker, t);
	worker(0);
	for (auto& w : workers)
		w.join();

	// цикл join
	// создаём потоки
	// ожидание
	return partial_sums[0].value / n;
}


