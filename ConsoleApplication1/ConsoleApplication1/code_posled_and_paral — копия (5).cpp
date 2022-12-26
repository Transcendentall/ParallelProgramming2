#include <chrono>
#include <thread>
#include <vector>
#include <iostream>
#include <omp.h>
#include <algorithm>
#include <memory>
#include <mutex>
#include <immintrin.h>
#include <cstdio>
#include <cstdlib>
#include <iomanip>

#define CACHE_LINE 64u

#define N (10000000)

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

double average_par_1(const double* v, size_t m) {
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
		for (size_t i = t; i < m; i += T) {
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

double average_par_2(const double* v, size_t n) {
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

		size_t n_t, i_0;

		if (t < n % T) {
			n_t = n / T + 1;
			i_0 = n_t * t;
		}
		else {
			n_t = n / T;
			i_0 = t * (n / T) + (n % T);
		}

		for (size_t i = i_0; i < n_t; ++i) {
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









class lambda_1
{
public:
	std::mutex& mtx1;
	std::mutex& mtx2;
	int& x1, x2;
	void operator()()
	{
		x1 += 1;
		x2 += 1;
	}
	lambda_1(std::mutex& m1, int n1, std::mutex& m2, int n2)
		: mtx1(m1), mtx2(m2), x1(n1), x2(n2)
	{}
};


double average_cs_cpp(double* V, size_t n)
{
	// Взаимоблокировка (deadlock)
	std::mutex mtx1, mtx2;
	int x1 = 0, x2 = 0;
	auto tp1 = [&mtx1, &x1, &mtx2, &x2]()
	{
		x1 += 1;
		x2 += 1;
	};
	auto tp2 = [&mtx2, &x2, &mtx1, &x1]()
	{
		x2 -= 1;
		x1 -= 1;
	};
	auto th1 = std::thread(tp1), th2 = std::thread(tp2);
	th1.join();
	th2.join();
	std::cout << "x1 = " << x1 << ", x2 = " << x2 << "\n"; // "x1 = 0, x2 = 0"
	return 0;
}

/*
int main() {
	std::cout << "AveragePar1:" << std::endl;
	measure_scalability(average_par_1);
	std::cout << "AveragePar2:" << std::endl;
	measure_scalability(average_par_2);



	auto P = omp_get_num_procs();
	auto partial_res = std::make_unique<result_t[]>(P);
	auto v = std::make_unique<double[]>(N);
	fillVector(v.get(), N);
	std::cout << "AverageCsCpp:" << std::endl;
	average_cs_cpp(v.get(), N);
}
*/





/*
for (i = i_0; i < n_t; ++i)
	local_sum += V[i];
#pragma omp critical
{
	result += local_sum;
}
}
return result / N;
}
*/



double average(const double* V, size_t m) {
	double result = 0;
	int n = 100;

#pragma omp parallel shared(result)
	{
		unsigned T = omp_get_num_threads();
		unsigned t = omp_get_thread_num();
		size_t n_t = n / T, i_0 = n % T;
		if (t < i_0)
			i_0 = ++n_t * t;
		else
			i_0 += t * n_t;
	}
}



int main()
{
	//	std::size_t N = 1000000;
	auto m = std::make_unique<double[]>(N);
	std::generate_n(&m[0], N, []() {
		static int i;
		return i++;  });
	std::cout << " Average value: "
		<< average(m.get(), N) << "\n";
	return 0;

}







// 07.11.2022



extern int global;
std::mutex mtx;
void thread_proc(int x)
{
	int x_sq = x * x;
	x_sq = x_sq + x_sq + x * x;
	{
		std::scoped_lock lck{ mtx };
		++global;
		// memory order acquire и release
	}
	return;
}



queue q;
auto reader()
{
	while (q_empty)
		continue;
	int val = q.pop();
	// memory order (fence?) acquire
	q.empty = true;

	return val;
}

auto writer(int x)
{
	q.push(x);
	// memory order release
	q_empty = false;

	q.push(x);
}


// xor в C/C++ называется ^

// 1) используя мьютекс



// 2) используя атомарность (в OpenMP : pragma omp atomic { атомарная операция })

// 2.1) в OpenMP:
/*
#pragma omp atomic
{
	r^ = vi;
	// атомарная операция
}
*/

// 2.2) в C++:
std::atomic<std::size_t> r;
//fetch_add(10);
//memory_order_seq_cst(memory_order-acq_rel


//


unsigned chechsum(const unsigned* v, std::size_t n) {
	unsigned total_sum = 0;
	size_t nt, i0;
#pragma omp parallel {


	unsigned T = omp_get_num_threads();
	unsigned t = omp_get_thread_num();
	if (t < n % T) {
		nt = n / T + 1;
		i0 = nt * t;
	}
	else {
		nt = n / T;
		i0 = t - (n / T) + (n % T);
	}

	unsigned p_sum = 0;
	for (size_t i = i0; i < nt; i++) {
		p_sum ^= v[i];
	}
#pragma omp atomic {
	total_sum ^= p_sum;
	}
	}
	return total_sum;
}





// через mutex (без atomic)
unsigned checkSumMutex(const unsigned* v, size_t n) {
	unsigned total_sum = 0;
	std::mutex mtx;
	std::vector<std::thread> workers;
	auto worker = [&total_sum, &mtx, v, n](unsigned t) {
		unsigned local_sum = 0;
		unsigned T = omp_get_num_threads();
		size_t nt, i0;
		if (t < n % T) {
			nt = n / T + 1;
			i0 = nt * t;
		}
		else {
			nt = n / T;
			i0 = t - (n / T) + (n % T);
		}

		for (size_t i = i0; i < nt; i++) {
			local_sum ^= v[i];
		}

		std::scoped_lock lock{ mtx };
		total_sum ^= local_sum;
	};

	for (unsigned t = 0; t < get_num_threads(); ++t)
		workers.emplace_back(worker, t);
	worker(0);
	for (auto& w : workers)
		w.join();
	return total_sum;

}





// через atomic (без mutex)


unsigned checkSumAtomic(const unsigned* v, size_t n) {
	std::atomic<unsigned> total_sum{ 0 };
	std::vector<std::thread> workers;
	auto worker = [&total_sum, v, n](unsigned t) {
		unsigned local_sum = 0;
		unsigned T = omp_get_num_threads();
		size_t nt, i0;
		if (t < n % T) {
			nt = n / T + 1;
			i0 = nt * t;
		}
		else {
			nt = n / T;
			i0 = t - (n / T) + (n % T);
		}

		for (size_t i = i0; i < nt; i++) {
			local_sum ^= v[i];
		}

		total_sum.fetch_xor(local_sum, std::memory_order_relaxed);
		total_sum ^= local_sum;
	};

	for (unsigned t = 0; t < get_num_threads(); ++t)
		workers.emplace_back(worker, t);
	worker(0);
	for (auto& w : workers)
		w.join();
	return total_sum;
	// total_sum.load(std::memory_order_release);
	// total_sum.load(std::memory_order_relaxed); - можно использовать в некоторых случаях вместо release
}





/*
// точка входа
for (size_t i = 0; i < n; i++)
	v[i] = 0x12345678;
unsigned cs = checksum(v, n);
std::cout << std::hex << cs << '\n';
*/


class my_barrier
{
	std::condition_variable cv; std::mutex mtx;
	const unsigned T; bool gen;
	unsigned wT = 0;
public:
//	my_barrier(unsigned threads) : T(threads) { }
//	void arrive_and_wait();
};


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

/*
#include <mutex>
class my_barrier { };
# if !defined(__GNUC__) || __GNUC__ >= 11
# include <barrier>

typedef std::barrier barrier;
# else
typedef my_barrier barrier;
# endif

# include <atomic>
int main(int argc, char** argv)
{
	unsigned T = 100;
	std::atomic<int> v(0);
	std::vector<std::thread> thr;

//	for (unsigned t = 0; t<T; ++t)
//		thr.emplace_back([&v]() { ++v; });
//		for (auto&th:thr) th.join();

	barrier bar { T };

	auto proc = [&v, &bar, T] {
	++v;
	// bar.arrive_and_wait();
	if (v == T)
		std::cout << "V==T\n";
	--v;
	};
	for (unsigned t=0; t<T; ++t)
		thr.emplace_back(proc);
	for (auto& th:thr)
		th.join();
	return 0;

}
*/


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