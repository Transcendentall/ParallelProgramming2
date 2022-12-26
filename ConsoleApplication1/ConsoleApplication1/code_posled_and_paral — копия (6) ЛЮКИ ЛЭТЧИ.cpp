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


#include <mutex>
class Latch
{
	std::condition_variable cv; std::mutex mtx;
	const unsigned T; bool gen;
	unsigned wT = 0;
public:
	void arrive_and_wait();
};



# if !defined(__GNUC__) || __GNUC__ >= 11
# include <barrier>

typedef std::barrier<> barrier;
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

	barrier bar (T);

	auto proc = [&v, &bar, T] {
	++v;
	 bar.arrive_and_wait();
	if (v-- == T)
		std::cout << "V==T\n";
	};
	for (unsigned t=0; t<T; ++t)
		thr.emplace_back(proc);
	for (auto& th:thr)
		th.join();
	return 0;

}




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


/*
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
*/