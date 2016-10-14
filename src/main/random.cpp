#include <iostream>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include "random.h"
#include <sched.h>
#include <sys/types.h>
#include <unistd.h>

long Random::m_seed = 0;

void Random::set_seed(long seed)
{
	//m_seed = seed*67+43*sched_getcpu()+37*(int)getpid()+time(NULL);
	m_seed = seed;
	srand((unsigned) m_seed);
	srand48(m_seed);
}

long Random::get_seed()
{
	return m_seed;
}

double Random::rnd(double max_val)
{
	return drand48() * max_val;
}

int Random::rnd(int max_val)
{
	assert(max_val > 0);
	int limit = RAND_MAX / max_val;

	for ( ; ; )
	{
		int val = rand();

		// (make sure the random numbers are not biased)

		if (val / max_val < limit)
		{
			return val % max_val;
		}
	}
}

