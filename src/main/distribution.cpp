
#include <cassert>
#include <iostream>
#include <algorithm>
#include "random.h"
#include "distribution.h"

Distribution::Distribution()
	: m_prob_total(0.0), m_sorted(false)
{
}

void Distribution::clear()
{
	m_val.clear();
	m_sorted = false;
}

void Distribution::add(double probability, int value)
{
	assert(probability > 0.0);
	Dist_Element e(probability, value);
	m_val.push_back(e);
	m_sorted = false;
}

int Distribution::select()
{
	assert(m_val.size() != 0);

	if (!m_sorted)
	{
		sort_values();
	}

	int n = binary_search(Random::rnd(m_prob_total));
	return m_val[n].val;
}

void Distribution::sort_values()
{
	std::sort(m_val.begin(), m_val.end());
	m_prob_total = 0.0;

	for (unsigned int n = 0;n < m_val.size();n++)
	{
		m_val[n].sum_to_here = (m_prob_total += m_val[n].prob);
	}

	m_sorted = true;
}

int Distribution::binary_search(double val) const
{
	int low = -1;
	int high = (int) m_val.size() - 1;

	while (low < high)
	{
		// assert(low >= 0);
		assert(high < (int) m_val.size());

		int m = (low + high) / 2;
		assert(m >= 0);

		if (val < m_val[m].sum_to_here)
		{
			high = m - 1;
		}
		else
		{
			low = m + 1;
		}
	}

	if (low < 0) { return 0; }
	if (val > m_val[low].sum_to_here) { return low + 1; }
	return low;
}

void Distribution::dump(std::ostream &out /*= std::cout*/)
{
	out << "{";

	for (unsigned int n = 0;n < m_val.size();n++)
	{
		out << " (" << m_val[n].val
			<< ' ' << m_val[n].prob
			<< ' ' << m_val[n].sum_to_here
			<< ')';
	}

	out << " }\n";
}

void Distribution::self_test()
{
	int num;

	for (num = 1;num < 1000;num++)
	{
		Distribution d;

		for (int n = 0;n < num;n++)
		{
			d.add(1.0, n);
		}

		// need to sort the values, since we are calling binary_search()
		// not select() (which ensures that they are sorted first)

		d.sort_values();
		// d.dump();

		for (double v = 0.5;v < (double) num;v++)
		{
			int i = d.binary_search(v);
			
			if (i != (int) v)
			{
				std::cout << "Error!\n";
				assert(false);
			}
		}

		std::cout << num << " ok\n";

		// select some random values and make sure they are all present

		d.clear();

		bool found[1100];

		for (int n = 0;n < num;n++)
		{
			d.add(Random::rnd(10.0) + 5.0, n + 100);
			found[n + 100] = false;
		}

		for (int k = 0;k < num * 50;k++)
		{
			int x = d.select();
			assert(x >= 100 && x < 100 + num);
			found[x] = true;
		}
		
		for (int j = 0;j < num;j++)
		{
			assert(found[j + 100]);
		}
	}
}

#ifdef SELF_TEST

int main(int argc, char **argv)
{
	Distribution::self_test();
	return 0;
}

#endif // SELF_TEST

