
#include <cassert>
#include <cmath>
#include <iostream>
#include "common.h"
#include "peptide.h"
#include "ribosome.h"

Ribosome::Ribosome()
{
}

Ribosome::~Ribosome()
{
}

double Ribosome::score(const Peptide& p, bool verbose /*= false*/)
{
	if (p.full_grown() || p.length() < 2)
	{
		return 0.0;
	}

	// most recently extruded residue
	int most_recent = (reverseSaint ? p.start() : p.end());
	int n_start = p.start();
	int n_end = p.end();

	if (reverseSaint)
	{
		n_start++;
	}
	else
	{
		n_end--;
	}

	double x_limit = p.atom_pos(most_recent, Atom_CA).x - 5.0;
	double total = 0.0;

	for (int n = n_start;n <= n_end;n++)
	{
		double x = p.atom_pos(n, Atom_CA).x;

		if (x < x_limit)
		{
			double diff = (double) (x_limit - x);
			total += diff * diff;

			if (verbose)
			{
				std::cout << "Ribosome scorer: x diff = " << diff << "\n";
			}
		}
	}

	if (verbose)
	{
		std::cout << "Ribosome scorer: total = " << total << "\n";
	}

	// no need to normalise for length -- same penalty for
	// the same number of atoms beyond the limit for all lengths

	return total;
}

