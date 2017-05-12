#include <cstdlib> // PG added this
#include <cassert> // PG changed from assert.h to cassert
#include <cmath>
#include <iostream>
#include <fstream>
#include "peptide.h"
#include "common.h"
#include "atom_id.h"
#include "scorer_combined.h"
#include "lennard_jones.h"

//#define RAW_SCORE

Lennard_Jones::LJ_Params Lennard_Jones::m_lj[Num_Backbone][Num_Backbone];

Lennard_Jones::LJ_Params::LJ_Params()
{
}

Lennard_Jones::LJ_Params::LJ_Params(double c12_val, double c6_val)
	: c12(c12_val), c6(c6_val)
{
	// sigma is the distance at which the potential equals zero
	double sigma = pow(c12 / c6, 1.0 / 6.0);

	// recommended cutoff distance
	// (http://en.wikipedia.org/wiki/Lennard-Jones_potential)
	//
	// This causes a small discontinuity that is probably not worth
	// fixing for our purposes (the jump is about 1.5 percent of the
	// well depth).
	this->max_dist = 2.5 * sigma;
}

// return a single letter representing an atom type
// (backbone atoms only)

inline char atom_char1(Atom_Id a)
{
	switch (a)
	{
		case Atom_C:
		case Atom_CA:
		case Atom_CB:
			return 'C';

		case Atom_N:
			return 'N';

		case Atom_O:
			return 'O';

		default:
			std::cerr << "Error: Illegal value in atom_char1() in "
				<< __FILE__ << "\n";
			exit(1);
	}
}

class LJ_Static_Init
{
public:
	// Values taken from:
	//
	// http://web.njit.edu/all_topics/Prog_Lang_Docs/html/autodock/
	//		AD3.a.0UserGuide.html#29613
	//
	// (Table 1, "Self-consistent Lennard-Jones 12-6 parameters")

	LJ_Static_Init()
	{
		// parameters taken from
		// http://dasher.wustl.edu/tinker/distribution/params/oplsaa.prm

		const double cc_sigma = 3.75;
		const double cc_epsilon = 0.105;

		const double oo_sigma = 2.96;
		const double oo_epsilon = 0.21;

		const double nn_sigma = 3.25;
		const double nn_epsilon = 0.17;

		for (int a1 = 0;a1 < Num_Backbone;a1++)
		{
			for (int a2 = 0;a2 < Num_Backbone;a2++)
			{
				char c1 = atom_char1((Atom_Id) a1);
				char c2 = atom_char1((Atom_Id) a2);
				double e, s;

				if (c1 > c2)
				{
					swap(c1, c2);
				}

				if (c1 == 'C' && c2 == 'C')
				{
					e = cc_epsilon;
					s = cc_sigma;
				}
				else
				if (c1 == 'C' && c2 == 'N')
				{
					e = sqrt(cc_epsilon * nn_epsilon);
					s = 0.5 * (cc_sigma + nn_sigma);
				}
				else
				if (c1 == 'C' && c2 == 'O')
				{
					e = sqrt(cc_epsilon * oo_epsilon);
					s = 0.5 * (cc_sigma + oo_sigma);
				}
				else
				if (c1 == 'N' && c2 == 'N')
				{
					e = nn_epsilon;
					s = nn_sigma;
				}
				else
				if (c1 == 'N' && c2 == 'O')
				{
					e = sqrt(nn_epsilon * oo_epsilon);
					s = 0.5 * (nn_sigma + oo_sigma);
				}
				else
				if (c1 == 'O' && c2 == 'O')
				{
					e = oo_epsilon;
					s = oo_sigma;
				}
				else
				{
					std::cerr << "Error in LJ_Static_Init()\n";
					exit(1);
				}

				double c12 = e * pow(s, 12.0);
				double c6 = 2.0 * e * pow(s, 6.0);

				Lennard_Jones::m_lj[a1][a2] =
					Lennard_Jones::LJ_Params(c12, c6);

				// check if minimum is where it should be: at (s, -e)
				if (fabs(c12 / pow(s, 12.0) - c6 / pow(s, 6.0) + e) > 0.001)
				{
					std::cerr << "Incorrect values in LJ_Static_Init()\n";
					exit(1);
				}
			}
		}
	}
};

namespace
{
	LJ_Static_Init lj_static_init;
}

Lennard_Jones::Lennard_Jones()
{
}

Lennard_Jones::~Lennard_Jones()
{
}

double Lennard_Jones::score(const Peptide& p, bool verbose /*= false*/)
{
	double total = 0.0;
	
	for (int n1 = p.start() + 2;n1 <= p.end();n1++)
	{
		for (int a1 = 0;a1 < Num_Backbone;a1++)
		{
			if (!p.atom_exists(n1, (Atom_Id) a1))
			{
				continue;
			}

			Point p1 = p.atom_pos(n1, (Atom_Id) a1);

			
			// (n1 & n2 not in the same or an adjacent residue)

			for (int n2 = p.start();n2 < n1 - 1;n2++)
			{
				for (int a2 = 0;a2 < Num_Backbone;a2++)
				{
					if (!p.atom_exists(n2, (Atom_Id) a2))
					{
						continue;
					}

					const LJ_Params &lj = m_lj[a1][a2];
					Point p2 = p.atom_pos(n2, (Atom_Id) a2);

					if (fabs(p1.x - p2.x) < lj.max_dist &&
						fabs(p1.y - p2.y) < lj.max_dist &&
						fabs(p1.z - p2.z) < lj.max_dist)
					{
						double d = p1.distance(p2);

						if (d < lj.max_dist)
						{
							double old_total = total;

							if (d < 1.0)
							{
								// treat distance as 1.0
								// (to avoid ridiculous values)
								total += lj.c12 - lj.c6;
							}
							else
							{
								double d6 = square(d * d * d);
								total += (lj.c12 / square(d6)) - (lj.c6 / d6);
							}

							if (verbose)
							{
								std::cout << "LJ "
									<< n2 << " " << n1
									<< " = " << total - old_total
									<< " dist "
									<< d
									<< " " << n2 << " "
									<< Atom_Type((Atom_Id) a2).name()
									<< " at " << p2
									<< " & " << n1 << " "
									<< Atom_Type((Atom_Id) a1).name()
									<< " at " << p1
									<< "\n";
							}
						}
					}
				}
			}
		}
	}

    if(1)
    {   
	    int len = p.length();
        //std::cout << "!!! RAW LJ SCORE = " << total << ", len = " << len << "\n";

    	// normalise so that all score types have approximately the same range
	    total /= (len + 200.0);

    	total *= 350.0;
    }

	return total;
}

bool Lennard_Jones::steric_clash(const Peptide &p, double max_total /* -1.0*/)
{
	double total = 0.0;

	for (int n1 = p.start() + 2;n1 <= p.end();n1++)
	{
		for (int a1 = 0;a1 < Num_Backbone;a1++)
		{
			if (!p.atom_exists(n1, (Atom_Id) a1))
			{
				continue;
			}

			Point p1 = p.atom_pos(n1, (Atom_Id) a1);

			// (n1 & n2 not in the same or an adjacent residue)

			for (int n2 = p.start();n2 < n1 - 1;n2++)
			{
				for (int a2 = 0;a2 < Num_Backbone;a2++)
				{
					if (!p.atom_exists(n2, (Atom_Id) a2))
					{
						continue;
					}

					const LJ_Params &lj = m_lj[a1][a2];
					Point p2 = p.atom_pos(n2, (Atom_Id) a2);

					if (fabs(p1.x - p2.x) < lj.max_dist &&
						fabs(p1.y - p2.y) < lj.max_dist &&
						fabs(p1.z - p2.z) < lj.max_dist)
					{
						double d = p1.distance(p2);

						if (d < lj.max_dist)
						{
							double t;

							if (d < 1.0)
							{
								// treat distance as 1.0
								// (to avoid ridiculous values)
								t = lj.c12 - lj.c6;
							}
							else
							{
								double d6 = square(d * d * d);
								t = (lj.c12 / square(d6)) - (lj.c6 / d6);
							}

							// for purposes of creating decoys, Threshold
							// must be more than the value for any of the
							// native structures

							static const double Threshold = 20.0;

							if (t > Threshold)
							{
								return true;
							}

							total += t;
						}
					}
				}
			}
		}
	}

	if (max_total != -1.0 && total > max_total)
	{
		return true;
	}

	return false;
}

