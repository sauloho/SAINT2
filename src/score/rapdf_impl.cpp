
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cmath>
#include "atom.h"
#include "amino.h"
#include "peptide.h"
#include "c_file.h"
#include "scorer_combined.h"
#include "rapdf_impl.h"

//#define RAW_SCORE

RAPDF_impl::RAPDF_impl()
	: m_data_loaded(false)
{
	m_data = NULL;
}

RAPDF_impl::~RAPDF_impl()
{
	if (m_data != NULL)
	{
		for (int b1 = 0;b1 < NUM_RAPDF_IDS;b1++)
		{
			for (int b2 = 0;b2 <= b1;b2++)
			{
				delete [] m_data[b1][b2];
			}

			delete [] m_data[b1];
		}

		delete [] m_data;
	}
}

void RAPDF_impl::set_data_file(const std::string &filename)
{
	m_filename = filename;
}

void RAPDF_impl::load_data()
{
	if (m_data_loaded)
	{
		return;
	}

	if (m_filename.empty())
	{
		std::cerr << "Error: RAPDF data file undefined\n";
		exit(1);
	}

//std::cout << "LOADING RAPDF DATA\n";

	// exits with error message if file does not exist
	C_File file(m_filename, "r", "RAPDF data file");

	// std::cout << "Loading RAPDF data" << std::endl;

	static const int Max_Len = 1000;
	char buffer[Max_Len];

	if (!file.next_line(buffer, Max_Len) ||
		sscanf(buffer, "%d %d", &m_first_bin, &m_top_bin) != 2)
	{
		std::cerr << "Error: expected RAPDF data file "
			<< m_filename
			<< " to start with two numbers (first bin, top bin)\n";
		exit(1);
	}

	int b1, b2, d;

	m_data = new double** [NUM_RAPDF_IDS];

	for (b1 = 0;b1 < NUM_RAPDF_IDS;b1++)
	{
		m_data[b1] = new double* [NUM_RAPDF_IDS];

		for (b2 = 0;b2 < NUM_RAPDF_IDS;b2++)
		{
			m_data[b1][b2] = new double[m_top_bin];
		}
	}

	for (b1 = 0;b1 < NUM_RAPDF_IDS;b1++)
	{
		std::string b1_aa_str = Amino::rapdf_amino(b1).abbr();
		std::string b1_a_str = Amino::rapdf_atom(b1).name();

		for (b2 = 0;b2 <= b1;b2++)
		{
			std::string b2_aa_str = Amino::rapdf_amino(b2).abbr();
			std::string b2_a_str = Amino::rapdf_atom(b2).name();

			char aa1_str[100], a1_str[100];
			char aa2_str[100], a2_str[100];

			if (!file.next_line(buffer, Max_Len) ||
				sscanf(buffer, "%s %s %s %s",
				aa1_str, a1_str, aa2_str, a2_str) != 4 ||
				std::string(aa1_str) != b1_aa_str ||
				std::string(a1_str) != b1_a_str ||
				std::string(aa2_str) != b2_aa_str ||
				std::string(a2_str) != b2_a_str)
			{
				std::cerr << "Error on line " << file.line_num()
					<< " of " << m_filename
					<< ": expected \""
					<< b1_aa_str << " " << b1_a_str << "  "
					<< b2_aa_str << " " << b2_a_str << "\"\n";
				exit(1);
			}

			int dist;
			double val;

			for (d = m_first_bin;d < m_top_bin;d++)
			{
				if (!file.next_line(buffer, Max_Len) ||
					sscanf(buffer, "%d %lf", &dist, &val) != 2 ||
					dist != d)
				{
					std::cerr << "Error on line " << file.line_num()
						<< " of " << m_filename
						<< ": expected \""
						<< d
						<< "\" followed by value\n";
					exit(1);
				}

				m_data[b1][b2][d] = m_data[b2][b1][d] = val;
			}
		}
	}

	if (file.next_line(buffer, Max_Len))
	{
		std::cerr << "Warning: ignoring extra data on line "
			<< file.line_num()
			<< " of "
			<< m_filename
			<< "\n";
	}

	// std::cout << "Loaded RAPDF data" << std::endl;
	m_data_loaded = true;
}

double RAPDF_impl::score(const Peptide &p, bool verbose)
{
	// (does nothing if already loaded)
	load_data();

	double total = 0.0;

	for (int n1 = p.start() + 2;n1 <= p.end();n1++)
	{
		const Residue &res1 = p.res(n1);

		for (int a1 = 0;a1 < res1.num_atoms();a1++)
		{
			Atom_Id t1 = res1.atom(a1).type().type();

			if (t1 == Atom_Undef)
			{
				continue;
			}

			int id_1 = res1.amino().rapdf_id(t1);

			/*
			if (a1 >= Num_Backbone)
			{
				std::cerr << "HERE !!!\n";
				exit(1);
			}
			*/

			Point pos1 = p.atom_pos(n1, t1);

			for (int n2 = p.start();n2 < n1 - 1;n2++)
			{
				const Residue &res2 = p.res(n2);

				for (int a2 = 0;a2 < res2.num_atoms();a2++)
				{
					Atom_Id t2 = res2.atom(a2).type().type();

					if (t2 == Atom_Undef)
					{
						continue;
					}

					int id_2 = res2.amino().rapdf_id(t2);

					Point pos2 = p.atom_pos(n2, t2);

					if (pos1.closer_than(m_top_bin, pos2))
					{
						double d_dist = pos1.distance(pos2);
						int dist = (int) d_dist;
						double frac = d_dist - (double) dist;

						assert(dist < m_top_bin);

						int other_bin = (frac > 0.5 ? dist + 1 : dist - 1);

						if (dist < m_first_bin)
						{
							dist = m_first_bin;
						}

						if (other_bin < m_first_bin ||
							other_bin >= m_top_bin)
						{
							total += m_data[id_1][id_2][dist];
						}
						else
						{
							double weight =
								(frac > 0.5 ? frac - 0.5 : 0.5 - frac);
							double total_weight = weight + 0.5;

							double amount = (m_data[id_1][id_2][dist] * 0.5 +
							 	m_data[id_1][id_2][other_bin] * weight) /
							 	total_weight;

							total += amount;

							if (verbose)
							{
							}
						}
					}
					else
					{
						double amount = m_data[id_1][id_2][m_top_bin - 1];
						total += amount;
					}

					/*
					int dist = pos1.int_distance(pos2);

					if (dist < m_top_bin)
					{
						if (dist < m_first_bin)
						{
							dist = m_first_bin;
						}

						total += m_data[id_1][id_2][dist];
					}
					*/
				}
			}
		}
	}

	// Alter the score so that it has about the same distribution
	// for all lengths

#ifndef RAW_SCORE
	int len = p.length();

	if (len <= SHORT_PEPTIDE)
	{
		total /= (double) len;
		// put on the same scale as the other scoring terms
		total *= 7.0;
	}
	else
	{
		total /= sqrt((double) len);
		// put on the same scale as the other scoring terms
		total *= 0.75;
	}
#endif // RAW_SCORE

	return total;
}

