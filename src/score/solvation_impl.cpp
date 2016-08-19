
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>

#include "peptide.h"
#include "residue.h"
#include "atom.h"
#include "amino.h"
#include "c_file.h"
#include "scorer_combined.h"
#include "solvation.h"
#include "solvation_impl.h"

Solvation_impl::Solvation_impl() :
	m_data_loaded(false)
{
}

Solvation_impl::~Solvation_impl()
{
}

void Solvation_impl::set_data_file(const std::string &filename)
{
	m_filename = filename;
}

void Solvation_impl::load_data()
{
	if (m_data_loaded)
	{
		return;
	}

	if (m_filename.empty())
	{
		std::cerr << "Error: Solvation data file undefined\n";
		exit(1);
	}

//std::cout << "LOADING SOLVATION DATA\n";

	// exits with error message if file does not exist
	C_File file(m_filename, "r", "Solvation data file");

	static const int Max_Len = 1000;
	char buffer[Max_Len];
	int last_bin;

	if (!file.next_line(buffer, Max_Len) ||
		sscanf(buffer, "%d %d %lf",
			   &m_first_bin, &last_bin, &m_solv_dist) != 3)
	{
		std::cerr << "Error: expected solvation data file "
			<< m_filename
			<< " to start with three numbers (min count, max count, radius)\n";
		exit(1);
	}

	m_top_bin = last_bin + 1;
	m_bin.resize(Amino::Num);

	for (int a = 0;a < Amino::Num;a++)
	{
		if (!file.next_line(buffer, Max_Len) ||
			std::string(buffer, 3) != std::string(Amino(a).abbr()))
		{
			std::cerr << "Error: expected \""
				<< Amino(a).abbr() << "\" on line "
				<< file.line_num()
				<< " of solvation data file "
				<< m_filename
				<< "\n";
			exit(1);
		}

		m_bin[a].resize(m_top_bin);

		for (int c = m_first_bin;c < m_top_bin;c++)
		{
			int bin;

			if (!file.next_line(buffer, Max_Len) ||
				sscanf(buffer, "%d %lf", &bin, &m_bin[a][c]) != 2 ||
				bin != c)
			{
				std::cerr << "Error: expected \""
					<< c << "\" followed by value on line "
					<< file.line_num()
					<< " of "
					<< m_filename
					<< "\n";
				exit(1);
			}
		}
	}

	if (file.next_line(buffer, Max_Len))
	{
		std::cerr << "Warning: ignoring extra data on line "
			<< file.line_num()
			<< " of solvation data file "
			<< m_filename
			<< "\n";
	}

	m_data_loaded = true;
}

inline Point cbeta_pos(const Peptide &p, int n, bool *failed)
{
	Atom_Id a = (p.res(n).amino().is_glycine() ? Atom_CA : Atom_CB);

	*failed = !(p.atom_exists(n, a));

	if (*failed)
	{
		Point pp;
		return pp;
	}
	else
	{
		return p.atom_pos(n, a);
	}
}

double Solvation_impl::score(const Peptide &p, bool verbose)
{
	// (does nothing if already loaded)
	load_data();

	std::vector<int> count;
	count.resize(p.full_length());

	int n;
	for (n = p.start();n <= p.end();n++)
	{
		count[n] = 0;
	}

	for (n = p.start() + 1;n <= p.end();n++)
	{
		bool failed = false;
		Point pos = cbeta_pos(p, n, &failed);
		if (failed) { continue; }
		
		if (failed)
		{
			continue;
		}

		for (int m = p.start();m < n;m++)
		{
			Point p2 = cbeta_pos(p, m, &failed);
			if (failed) { continue; }

			if (pos.closer_than(m_solv_dist, p2))
			{
				count[n]++;
				count[m]++;
			}
		}
	}

	double total = 0.0;

	for (n = p.start();n <= p.end();n++)
	{
		int a = p.res(n).amino().num();
		int c = count[n];

		if (c < m_first_bin) { c = m_first_bin; }
		else
		if (c >= m_top_bin) { c = m_top_bin - 1; }

/*
std::cout
	<< m_bin[a][c]
	<< "  "
	<< n << "/" << p.res(n).res_seq_str()
	<< " " << p.res(n).amino().name()
	<< " (line " << p.res(n).pdb_line()
	<< ")  count = " << count[n]
	<< "\n";
*/

		total += m_bin[a][c];

		if (verbose)
		{
		}
	}

#ifndef RAW_SCORE
	// Alter the score so that it has about the same distribution
	// for all lengths

	int len = p.length();

	if (len <= SHORT_PEPTIDE)
	{
		total /= sqrt((double) len);

		// normalise so that all score types have approximately the same range
		total *= 60.0;
	}
	else
	{
		// normalise so that all score types have approximately the same range
		total *= 2.7;
	}

#endif // RAW_SCORE

	return total;
}

