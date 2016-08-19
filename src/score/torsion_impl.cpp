
#include <iostream>
#include <string>
#include <cstdio>
#include <cassert>
#include "peptide.h"
#include "geom.h"
#include "c_file.h"
#include "stream_printf.h"
#include "scorer_combined.h"
#include "torsion.h"
#include "torsion_impl.h"

Torsion_impl::Torsion_impl() :
	m_data_loaded(false)
{
}

Torsion_impl::~Torsion_impl()
{
}

void Torsion_impl::set_data_file(const std::string &filename)
{
	m_filename = filename;
}

void Torsion_impl::load_data()
{
	if (m_data_loaded)
	{
		return;
	}

	if (m_filename.empty())
	{
		std::cerr << "Error: Torsion data file undefined\n";
		exit(1);
	}

	// exits with error message if the file does not exist
	C_File file(m_filename, "r", "Torsion data file");

	static const int Max_Len = 1000;
	char buffer[Max_Len];

	for (int a = 0;a < Amino::Num;a++)
	{
		if (!file.next_line(buffer, Max_Len) ||
			std::string(buffer, 3) != std::string(Amino(a).abbr()))
		{
			std::cerr << "Error: expected \""
				<< Amino(a).abbr() << "\" on line "
				<< file.line_num()
				<< " of torsion data file "
				<< m_filename
				<< "\n";
			exit(1);
		}

		int ph, ps;
		double val;

		for (int phi_bin = 0;phi_bin < TORSION_BINS;phi_bin++)
		{
			for (int psi_bin = 0;psi_bin < TORSION_BINS;psi_bin++)
			{
				if (!file.next_line(buffer, Max_Len) ||
					sscanf(buffer, "%d %d %lf", &ph, &ps, &val) != 3 ||
					ph != phi_bin || ps != psi_bin)
				{
					std::cerr << "Error: expected \""
						<< phi_bin << ' ' << psi_bin
						<< "\" followed by value on line "
						<< file.line_num()
						<< " of "
						<< m_filename
						<< "\n"
						<< "(TORSION_BINS = " << TORSION_BINS << ")\n";
					exit(1);
				}
				
				m_data[phi_bin][psi_bin][a] = val;
			}
		}
	}

	if (file.next_line(buffer, Max_Len))
	{
		std::cerr << "Warning: ignoring extra data on line "
			<< file.line_num()
			<< " of torsion data file "
			<< m_filename
			<< "\n";
	}

	m_data_loaded = true;
}

double Torsion_impl::score(const Peptide& p, bool verbose)
{
	// (does nothing if data already loaded)
	load_data();

	double total = 0.0;
	int phi_bin, psi_bin;

	// (first and last residue are missing one torsion angle, so
	// exclude them)

	for (int n = p.start() + 1;n < p.end();n++)
	{
		if (p.res(n).amino().is_glycine())
		{
			// glycine torsion propensities cannot be calculated reliably,
			// since it occupies regions of the Ramachandran plot that
			// no other amino acid does
			continue;
		}

		if (Torsion::get_phi_psi_bin(p, n, &phi_bin, &psi_bin))
		{
			int a = p.res(n).amino().num();
			total += m_data[phi_bin][psi_bin][a];

			/*
			double phi = rad2deg(p.conf().phi(n));
			double psi = rad2deg(p.conf().psi(n));

			std::cout 
				<<  m_data[phi_bin][psi_bin][a] << " "
				<< n << " " << p.res(n).amino().abbr()
				<< Printf(" %.1f", phi)
				<< Printf(" %.1f", psi)
				<< "\n";
			*/

			if (verbose)
			{
				std::cout
					<< m_data[phi_bin][psi_bin][a]
					<< " tor ("
					<< n << ") "
					<< p.res(n).amino().abbr()
					<< ' '
					<< Printf("%6.1f", range_m180_180(rad2deg(p.conf().phi(n))))
					<< ' '
					<< Printf("%6.1f", range_m180_180(rad2deg(p.conf().psi(n))))
					<< " [" << phi_bin << ' ' << psi_bin << "]\n";
			}
		}
	}

#ifndef RAW_SCORE
	int len = p.length();

	if (len <= SHORT_PEPTIDE)
	{
		// normalise so that all score types have approximately the same range
		total *= 8.0;
	}
	else
	{
		// make length independent
		total /= log((double) len);

		// normalise so that all score types have approximately the same range
		total *= 35.0;
	}
#endif // RAW_SCORE

	return total;
}

