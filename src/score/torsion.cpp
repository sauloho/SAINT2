
#include "torsion.h"
#include "torsion_impl.h"
#include "peptide.h"
#include "geom.h"
#include "scorer_combined.h"

Torsion::Torsion()
{
	m_short = new Torsion_impl;
	m_long = new Torsion_impl;
}

Torsion::~Torsion()
{
	delete m_short;
	delete m_long;
}

void Torsion::set_short_data_file(const std::string &filename)
{
	m_short->set_data_file(filename);
}

void Torsion::set_long_data_file(const std::string &filename)
{
	m_long->set_data_file(filename);
}

double Torsion::score(const Peptide &p, bool verbose)
{
	if (p.length() <= SHORT_PEPTIDE)
	{
		return m_short->score(p, verbose);
	}
	else
	{
		return m_long->score(p, verbose);
	}
}

bool Torsion::get_phi_psi_bin(const Peptide &p, int n,
	int *phi_bin, int *psi_bin)
{
	double phi = p.conf().phi(n);
	double psi = p.conf().psi(n);

	if (phi == TORSION_UNKNOWN || psi == TORSION_UNKNOWN)
	{
		return false;
	}

	*phi_bin = (int) (rad2deg(phi) / (double) TORSION_SIZE);
	*psi_bin = (int) (rad2deg(psi) / (double) TORSION_SIZE);

	assert(*phi_bin >= 0 && *phi_bin < TORSION_BINS &&
		   *psi_bin >= 0 && *psi_bin < TORSION_BINS);

	return true;
}

