
#include "solvation.h"
#include "solvation_impl.h"
#include "peptide.h"
#include "scorer_combined.h"

Solvation::Solvation()
{
	m_short = new Solvation_impl;
	m_long = new Solvation_impl;
}

Solvation::~Solvation()
{
	delete m_short;
	delete m_long;
}

void Solvation::set_short_data_file(const std::string &filename)
{
	m_short->set_data_file(filename);
}

void Solvation::set_long_data_file(const std::string &filename)
{
	m_long->set_data_file(filename);
}

double Solvation::score(const Peptide &p, bool verbose, bool continuous)
{
	if (p.length() <= SHORT_PEPTIDE || continuous)
	{
		return m_short->score(p, verbose, continuous);
	}
	else
	{
		return m_long->score(p, verbose);
	}
}

