
#include "rapdf.h"
#include "rapdf_impl.h"
#include "peptide.h"
#include "scorer_combined.h"

RAPDF::RAPDF()
{
	m_short = new RAPDF_impl;
	m_long = new RAPDF_impl;
}

RAPDF::~RAPDF()
{
	delete m_short;
	delete m_long;
}

void RAPDF::set_short_data_file(const std::string &filename)
{
	m_short->set_data_file(filename);
}

void RAPDF::set_long_data_file(const std::string &filename)
{
	m_long->set_data_file(filename);
}

double RAPDF::score(const Peptide &p, bool verbose, bool continuous)
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

