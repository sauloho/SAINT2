#ifndef SOLVATION_H_INCLUDED
#define SOLVATION_H_INCLUDED

#include <iostream>
#include <vector>
class Peptide;
class Solvation_impl;

/**
 * @brief Solvation_impl Potential.
 *
 * The relative solvant accessibility is estimated as the number of
 * other C_beta atoms within a sphere of radius 10A centred on the
 * residues's C_beta atom. [Jones, (1999) J. Mol. Biol. 287:797-815].
 */
class Solvation
{
public:
	Solvation();
	~Solvation();

	// score a peptide (low scores are better)
	double score(const Peptide& peptide, bool verbose = false, bool continuous = false);

	// set the name of the solvation data file
	void set_short_data_file(const std::string &filename);
	void set_long_data_file(const std::string &filename);

private:
	Solvation_impl *m_short;    // for short proteins
	Solvation_impl *m_long;     // for long proteins
};

#endif // SOLVATION_H_INCLUDED
