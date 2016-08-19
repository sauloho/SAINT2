#ifndef TORSION_H_INCLUDED
#define TORSION_H_INCLUDED

#include <string>
#include <iostream>
#include <cstdlib> // PG added this

class Peptide;
class Torsion_impl;

// size of each torsion bin (degrees)
#define TORSION_SIZE 3

#define TORSION_BINS (360 / TORSION_SIZE)

/**
 * @brief Torsion angle potential.
 *
 * The torsion angle potential is derived from the porpensitites of each
 * amino acid to adopt different ocmbinations of (phi, psi) torsion
 * angles [Ramachandran et al. 1963].  The backbone (phi, psi) angle
 * space for each of the 20 amino acids was divided into N x N degree
 * bins.  No attempt was made to include side chain propensities.
 * The total energy is again summed over all individual contributions
 * in a protein.
 */
class Torsion
{
public:
	Torsion();
	~Torsion();

	// score a peptide (low scores are better)
	double score(const Peptide& peptide, bool verbose = false);

	// set the name of the torsion data file
	void set_short_data_file(const std::string &filename);
	void set_long_data_file(const std::string &filename);

	// get the bins to use for residue n
	static bool get_phi_psi_bin(const Peptide &p, int n,
		int *phi_bin, int *psi_bin);

private:
	void load_data();

private:
	Torsion_impl *m_short;    // for short proteins
	Torsion_impl *m_long;     // for long proteins
};

#endif	// TORSION_H_INCLUDED
