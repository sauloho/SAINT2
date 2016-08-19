#ifndef ORIENTATION_H_INCLUDED
#define ORIENTATION_H_INCLUDED

#include <iostream>
#include "amino.h"

class Peptide;
class Orientation_impl;

//#define ORIENT_DISTS	8
#define ORIENT_DISTS	18
#define ORIENT_ANGLES	10

/**
 * @brief Orientation Potential.
 */
class Orientation
{
public:
	Orientation();
	~Orientation();

	// score a peptide (low scores are better)
	double score(const Peptide& peptide, bool verbose = false);

	// set the name of the orientation data file
	void set_short_data_file(const std::string &filename);
	void set_long_data_file(const std::string &filename);

	// get the distance and angle bin between two residues (index n and m)
	// Returns false if there is no bin (too far apart or missing atoms)
	static bool get_bins(const Peptide &p, int n, int m,
		int *dist_bin, int *angle_bin);

	// dump the values used for scoring
	void dump(std::ostream &out = std::cout);

private:
	Orientation_impl *m_short;    // for short proteins
	Orientation_impl *m_long;     // for long proteins
};

#endif // ORIENTATION_H_INCLUDED

