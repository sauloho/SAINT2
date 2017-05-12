#ifndef ORIENTATION_IMPL_H_INCLUDED
#define ORIENTATION_IMPL_H_INCLUDED

#include <iostream>
#include "amino.h"

class Peptide;
class Residue;
class Atom;

/**
 * @brief Orientation Potential.
 */
class Orientation_impl
{
public:
	Orientation_impl();
	~Orientation_impl();

	// score a peptide (low scores are better)
	double score(const Peptide& peptide, bool verbose = false, bool continuous = false);

	// set the name of the orientation data file
	void set_data_file(const std::string &filename);

	// dump the values used for scoring
	void dump(std::ostream &out = std::cout);

private:
	void load_data();

private:
	std::string m_filename;		// name of orientation data file
	bool m_data_loaded;					// whether data has been loaded

	double m_data[ORIENT_DISTS][ORIENT_ANGLES][Amino::Num][Amino::Num];
};

#endif // ORIENTATION_IMPL_H_INCLUDED

