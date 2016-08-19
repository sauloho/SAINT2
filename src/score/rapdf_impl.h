#ifndef RAPDF_IMPL_H_INCLUDED
#define RAPDF_IMPL_H_INCLUDED

#include "rapdf.h"
#include "amino.h"
#include <string>
#include <vector>
#include <iostream>

// forward declarations
class Peptide;
class Atom;

class RAPDF_impl
{
public:
	// constructor
	RAPDF_impl();

	// destructor
	~RAPDF_impl();

	// score a peptide (low scores ate better)
	double score(const Peptide &p, bool verbose = false);

	// set the name of the data file
	void set_data_file(const std::string &filename);

private:
    // disable copy and assignment by making them private
	RAPDF_impl(const RAPDF_impl&);
	RAPDF_impl &operator = (const RAPDF_impl&);

	// read data file (if it has not already been read)
	void load_data();

private:
	std::string m_filename;		// data file
	bool m_data_loaded;			// whether data has been read

	int m_first_bin;			// lowest distance bin
	int m_top_bin;				// highest index + 1 of m_bin[n1][n2][*]

	// (dimensions: [NUM_RAPDF_IDS][NUM_RAPDF_IDS][RAPDF_TOP_BIN])
	double ***m_data;
};

#endif // RAPDF_IMPL_INCLUDED
