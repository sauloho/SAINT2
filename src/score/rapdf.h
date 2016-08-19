#ifndef RAPDF_H_INCLUDED
#define RAPDF_H_INCLUDED

#include "rapdf.h"
#include "amino.h"
#include <string>
#include <vector>
#include <iostream>

// forward declarations
class Peptide;
class RAPDF_impl;

class RAPDF
{
public:
	// constructor
	RAPDF();

	// destructor
	~RAPDF();

	// score a peptide (low scores ate better)
	double score(const Peptide &p, bool verbose = false);

	// set the name of the RAPDF data file
	void set_short_data_file(const std::string &filename);
	void set_long_data_file(const std::string &filename);

private:
    // disable copy and assignment by making them private
	RAPDF(const RAPDF&);
	RAPDF &operator = (const RAPDF&);

private:
	RAPDF_impl *m_short;	// for short proteins
	RAPDF_impl *m_long;		// for long proteins
};

#endif // RAPDF_INCLUDED
