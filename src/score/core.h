#ifndef CORE_H_INCLUDED
#define CORE_H_INCLUDED

#include "core.h"
#include "amino.h"
#include <string>
#include <vector>
#include <iostream>


//#define ORIENT_DISTS	8
#define ORIENT_DISTS	18
#define ORIENT_ANGLES	10

// forward declarations
class Peptide;
class CORE_impl;

class CORE
{
public:
	// constructor
	CORE();

	// destructor
	~CORE();

	// score a peptide (low scores ate better)
	double score(const Peptide &p, double weight1, double weight2, bool verbose = false);

	// set the name of the RAPDF data file
	void set_short_data_file(const std::string &filename);
	void set_long_data_file(const std::string &filename);
	// set the name of the Orientation data file
	void set_short_data_file_ori(const std::string &filename);
	void set_long_data_file_ori(const std::string &filename);

private:
    // disable copy and assignment by making them private
	CORE(const CORE&);
	CORE &operator = (const CORE&);

private:
	CORE_impl *m_short;	// for short proteins
	CORE_impl *m_long;	// for long proteins
};

#endif // CORE_H_INCLUDED
