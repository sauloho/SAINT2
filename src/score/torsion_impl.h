#ifndef TORSION_IMPL_H_INCLUDED
#define TORSION_IMPL_H_INCLUDED

#include <string>
#include <iostream>
#include "torsion.h"
#include "amino.h"

class Peptide;

class Torsion_impl
{
public:
	Torsion_impl();
	~Torsion_impl();

	// set the name of the torsion data file
	void set_data_file(const std::string &filename);

	// score a peptide (low scores are better)
	double score(const Peptide& peptide, bool verbose = false);

private:
	void load_data();

private:
	std::string m_filename;
	bool m_data_loaded;

	// (first two dimensions are phi and psi / TORSION_SIZE)
	double m_data[TORSION_BINS][TORSION_BINS][Amino::Num];
};

#endif	// TORSION_IMPL_H_INCLUDED
