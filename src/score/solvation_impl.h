#ifndef SOLVATION_IMPL_H_INCLUDED
#define SOLVATION_IMPL_H_INCLUDED

#include <iostream>
#include <vector>
class Peptide;
class Residue;
class C_File;
class Atom;

class Solvation_impl
{
public:
	Solvation_impl();
	~Solvation_impl();

	// score a peptide (low scores are better)
	double score(const Peptide& peptide, bool verbose = false, bool continuous = false);

	// set the name of the solvation data file
	void set_data_file(const std::string &filename);

private:
	void load_data();

private:
	typedef std::vector<double> Double_Vec;
	typedef std::vector<Double_Vec> Double_Vec_Vec;

	Double_Vec_Vec m_bin;

	std::string m_filename;		// data file
	bool m_data_loaded;			// whether data has been loaded
	int m_first_bin;			// lowest valid index of m_bin[amino][*]
	int m_top_bin;				// hightest valid index + 1 of m_bin[amino][*]
	double m_solv_dist;			// "near" distance between CB atoms (eg. 10.0)
};

#endif // SOLVATION_IMPL_H_INCLUDED
