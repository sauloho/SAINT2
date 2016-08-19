#ifndef CORE_IMPL_H_INCLUDED
#define CORE_IMPL_H_INCLUDED

#include "core.h"
#include "amino.h"
#include <string>
#include <vector>
#include <iostream>

class Peptide;
class Residue;
class Atom;

class CORE_impl
{
	friend class CORE_Static_Init;
public:
	// constructor
	CORE_impl();
	// destructor
	~CORE_impl();

	// score a peptide (low scores are better)
	double score(const Peptide& peptide,double w_LJ, double w_RAPDF, bool verbose = false);
	
	// set the name of the data file for RAPDF
	void set_data_file(const std::string &filename);
	// set the name of the data file for Orientation
	void set_data_file_ori(const std::string &filename);

private:
    // disable copy and assignment by making them private
	CORE_impl(const CORE_impl&);
	CORE_impl &operator = (const CORE_impl&);

	// read data file (if it has not already been read)
	void load_data();
	void load_rapdf_ids();

protected:
	struct CORE_Params
	{
		double c12;
		double c6;
		double max_dist;
		double sigma;
	
	
		CORE_Params();
		CORE_Params(double c12_val, double c6_val);
	};

	std::string m_filename;		// data file for RAPDF
	std::string m_filename_ori;	// data file for Orientation
	bool m_data_loaded;			// whether data has been read

	int m_first_bin;			// lowest distance bin
	int m_top_bin;				// highest index + 1 of m_bin[n1][n2][*]

	// (dimensions: [NUM_RAPDF_IDS][NUM_RAPDF_IDS][RAPDF_TOP_BIN])
	double m_data_ori[ORIENT_DISTS][ORIENT_ANGLES][Amino::Num][Amino::Num];
	double ***m_data;
	int **m_rapdf_ids;

	static CORE_Params m_lj[Num_Backbone][Num_Backbone];
};

#endif // CORE_IMPL_H_INCLUDED

