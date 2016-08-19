// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2012 Saulo H.P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2015-03-23
//

#ifndef PREDTOR_INCLUDED
#define PREDTOR_INCLUDED

class Peptide;
/**
 * 
 * Predicted Torsion Angle Score: 
 *
 * Compute the difference between predicted and observed torsion angles. 
 * 
 */
class PredTor
{
public:
	// Constructor
	PredTor(); 
	// Destructor
	~PredTor();

	/* This method returns the random Score for the Peptide! */
	double score(const Peptide& peptide, bool verbose = false);

	// Set the name of the Predicted Torsion Angle data file
	void set_short_data_file(const std::string &filename);
	void set_long_data_file(const std::string &filename);

private:
	void load_data(const Peptide& p);
	//
private:
	std::string m_filename;
	bool m_data_loaded;					// whether data has been loaded
	int m_maxtor;
	double *m_phipsi[2];
};

#endif // PREDTOR_INCLUDED
