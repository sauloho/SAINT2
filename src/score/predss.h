// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2014 Saulo H.P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2014-27-04
//

#ifndef PREDSS_INCLUDED
#define PREDSS_INCLUDED

class Peptide;
/**
 * 
 * Predicted Secondary Structure's Scoring Class:
 *
 * Using an input predicted SS file, this function computes how much
 * the growing peptide deviates from the predicted secondary structure.
 * 
 */
class PredSS
{
public:
	// Constructor
	PredSS(); 
	// Destructor
	~PredSS();

	/* This method returns the random Score for the Peptide! */
	double score(const Peptide& peptide, bool verbose = false);

	// Set the name of the Predicted Secondary Structure data file
	void set_short_data_file(const std::string &filename);
	void set_long_data_file(const std::string &filename);

private:
	void load_data();
	//
private:
	std::string m_filename;
	bool m_data_loaded;					// whether data has been loaded
	char m_data[1000];
};

#endif // PREDSS_INCLUDED
