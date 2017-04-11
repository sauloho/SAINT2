// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2012 Saulo H.P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2012-08-06
//

#ifndef ELEANOR_INCLUDED
#define ELEANOR_INCLUDED

class Peptide;
/**
 * 
 * Saulo's Scoring Class is a very intuitive class:
 *
 * Add a random score component to represent noise.
 * 
 */
class Eleanor
{
public:
	// Constructor
	Eleanor(); 
	// Destructor
	~Eleanor();

	/* This method returns the random Score for the Peptide! */
	double score(const Peptide& peptide, bool verbose = false);

	// Set the name of the Span data file
	void set_short_data_file(const std::string &filename);
	void set_long_data_file(const std::string &filename);

private:
	void load_data(const Peptide& p);
	//
private:
	// not in use at the moment
	std::string m_filename;

	bool m_data_loaded;					// whether data has been loaded
	// an array of indeces for looking up each amino acid in the pot table
	int m_aaindeces[1000];
	int m_previous_len;

	// the membrane potential table
	double m_mempot[20][34];

	// the order that amino acids appear as rows in the table
	std::string zpot_residues_list[20];

};

#endif // ELEANOR_INCLUDED
