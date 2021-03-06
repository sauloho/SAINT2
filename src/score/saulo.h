// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2012 Saulo H.P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2012-08-06
//

#ifndef SAULO_INCLUDED
#define SAULO_INCLUDED

class Peptide;
/**
 * 
 * Saulo's Scoring Class is a very intuitive class:
 *
 * Add a random score component to represent noise.
 * 
 */
class Saulo
{
public:
	// Constructor
	Saulo(); 
	// Destructor
	~Saulo();

	/* This method returns the random Score for the Peptide! */
	double score(const Peptide& peptide, bool verbose = false);

	// Set the name of the Contact Map data file
	void set_short_data_file(const std::string &filename);
	void set_long_data_file(const std::string &filename);

private:
	void load_data(const Peptide& p);
	//
private:
	std::string m_filename;
	bool m_data_loaded;					// whether data has been loaded
	int *m_con[2],num_con;
    int m_previous_len;
};

#endif // SAULO_INCLUDED
