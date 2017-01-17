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
	std::string m_filename;
	bool m_data_loaded;					// whether data has been loaded
	int *m_con[2],num_con;
	// number of spans, read from span file
	int m_num_spans;
	// a vector of span starts and a vector of span ends
	//int *m_spans[2];
	// a vector of whether each residue is in centre of a span
	int m_layer[1000];
	// a vector of whether each residue is not in a span
	int m_loop[1000];
	int m_previous_len;
};

#endif // ELEANOR_INCLUDED
