// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2012 Saulo H.P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2012-08-22
//

#ifndef CONTACT_INCLUDED
#define CONTACT_INCLUDED

class Peptide;

/**
 * 
 * Contact Information Scoring Class:
 *
 * Penalizes contacts and anti-contacts that are not satisfied by the model.
 * 
 */

class Contact
{
public:
	// Constructor
	Contact(); 
	// Destructor
	~Contact();

	/* This method returns the random Score for the Peptide! */
	double score(const Peptide& peptide, bool verbose = false) const;

	// Set the name of the Contact Map data file
	void set_short_data_file(const std::string &filename);
	void set_long_data_file(const std::string &filename);

	//
private:
	std::string m_filename;
};

#endif // CONTACT_INCLUDED
