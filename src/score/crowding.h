// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2012 Saulo H.P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2012-08-22
//

#ifndef CROWDING_INCLUDED
#define CROWDING_INCLUDED

class Peptide;
/**
 * 
 * Crowding scoring class is an intuitive class:
 *
 * Adds a score penalty for every pair of C-alphas more than
 * 30 angstroms apart. 
 *
 */
class Crowding
{
public:
	// Constructor
	Crowding(); 
	// Destructor
	~Crowding();

	/* This method returns the random Score for the Peptide! */
	double score(const Peptide& peptide, bool verbose = false) const;
};

#endif // CROWDING_INCLUDED
