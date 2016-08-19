// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2012 Saulo H.P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2012-08-23
//

#ifndef RGYR_INCLUDED
#define RGYR_INCLUDED

class Peptide;
/**
 * 
 * A scoring class that calculates the Radius of gyration (Rg) of the growing peptide 
 * and adds a score penalty proportional to the difference between the expected
 * Rg.
 *
 */
class Rgyr
{
public:
	// Constructor
	Rgyr(); 
	// Destructor
	~Rgyr();

	/* This method returns the random Score for the Peptide! */
	double score(const Peptide& peptide, bool verbose = false) const;
};

#endif // Rgyr_INCLUDED
