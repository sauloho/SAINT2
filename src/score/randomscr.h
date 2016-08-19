// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2012 Saulo H.P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2012-08-06
//

#ifndef RANDSCR_INCLUDED
#define RANDSCR_INCLUDED

class Peptide;
/**
 * 
 * Saulo's Scoring Class is a very intuitive class:
 *
 * Add a random score component to represent noise.
 * 
 */
class Randomscr
{
public:
	Randomscr(); 
	~Randomscr();

	/* This method returns the random Score for the Peptide! */
	double score(const Peptide& peptide, bool verbose = false) const;
};

#endif // RANDSCR_INCLUDED
