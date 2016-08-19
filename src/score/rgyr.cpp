// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2012 Saulo H. P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2012-08-23
//
// Radius of Gyration scoring method.
//
// A scoring class that calculates the Radius of gyration (Rg) of the growing peptide 
// and adds a score penalty proportional to the difference between the expected
// Rg.

#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

#include "peptide.h"
#include "atom.h"
#include "atom_id.h"
#include "scorer_combined.h"
#include "rgyr.h"

Rgyr::Rgyr()
{
}

Rgyr::~Rgyr()
{
}

double Rgyr::score(const Peptide& p, bool verbose) const
{
	int len = p.length();	
	double total=0.0;

//	total = pow(p.radius_of_gyr() - (1.484*pow(len,0.4)+5.33333), 2 )  ;
	if(p.length() > 50)
		total = p.radius_of_gyr();

	return total;
}

