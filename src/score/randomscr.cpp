// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2012 Saulo H. P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2012-08-06
//
// Saulo's scoring method.
//
// This "very intrincated" scoring system gives a random score component to the peptides.
// The idea behind this unusual scoring method is to add noise to the scoring function.
// 
// 

#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iostream>

#include "peptide.h"
#include "atom.h"
#include "scorer_combined.h"
#include "randomscr.h"

Randomscr::Randomscr()
{
}

Randomscr::~Randomscr()
{
}

double Randomscr::score(const Peptide& p, bool verbose) const
{
	int len = p.length();
	double total=0.0;
	std::srand(time(0));
	
	total = ((double)std::rand()/(double)RAND_MAX)*100.0;
	total /= sqrt((double) len);
	total = total * 28.0 + 50.0;

	std::cout << "Total = " << total << " !!\n";	
	return total;
}

