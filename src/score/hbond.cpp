// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2009 Jonathan Ellis
//
// Jonathan Ellis <jonathan.ellis@mq.edu.au>
// 2009-06-24
//
// Hydrogen bond potential
//
// A simple hydrogen bond potential is implemented as summing over minus
// the number of backbone hydrogen bonds.  Rather than adding explicit
// hydrogen atoms and calculating the relative orientation of the N-H
// and C-O bonds, three distinct distances are used to estimate the
// presence of a backbone hydrogen bond.  These are
//
//   2 < d(N_i, O_j) < 4, 
//   d(N_i, O_j) < d(N_i, C_j), and 
//   d(N_i, O_j) < d(C_alpha_i, O_j),
//  
// where d(A_i, B_j) is the distance, in A, between atom A of residue i
// and atom B of residue j (|i - j| >= 1).  This simple scheme captures
// the intuitive concept of a hydrogen bond being a contact between a
// backbone N and O atom in an angle where both are facing each other,
// i.e., the C_alpha and C atoms are more distant.

#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "peptide.h"
#include "atom.h"
#include "scorer_combined.h"
#include "hbond.h"

//#define RAW_SCORE

HBond::HBond()
{
}

HBond::~HBond()
{
}

double HBond::score(const Peptide& p, bool verbose) const
{
	int num_hbonds = 0;

	for (int i = p.start();i <= p.end();i++)
	{
		if (!(p.atom_exists(i, Atom_CA) && p.atom_exists(i, Atom_N)))
		{
			continue;
		}

		Point ca_i = p.atom_pos(i, Atom_CA);
		Point n_i = p.atom_pos(i, Atom_N);
	
		for (int j = p.start();j <= p.end();j++)
		{
			if (i != j)
			{
				if (!(p.atom_exists(j, Atom_C) && p.atom_exists(j, Atom_O)))
				{
					continue;
				}

				Point c_j = p.atom_pos(j, Atom_C);
				Point o_j = p.atom_pos(j, Atom_O);

				double dist  = n_i.distance(o_j);

				if (dist > 2.0 &&
					dist < 4.0 &&
					dist < n_i.distance(c_j) &&
					dist < ca_i.distance(o_j))
				{
					if (verbose)
					{
					}

					num_hbonds++;

					// cannot form a hydrogen bond with more than one
					// other residue at the same time, so stop inner loop
					break;
				}
			}
		}
	}

	double total = ((double) -num_hbonds);

#ifndef RAW_SCORE

	int len = p.length();

	// Alter the score so that it has about the same distribution
	// for all lengths

	total /= sqrt((double) len);

	// normalise so that all score types have approximately the same range
	total = total * 28.0 + 50.0;

#endif // RAW_SCORE

	return total;
}

