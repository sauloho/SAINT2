// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2009 Jonathan Ellis
//
// Jonathan Ellis <jonathan.ellis@mq.edu.au>
// 2009-06-24
//

#ifndef HBOND_H_INCLUDED
#define HBOND_H_INCLUDED

class Peptide;

/**
 * @brief Hydrogen bond potential.
 *
 * A simple hydrogen bond potential is implemented as summing over minus
 * the number of backbone hydrogen bonds.  Rather than adding explicit
 * hydrogen atoms and calculating the relative orientation of the N-H
 * and C-O bonds, three distinct distances are used to estimate the
 * presence of a backbone hydrogen bond.  These are,
 *
 * <ol>
 *   <li> 2 <= d(N_i, O_j) <= 4,
 *
 *   <li> d(N_i, O_j) < d(N_i, C_j), and
 *
 *   <li> d(N_i, O_j) < d(C_alpha_i, O_j),
 * </ol>
 *
 * where d(A_i, B_j) is the distance, in A, between atom A of residue i
 * and atom B of residue j (|i - j| >= 1).  This simple scheme captures
 * the intuitive concept of a hydrogen bond being a contact between a
 * backbone N and O atom in an angle where both are facing each other,
 * i.e., the C_alpha and C atoms are more distant.
 */
class HBond
{
public:
	HBond();
	~HBond();

	double score(const Peptide& peptide, bool verbose = false) const;
};

#endif // HBOND_H_INCLUDED

