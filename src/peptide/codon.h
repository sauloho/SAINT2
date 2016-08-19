
#ifndef CODON_H_INCLUDED
#define CODON_H_INCLUDED

// Class representing a single codon.
// 
// Each codon corresponds to a number from 0 to 63.
//
// It is not necessary for other classes to know how to convert codon strings
// to/from these numbers (since they can just use the conversion functions in
// this class), but for reference, the number equals
//
// nuc0 * 16 + nuc1 * 4 + nuc2
//
// where nuc0, nuc1 and nuc2 are the three nucleotides in 5' to 3' order,
// with the following values for each type:
//
// A   = 0
// C   = 1
// G   = 2
// T/U = 3
//
// eg. AAA = 0, TTT = 63, CTG = 30 (1*16 + 3*4 + 2)
//
// For convenience, a complete list of all codons in ascending
// order is included at the end of codon.cpp.

#include <string>
#include <vector>
#include "amino.h"

class Codon_Helper;
class Static_Init;

class Codon
{
	friend class Codon_Helper;
	friend class Static_Init;

public:
	// number of codon types
	static const int Num = 64;

public:
	Codon();

	// create the nth type of code (0 to Num - 1)
	Codon(int n);

	// construct a codon from a three-nucleotide sequence (case insensitive)
	// (5' to 3' mRNA order; success can be checked with the illegal() function)
	Codon(char nuc0, char nuc1, char nuc2);

	// construct a codon from the first three characters of a string
	Codon(const char *str);

	// get the amino acid for this codon
	// (returns an "illegal" amino acid if it is a stop codon)
	Amino to_amino() const;

	// get the nth nucleotide for the codon (n = 0, 1 or 2)
	// The optional argument is the letter to use for thymine/uracil
	char nucleotide(int n, char t_u_value = 'T') const;

	// convert to a 3-letter string.
	// The optional argument is the letter to use for thymine/uracil
	std::string str(char t_u_value = 'T') const;

	// convert to an integer from 0 to Num - 1
	int num() const
	{ return m_val; }

	// check if the codon is illegal (eg. illegal values used in constructor)
	bool illegal() const
	{ return (m_val == c_illegal); }

	// alias for illegal()
	bool unknown() const
	{ return illegal(); }

	// reset to unknown/illegal state
	void make_unknown()
	{ m_val = c_illegal; }

	// compare for equality
	bool operator == (const Codon &other) const
	{ return (m_val == other.m_val); }

	// compare for inequality
	bool operator != (const Codon &other) const
	{ return (m_val != other.m_val); }

	// check if this is a stop codon
	bool is_stop_codon() const;

	// get the start codon (same codon as Methionine)
	static Codon start_codon();

protected:
	// do required static initialisations (called by Static_Init)
	static void do_static_init();

	// call destructors for objects created in do_static_init();
	static void do_static_delete();

private:
	// convert three nucleotides to a codon number (case insensitive)
	// Returns c_illegal on error
	int nucs_to_num(char nuc0, char nuc1, char nuc2);

	// convert a nucleotide character to a number from 0 to 3
	// (case insensitive)
	// Returns -1 on error
	int nuc_to_num(char nuc);

private:
	// helper object for performing conversions
	//
	// (This is a pointer so that its time of creation can be managed
	// by the Static_Init class, since there is a dependency order:
	// Amino::Amino_Helper must be created first)
	static Codon_Helper *m_helper;

	// illegal value
	static const int c_illegal = Num;

	// codon type
	unsigned char m_val;
};

// a vector of codons
typedef std::vector<Codon> Codon_Vec;

#endif // CODON_H_INCLUDED

