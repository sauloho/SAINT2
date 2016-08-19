#ifndef AMINO_H_INCLUDED
#define AMINO_H_INCLUDED

#include <vector>
#include "atom_id.h"

#define NUM_RAPDF_IDS   167

class Atom_Type;
class Amino_Helper;
class Static_Init;

// Class representing a single amino acid

class Amino
{
	friend class Amino_Helper;
	friend class Atom_Type;
	friend class Static_Init;
	friend class CORE_impl;

public:
	// number of amino acid types
	static const int Num = 20;

	// number of amino acid types (including ambiguous, unknown and undefined
	// types)
	static const int Full_Num = 27;

	// ids of amino acids with special properties
	static const int AA_Cysteine = 4;
	static const int AA_Glycine = 7;
	static const int AA_Proline = 14;

public:
	// default constructor
	Amino() : m_val(c_undefined)
	{ }

	// constructor from a single character (case insensitive)
	// (success or failure can be checked with the undefined() function)
	Amino(char ch);

	// constructor from an integer
	// (standard amino acids are numbered 0 to Num-1)
	Amino(int val);

	// constructor from three-letter abbreviation
	Amino(const char *abbr_str);

	// get the unique RAPDF id for a particular atom within this amino acid
	// (returns -1 if there is no such atom within the amino acid)
	int rapdf_id(Atom_Type t);

	// get the atom type for nth RAPDF id
	static Atom_Type rapdf_atom(int n);

	// get the amino acid type for nth RAPDF id
	static Amino rapdf_amino(int n);

	// one-letter code
	char code() const;

	// three-letter abbreviation
	const char *abbr() const;

	// full name
	const char *name() const;

	// numerical index (opposite of constructor from int)
	int num() const
	{ return m_val; }

	// convenience functions for checking the amino acid type

	bool is_glycine() const
	{ return (m_val == AA_Glycine); }

	bool is_proline() const
	{ return (m_val == AA_Proline); }

	bool is_cysteine() const
	{ return (m_val == AA_Cysteine); }

	// check if the amino acid type is one of the standard twenty
	// (ie. not ambiguous, rare or unknown)
	bool standard() const;

	// check if the amino acid type is ambiguous
	bool ambiguous() const;

	// check if the amino acid type is rare (not one of the standard
	// twenty, but not ambiguous or unknown)
	bool rare() const;

	// check if the amino acid type is unknown ("X")
	bool unknown() const;

	// check if this is an undefined amino acid
	// (eg. failed to convert from char)
	bool undefined() const
	{ return (m_val == c_undefined); }

	// make this amino acid undefined
	void make_undefined()
	{ m_val = c_undefined; }

	// equality
	bool operator == (const Amino &other) const
	{ return (m_val == other.m_val); }

	// inequality
	bool operator != (const Amino &other) const
	{ return (m_val != other.m_val); }

	// number of atoms expected in this kind of amino acid
	int expected_atoms() const;

	// nth expexted atom type in this kind of amino acid
	Atom_Id expected_atom(int n) const;

protected:
	// do required static initialisations (called by Static_Init)
	static void do_static_init();

	// call destructors for objects created in do_static_init();
	static void do_static_delete();

	// get the unique id for a particular type of atom within an
	// amino acid
	static int lookup_atom_type(const Amino &a, Atom_Type t);

private:
	// helper object for performing conversions
	//
	// (This is a pointer so that its time of creation can be managed
	// by the Static_Init class, since there is a dependency order:
	// Codon::Codon_Helper must not be before this object)
	static Amino_Helper *m_helper;

	// undefined value
	static const int c_undefined = Full_Num - 1;

	// the type of amino acid
	unsigned char m_val;
};

// a vector of amino acids
typedef std::vector<Amino> Amino_Vec;

#endif // AMINO_H_INCLUDED

