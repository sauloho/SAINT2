#ifndef RESIDUE_H_INCLUDED
#define RESIDUE_H_INCLUDED

#include <vector>
#include "amino.h"
#include "codon.h"
#include "atom.h"
#include "atom_type.h"

// A single residue.

class Residue
{
	friend class Peptide;
	friend class CORE_impl;
public:
	Residue();
	Residue(Amino a, int res_seq = -1, char i_code = ' ');

	Residue(const Residue &other);
	Residue& operator = (const Residue &other);

	// set the amino acid value
	void set_amino(Amino a);

	// set the codon value
	void set_codon(Codon c);

	// get the amino acid for this residue
	Amino amino() const
	{ return m_amino; }

	// get the codon for this residue (may be unknown)
	Codon codon() const
	{ return m_codon; }

	// reset all values
	void clear();

	// set the residue sequence number
	void set_res_seq(int res_seq, char i_code = ' ');

	int res_seq() const
	{ return m_res_seq; }

	char i_code() const
	{ return m_i_code; }

	// get a string made up of the sequence number and i-code
	std::string res_seq_str() const;

	// set the line number in the PDB file of the start of the residue
	void set_pdb_line(int line_num)
	{ m_pdb_line = line_num; }

	int pdb_line() const
	{ return m_pdb_line; }

	// number of atoms in this residue
	int num_atoms() const
	{ return (int) m_atom.size(); }

	// add all backbone atoms (if they do not exist already)
	void allocate_backbone_atoms();

	// add a new atom to the residue (returns a pointer to it)
	// If ok_if_exists is false, it is an error if an atom of that type
	// already exists in the residue (use atom_exists() to check first).
	// The amino acid type for the residue must have been set before
	// calling this function.
	Atom *add_atom(Atom_Type t, bool ok_if_exists = false);

	// get a pointer to a particular atom type in the residue;
	// the atom is created if it does not already exist.
	// The amino acid type for the residue must have been set before
	// calling this function.
	Atom *get_or_add_atom(Atom_Type t);

	// remove an atom from the residue
	void remove_atom(Atom_Type t);

	// get the nth atom
	const Atom &atom(int n) const;
	Atom &atom(int n);

	// find the index of a particular atom type.
	// Returns -1 if the atom does not exist (always
	// -1 if the amino acid type has not yet been set)
	int atom_index(Atom_Type t) const;

	// check if a particular atom type is present
	bool atom_exists(Atom_Type t) const
	{ return (atom_index(t) != -1); }

	// get the number of missing backbone atoms in this residue
	int num_missing_backbone() const;

	// get the atom with a particular type
	const Atom &atom(Atom_Type t) const;
	Atom &atom(Atom_Type t);

	void remove_non_backbone_atoms();

	// set the "missing residues(s) between this one and the next one" flag
	void set_missing_after(bool val)
	{ m_missing_after = val; }

	// get the value of the "missing residues(s) between this one and the
	// next one" flag
	bool missing_after() const
	{ return m_missing_after; }

private:
	// amino acid for this residue (never unknown)
	Amino m_amino;

	// codon for this residue (if unknown, m_codon.unknown() returns true)
	Codon m_codon;

	// all atoms for this residue type
	// (the first five atoms are always N, CA, C, O and CB (except for Glycine
	// which has no CB); other atoms (if present) may be in any order
	Atom_Vec m_atom;

	// sequence number in PDB file (or -1 if not from a PDB file)
	int m_res_seq;

	// i-code (letter after sequence number) (or ' ' if not from a PDB file)
	int m_i_code;

	// line number of start of residue in PDB file (0 if not from a PDB file)
	int m_pdb_line;

	// whether there are residues missing between this one and the next one
	bool m_missing_after;
};

typedef std::vector<Residue> Residue_Vec;

#endif // RESIDUE_H_INCLUDED

