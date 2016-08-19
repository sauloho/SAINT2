#ifndef ATOM_H_INCLUDED
#define ATOM_H_INCLUDED

#include <vector>
#include <string>
#include "amino.h"
#include "atom_type.h"
#include "pdb_atom_rec.h"
#include "point.h"

// A single, physical atom.

class Atom
{
	friend class CORE_impl;
public:
	Atom();
	Atom(Amino a, Atom_Type t);

	~Atom();

	// copy constructor and assignment operator
	Atom(const Atom &other);
	Atom &operator = (const Atom &other);

	Atom_Type type() const
	{ return m_type; }

	Atom_Id id() const
	{ return m_type.type(); }

	// returns false if the atom type does not exist in the
	// amino acid
	bool set_type(Amino a, Atom_Type t);

	// returns -1 if either the atom is undefined
	// or the atom type does not exist in the amino acid
	int rapdf_id() const
	{ return m_rapdf_id; }

	Amino amino() const
	{ return m_amino; }

	bool undefined() const
	{ return m_type.undefined(); }

	void make_undefined();

	// set the PDB file line number for the atom record
	void set_pdb_line(int line_num)
	{ m_pdb_line = line_num; }

	// get the line number in the PDB file for this atom
	int pdb_line() const
	{ return m_pdb_line; }

	// void set_record(const PDB_Atom_Rec &rec);

	// get the data from the original ATOM line in the PDB file this
	// atom was read from; returns NULL if it was not read from a PDB file
	PDB_Atom_Rec *pdb_rec() const
	{ return m_pdb_rec; }

	// set the atom's "bad position" status
	void bad_pos(bool val)
	{ m_bad_pos = val; }

	// get the atom's "bad position" status
	bool bad_pos() const
	{ return m_bad_pos; }

private:
	Amino m_amino;
	Atom_Type m_type;
	int m_rapdf_id;

	// whether the atom's position is suspicious
	bool m_bad_pos;

	// line number in PDB file (0 if not from a PDB file)
	int m_pdb_line;

	PDB_Atom_Rec *m_pdb_rec;	// data from PDB file (or NULL)
};

typedef std::vector<Atom> Atom_Vec;

#endif // ATOM_H_INCLUDED

