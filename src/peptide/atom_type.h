#ifndef ATOM_TYPE_H_INCLUDED
#define ATOM_TYPE_H_INCLUDED

#include <iostream>
#include <string>
#include "atom_id.h"

class Atom_Type_Helper;

// Wrapper for enum Atom_Id

class Atom_Type
{
	friend class Atom_Type_Helper;
	friend class CORE_impl;
public:
	Atom_Type();
	Atom_Type(Atom_Id id);

	// create from a 4 character PDB atom name (" CA " for alpha carbon, etc.)
	// If the atom type is not recognised (eg. hydrogen), the result will
	// be type Atom_Undef (use undefined() function to test for this)
	Atom_Type(const char *pdb_atom_name);

	void set_type(Atom_Id t)
	{ m_type = t; }

	Atom_Id type() const
	{ return m_type; }

	bool undefined() const
	{ return (m_type == Atom_Undef); }

	// check if this is a backbone atom (including C alpha & beta)
	bool is_backbone() const
	{ return (m_type < Num_Backbone); }

	// name (not including amino acid), eg. "C", "CB", "OG2"
	const char *name() const;

	// 4-character PDB name, eg. " CA "
	const char *pdb_name() const;

	/*
	// check if the atom type in a pdb format ATOM line is a hydrogen
	static bool pdb_atom_type_is_H(const char *s);
	*/

	// equality
	bool operator == (const Atom_Type &other)
	{ return (m_type == other.m_type); }

	bool operator == (Atom_Id other)
	{ return (m_type == other); }

	// inequality
	bool operator != (const Atom_Type &other)
	{ return (m_type != other.m_type); }

	bool operator != (Atom_Id other)
	{ return (m_type != other); }

	// debugging function
	void dump(std::ostream &out = std::cout) const;

	// internal initialisation/cleanup for static data members of this class
	static void do_static_init();
	static void do_static_delete();

private:
	static Atom_Type_Helper *m_helper;
	Atom_Id m_type;
};

#endif // ATOM_TYPE_H_INCLUDED
