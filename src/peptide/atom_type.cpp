
#include <cassert>
#include <iostream>
#include <cstring>
#include <string>
#include <map>
#include "atom_type.h"

class Atom_Type_Helper
{
public:
	Atom_Type_Helper();

	const char *name(Atom_Id id);
	const char *pdb_name(Atom_Id id);
	Atom_Id lookup(const std::string &atom_name);

private:
	// atom names (variable length)
	static const char *m_atom_name[Num_Atom_Ids];

	// 4-letter (PDB) atom names
	std::string m_atom_pdb_name[Num_Atom_Ids];

	// map from 4-letter (PDB) atom names to Atom_Ids
	typedef std::map<std::string,Atom_Id> Atom_Id_Map;
	Atom_Id_Map m_id_map;
};

// all atom names (same order as enum Atom_Id in atom_type.h)
const char *Atom_Type_Helper::m_atom_name[Num_Atom_Ids] =
{
	"N", "CA", "C", "O", "CB", "CG", "CG1", "CG2", "CD", "CD1",
	"CD2", "CE", "CE1", "CE2", "CE3", "CZ", "CZ2", "CZ3", "CH2",
	"ND1", "ND2", "NE", "NE1", "NE2", "NZ", "NH1", "NH2", "OG",
	"OG1", "OD1", "OD2", "OE1", "OE2", "OH", "SG", "SD", "???"
};

Atom_Type_Helper::Atom_Type_Helper()
{
	for (int n = 0;n < Num_Atom_Ids;n++)
	{
		// create the 4-letter atom name
		char name[10];

		int len = strlen(m_atom_name[n]);
		assert(len >= 1 && len <= 4);

		if (len == 4)
		{
			strcpy(name, m_atom_name[n]);
		}
		else
		{
			strcpy(name, "    ");
			strncpy(name + 1, m_atom_name[n], strlen(m_atom_name[n]));
		}

		m_atom_pdb_name[n] = name;
		m_id_map[name] = (Atom_Id) n;
	}
}

const char *Atom_Type_Helper::name(Atom_Id id)
{
	assert(id >= 0 && id < Num_Atom_Ids);
	return m_atom_name[(int) id];
}

const char *Atom_Type_Helper::pdb_name(Atom_Id id)
{
	assert(id >= 0 && id < Num_Atom_Ids);
	return m_atom_pdb_name[(int) id].c_str();
}

Atom_Id Atom_Type_Helper::lookup(const std::string &atom_name)
{
	Atom_Id_Map::const_iterator i = m_id_map.find(atom_name);

	if (i == m_id_map.end())
	{
		return Atom_Undef;
	}
	else
	{
		return i->second;
	}
}

////////////////////////////////////////////////////

// static member variables
Atom_Type_Helper *Atom_Type::m_helper = NULL;

void Atom_Type::do_static_init()
{
	m_helper = new Atom_Type_Helper;
}

void Atom_Type::do_static_delete()
{
	delete m_helper;
}

Atom_Type::Atom_Type()
	: m_type(Atom_Undef)
{
}

Atom_Type::Atom_Type(Atom_Id id)
	: m_type(id)
{
}

Atom_Type::Atom_Type(const char *pdb_atom_name)
{
	m_type = m_helper->lookup(pdb_atom_name);
}

const char *Atom_Type::name() const
{
	return m_helper->name(m_type);
}

const char *Atom_Type::pdb_name() const
{
	return m_helper->pdb_name(m_type);
}

void Atom_Type::dump(std::ostream &out /*= std::cout*/) const
{
	out << "<" << name() << " " << (int) m_type << ">";
}

