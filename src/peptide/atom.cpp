
#include <cmath>
#include <string>
#include "atom.h"

Atom::Atom()
	: m_rapdf_id(-1), m_pdb_line(0), m_pdb_rec(NULL)
{
}

Atom::Atom(Amino a, Atom_Type t)
	: m_amino(a), m_type(t), m_pdb_line(0), m_pdb_rec(NULL)
{
	m_rapdf_id = m_amino.rapdf_id(m_type);
}

Atom::Atom(const Atom &other)
{
	// call assignment operator
	*this = other;
}

Atom &Atom::operator = (const Atom &other)
{
	m_amino = other.m_amino;
	m_type = other.m_type;
	m_rapdf_id = other.m_rapdf_id;
	m_pdb_line = other.m_pdb_line;

	if (other.m_pdb_rec == NULL)
	{
		m_pdb_rec = NULL;
	}
	else
	{
		m_pdb_rec = new PDB_Atom_Rec(*other.m_pdb_rec);
	}

	return *this;
}

void Atom::make_undefined()
{
	m_type.set_type(Atom_Undef);
}

Atom::~Atom()
{
	delete m_pdb_rec;
}

bool Atom::set_type(Amino a, Atom_Type t)
{
	m_amino = a;
	m_type = t;
	m_rapdf_id = m_amino.rapdf_id(m_type);
	return (m_rapdf_id != -1);
}

