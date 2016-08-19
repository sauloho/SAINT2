
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <string>
#include "residue.h"

Residue::Residue()
	: m_res_seq(-1), m_i_code(' '), m_pdb_line(0), m_missing_after(false)
{
}

Residue::Residue(Amino a, int res_seq /*= -1*/, char i_code /*= ' '*/)
	: m_res_seq(res_seq), m_i_code(i_code), m_pdb_line(0),
	  m_missing_after(false)
{
	set_amino(a);
}

Residue::Residue(const Residue &other)
	: m_amino(other.m_amino), m_codon(other.m_codon), m_atom(other.m_atom),
		m_res_seq(other.m_res_seq), m_i_code(other.m_i_code),
		m_pdb_line(other.m_pdb_line),
		m_missing_after(other.m_missing_after)
{
}

Residue& Residue::operator = (const Residue &other)
{
	m_amino = other.m_amino;
	m_codon = other.m_codon;
	m_atom = other.m_atom;
	m_res_seq = other.m_res_seq;
	m_i_code = other.m_i_code;
	m_pdb_line = other.m_pdb_line;
	m_missing_after = other.m_missing_after;
	return *this;
}

void Residue::clear()
{
	m_amino.make_undefined();
	m_codon.make_unknown();
	m_atom.resize(0);
	m_res_seq = -1;
	m_i_code = ' ';
	m_pdb_line = 0;
	m_missing_after = false;
}

void Residue::set_res_seq(int res_seq, char i_code /*= ' '*/)
{
	m_res_seq = res_seq;
	m_i_code = i_code;
}

std::string Residue::res_seq_str() const
{
	char s[100];

	if (m_i_code == ' ')
	{
		sprintf(s, "%d", m_res_seq);
	}
	else
	{
		sprintf(s, "%d%c", m_res_seq, m_i_code);
	}

	return std::string(s);
}

void Residue::set_amino(Amino a)
{
	// must be one of the twenty standard amino acids
	// (not an ambiguous, rare or "unknown" amino acid)
	assert(a.standard());
	m_amino = a;

	// backbone atoms always exist, and are in a fixed order
	// (the atom types will initially be Atom_Undef)

	if (m_amino.is_glycine())
	{
		// exclude beta carbon
		m_atom.resize(Num_Backbone - 1);
	}
	else
	{
		m_atom.resize(Num_Backbone);
	}
}

void Residue::set_codon(Codon c)
{
	m_codon = c;
}

void Residue::allocate_backbone_atoms()
{
	// amino acid type needs to have been defined
	assert(m_amino.standard());

	int num = Num_Backbone;

	if (m_amino.is_glycine())
	{
		num--;
	}

	// !!!
	// if ((int) m_atom.size() < num)
	m_atom.resize(num);

	for (int n = 0;n < num;n++)
	{
		m_atom[n].set_type(m_amino, (Atom_Id) n);
	}
}

Atom *Residue::add_atom(Atom_Type t, bool ok_if_exists /*= false*/)
{
	if (m_amino.undefined())
	{
		std::cerr << "Error: called Residue::add_atom() without setting "
			"amino acid type first\n";
		exit(1);
	}

	assert(m_atom.size() != 0);
	int n = atom_index(t);

	if (n != -1)
	{
		if (!ok_if_exists)
		{
			std::cerr << "Error: tried to add duplicate "
				<< t.name() << " atom to "
				<< m_amino.name() << " residue\n";
			exit(1);
		}

		assert(m_atom[n].type() == t);
	}
	else
	{
		n = (int) t.type();

		if (n < Num_Backbone)
		{
			if (t == Atom_CB && m_amino.is_glycine())
			{
				std::cerr << "Error: tried to add CB to glycine\n";
				exit(1);
			}
	
			assert(n < (int) m_atom.size());
		}
		else
		{
			n = (int) m_atom.size();
			m_atom.resize((unsigned int) (n + 1));
		}

		m_atom[n].set_type(m_amino, t);
	}

	return &m_atom[n];
}

Atom *Residue::get_or_add_atom(Atom_Type t)
{
	return add_atom(t, true);
}

void Residue::remove_atom(Atom_Type t)
{
	assert(m_atom.size() != 0);
	int n = atom_index(t);

	if (n != -1)
	{
		m_atom[n].make_undefined();
	}
}

Atom &Residue::atom(int n)
{
	assert(n >= 0 && n < num_atoms());
	return m_atom[n];
}

const Atom &Residue::atom(int n) const
{
	assert(n >= 0 && n < num_atoms());
	return m_atom[n];
}

Atom &Residue::atom(Atom_Type t)
{
	int n = atom_index(t);

	if (n == -1)
	{
		std::cerr
			<< "Fatal error: atom not found in Residue::atom(Atom_Type)\n";
		exit(1);
	}

	return m_atom[n];
}

const Atom &Residue::atom(Atom_Type t) const
{
	// call non-const version of this function
	return (const_cast<Residue*>(this))->atom(t);
}

int Residue::atom_index(Atom_Type t) const
{
	int tt = (int) t.type();

	// backbone atoms (0 .. Num_Backbone - 1) always have the same
	// index in m_atom

	if (tt < Num_Backbone)
	{
		// The following expression checks several things at once:
		//
		// - if m_atom.size() == 0 (which will happen if the amino acid
		//   for this residue has not been set)
		// - if this is Glycine and the atom is CB (in which case tt == 5
		//   and m_atom.size() == 4)
		// - if the atom has not been added with add_atom (in which case
		//   tt < size, but the atom will be undefined)

		if (tt < (int) m_atom.size() && !m_atom[tt].undefined())
		{
			return tt;
		}

		return -1;
	}

	// search for the atom
	
	for (int n = Num_Backbone;n < (int) m_atom.size();n++)
	{
		if (m_atom[n].type() == t)
		{
			return n;
		}
	}

	return -1;
}

int Residue::num_missing_backbone() const
{
	int bb = Num_Backbone;

	if (m_amino.is_glycine())
	{
		bb--;
	}

	int total = 0;

	for (int i = 0;i < bb;i++)
	{
		Atom_Id t = (Atom_Id) i;

		if (!atom_exists(t))
		{
			total++;
		}
	}

	return total;
}

void Residue::remove_non_backbone_atoms()
{
	int bb = Num_Backbone;

	if (m_amino.is_glycine())
	{
		bb--;
	}
	
	if ((int) m_atom.size() > bb)
	{
		m_atom.resize(bb);

		// Note: resize() to a smaller size doesn't deallocate memory,
		// so copy the atoms to a temporary vector and swap their buffers

		Atom_Vec temp(m_atom);
		m_atom.swap(temp);
	}
}

