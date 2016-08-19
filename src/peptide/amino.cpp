#include <strings.h> // PG added this
#include <iostream>
#include <cstdlib>
#include <cctype>
#include <cassert>
#include "amino.h"
#include "atom_type.h"

// Helper for class Amino (performs conversions between types).

class Amino_Helper
{
	friend class Amino;

public:
	// amino acid classifications
	enum Classif { Standard, Ambig, Rare, Unknown, Undefined };

public:
	Amino_Helper();

	int char_to_val(char ch)
	{
		return m_lookup[(unsigned char) ch];
	}

	int abbr_to_val(const char *abbr_str);

	char val_to_char(int val)
	{
		assert(val >= 0);
		return c_desc[val].code;
	}

	const char *val_to_abbr(int val)
	{
		assert(val >= 0);
		return c_desc[val].abbr;
	}

	const char *val_to_name(int val)
	{
		assert(val >= 0);
		return c_desc[val].name;
	}

	// check if amino acid type val has a particular classification
	bool classification(int val, Classif c)
	{
		assert(val >= 0);
		return (c_desc[val].classif == c);
	}

	int lookup_atom_type(Amino a, Atom_Type t);

	int get_num_atoms(int val);
	Atom_Id get_expected_atom(int val, int n);

private:
	// maximum value for a char
	static const int c_char_max = 256;

	// lookup table (character value => amino acid number)
	static int m_lookup[c_char_max];

	// descriptor for an amino acid type
	struct Desc
	{
		const char *name;
		const char *abbr;
		char code;
		Classif classif;
		int aa_index;		// index of start of amino acid in aa_atom[] array
		int num_atoms;		// number of atoms (including backbone)
	};

	static Desc c_desc[Amino::Full_Num];

protected:
	// the order of this array is the order required by
	// the RAPDF scoring function

	// amino acid atoms

	struct AA_Atom
	{
		char aa;			// amino acid (code)
		Atom_Id type;		// C, CA, etc.
	};

	static AA_Atom aa_atom[];
	static int c_num_rapdf_ids;
};

// static data members
int Amino_Helper::m_lookup[c_char_max];

Amino_Helper::Desc Amino_Helper::c_desc[Amino::Full_Num] =
{
	// standard amino acid types
	// (the last two values, aa_index and num_atoms, are calculated
	// on initialisation from the aa_atom[] array)

	{ "Alanine",		"ALA", 'A', Standard, 0, 0 },
	{ "Arginine",		"ARG", 'R', Standard, 0, 0 },
	{ "Asparagine",		"ASN", 'N', Standard, 0, 0 },
	{ "Aspartic acid",	"ASP", 'D', Standard, 0, 0 },
	// constant AA_Cysteine must be equal to this index
	{ "Cysteine",		"CYS", 'C', Standard, 0, 0 },
	{ "Glutamic acid",	"GLU", 'E', Standard, 0, 0 },
	{ "Glutamine",		"GLN", 'Q', Standard, 0, 0 },
	// constant AA_Glycine must be equal to this index
	{ "Glycine",		"GLY", 'G', Standard, 0, 0 },
	{ "Histidine",		"HIS", 'H', Standard, 0, 0 },
	{ "Isoleucine",		"ILE", 'I', Standard, 0, 0 },
	{ "Leucine",		"LEU", 'L', Standard, 0, 0 },
	{ "Lysine",			"LYS", 'K', Standard, 0, 0 },
	{ "Methionine",		"MET", 'M', Standard, 0, 0 },
	{ "Phenylalanine",	"PHE", 'F', Standard, 0, 0 },
	// constant AA_Proline must be equal to this index
	{ "Proline",		"PRO", 'P', Standard, 0, 0 },
	{ "Serine",			"SER", 'S', Standard, 0, 0 },
	{ "Threonine",		"THR", 'T', Standard, 0, 0 },
	{ "Tryptophan",		"TRP", 'W', Standard, 0, 0 },
	{ "Tyrosine",		"TYR", 'Y', Standard, 0, 0 },
	{ "Valine",			"VAL", 'V', Standard, 0, 0 },

	// ambiguous types

	{ "Asparagine or Aspartic acid",	"Asx", 'B', Ambig, 0, 0 },
	{ "Glutamine or Glutamic acid",		"Glx", 'Z', Ambig, 0, 0 },
	{ "Leucine or Isoleucine",			"Xle", 'J', Ambig, 0, 0 },

	// unknown type

	{ "Unknown",		"Xaa", 'X', Unknown, 0, 0 },

	// rare amino acids

	{ "Pyrrolysine",	"Pyl", 'O', Rare, 0, 0 },
	{ "Selenocysteine",	"Sec", 'U', Rare, 0, 0},

	// undefined type (must come last, ie. have index value Full_Num)

	{ "Undefined",		"???", '?', Undefined, 0, 0 }
};

// Note: the order of aa_atom[] entries corresponds to the "atom id"
// required by the RAPDF scoring function

Amino_Helper::AA_Atom Amino_Helper::aa_atom[] =
{
	// Cysteine
	
	{ 'C', Atom_N },	// 0
	{ 'C', Atom_CA },
	{ 'C', Atom_C },
	{ 'C', Atom_O },
	{ 'C', Atom_CB },
	{ 'C', Atom_SG },
	
	// Glutamine
	
	{ 'Q', Atom_N },	// 6
	{ 'Q', Atom_CA },
	{ 'Q', Atom_C },
	{ 'Q', Atom_O },
	{ 'Q', Atom_CB },
	{ 'Q', Atom_CG },
	{ 'Q', Atom_CD },
	{ 'Q', Atom_OE1 },
	{ 'Q', Atom_NE2 },
	
	// Aspartic acid
	
	{ 'D', Atom_N },	// 15
	{ 'D', Atom_CA },
	{ 'D', Atom_C },
	{ 'D', Atom_O },
	{ 'D', Atom_CB },
	{ 'D', Atom_CG },
	{ 'D', Atom_OD1 },
	{ 'D', Atom_OD2 },
	
	// Serine
	
	{ 'S', Atom_N },	// 23
	{ 'S', Atom_CA },
	{ 'S', Atom_C },
	{ 'S', Atom_O },
	{ 'S', Atom_CB },
	{ 'S', Atom_OG },
	
	// Valine
	
	{ 'V', Atom_N },	// 29
	{ 'V', Atom_CA },
	{ 'V', Atom_C },
	{ 'V', Atom_O },
	{ 'V', Atom_CB },
	{ 'V', Atom_CG1 },
	{ 'V', Atom_CG2 },
	
	// Methionine
	
	{ 'M', Atom_N },	// 36
	{ 'M', Atom_CA },
	{ 'M', Atom_C },
	{ 'M', Atom_O },
	{ 'M', Atom_CB },
	{ 'M', Atom_CG },
	{ 'M', Atom_SD },
	{ 'M', Atom_CE },
	
	// Proline
	
	{ 'P', Atom_N },	// 44
	{ 'P', Atom_CA },
	{ 'P', Atom_C },
	{ 'P', Atom_O },
	{ 'P', Atom_CB },
	{ 'P', Atom_CG },
	{ 'P', Atom_CD },
	
	// Lysine
	
	{ 'K', Atom_N },	// 51
	{ 'K', Atom_CA },
	{ 'K', Atom_C },
	{ 'K', Atom_O },
	{ 'K', Atom_CB },
	{ 'K', Atom_CG },
	{ 'K', Atom_CD },
	{ 'K', Atom_CE },
	{ 'K', Atom_NZ },
	
	// Threonine
	
	{ 'T', Atom_N },	// 60
	{ 'T', Atom_CA },
	{ 'T', Atom_C },
	{ 'T', Atom_O },
	{ 'T', Atom_CB },
	{ 'T', Atom_OG1 },
	{ 'T', Atom_CG2 },
	
	// Phenylalanine
	
	{ 'F', Atom_N },	// 67
	{ 'F', Atom_CA },
	{ 'F', Atom_C },
	{ 'F', Atom_O },
	{ 'F', Atom_CB },
	{ 'F', Atom_CG },
	{ 'F', Atom_CD1 },
	{ 'F', Atom_CD2 },
	{ 'F', Atom_CE1 },
	{ 'F', Atom_CE2 },
	{ 'F', Atom_CZ },
	
	// Alanine
	
	{ 'A', Atom_N },	// 78
	{ 'A', Atom_CA },
	{ 'A', Atom_C },
	{ 'A', Atom_O },
	{ 'A', Atom_CB },
	
	// Histidine
	
	{ 'H', Atom_N },	// 83
	{ 'H', Atom_CA },
	{ 'H', Atom_C },
	{ 'H', Atom_O },
	{ 'H', Atom_CB },
	{ 'H', Atom_CG },
	{ 'H', Atom_ND1 },
	{ 'H', Atom_CD2 },
	{ 'H', Atom_CE1 },
	{ 'H', Atom_NE2 },
	
	// Glycine
	
	{ 'G', Atom_N },	// 93
	{ 'G', Atom_CA },
	{ 'G', Atom_C },
	{ 'G', Atom_O },
	
	// Isoleucine
	
	{ 'I', Atom_N },	// 97
	{ 'I', Atom_CA },
	{ 'I', Atom_C },
	{ 'I', Atom_O },
	{ 'I', Atom_CB },
	{ 'I', Atom_CG1 },
	{ 'I', Atom_CG2 },
	{ 'I', Atom_CD1 },
	
	// Glutamic acid
	
	{ 'E', Atom_N },	// 105
	{ 'E', Atom_CA },
	{ 'E', Atom_C },
	{ 'E', Atom_O },
	{ 'E', Atom_CB },
	{ 'E', Atom_CG },
	{ 'E', Atom_CD },
	{ 'E', Atom_OE1 },
	{ 'E', Atom_OE2 },
	
	// Leucine
	
	{ 'L', Atom_N },	// 114
	{ 'L', Atom_CA },
	{ 'L', Atom_C },
	{ 'L', Atom_O },
	{ 'L', Atom_CB },
	{ 'L', Atom_CG },
	{ 'L', Atom_CD1 },
	{ 'L', Atom_CD2 },
	
	// Arginine
	
	{ 'R', Atom_N },	// 122
	{ 'R', Atom_CA },
	{ 'R', Atom_C },
	{ 'R', Atom_O },
	{ 'R', Atom_CB },
	{ 'R', Atom_CG },
	{ 'R', Atom_CD },
	{ 'R', Atom_NE },
	{ 'R', Atom_CZ },
	{ 'R', Atom_NH1 },
	{ 'R', Atom_NH2 },
	
	// Tryptophan
	
	{ 'W', Atom_N },	// 133
	{ 'W', Atom_CA },
	{ 'W', Atom_C },
	{ 'W', Atom_O },
	{ 'W', Atom_CB },
	{ 'W', Atom_CG },
	{ 'W', Atom_CD1 },
	{ 'W', Atom_CD2 },
	{ 'W', Atom_NE1 },
	{ 'W', Atom_CE2 },
	{ 'W', Atom_CE3 },
	{ 'W', Atom_CZ2 },
	{ 'W', Atom_CZ3 },
	{ 'W', Atom_CH2 },
	
	// Asparagine
	
	{ 'N', Atom_N },	// 147
	{ 'N', Atom_CA },
	{ 'N', Atom_C },
	{ 'N', Atom_O },
	{ 'N', Atom_CB },
	{ 'N', Atom_CG },
	{ 'N', Atom_OD1 },
	{ 'N', Atom_ND2 },
	
	// Tyrosine
	
	{ 'Y', Atom_N },	// 155
	{ 'Y', Atom_CA },
	{ 'Y', Atom_C },
	{ 'Y', Atom_O },
	{ 'Y', Atom_CB },
	{ 'Y', Atom_CG },
	{ 'Y', Atom_CD1 },
	{ 'Y', Atom_CD2 },
	{ 'Y', Atom_CE1 },
	{ 'Y', Atom_CE2 },
	{ 'Y', Atom_CZ },
	{ 'Y', Atom_OH },

	// end marker

	{ '\0', Atom_C }	// 167
};

Atom_Type Amino::rapdf_atom(int n)
{
	return Atom_Type(Amino_Helper::aa_atom[n].type);
}

Amino Amino::rapdf_amino(int n)
{
	return Amino(Amino_Helper::aa_atom[n].aa);
}

Amino_Helper::Amino_Helper()
{
	// make sure "unknown" amino acid type has the right value

	if (c_desc[Amino::c_undefined].classif != Undefined)
	{
		std::cerr << "Internal error in Amino_Helper::Amino_Helper(): "
			"\"Undefined\" amino acid type has incorrect index\n";
		exit(1);
	}

	// initialise lookup table

	for (int n = 0;n < c_char_max;n++)
	{
		m_lookup[n] = Amino::c_undefined;
	}

	for (int a = 0;a < Amino::Full_Num;a++)
	{
		if (!isalpha(c_desc[a].code) && c_desc[a].code != '?')
		{
			std::cerr << "Internal error in Amino_Helper::Amino_Helper(): "
				"non-alphabetical amino acid code '" << c_desc[a].code << "'\n";
			exit(1);
		}

		unsigned char ch = c_desc[a].code;

		if (m_lookup[ch] != Amino::c_undefined)
		{
			std::cerr << "Internal error in Amino_Helper::Amino_Helper(): "
				"duplicate amino acid code '" << c_desc[a].code << "'\n";
			exit(1);
		}

		m_lookup[tolower(ch)] = a;
		m_lookup[toupper(ch)] = a;
	}

	// set aa_index and num_atoms for each amino acid type

	unsigned char prev_aa = '\0';
	int val = -1;
	int prev_val = -1;

	int m;
	for (m = 0;aa_atom[m].aa != '\0';m++)
	{
		if (aa_atom[m].aa != prev_aa)
		{
			if (prev_val != -1)
			{
				c_desc[val].num_atoms = m - c_desc[prev_val].aa_index;
			}

			prev_aa = aa_atom[m].aa;
			val = m_lookup[prev_aa];
			c_desc[val].aa_index = m;

			prev_val = val;
		}
	}

	assert(val != -1 && prev_val != -1);
	c_desc[val].num_atoms = m - c_desc[prev_val].aa_index;

	assert(m == NUM_RAPDF_IDS);

	// make sure all values are set

	for (int i = 0;i < Amino::Num;i++)
	{
		assert(c_desc[i].num_atoms > 0);
	}
}

int Amino_Helper::lookup_atom_type(Amino a, Atom_Type t)
{
	int n = c_desc[a.m_val].aa_index;
	int end = n + c_desc[a.m_val].num_atoms;

// std::cout << "[lookup " << a.code() << ", pos "
// << (int) p << ", element " << (int) e
// << ", sub " << subscript << "]\n";

	for ( ;n < end;n++)
	{
		if (aa_atom[n].type == t.type())
		{
			return n;
		}
	}

// std::cout << "[undefined]\n";
	return -1;
}

int Amino_Helper::get_num_atoms(int val)
{
	return c_desc[val].num_atoms;
}

Atom_Id Amino_Helper::get_expected_atom(int val, int n)
{
	int i = c_desc[val].aa_index;
	return aa_atom[i + n].type;
}

int Amino_Helper::abbr_to_val(const char *abbr_str)
{
	// (this could be speeded up by replacing the linear search with
	// a std::map, but speed probably isn't important for this function)

	for (int n = 0;n < Amino::Full_Num;n++)
	{
		if (strcasecmp(abbr_str, c_desc[n].abbr) == 0)
		{
			return n;
		}
	}

	return Amino::c_undefined;
}

/////////////////////////////////////////////////////////////////////

// static data members
Amino_Helper *Amino::m_helper = NULL;

Amino::Amino(char ch)
	: m_val(m_helper->char_to_val(ch))
{
}

Amino::Amino(int val)
	: m_val(val)
{
	assert(val >= 0 && val < Full_Num);
}

Amino::Amino(const char *abbr_str)
{
	m_val = m_helper->abbr_to_val(abbr_str);
}

int Amino::rapdf_id(Atom_Type t)
{
	return m_helper->lookup_atom_type(m_val, t);
}

char Amino::code() const
{
	return m_helper->val_to_char(m_val);
}

const char *Amino::abbr() const
{
	return m_helper->val_to_abbr(m_val);
}

const char *Amino::name() const
{
	return m_helper->val_to_name(m_val);
}

bool Amino::standard() const
{
	return m_helper->classification(m_val, Amino_Helper::Standard);
}

bool Amino::ambiguous() const
{
	return m_helper->classification(m_val, Amino_Helper::Ambig);
}

bool Amino::rare() const
{
	return m_helper->classification(m_val, Amino_Helper::Rare);
}

bool Amino::unknown() const
{
	return m_helper->classification(m_val, Amino_Helper::Unknown);
}

int Amino::expected_atoms() const
{
	return m_helper->get_num_atoms(m_val);
}

Atom_Id Amino::expected_atom(int n) const
{
	assert(n >= 0 && n < expected_atoms());
	return m_helper->get_expected_atom(m_val, n);
}

void Amino::do_static_init()
{
	m_helper = new Amino_Helper;
}

void Amino::do_static_delete()
{
	delete m_helper;
}

