
#include <cstdlib>
#include <iostream>
#include <cctype>
#include <cassert>
#include "amino.h"
#include "codon.h"

class Codon_Helper
{
public:
	Codon_Helper();

	// get the amino acid corresponding to a codon.
	//
	// Returns a number from 0 to (Amino::Num - 1), or -1 for a stop codon
	// (the start codon is the same as Methionine)

	int get_amino_num(Codon c)
	{ return m_amino_num[c.num()]; }

	bool is_stop_codon(Codon c)
	{ return (m_amino_num[c.num()] == -1); }

	// get the start codon
	Codon start_codon() const
	{ return m_start_codon; }

private:
	static const int Groups = Amino::Num + 1;

	static const char *c_group[Groups];

	// the start codon is the same as Methionine
	static const char c_start_codon_res = 'M';

	// number of the start codon
	int m_start_codon;

	// amino acid number for each codon (-1 means stop codon)
	int m_amino_num[Codon::Num];
};

// static data members

const char *Codon_Helper::c_group[Codon_Helper::Groups] =
{
	"A GCT GCC GCA GCG",			// Alanine
	"R CGT CGC CGA CGG AGA AGG",	// Arginine
	"N AAT AAC",					// Asparagine
	"D GAT GAC",					// Aspartic acid
	"C TGT TGC",					// Cysteine
	"E GAA GAG",					// Glutamic acid
	"Q CAA CAG",					// Glutamine
	"G GGT GGC GGA GGG",			// Glycine
	"H CAT CAC",					// Histidine
	"I ATT ATC ATA",				// Isoleucine
	"L TTA TTG CTT CTC CTA CTG",	// Leucine
	"K AAA AAG",					// Lysine
	"M ATG",						// Methionine
	"F TTT TTC",					// Phenylalanine
	"P CCT CCC CCA CCG",			// Proline
	"S TCT TCC TCA TCG AGT AGC",	// Serine
	"T ACT ACC ACA ACG",			// Threonine
	"W TGG",						// Tryptophan
	"Y TAT TAC",					// Tyrosine
	"V GTT GTC GTA GTG",			// Valine
	"- TAG TGA TAA"					// stop codons

	// the start codon "ATG" is not listed because it is also
	// used for Methionine
};

Codon_Helper::Codon_Helper()
{
	static const int undef = -3;
	m_start_codon = -1;
	bool seen_aa[Amino::Num];

	int n;
	for (n = 0;n < Codon::Num;n++)
	{
		m_amino_num[n] = undef;
	}

	for (n = 0;n < Amino::Num;n++)
	{
		seen_aa[n] = false;
	}

	int count = 0;

	for (int g = 0;g < Groups;g++)
	{
		const char *s = c_group[g];
		bool is_start_codon = (*s == c_start_codon_res);
		int num;

		if (*s == '-')
		{
			// stop codon
			num = -1;
		}
		else
		{
			Amino a(*s);
			if (!a.standard())
			{
				std::cerr << "Internal error in Codon_Helper::Codon_Helper(): "
					"illegal residue type '" << *s << "'\n";
				exit(1);
			}

			num = a.num();

			if (seen_aa[num])
			{
				std::cerr << "Internal error in Codon_Helper::Codon_Helper(): "
					"duplicate residue type '" << *s << "' ["
					<< num << "]\n";
				exit(1);
			}

			seen_aa[num] = true;
		}

		for (s++;*s;s += 3)
		{
			assert(*s == ' ');
			s++;

			Codon c(s);

			if (c.illegal())
			{
				std::cerr << "Internal error in Codon_Helper::Codon_Helper(): "
					"illegal codon '" << s[0] << s[1] << s[2] << "'\n";
				exit(1);
			}

			int c_num = c.num();

			if (m_amino_num[c_num] != undef)
			{
				std::cerr << "Internal error in Codon_Helper::Codon_Helper(): "
					"duplicate codon '" << s[0] << s[1] << s[2] << "'\n";
				exit(1);
			}

			m_amino_num[c_num] = num;
			count++;

			if (is_start_codon)
			{
				m_start_codon = c_num;
			}
		}
	}

	if (count != Codon::Num)
	{
		std::cerr << "Internal error in Codon_Helper::Codon_Helper(): "
			"only " << count << " codons defined (<" << Codon::Num << ")\n";
		exit(1);
	}

	for (n = 0;n < Amino::Num;n++)
	{
		if (!seen_aa[n])
		{
			std::cerr << "Internal error in Codon_Helper::Codon_Helper(): "
				"no codons for residue type '" << Amino(n).code() << "'\n";
			exit(1);
		}
	}

	if (m_start_codon == -1)
	{
		std::cerr << "Internal error in Codon_Helper::Codon_Helper(): "
			"start codon undefined\n";
		exit(1);
	}
}

//////////////////////////////////////////////////////////////////

// static data members
Codon_Helper *Codon::m_helper = NULL;

Codon::Codon()
	: m_val(c_illegal)
{
}

Codon::Codon(int n)
	: m_val(n)
{
	assert(n >= 0 && n < Num);
}

Codon::Codon(char nuc0, char nuc1, char nuc2)
{
	m_val = nucs_to_num(nuc0, nuc1, nuc2);
}

Codon::Codon(const char *str)
{
	assert(str[0] != '\0' && str[1] != '\0' && str[2] != '\0');
	m_val = nucs_to_num(str[0], str[1], str[2]);
}

int Codon::nucs_to_num(char nuc0, char nuc1, char nuc2)
{
	int n0 = nuc_to_num(nuc0);
	int n1 = nuc_to_num(nuc1);
	int n2 = nuc_to_num(nuc2);

	if (n0 == -1 || n1 == -1 || n2 == -1)
	{
		return c_illegal;
	}
	else
	{
		return n0 * 16 + n1 * 4 + n2;
	}
}

int Codon::nuc_to_num(char nuc)
{
	switch (tolower(nuc))
	{
		case 'a':
			return 0;
		case 'c':
			return 1;
		case 'g':
			return 2;
		case 't':
		case 'u':
			return 3;
		default:
			return -1;
	}
}

Amino Codon::to_amino() const
{
	Amino undef;				// ("undefined" amino acid by default)
	assert(undef.undefined());

	if (m_val == c_illegal)
	{
		return undef;
	}

	int n = m_helper->get_amino_num(m_val);

	if (n == -1)
	{
		return undef;
	}

	return Amino(n);
}

char Codon::nucleotide(int n, char t_u_value /*= 'T'*/) const
{
	assert(n >= 0 && n <= 2);
	assert(t_u_value == 'T' || t_u_value == 'U');

	if (m_val == c_illegal || n < 0 || n > 2)
	{
		return '?';
	}
	else
	{
		int v;

		// extract the appropriate two bits
		switch (n)
		{
			case 0:
				v = (m_val >> 4) & 0x03;
				break;

			case 1:
				v = (m_val >> 2) & 0x03;
				break;

			default:
				v = m_val & 0x03;
		}

		if (v == 3)
		{
			return t_u_value;
		}
		else
		{
			return "ACG"[v];
		}
	}
}

std::string Codon::str(char t_u_value /*= 'T'*/) const
{
	if (m_val == c_illegal)
	{
		return "???";
	}

	char buff[3];

	for (int n = 0;n < 3;n++)
	{
		buff[n] = nucleotide(n, t_u_value);
	}

	return std::string(buff, 3);
}

bool Codon::is_stop_codon() const
{
	return m_helper->is_stop_codon(m_val);
}

Codon Codon::start_codon()
{
	return m_helper->start_codon();
}

void Codon::do_static_init()
{
	m_helper = new Codon_Helper;
}

void Codon::do_static_delete()
{
	delete m_helper;
}

///////////////////////////////////////////////////////////////////
// For convenience -- a complete list of codons in ascending order:
//
// AAA - K
// AAC - N
// AAG - K
// AAT - N
// ACA - T
// ACC - T
// ACG - T
// ACT - T
// AGA - R
// AGC - S
// AGG - R
// AGT - S
// ATA - I
// ATC - I
// ATG - M
// ATT - I
// CAA - Q
// CAC - H
// CAG - Q
// CAT - H
// CCA - P
// CCC - P
// CCG - P
// CCT - P
// CGA - R
// CGC - R
// CGG - R
// CGT - R
// CTA - L
// CTC - L
// CTG - L
// CTT - L
// GAA - E
// GAC - D
// GAG - E
// GAT - D
// GCA - A
// GCC - A
// GCG - A
// GCT - A
// GGA - G
// GGC - G
// GGG - G
// GGT - G
// GTA - V
// GTC - V
// GTG - V
// GTT - V
// TAA - stop
// TAC - Y
// TAG - stop
// TAT - Y
// TCA - S
// TCC - S
// TCG - S
// TCT - S
// TGA - stop
// TGC - C
// TGG - W
// TGT - C
// TTA - L
// TTC - F
// TTG - L
// TTT - F

