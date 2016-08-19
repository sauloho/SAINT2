
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cstring>
#include <cctype>
#include "common.h"
#include "c_file.h"
#include "peptide.h"
#include "sequence.h"

// static data members

const char *Sequence::c_config_section = "sequence";

const char *Sequence::c_type_amino = "amino";
const char *Sequence::c_type_nuc = "nucleotide";

const char *Sequence::c_param_type = "type";
const char *Sequence::c_param_file = "file";
const char *Sequence::c_param_seq = "seq";
const char *Sequence::c_param_fasta_index = "index";
const char *Sequence::c_param_fasta_id = "id";

const char *Sequence::c_default_type = c_type_amino;
const int Sequence::c_default_fasta_index = 0;

Sequence::Sequence()
{
}

Sequence::Sequence(const Config &config)
{
	config.create_sequence(this);
}

Sequence::~Sequence()
{
}

void Sequence::create_from_peptide(const Peptide &p)
{
	bool codons_known = true;   // set to false if any codons are unknown
	clear_seq();

	for (int n = 0;n < p.full_length();n++)
	{
		m_amino.push_back(p.res(n).amino());
		Codon c = p.res(n).codon();

		if (codons_known)
		{
			if (c.unknown())
			{
				codons_known = false;
			}
			else
			{
				m_codon.push_back(c);
			}
		}
	}

	if (!codons_known)
	{
		m_codon.clear();
	}
}

void Sequence::print(std::ostream &out /*= std::cout*/,
	int aminos_per_line /*= 0*/)
{
	for (unsigned int n = 0;n < m_amino.size();n++)
	{
		if (n != 0 && aminos_per_line != 0 &&
			n % aminos_per_line == 0)
		{
			out << '\n';
		}

		out << m_amino[n].code();
	}

	out << '\n';
}

void Sequence::print_codons(std::ostream &out /*= std::cout*/,
	int codons_per_line /*= 0*/, char separator /*= '\0'*/,
	bool sep_at_end_of_line /*= false*/)
{
	if (!codons_known())
	{
		out << "[codon sequence unknown]\n";
		return;
	}

	for (unsigned int n = 0;n < m_codon.size();n++)
	{
		if (n != 0)
		{
			if (codons_per_line != 0 &&
				n % codons_per_line == 0)
			{
				if (sep_at_end_of_line)
				{
					out << separator;
				}

				out << '\n';
			}
			else
			if (separator != '\0')
			{
				out << separator;
			}
		}

		out << m_codon[n].str();
	}

	out << '\n';
}

void Sequence::create_from_params(const Param_List &params)
{
	std::string type = c_default_type;
	std::string filename;
	std::string seq;
	int fasta_index = 0;
	bool fasta_index_found = false;
	std::string fasta_id;

	Param_List::const_iterator i = params.begin();
	Param_List::const_iterator end = params.end();

	if (i == end)
	{
		std::cerr << Config::cmd() << ": sequence required\n";
		exit(1);
	}

	// parse the parameters

	for ( ;i != end;++i)
	{
        std::string value = tolower(i->value);
        std::string value_case_unchanged = i->value;

		std::string full_name = c_config_section;
		full_name += " ";
		full_name += i->name;

		if (i->name == c_param_type)
		{
			type = value;
		}
		else
		if (i->name == c_param_file)
		{
			if (!seq.empty())
			{
				std::cerr
					<< Config::cmd()
					<< ": warning: literal sequence being ignored, "
						"since sequence filename also specified\n";
				seq.clear();
			}

			filename = value_case_unchanged;
		}
		else
		if (i->name == c_param_seq)
		{
			if (!filename.empty())
			{
				std::cerr
					<< Config::cmd()
					<< ": warning: sequence filename being ignored, "
						"since literal sequence also specified\n";
				filename.clear();
			}
			seq = value;
		}
		else
		if (i->name == c_param_fasta_index)
		{
			fasta_index = parse_integer(i->value, full_name, 0);
			fasta_index_found = true;

			if (!fasta_id.empty())
			{
				fasta_id.clear();
				std::cerr
					<< Config::cmd()
					<< ": warning: FASTA id being ignored, "
						"since FASTA index also specified\n";
			}
		}
		else
		if (i->name == c_param_fasta_id)
		{
			fasta_id = i->value;

			if (fasta_index_found)
			{
				fasta_index_found = false;
				std::cerr
					<< Config::cmd()
					<< ": warning: FASTA index being ignored, "
						"since FASTA id also specified\n";
			}
		}
		else
		{
			std::cerr << Config::cmd()
				<< ": unknown " << c_config_section << " parameter \""
				<< i->name << "\"\n";
			exit(1);
		}
	}

	// don't bother checking if a FASTA index or id was specified
	// but there is no filename, since it isn't a serious error

	create_sequence(type, filename, seq, fasta_index, fasta_id);
}

void Sequence::create_sequence(const std::string &type,
	const std::string &filename, const std::string &seq,
	int fasta_index, const std::string &fasta_id)
{
	// allow abbreviations for type

	if (strncmp(type.c_str(), c_type_amino, type.length()) == 0)
	{
		if (!filename.empty())
		{
			read_amino_seq(filename, fasta_index, fasta_id);
		}
		else
		if (!seq.empty())
		{
			create_amino_seq(seq);
		}
		else
		{
			// exit with error message
			missing_sequence_err();
		}
	}
	else
	if (strncmp(type.c_str(), c_type_nuc, type.length()) == 0)
	{
		if (!filename.empty())
		{
			read_codon_seq(filename, fasta_index, fasta_id);
		}
		else
		if (!seq.empty())
		{
			create_codon_seq(seq);
		}
		else
		{
			// exit with error message
			missing_sequence_err();
		}
	}
	else
	{
		std::cerr << Config::cmd()
			<< ": unknown " << c_config_section << " type \""
			<< type << "\"\n";
		exit(1);
	}
}

void Sequence::missing_sequence_err()
{
	std::cerr << Config::cmd()
		<< ": " << c_config_section
		<< " parameters must include either a filename (\""
		<< c_param_file << " = ...\") or the sequence itself (\""
		<< c_param_seq << " = ...\")\n";
	exit(1);
}

void Sequence::read_amino_seq(const std::string &filename,
	int fasta_index, const std::string &fasta_id)
{
	std::string seq;

	C_File file(filename, "r", "amino acid sequence file");
	read_fasta_seq(file, fasta_index, fasta_id, &seq);
	create_amino_seq(seq);
}

void Sequence::read_amino_seq(const std::string &filename)
{
	read_amino_seq(filename, 0, "");
}

void Sequence::read_codon_seq(const std::string &filename,
	int fasta_index, const std::string &fasta_id)
{
	std::string seq;

	C_File file(filename, "r", "codon sequence file");
	read_fasta_seq(file, fasta_index, fasta_id, &seq);
	create_codon_seq(seq);
}

void Sequence::read_codon_seq(const std::string &filename)
{
	read_codon_seq(filename, 0, "");
}

void Sequence::read_fasta_seq(C_File &file, int index,
	const std::string &id, std::string *seq)
{
	seq->clear();
	seq->reserve(1000);

	int curr_index = -1;
	bool found = false;
	bool reading_seq = false;
	int ch;

	while ((ch = getc(file)) != EOF)
	{
		if (ch == '>')
		{
			if (reading_seq)
			{
				// previous sequence has finished being read
				return;
			}

			curr_index++;

			if (!id.empty())
			{
				// read the sequence id (immediately follows the ">")

				std::string curr_id;

				while (!isspace(ch = getc(file)) && ch != EOF)
				{
					curr_id += ch;
				}

				found = (id == curr_id);
			}
			else
			{
				found = (index == curr_index);
			}
		}

		if (reading_seq)
		{
			do
			{
				// avoid allocating very small amounts of extra memory at a time

				if (seq->capacity() == seq->length())
				{
					seq->reserve(seq->capacity() + 1000);
				}

				(*seq) += ch;
				ch = getc(file);

			} while (ch != '\n' && ch != EOF);
		}
		else
		{
			// skip to the end of the line

			while (ch != '\n' && ch != EOF)
			{
				ch = getc(file);
			}

			if (found)
			{
				reading_seq = true;
			}
		}
	}

	if (!found)
	{
		std::cerr << Config::cmd()
			<< ": FASTA file " << file.name()
			<< " does not contain a sequence with ";

		if (!id.empty())
		{
			std::cerr << "id " << id; 
		}
		else
		{
			std::cerr << "index " << index;
		}

		std::cerr << "\n";
		exit(1);
	}
}

Amino Sequence::get_check_amino(char ch, unsigned int position)
{
	if (!isalpha(ch))
	{
		std::cerr << Config::cmd()
			<< ": illegal character '"
			<< ch << "' in amino sequence [pos "
			<< position << "]\n";
		exit(1);
	}

	Amino a(ch);

	if (!a.standard())
	{
		std::cerr << Config::cmd() << ": warning: ";

		if (a.ambiguous())
		{
			std::cerr << "ambiguous amino";
		}
		else
		{
			std::cerr << "rare/unknown amino";
		}

		std::cerr << " '" << (char) toupper(ch) << "'";

		if (a.ambiguous())
		{
			std::cerr << " (" << a.name() << ")";
		}
			
		std::cerr << " [pos " << position << "]\n";
	}

	return a;
}

Codon Sequence::get_check_codon(char *nuc, unsigned int position)
{
	for (int n = 0;n < 3;n++)
	{
		char ch = toupper(nuc[n]);

		if (strchr("ACGTU", ch) == NULL)
		{
			std::cerr << Config::cmd()
				<< ": illegal character '"
				<< nuc[n] << "' in codon sequence [pos "
				<< position << "]\n";
			exit(1);
		}
	}

	Codon c(nuc);
	assert(!c.illegal());

	if (c.is_stop_codon())
	{
		std::cerr << Config::cmd()
			<< ": error: codon sequence includes stop codon \""
			<< c.str() << "\" [codon "
			<< position << "]\n";
		exit(1);
	}

	return c;
}

void Sequence::create_amino_seq(const std::string &seq)
{
	clear_seq();
	m_amino.reserve(seq.length());

	for (unsigned int n = 0;n < seq.length();n++)
	{
		char ch = seq[n];

		// allow spaces inside the sequence
		if (!isspace(ch))
		{
			// append to the end of the vector
			m_amino.push_back(get_check_amino(ch, m_amino.size()));
		}
	}
}

void Sequence::create_codon_seq(const std::string &seq)
{
	clear_seq();
	
	// (avoid calling reserve() with a value of 0)
	m_codon.reserve(max_val(1, (int) seq.length() / 3));

	char nuc[3];	// the last three nucleotides
	int i = 0;		// index in nuc[]

	for (unsigned int n = 0;n < seq.length();n++)
	{
		char ch = seq[n];

		if (!isspace(ch))
		{
			nuc[i++] = ch;

			if (i == 3)
			{
				m_codon.push_back(get_check_codon(nuc, m_codon.size()));
				i = 0;
			}
		}
	}

	if (i != 0)
	{
		std::cerr << Config::cmd()
			<< ": error: codon sequence length is not a multiple of three\n";
		exit(1);
	}

	create_amino_from_codons();
}

void Sequence::create_amino_from_codons()
{
	assert(m_codon.size() != 0);
	assert(m_amino.size() == 0);

	m_amino.reserve(m_codon.size());

	for (unsigned int n = 0;n < m_codon.size();n++)
	{
		m_amino.push_back(m_codon[n].to_amino());
	}
}

void Sequence::clear_seq()
{
	m_amino.clear();
	m_codon.clear();

	// (clear() doesn't work on some compilers)
	assert(m_amino.size() == 0);
	assert(m_codon.size() == 0);
}

void Sequence::print_template(std::ostream &out)
{
	out << "[" << capitalise(c_config_section) << "]\n"
		"\n"
		<< c_param_type << " = " << c_type_amino
			<< "\t\t# \"" << c_type_amino
			<< "\" or \"" << c_type_nuc << "\"\n"
		<< c_param_file << " = ...\t\t# filename (FASTA format)\n"
		<< "#" << c_param_seq << " = ...\t\t# sequence\n"
		<< "#" << c_param_fasta_index << " = " << c_default_fasta_index
			<< "\t\t# index of sequence in FASTA file (starting from 0)\n"
		<< "#" << c_param_fasta_id << " = ...\t\t# id of sequence "
			"in FASTA file\n"
		"\n";
}

