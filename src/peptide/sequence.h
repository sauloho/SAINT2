
#ifndef SEQUENCE_H_INCLUDED
#define SEQUENCE_H_INCLUDED

// Class representing an amino acid sequence (with an optional
// corresponding codon sequence).
// Sequences are stored in N to C order, so the N terminus has index 0.

#include <string>
#include <cassert>
#include <iostream>
#include "amino.h"
#include "codon.h"
#include "param_list.h"
#include "c_file.h"
#include "config.h"

class Peptide;

class Sequence
{
public:
	// constructor
	Sequence(const Config &config);

	// default constructor (required for use with vector)
	Sequence();

	// destructor
	~Sequence();

	// extract the sequence from a peptide
	void create_from_peptide(const Peptide &p);

	// get the length of the sequence
	int length() const
	{ return (int) m_amino.size(); }

	// check if the codon sequence is known
	bool codons_known() const
	{ return (m_codon.size() != 0); }

	// get the nth amino acid (N terminus is at index 0)
	Amino amino(int n) const
	{
		assert(n >= 0 && n < length());
		return m_amino[n];
	}

	// get the nth codon (N terminus is at index 0)
	Codon codon(int n) const
	{
		assert(codons_known());
		assert(n >= 0 && n < length());
		return m_codon[n];
	}

	// print the amino acid sequence, splitting at aminos_per_line
	// (single line if the value is 0).
	// Includes a terminating newline
	void print(std::ostream &out = std::cout, int aminos_per_line = 0);

	// print the codon sequence, splitting at codons_per_line
	// (single line if the value is 0).
	// If a separator is specified, it is printed between each codon.
	// Includes a terminating newline
	void print_codons(std::ostream &out = std::cout,
		int codons_per_line = 0, char separator = '\0',
		bool sep_at_end_of_line = false);

	// create sequence from configuration parameters
	void create_from_params(const Param_List &params);

	// create from an amino acid sequence
	void create_amino_seq(const std::string &seq);

	// create from a codon sequence
	void create_codon_seq(const std::string &seq);

	// read amino acid sequence (FASTA file format)
	// If fasta_id is not empty, fasta_index is ignored
	void read_amino_seq(
		const std::string &filename,	// FASTA file
		int fasta_index,				// index of sequence (from 0)
		const std::string &fasta_id		// id of sequence
	);

	// (shorthand for the above to read the first sequence in the file)
	void read_amino_seq(const std::string &filename);

	// read codon sequence (FASTA file format)
	// If fasta_id is not empty, fasta_index is ignored
	void read_codon_seq(
		const std::string &filename,	// FASTA file
		int fasta_index,				// index of sequence (from 0)
		const std::string &fasta_id		// id of sequence
	);

	// (shorthand for the above to read the first sequence in the file)
	void read_codon_seq(const std::string &filename);

	// print template parameters (with comments)
	static void print_template(std::ostream &out);

	// get the name of the configuration file section for this class
	static const char *config_section()
	{ return c_config_section; }

private:
    // disable implicit copy and assignment (for safety) by making them private
    Sequence(const Sequence&);
    Sequence& operator = (const Sequence&);

	// helper function for create_from_params()
	void create_sequence(const std::string &type,
		const std::string &filename, const std::string &seq,
		int fasta_index, const std::string &fasta_id);

	// create the amino acid sequence (assumes codon sequence already created)
	void create_amino_from_codons();

	// convert a char to an amino acid, reporting any errors
	// (exits on illegal character)
	Amino get_check_amino(char ch, unsigned int position);

	// convert three characters to a codon, reporting any errors
	// (exits on illegal character)
	Codon get_check_codon(char *nuc, unsigned int position);

	// exit with an error message about the sequence not being defined
	void missing_sequence_err();

	// extract a sequence from a FASTA file (amino acid or codons)
	// Searches for id if it is non-empty, otherwise reaths the nth
	// sequence in the file (ie. value of "index"; the first one is 0)
	void read_fasta_seq(C_File &file, int index, const std::string &id,
		std::string *seq);

	// clear all existing values
	void clear_seq();

private:
	// name of config file section
	static const char *c_config_section;

	// "type" values
	static const char *c_type_amino;
	static const char *c_type_nuc;

	// parameter names
	static const char *c_param_type;
	static const char *c_param_file;
	static const char *c_param_seq;
	static const char *c_param_fasta_index;
	static const char *c_param_fasta_id;

	// default values
	static const char *c_default_type;
	static const int c_default_fasta_index;

	Amino_Vec m_amino;		// amino acids
	Codon_Vec m_codon;		// corresponding codons (if known)
};

#endif // SEQUENCE_H_INCLUDED

