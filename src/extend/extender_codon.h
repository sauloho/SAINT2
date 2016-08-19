#ifndef EXTENDER_CODON_H_INCLUDED
#define EXTENDER_CODON_H_INCLUDED

// Type of Extender with different extrusion speeds for different codon
// types.

#include <iostream>
#include <string>
#include "codon.h"

// forward declarations
class Runner;
class Sequence;

class Extender_Codon : public Extender
{
public:
	// constructor
	Extender_Codon();

	// destructor
	virtual ~Extender_Codon();

	// set the name of the codon speed file
	void set_filename(const std::string &filename);

	/// Parse a config file parameter. Returns false if the parameter
	/// is not recognised.
	virtual bool parse_parameter(const std::string &name,
		const std::string &value);

	// verify that the parameters are consistent and complete
	// (otherwise exits with an error message)
	virtual void verify_parameters();

	// print sample config file parameters
	static void print_template(std::ostream &out, bool commented = true);

	// "type" value in config file
	static const char *type()
	{ return c_type; }

private:
	// read the codon speed file
	void read_file();

protected:
	virtual void calculate_num_moves(const Sequence &seq);

private:
	// config file parameters
	static const char *c_type;
	static const char *c_param_file;
	static const char *c_param_tunnel_length;

	// default parameter values

	std::string m_filename;				// codon speed file
	double m_codon_factor[Codon::Num];	// relative number of moves for codons
	int m_tunnel_length;				// number of residues in ribosome tunnel
};

#endif // EXTENDER_CODON_H_INCLUDED

