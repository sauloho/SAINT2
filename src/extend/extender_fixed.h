#ifndef EXTENDER_FIXED_H_INCLUDED
#define EXTENDER_FIXED_H_INCLUDED

// Basic type of extender (extrude next residue after a fixed time).
//
// Most of the functionality is handled by the base class, Extender.

#include <iostream>
#include <string>

class Sequence;

class Extender_Fixed : public Extender
{
public:
	// constructor
	Extender_Fixed();

	// destructor
	virtual ~Extender_Fixed();

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

protected:
	virtual void calculate_num_moves(const Sequence &seq);

private:
	// config file parameters
	static const char *c_type;
    static const char *c_param_move_distrib;
    static const char *c_default_move_distrib;

	// types

	enum Move_Distribution { Distrib_Fixed, Distrib_Linear };

	Move_Distribution m_move_distrib;   // "move distribution" parameter
};

#endif // EXTENDER_FIXED_H_INCLUDED

