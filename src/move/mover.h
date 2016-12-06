
#ifndef MOVER_H_INCLUDED
#define MOVER_H_INCLUDED

// Abstract base class for classes that handle things to do with
// movement (random moves, initial structure, extrusion, etc.)
//
// When a new Mover subclass is created, add it to the list in the
// for_each_Mover_subclass() macro in mover.cpp.
//
// In addition to the virtual functions in the base class, the new
// class should contain static functions type() and print_template()

#include <string>
#include <iostream>
#include "peptide.h"
#include "param_list.h"
#include "conformation.h"

// forward declarations
class Peptide;
class Runner;
class Run_Observer;

class Mover
{
public:
	// constructor
	Mover();

	// destructor
	virtual ~Mover();

	// create an object of the appropriate subclass
	static Mover *create(const Param_List &params);

	/// Parse a config file parameter. Returns false if the parameter
	/// is not recognised.
    virtual bool parse_parameter(const std::string &name,
        const std::string &value) = 0;

    // verify that the parameters are consistent and complete
    // (otherwise exits with an error message)
    virtual void verify_parameters() = 0;

	// set a peptide to its initial state (eg. single residue extruded)
	virtual void init_sequential(Peptide &s, int initial_length,
		Run_Observer *observer) = 0;

	// set a peptide to its initial state if building from segment
	virtual void init_sequential_from_segment(Peptide &s, int initial_length,
		Run_Observer *observer) = 0;

	// set a peptide to its initial state (eg. fully extended)
	virtual void init_non_sequential(Peptide &s, bool random_coil,
		Run_Observer *observer) = 0;

	// initialise the Mover from an existing peptide
	virtual void init_from_peptide(Peptide &p, Run_Observer *observer) = 0;

	// create a set of structures from a peptide; each one is a random
	// move away from the original structure
	virtual void do_random_move(Peptide &p, int num, bool exhaustive_for_pos,
		Conf_Vec &result, Run_Observer *observer) = 0;

	// extend the peptide by the requested number of residues
	virtual void extend(Peptide &s, int num_res,
		bool ribosome_wall, Run_Observer *observer) = 0;

	// print template config file section
	static void print_template(std::ostream &out);

	// get the name of the configuration file section for this class
	static const char *config_section()
	{ return c_config_section; }

private:
	// name of config file section
	static const char *c_config_section;

	// "type" parameter name
	static const char *c_param_type;

	// default values
	static const char *c_default_type;
};

#endif // MOVER_H_INCLUDED

