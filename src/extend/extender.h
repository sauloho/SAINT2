#ifndef EXTENDER_H_INCLUDED
#define EXTENDER_H_INCLUDED

#include <iostream>
#include <vector>
#include <cassert>
#include "param_list.h"

// forward declarations
class Runner;
class Peptide;
class Sequence;

/// Abstract base class for deciding when to extrude the next residue
/// (or group of residues). This class does not actually perform the
/// extrusion, it just decides when to do it.
///
/// When a new Extender subclass is created, add it to the list in
/// the for_each_Extender_subclass() macro in extender.cpp.
///
/// In addition to the virtual functions in the base class, the new
/// class should contain static functions type() and print_template().

class Extender
{
	friend class Extender_Friend;

public:
	/// Constructor.
	Extender();

	/// Destructor.
	virtual ~Extender();

	/// Create an object of the appropriate subclass.
	static Extender *create(const Param_List &params);

	/// Parse a config file parameter. Returns false if the parameter
	/// is not recognised.
	virtual bool parse_parameter(const std::string &name,
		const std::string &value) = 0;

	/// Verify that the parameters are consistent and complete
	/// (otherwise exits with an error message).
	virtual void verify_parameters() = 0;

	/// Get the number of residues to extrude (0 if it is not time
	/// for the next extrusion yet).
	int check_extend(Runner *runner);

	/// Called when a new run is about to start (after the initial extrusion).
	virtual void start_run(const Sequence &seq);

	// called after the peptide has just been extended
	virtual void after_extend(const Peptide &p, Runner *runner);

	// check if it is time to extrude the next residue(s)
	virtual bool must_extend(const Peptide &p, Runner *runner);

	// initial number of residues to extrude
	int initial_residues() const
	{ return m_initial_res; }

	// number of residues extruded at a time
	int residues_to_extrude() const
	{ return m_extrude_res; }

	// set initial number of residues to extrude
	void set_initial_res(int val);

	// set number of residues extruded at a time
	void set_extrude_res(int val);

	// set total number of moves to perform during the growth phase
	void set_growth_moves(long val);

	// get the number of moves to perform for the current length
	int curr_length_move_limit(const Peptide &p);

	// print config file template
	static void print_template(std::ostream &out);

	// get the name of the configuration file section for this class
	static const char *config_section()
	{ return c_config_section; }

	// get the name of the "residues to extrude at a time" parameter
	static const char *extrude_res_param()
	{ return c_param_extrude; }

protected:
	// work out the number of moves to perform after each extrusion
	// (fills in values in m_moves[] vector; undefined values should
	// be set to -1)
	virtual void calculate_num_moves(const Sequence &seq) = 0;

protected:
	std::vector<int> m_moves;	// number of moves at each length
								// (-1 means undefined)

	// name of config file section
	static const char *c_config_section;

	// config file parameters
    static const char *c_param_type;
    static const char *c_param_initial;
    static const char *c_param_extrude;
    static const char *c_param_growth_moves;

	// default parameter values
    static const char *c_default_type;
    static const int c_default_initial;
    static const int c_default_extrude;
    static const long c_default_growth_moves;

	// config file values

    int m_initial_res;		// initial number of residues to extrude
    int m_extrude_res;		// number of residues to extrude at a time
    long m_growth_moves;	// total number of moves during growth phase
};

#endif // EXTENDER_H_INCLUDED

