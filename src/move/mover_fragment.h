
#ifndef MOVER_FRAGMENT_H_INCLUDED
#define MOVER_FRAGMENT_H_INCLUDED

// A Mover that is based on fragment replacement.

#include <string>
#include "mover.h"
#include "fragment.h"

// forward declarations
class Peptide;

class Mover_Fragment : public Mover
{
public:
	// constructor
	Mover_Fragment();

	// destructor
	virtual ~Mover_Fragment();

	// set the fragment library location
	void set_library(const std::string &lib);

	/// Parse a config file parameter. Returns false if the parameter
	/// is not recognised.
    virtual bool parse_parameter(const std::string &name,
        const std::string &value);

    // verify that the parameters are consistent and complete
    // (otherwise exits with an error message)
    virtual void verify_parameters();

	// set a peptide to its initial state (non-sequential)
	virtual void init_non_sequential(Peptide &s, bool random_coil,
		Run_Observer *observer);

	// initialise the Mover from an existing peptide
	virtual void init_from_peptide(Peptide &p, Run_Observer *observer);

	// print sample config file parameters
	static void print_template(std::ostream &out, bool commented /*= true*/);

    // "type" value in config file
    static const char *type()
    { return c_type; }

	// dump internal state (debugging function)
	virtual void dump() = 0;

private:
	// disable copy and assignment by making them private
	Mover_Fragment(const Mover_Fragment&);
	Mover_Fragment &operator = (const Mover_Fragment&);

protected:
	// read the fragment library (if it has not been read already)
	void load_fragments(Run_Observer *observer);

	// add a new fragment
	virtual Fragment *add_fragment(int start_pos, int length) = 0;

	// called at end of load_fragments()
	virtual void after_fragments_loaded(int c_terminus) = 0;

	// transform the chain so that the centre of gravity of the CA atoms is
	// on the -x axis (assumes the most recently extruded residue is at
	// (0, 0, 0))
	void reorient_for_ribosome(Peptide &p);

protected:
	typedef std::vector<Fragment_Vec> Fragment_Vec_Vec;
	typedef std::vector<Fragment*> Fragment_Ptr_Vec;

	// "type" parameter name
	static const char *c_type;
	static const char *c_param_lib;
	static const char *c_param_double_replacement_prob;

	std::string m_lib;			// fragment library
	double m_double_replacement_prob;	// probability of doing two in a row
	bool m_fragments_loaded;	// whether fragment library has been read

private:
	int m_build_from_pos;	// first end position we're using fragments for
protected:
	// pointer to first end position we're using fragments for
	int* p_m_build_from_pos = &m_build_from_pos;
};

#endif // MOVER_FRAGMENT_H_INCLUDED

