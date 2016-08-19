
#ifndef SCORER_H_INCLUDED
#define SCORER_H_INCLUDED

// Abstract base class to perform scoring of protein structures.
//
// Note that a low score is better than a high score.
//
// When a new Scorer subclass is created, add it to the list in the
// for_each_Scorer_subclass() macro in scorer.cpp.
//
// In addition to the virtual functions in the base class, the new
// class should contain static functions type() and print_template()

#include <iostream>
#include "param_list.h"

// forward declarations
class Peptide;

class Scorer
{
public:
	// constructor
	Scorer();

	// destructor (must be virtual because this is a base class)
	virtual ~Scorer();

	// create the type of scorer specified by the parameters
	static Scorer *create(const Param_List &params);

	/// Parse a config file parameter. Returns false if the parameter
	/// is not recognised.
	virtual bool parse_parameter(const std::string &name,
		const std::string &value) = 0;

	// verify that the parameters are consistent and complete
	// (otherwise exits with an error message)
	virtual void verify_parameters() = 0;

	// score a peptide (low scores ate better)
	virtual double score(const Peptide &p, double progress = 1.0,
		double *progress1_score = NULL) = 0;

	// print a brief description of the type of scoring
	virtual void print_desc(std::ostream &out = std::cout) = 0;

	// print template config file section
	static void print_template(std::ostream &out);

	// get the name of the configuration file section for this class
	static const char *config_section()
	{ return c_config_section; }

	// whether to print individual scoring terms when score() is called
	void print_info_when_scoring(bool info_on);

	// set verbose flag (print full details of scoring)
	void set_verbose(bool val);

protected:
	bool print_info_when_scoring() const
	{ return m_score_info_on; }

	bool verbose() const
	{ return m_verbose; }

private:
	// name of config file section
	static const char *c_config_section;

    // "type parameter name
    static const char *c_param_type;

	// default values
	static const char *c_default_type;

	// set by print_info_when_scoring()
	bool m_score_info_on;

	// set by set_verbose()
	bool m_verbose;
};

#endif // SCORER_H_INCLUDED

