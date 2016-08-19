
#ifndef STRATEGY_H_INCLUDED
#define STRATEGY_H_INCLUDED

// Abstract base class representing the strategy used for determining
// whether to accept or reject a move (based on the score).
//
// When a new Strategy subclass is created, add it to the list in the
// for_each_Strategy_subclass() macro in strategy.cpp.
//
// In addition to the virtual functions in the base class, the new
// class should contain static functions type() and print_template()

#include "common.h"
#include "param_list.h"

// forward declarations
class Config;
class Runner;
class Run_Observer;
class Sequence;

class Strategy
{
public:
	// default constructor
	Strategy();
	
	// destructor
	virtual ~Strategy();

	// create an object of the appropriate subclass
	static Strategy *create(const Param_List &params);

	/// Parse a config file parameter. Returns false if the parameter
	/// is not recognised.
	virtual bool parse_parameter(const std::string &name,
		const std::string &value) = 0;

	// verify that the parameters are consistent and complete
	// (otherwise exits with an error message)
	virtual void verify_parameters() = 0;

	// get the number of candidate peptides required for this strategy
	// (eg. Strict_Improvement requires only one; Boltzmann requires
	// a set of candidates to choose from)
	virtual int num_candidates() = 0;

	// returns the index of the candidate to select
	// or -1 to select none of them
	virtual int select(
		double old_score,			// in - previous peptide's score
		const Double_Vec &new_score	// in - new peptide scores
									//      (the length of this vector
									//      is num_candidates())
	) = 0;

	// called when a new run is about to start (before the first move)
	virtual void start_run(Runner *runner) = 0;

	// called after run is ended
	virtual void end_run(Runner *runner) = 0;

	// check whether it is time to stop making moves
	virtual bool stop() = 0;

	// for strategies which have a specialised way of performing runs
	// (eg. strategies that maintain an ensemble of structures, such
	// as a genetic algorithm).
	// If this function is overridden, it should always return true.
	virtual bool do_runs(Runner &runner, Sequence &seq, Run_Observer &observer)
	{ return false; }

	// print config file template
	static void print_template(std::ostream &out);

	// get the name of the configuration file section for this class
	static const char *config_section()
	{ return c_config_section; }

private:
	// name of config file section
    static const char *c_config_section;

	// parameter names
	static const char *c_param_type;

	// default values
	static const char *c_default_type;
};

#endif // STRATEGY_H_INCLUDED

