
#ifndef STRATEGY_BOLTZ_H_INCLUDED
#define STRATEGY_BOLTZ_H_INCLUDED

// A score improvement strategy that involves generating several
// different structures, then selecting one of them based on the
// Boltzmann distribution (higher probabilities for structures
// that improve the score more).
//
// The "temperature" value affects the probability distribution
// (a temperature of 0.0 makes the distribution completely even,
// and a high temperature makes it more like a "greedy" algorithm
// that always chooses the structure with the maximum score
// improvement).

#include <string>
#include "common.h"
#include "strategy.h"

class Strategy_Boltz : public Strategy
{
public:
	// default constructor
	Strategy_Boltz();

	// destructor
	virtual ~Strategy_Boltz();

	// set the temperature (val >= 0.0)
	void set_temperature(double val);

	// set the number of candidates to select from
	void set_number(int num);

	/// Parse a config file parameter. Returns false if the parameter
	/// is not recognised.
	virtual bool parse_parameter(const std::string &name,
		const std::string &value);

	// verify that the parameters are consistent and complete
	// (otherwise exits with an error message)
	virtual void verify_parameters();

	// get the number of candidate peptides required
	virtual int num_candidates()
	{ return m_num; }

	// returns the index of the peptide to select,
	// or -1 to select none of them
	virtual int select(
        double old_score,            // in - previous peptide's score
        const Double_Vec &new_score   // in - new peptide scores
    );

	// called when a new run is about to start (before the first move)
	virtual void start_run(Runner *runner);

	// called after run is ended
	virtual void end_run(Runner *runner);

	// check whether it is time to stop making moves
	virtual bool stop();

	// print sample config file parameters
    static void print_template(std::ostream &out, bool commented = true);

    // "type" value in config file
    static const char *type()
    { return c_type; }


private:
	// config file parameters
    static const char *c_type;
    static const char *c_param_temp;
    static const char *c_param_num;
    static const double c_default_temp;
    static const int c_default_num;

	double m_temp;		// temperature
	int m_num;			// number of peptides to choose from

	Double_Vec m_prob;	// probability of selecting each peptide
						// (is a member variable to save re-creating it
						// each time it is needed)
};

#endif // STRATEGY_BOLTZ_H_INCLUDED

