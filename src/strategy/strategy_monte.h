
#ifndef STRATEGY_MONTE_H_INCLUDED
#define STRATEGY_MONTE_H_INCLUDED

#include <string>
#include "common.h"
#include "strategy.h"

// A score improvement strategy that uses the Metropolis Monte Carlo
// criterion: if the new structure improves the score (or leave it
// unchanged), accept it; otherwise the probability of acceptance
// depends on the score change and the temperature.
//
// A temperature of 0.0 means every structure is accepted, and a
// high temperature makes it more likely that moves that make the
// score worse will be rejected.

class Strategy_Monte : public Strategy
{
public:
	// constructor
	Strategy_Monte();

	// destructor
	virtual ~Strategy_Monte();

	// set the temperature (val >= 0.0)
	void set_temperature(double val);

	/// Parse a config file parameter. Returns false if the parameter
	/// is not recognised.
    virtual bool parse_parameter(const std::string &name,
        const std::string &value);

    // verify that the parameters are consistent and complete
    // (otherwise exits with an error message)
    virtual void verify_parameters();

	// get the number of candidate peptides required
	virtual int num_candidates()
	{ return 1; }

	// returns 0 to accept the peptide or -1 to reject it
	virtual int select(
        double old_score,            // in - previous peptide's score
        const Double_Vec &new_score   // in - new peptide score
                                    //      (vector is always length 1)
    );

	// called when a new run is about to start (before the first move)
	virtual void start_run(Runner *runner);

	// called after run is ended
	virtual void end_run(Runner *runner);

	// check whether it is time to stop making moves
	virtual bool stop();

	// set the number of candidates to try at once
	void set_num_candidates(int num);

	// print sample config file parameters
	static void print_template(std::ostream &out, bool commented = true);

	// "type" value in config file
	static const char *type()
	{ return c_type; }

private:
	// config file parameters
    static const char *c_type;
    static const char *c_param_temp;
    static const char *c_param_candidates;
	static const double c_default_temp;
    static const int c_default_candidates;

	double m_temp;				// temperature
	int m_candidates;			// number of candidate structures to try at once
};

#endif // STRATEGY_MONTE_H_INCLUDED

