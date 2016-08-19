
#ifndef STRATEGY_ALWAYS_H_INCLUDED
#define STRATEGY_ALWAYS_H_INCLUDED

#include "strategy.h"

// Always accept new structures.
// This class is for testing purposes only.

class Strategy_Always : public Strategy
{
public:
	// constructor
	Strategy_Always();

	// destructor
	virtual ~Strategy_Always();

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
	// (for this class, always returns 0)
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

	// set the number of candidate peptides to try at once
	void set_num_candidates(int num);

	// "type" value in config file
    static const char *type()
    { return c_type; }

	// print sample config file parameters
    static void print_template(std::ostream &out, bool commented = true);

private:
    // config file parameters
    static const char *c_type;
};

#endif // STRATEGY_ALWAYS_H_INCLUDED

