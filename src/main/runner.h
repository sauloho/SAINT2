#ifndef RUNNER_H_INCLUDED
#define RUNNER_H_INCLUDED

#include <string>
#include "param_list.h"
#include "peptide.h"

// forward declarations
class Config;
class Sequence;
class Scorer;
class Strategy;
class Extender;
class Mover;
class Run_Observer;

/// @brief Class that performs a series of "runs" on a sequence, producing
/// a Peptide structure at the end of each run.

class Runner
{
public:
	/// @brief Constructor.
	Runner(Config &config);

	/// @brief Destructor.
	~Runner();

	/// @brief Parse configuration file parameters related to the Runner
	/// (eg. the number of runs to perform).
	void parse_params(const Param_List &params);

	/// @brief Perform the runs. For each run, a Peptide is generated for
	/// the specified sequence, and the \a observer is told about it.
	void do_runs(Sequence &seq, Run_Observer &observer);

	/// @brief Set the number of runs to perform.
	void set_num_runs(int num_runs);
	
	/// @brief Set the output file.
	void set_output(std::string outfile);

	/// @brief The number of runs to perform.
	int num_runs() const
	{ return m_num_runs; }

	/// @brief The output file.
	std::string output() 
	{ return m_outfile; }

	/// @brief The current run number (starting from 0)
	int run_number() const
	{ return m_run; }

	/// @brief Get the current structure.
	Peptide &peptide()
	{ return m_peptide; }

	/// @brief Get the current structure (const version).
	const Peptide &peptide() const
	{ return m_peptide; }

	/// @brief Get the current structure's score.
	double score() const
	{ return m_curr_score; }

	/// @brief Get the score before the last move was made.
	double prev_score() const
	{ return m_prev_score; }

	/// @brief Check if the last attempted move failed.
	bool last_move_failed() const
	{ return m_move_failed; }

	/// @brief Get the number of moves at the current peptide length
	/// (for the current run).
	int curr_length_moves() const
	{ return m_curr_length_moves; }

	/// @brief Get the number of times in a row no new structure was selected.
	int no_selection_count() const
	{ return m_no_sel_count; }

	/// @brief Set the Scorer to use (deleted by the Runner on destruction).
	void set_scorer(Scorer *s);

	/// @brief Set the Strategy to use (deleted by the Runner on destruction).
	void set_strategy(Strategy *s);

	/// @brief Set the Extender to use (deleted by the Runner on destruction).
	void set_extender(Extender *e);

	/// @brief Set the Mover to use (deleted by the Runner on destruction).
	void set_mover(Mover *m);

	/// @brief The Scorer currently being used.
	Scorer *scorer()
	{ return m_scorer; }

	/// @brief The Strategy currently being used.
	Strategy *strategy()
	{ return m_strategy; }

	/// @brief The Extender currently being used.
	Extender *extender()
	{ return m_extender; }

	/// @brief The Mover currently being used.
	Mover *mover()
	{ return m_mover; }

	bool native_known()
	{ return m_native_known; }

	const Peptide &native()
	{ return m_native_peptide; }

	long move_limit() const
	{ return m_move_limit; }

	/// @brief Print a template for the Runner section of the config file.
	static void print_template(std::ostream &out);

	/// @brief Get config file section name.
	static const char *config_section();

private:
	// disable copy and assignment by making them private
	Runner(const Runner&);
	Runner &operator = (const Runner&);

private:
	/// name of config file section corresponding to the Runner class.
	static const char *m_config_section;
	std::string m_outfile;

	// config file parameters

    static const char *c_param_sequential;
    static const char *c_param_print_intermediates;
    //static const char *c_param_coil;
	static const char *c_param_codon_speed_file;
    static const char *c_param_initial_res;
	static const char *c_param_move_limit;
	static const char *c_param_no_sel_limit;
	//static const char *c_param_scwrl_executable;
	static const char *c_param_reverse;	// for global variable reverseSaint
	static const char *c_param_start_struct;
	static const char *c_param_native_struct;

	// default parameter values

    static const bool c_default_sequential;
    static const bool c_default_print_intermediates;
    //static const bool c_default_coil;
    static const int c_default_initial_res;
	static const long c_default_move_limit;
	static const long c_default_no_sel_limit;

	/// number of runs to perform
	int m_num_runs;

    /// boolean flag to indicate when to stop because of contacts

	/// scoring object
	Scorer *m_scorer;

	/// score strategy object
	Strategy *m_strategy;
	
	/// extrusion object
	Extender *m_extender;

	/// random movement object
	Mover *m_mover;

	/// current run (starting from 0)
	int m_run;

	/// current structure
	Peptide m_peptide;

	/// score for current structure
	double m_curr_score;

	/// previous structure's score
	double m_prev_score;

	/// @brief whether the last move failed
	///
	/// (whether the last call to m_strategy->select() returned -1)
	bool m_move_failed;

	/// number of moves at the current peptide length
	long m_curr_length_moves;

	/// @brief number of times in a row nothing has been selected
	///
	/// (reset to 0 after peptide is extended)
	long m_no_sel_count;

	// configuration file parameters

	/// sequential or non-sequential
	bool m_sequential;

    /// print intermediates??
    bool m_print_intermediates;

	/// whether non-sequential starts as a random coil or fully extended
	//bool m_coil;

	/// maximum number of moves in a run (excluding moves during growth)
	long m_move_limit;

	/// stop after nothing is selected this many times in a row (if full grown)
	long m_no_sel_limit;

	/// starting structure filename
	std::string m_start_struct;

	/// native structure filename
	std::string m_native_struct;

	//// native structure (if known)
	Peptide m_native_peptide;

	bool m_native_known;
};

#endif // RUNNER_H_INCLUDED
