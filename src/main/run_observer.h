#ifndef RUN_OBSERVER_H_INCLUDED
#define RUN_OBSERVER_H_INCLUDED

#include <string>
#include "peptide.h"
#include "conformation.h"
#include "common.h"

// forward declarations
class Runner;
class Config;

/// Class for monitoring a Runner (to print progress, statistics, etc.
/// and do something with the final structures, eg. write to a file).
///
/// Example:
/// <pre>
/// Class Reporter : public Run_Observer
/// {
/// public:
///    ...
///
///    virtual void start_run(Runner *r)
///    {
///        std::cout << "Starting run #" << r->run_number() << "\n";
///    }
/// };
///
/// Config config(argc, argv);
/// Runner runner(config);
/// Sequence seq(config);
///
/// Reporter rep;
/// runner.do_runs(seq, &rep);
/// </pre>

class Run_Observer
{
public:
	/// Constructor.
	Run_Observer(const Config &config)
		: m_config(config)
	{ }

	/// Destructor.
	virtual ~Run_Observer()
	{ }

	/// Called before the first run.
	virtual void before_start(Runner *r)
	{ }

	/// Called after the last run.
	virtual void after_end(Runner *r)
	{ }

	/// Called at the beginning of every run (after the peptide is initialised).
	virtual void start_run(Runner *r)
	{ }

	/// Called at the end of every run.
	virtual void end_run(Runner *r)
	{ }

	/// Called after every move (whether it succeeded or failed).
	/// @param r The Runner object.
	/// @param candidate Set of potential structures.
	/// @param score Score for each \a candidate structure.
	/// @param choice Which candidate was chosen (index in \a candidate;
	/// -1 means none).
	/// @param best_so_far Whether candidate chosen has the best score so far
	virtual void after_move(Runner *r, const Conf_Vec &candidate,
		const Double_Vec &score, int choice, bool best_so_far,
		double full_score)
	{ }

	/// Called after more residues are extruded (including initial residues)
	virtual void after_extend(Runner *r, int num_extruded)
	{ }

	/// General message to be printed.
	virtual void msg(const std::string &s)
	{ }

protected:
	const Config &m_config;
};

#endif // RUN_OBSERVER_H_INCLUDED
