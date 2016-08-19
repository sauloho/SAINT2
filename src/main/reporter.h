#ifndef REPORTER_H_INCLUDED
#define REPORTER_H_INCLUDED

#include <string>
#include "conformation.h"
#include "run_observer.h"

class Config;

/// @brief A Run_Observer printing the progress and results of each run
/// performed by a Runner.

class Reporter : public Run_Observer
{
public:
	/// @brief Constructor.
	Reporter(const Config &config);

	/// @brief Destructor.
	virtual ~Reporter();

	/// @brief Called before the first run (the Peptide object exists,
	/// but only its residue sequence has been initialised, not the atom
	/// coordinates).
	virtual void before_start(Runner *r);

	/// @brief Called after the last run has finished.
	///
	virtual void after_end(Runner *r);

	/// @brief Called at the beginning of every run (the Peptide has its
	/// initial structure); preceeds the first call to after_extend().
	virtual void start_run(Runner *r);

	/// @brief Called at the end of every run (the Peptide has its final
	/// structure).
	virtual void end_run(Runner *r);

	/// @brief Called after every move, whether the move succeeded or failed.
	///
	/// @param r The Runner being observed.
	/// @param candidate List of new structures to choose from.
	/// @param score Score for each candidate.
	/// @param choice Which candidate was chosen (-1 means none).
	/// @param best_so_far Whether candidate chosen has the best score so far
	virtual void after_move(Runner *r, const Conf_Vec &candidate,
		const Double_Vec &score, int choice, bool best_so_far,
		double full_score);

	/// @brief Called whenever the peptide is extended (including the initial
	/// extension).
	virtual void after_extend(Runner *r, int num_extruded);

	/// @brief A general message to be printed by the Reporter.
	virtual void msg(const std::string &s);

private:
	void compare_to_native(Runner *r);
	void compare_to_native_contacts(Runner *r);
	void compare_to_native_angles(Runner *r);

	// print a matlab statement to set pixel(x, y) to a colour
	// based on a torsion angle difference (in radians)
	void matlab_pixel(int x, int y, double angle_diff);
};

#endif // REPORTER_H_INCLUDED
