
#include <iostream>
#include <sys/times.h>
#include "config.h"
#include "runner.h"
#include "run_observer.h"
#include "sequence.h"
#include "peptide.h"

using namespace std;

class Time_Recorder : public Run_Observer
{
public:
    /// @brief Constructor.
    Time_Recorder(const Config &config)
		: Run_Observer(config),
		  m_first_extend(true),
		  m_move_count(0)
	{
	}

    /// @brief Destructor.
    virtual ~Time_Recorder()
	{
	}

    /// @brief Called before the first run (the Peptide object exists,
    /// but only its residue sequence has been initialised, not the atom
    /// coordinates).
    virtual void before_start(Runner *r)
	{
	}

    /// @brief Called after the last run has finished.
    ///
    virtual void after_end(Runner *r)
	{
	}

    /// @brief Called at the beginning of every run (the Peptide has its
    /// initial structure); preceeds the first call to after_extend().
    virtual void start_run(Runner *r)
	{
		cout << "Starting run\n";
		m_prev_time = user_time();
		m_move_count = 0;
		m_first_extend = true;
	}

    /// @brief Called at the end of every run (the Peptide has its final
    /// structure).
    virtual void end_run(Runner *r)
	{
		print_time_diff(r->peptide().length(), user_time());
		cout << "Ending run\n";
	}

    /// @brief Called after every move, whether the move succeeded or failed.
    ///
    /// @param r The Runner being observed.
    /// @param candidate List of new structures to choose from.
    /// @param score Score for each candidate.
    /// @param choice Which candidate was chosen (-1 means none).
    /// @param best_so_far Whether candidate chosen has the best score so far
    virtual void after_move(Runner *r, const Conf_Vec &candidate,
        const Double_Vec &score, int choice, bool best_so_far,
        double full_score)
	{
		m_move_count++;
		if (choice != 0) { cout << "Choice " << choice << "\n"; }
	}

    /// @brief Called whenever the peptide is extended (including the initial
    /// extension).
    virtual void after_extend(Runner *r, int num_extruded)
	{
		clock_t t = user_time();

		if (m_first_extend)
		{
			m_first_extend = false;
		}
		else
		{
			print_time_diff(r->peptide().length() - 1, t);
		}

		m_prev_time = t;
		m_move_count = 0;

		if (r->peptide().full_grown()) { cout << "(full grown)\n"; }
	}

    /// @brief A general message to be printed by the Time_Recorder.
    virtual void msg(const std::string &s)
	{
	}

private:
	clock_t user_time()
	{
		struct tms t;
		(void) times(&t);
		return t.tms_utime;
	}

	void print_time_diff(int length, clock_t t)
	{
		cout << "Length " << length
			<< " moves " << m_move_count
			<< " time " << t - m_prev_time << "\n";
	}

private:
	clock_t m_prev_time;
	bool m_first_extend;
	int m_move_count;
};

int main(int argc, const char **argv)
{
	if (argc != 2)
	{
		cerr << "Usage: " << argv[0] << " config_file\n";
		exit(1);
	}

	Config config(2, argv);
	Runner runner(config);
	Sequence seq(config);
	Time_Recorder time_rec(config);

	runner.set_num_runs(1);
	runner.do_runs(seq, time_rec);

	return 0;
}

