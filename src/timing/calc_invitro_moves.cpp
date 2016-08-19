#include <cstdlib> // PG added this
#include <iostream>
#include <string>
#include <cmath>
#include "config.h"
#include "runner.h"
#include "sequence.h"
#include "extender.h"

using namespace std;

// this class can access protected members of class Extender

class Extender_Friend
{
public:
	Extender_Friend(Extender *e) : m_extender(e)
	{
	}

	void do_stuff(Sequence &seq, int non_growth_moves, bool dump)
	{
		m_extender->calculate_num_moves(seq);
		double total_ctime = 0.0;
		int total_moves = 0;

		if (dump)
		{
			cout << "Length cotrans_moves time\n";
		}

		for (unsigned n = 0;n < m_extender->m_moves.size();n++)
		{
			if (m_extender->m_moves[n] != -1)
			{
				total_moves += m_extender->m_moves[n];
				double cotrans_time = pow((double) n, 1.83) / 2890.0 *
					m_extender->m_moves[n];
				total_ctime += cotrans_time;

				if (dump)
				{
					cout << n << " " << m_extender->m_moves[n]
						<< " " << cotrans_time << "\n";
				}
			}
		}

		double invitro_move_time = pow((double) seq.length(), 1.85) / 3500.0;
		double result = total_ctime / invitro_move_time;

		if (dump)
		{
			cout << "Cotrans moves:\n" << total_moves << "\n";
			cout << "Total cotrans time:\n" << total_ctime << "\n";
			cout << "In vitro time for one move (length " << seq.length()
				<< "):\n"
				<< invitro_move_time << "\n";
			cout << "Equivalent in vitro moves:\n"
				<< (int) (result + 0.5) << "\n";
			cout << "Including non-growth moves from config file:\n";
		}

		cout << (int) (result + non_growth_moves + 0.5) << "\n";
	}

private:
	Extender *m_extender;
};

void show_usage(const char **argv)
{
	cerr << "Usage: " << argv[0] << " config_file [-d]\n\n"
		<< "Calculate the equivalent number of in vitro moves for a "
		   "cotranslational run.\n\n"
		   "Use \"-d\" to dump all of the moves.\n\n";
}

int main(int argc, const char **argv)
{
	if (argc != 2 && argc != 3)
	{
		show_usage(argv);
		exit(1);
	}

	bool dump = false;

	if (argc == 3)
	{
		if (argv[2][0] != '-')
		{
			show_usage(argv);
			exit(1);
		}
		else
		if (string(argv[2]) == "-d")
		{
			dump = true;
		}
		else
		{
			cerr << argv[0] << ": unknown option " << argv[2] << "\n";
			exit(1);
		}
	}

	Config config(2, argv);
	Sequence seq(config);

	// TO DO: check if is equential

	Runner runner(config);

	Extender_Friend f(runner.extender());
	f.do_stuff(seq, (int) runner.move_limit(), dump);

	return 0;
}

