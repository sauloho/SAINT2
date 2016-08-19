
#include <iostream>
#include <cmath>
#include <vector>
#include "geom.h"
#include "runner.h"
#include "config.h"
#include "reporter.h"

// defined PRINT_ALL to print lots of debug output
//#define PRINT_ALL

Reporter::Reporter(const Config &config)
	: Run_Observer(config)
{
}

Reporter::~Reporter()
{
}

void Reporter::before_start(Runner *r)
{
#ifdef PRINT_ALL
	std::cout << "* Start of runs\n";
#endif // PRINT_ALL
}

void Reporter::after_end(Runner *r)
{
#ifdef PRINT_ALL
	std::cout << "* End of runs\n";
#endif // PRINT_ALL
}

void Reporter::start_run(Runner *r)
{
#ifdef PRINT_ALL
	std::cout << "* Starting run #"
		<< r->run_number()
		<< std::endl;
#endif // PRINT_ALL
}

void Reporter::end_run(Runner *r)
{
#ifdef PRINT_ALL
	std::cout << "* Ending run #"
		<< r->run_number()
		<< "; final score = " << r->score()
		<< std::endl;
#endif // PRINT_ALL

	//r->peptide().call_scwrl();
	r->peptide().write_pdb(m_config.outfilename(r->run_number()).c_str());
}

void Reporter::after_move(Runner *r, const Conf_Vec &candidate,
	const Double_Vec &score, int choice, bool best_so_far, double full_score)
{
#ifdef PRINT_ALL
//	std::cout << "RG " << r->curr_length_moves() << ' '
//		<< r->peptide().radius_of_gyr() << "\n";

// return;	// !!!
	std::cout << "* Move " << r->curr_length_moves()
		<< " L " << r->peptide().length()
		<< ": " << r->prev_score()
		<< ", [";

	for (unsigned int n = 0;n < score.size();n++)
	{
		std::cout << ' ' << score[n] << " (";
		double diff = score[n] - r->prev_score();
		std::cout << (diff >= 0.0 ? "+" : "") << diff << ")";
	}

	std::cout << " ] : ";

	if (choice == -1)
	{
		std::cout << "none (count = "
			<< r->no_selection_count()
			<< ")";
	}
	else
	{
		std::cout << score[choice] << " (";
		double diff = score[choice] - r->prev_score();
		std::cout << (diff >= 0.0 ? "+" : "") << diff << ")";
	}

	if (best_so_far)
	{
		std::cout << " = best {" << full_score << "}";
	}

	std::cout << std::endl;
	compare_to_native(r);
#endif // PRINT_ALL
}

void Reporter::compare_to_native(Runner *r)
{
	if (!r->native_known())
	{
		return;
	}

	// compare_to_native_angles(r);
	// compare_to_native_contacts(r);

	const Peptide &p = r->peptide();
	const Peptide &nat = r->native();

	std::cout << "DD " << r->curr_length_moves()
		<< ' ' << nat.calc_rmsd(p) << '\n';
}

void Reporter::compare_to_native_contacts(Runner *r)
{
	static const double ContactDist = 10.0;

	const Peptide &p = r->peptide();
	const Peptide &nat = r->native();

	int move = r->curr_length_moves();

	std::vector<int> missing, extra;
	missing.resize(p.full_length());
	extra.resize(p.full_length());

	int i;
	for (i = p.start();i <= p.end();i++)
	{
		missing[i] = extra[i] = 0;
	}

	for (int n = p.start() + 2;n <= p.end();n++)
	{
		for (int m = p.start();m < n - 1;m++)
		{
			double d = p.atom_dist(n, Atom_CA, m, Atom_CA);
			double d_nat = nat.atom_dist(n, Atom_CA, m, Atom_CA);

			if (d_nat < ContactDist && d > ContactDist)
			{
				missing[n]++;
				missing[m]++;
			}
			else
			if (d_nat > ContactDist && d < ContactDist)
			{
				extra[n]++;
				extra[m]++;
			}
		}
	}

	std::cout << "CC " << move;

	for (i = p.start();i <= p.end();i++)
	{
		std::cout << ' ' << missing[i] << ' ' << extra[i];
	}

	std::cout << '\n';
}

void Reporter::compare_to_native_angles(Runner *r)
{
	const Peptide &p = r->peptide();
	const Peptide &nat = r->native();

	int move = r->curr_length_moves();

	for (int n = p.start();n <= p.end();n++)
	{
		if (n != p.start())
		{
			double phi = p.conf().phi(n);
			double nat_phi = nat.conf().phi(n);
			double phi_diff = range_0_2pi(nat_phi - phi);
//std::cout << "PHI DIFF: " << nat_phi << " - " << phi
//	<< " = " << nat_phi - phi << " = " << phi_diff << "\n";
			matlab_pixel(move, n * 2, phi_diff);
		}

		if (n != p.end())
		{
			double psi = p.conf().psi(n);
			double nat_psi = nat.conf().psi(n);
			double psi_diff = range_0_2pi(nat_psi - psi);
			matlab_pixel(move, n * 2 + 1, psi_diff);
		}
	}
}

void Reporter::matlab_pixel(int x, int y, double angle_diff)
{
	static const double ThreePiOnTwo =  M_PI + M_PI_2;
	double r, g, b;

	if (angle_diff < M_PI_2)
	{
		// white to red
		r = 1.0;
		g = b = 1.0 - (angle_diff / M_PI_2);
	}
	else
	if (angle_diff < M_PI)
	{
		// red to black
		r = (M_PI - angle_diff) / M_PI_2;
		g = b = 0.0;
	}
	else
	if (angle_diff < ThreePiOnTwo)
	{
		// black to blue
		b = (angle_diff - M_PI) / M_PI_2;
		g = r = 0.0;
	}
	else
	{
		// blue to white
		b = 1.0;
		g = r = (angle_diff - ThreePiOnTwo) / M_PI_2;
	}

	r = floor(r * 100.0 + 0.5) * 0.01;
	g = floor(g * 100.0 + 0.5) * 0.01;
	b = floor(b * 100.0 + 0.5) * 0.01;

//	std::cout << "% MM " << angle_diff
//		<< " (" << range_m180_180(rad2deg(angle_diff)) << ")\n";

	std::cout << "MM(" << y << "," << x << ",:)=["
		<< r << " " << g << " " << b << "];\n";
}

void Reporter::after_extend(Runner *r, int num_extruded)
{
#ifdef PRINT_ALL
	std::cout << "* Extended by " << num_extruded
		<< " to length " << r->peptide().length();

	if (r->peptide().full_grown())
	{
		std::cout << " -- full grown";
	}

	std::cout << "; score = " << r->score()
		<< "\n";

	/*
	** check if Mover_Fragment::reorient_for_ribosome() worked
	*

	Peptide &p = r->peptide();
	Point a(0.0, 0.0, 0.0);

	for (int n = p.start();n <= p.end();n++)
	{
		a.add(p.atom_pos(n, Atom_CA));
	}

	a.scale(1.0 / (double) p.length());

	std::cout << "@@@ First CA = "
		<< r->peptide().atom_pos(0, Atom_CA)
		<< ", mid CA = " << a << "\n";
	*/
#endif // PRINT_ALL
}

void Reporter::msg(const std::string &s)
{
#ifdef PRINT_ALL
	std::cout << "* Msg: " << s << "\n";
#endif // PRINT_ALL
}

