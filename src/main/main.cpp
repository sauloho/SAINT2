#include <cstring> // PG added this
#include <cstdlib> // PG added this
#include <iostream>
#include <csignal>
#include "config.h"
#include "sequence.h"
#include "runner.h"
#include "scorer.h"
#include "lennard_jones.h"
#include "reporter.h"
#include "stream_printf.h"
#include "geom.h"
#include "static_init.h"

/// @file main() and related functions.

/// @brief Read a PDB file and print its residue sequence and score.
/// @param filename PDB file.
/// @param chain Chain id (' ' means the first chain in the file).
/// @param scorer Used to calculate score.

void print_pdb_info(const std::string &filename, char chain,
	Scorer *scorer, bool show_torsion, bool backbone_only)
{
	Peptide p;

	if (!p.read_pdb(filename.c_str(), chain))
	{
		std::cerr << "Errors found in PDB file\n";
		exit(1);
	}

	p.conf().calc_torsion_angles();

	std::cout
		<< "File " << filename
		<< "\nChain '" << p.chain()
		<< "'\nLength " << p.length()
		<< "\n";

	/*
	if (p.many_non_backbone())
	{
		std::cout << "Previous contacts: "
			<< p.previous_contact_proportion()
			<< "\nInterleave proportion: "
			<< p.interleaving_proportion()
			<< "\n";
	}
	*/

	std::cout << '\n';

	if (backbone_only)
	{
		std::cout << "Removing non-backbone atoms\n\n";
		p.remove_non_backbone_atoms();
	}

	Sequence seq;
	seq.create_from_peptide(p);
	seq.print(std::cout, 60);
	std::cout << '\n';

	if (show_torsion)
	{
		for (int n = 0;n < p.length();n++)
		{
			double phi = p.conf().phi(n);
			double psi = p.conf().psi(n);
			double omega = p.conf().omega(n);

			if (phi != TORSION_UNKNOWN)
			{
				phi = range_m180_180(rad2deg(phi));
			}

			if (psi != TORSION_UNKNOWN)
			{
				psi = range_m180_180(rad2deg(psi));
			}

			if (omega != TORSION_UNKNOWN)
			{
				omega = range_m180_180(rad2deg(omega));
			}

			double n_angle = TORSION_UNKNOWN;
			double ca_angle = TORSION_UNKNOWN;
			double c_angle = TORSION_UNKNOWN;

			if (n > 0 &&
				p.atom_exists(n - 1, Atom_C) &&
				p.atom_exists(n, Atom_N) &&
				p.atom_exists(n, Atom_CA))
			{
				n_angle = rad2deg(angle_formed(p.atom_pos(n - 1, Atom_C),
					p.atom_pos(n, Atom_N), p.atom_pos(n, Atom_CA)));
			}

			if (p.atom_exists(n, Atom_N) &&
				p.atom_exists(n, Atom_CA) &&
				p.atom_exists(n, Atom_C))
			{
				ca_angle = rad2deg(angle_formed(p.atom_pos(n, Atom_N),
					p.atom_pos(n, Atom_CA), p.atom_pos(n, Atom_C)));
			}

			if (n < p.length() - 1 &&
				p.atom_exists(n, Atom_CA) &&
				p.atom_exists(n, Atom_C) &&
				p.atom_exists(n + 1, Atom_N))
			{
				c_angle = rad2deg(angle_formed(p.atom_pos(n, Atom_CA),
					p.atom_pos(n, Atom_C), p.atom_pos(n + 1, Atom_N)));
			}

			std::cout << Printf("%-4s ", p.res(n).res_seq_str().c_str())
				<< p.res(n).amino().abbr()
				<< Printf("  %6.1f ", phi)
				<< Printf(" %6.1f ", psi)
				<< Printf(" %6.1f ", omega)
				<< Printf(" %6.1f ", n_angle)
				<< Printf(" %6.1f ", ca_angle)
				<< Printf(" %6.1f", c_angle)
				<< "\n";
		}
	}
	else
	{
		scorer->print_info_when_scoring(true);
		double sc = scorer->score(p);
		scorer->print_desc();
		std::cout << " score = " << sc << "\n";

		//p.dump_backbone();
	}

	std::cout
		<< "\nRadius of gyration = " << p.radius_of_gyr()
		<< "\nDiameter = " << p.diameter() << "\n";

	Lennard_Jones lj_scorer;

	if (lj_scorer.steric_clash(p))
	{
		std::cout << "Steric clashes found\n";
	}

	// std::cout << "\nWriting to \"out\"\n";
	// p.write_pdb("out");
}

void sig_handler(int s)
{
	std::cerr << "!! signal " << s
		<< ": " << strsignal(s) << std::endl;
	signal(s, SIG_DFL);
	raise(s);
}

void setup()
{
	signal(1, sig_handler);
	signal(2, sig_handler);
	signal(3, sig_handler);
	signal(4, sig_handler);
	signal(5, sig_handler);
	signal(6, sig_handler);
	signal(7, sig_handler);
	signal(8, sig_handler);
	signal(9, sig_handler);
	signal(10, sig_handler);
	signal(11, sig_handler);
	signal(12, sig_handler);
	signal(13, sig_handler);
	signal(14, sig_handler);
	signal(15, sig_handler);
	signal(16, sig_handler);
	signal(17, sig_handler);
	signal(18, sig_handler);
	signal(19, sig_handler);
	signal(20, sig_handler);
	signal(21, sig_handler);
	signal(22, sig_handler);
	signal(23, sig_handler);
	signal(24, sig_handler);
	signal(25, sig_handler);
	signal(26, sig_handler);
	signal(27, sig_handler);
	signal(28, sig_handler);
	signal(29, sig_handler);
	signal(30, sig_handler);
	signal(31, sig_handler);
	signal(32, sig_handler);
}

int main(int argc, const char *argv[])
{
	setup();

	// parse command line and configuration file
	Config config(argc, argv);
	Runner runner(config);
	if (!config.pdb_filename().empty())
	{
		print_pdb_info(config.pdb_filename(), config.pdb_chain(),
			runner.scorer(), config.show_torsion_angles(),
			config.backbone_only());
		return 0;
	}
	Sequence seq(config);
	Reporter reporter(config);
	runner.do_runs(seq, reporter);

	// runner.peptide().call_scwrl();
	// runner.peptide().write_pdb(config.output_file().c_str());

	return 0;	// success
}

