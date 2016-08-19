
#include <iostream>
#include <cstdio>
#include "peptide.h"
#include "torsion.h"
#include "common.h"
#include "c_file.h"

#define MIN_RELIABLE_VAL 5

int bin_total[TORSION_BINS][TORSION_BINS][Amino::Num];

void show_usage(const char **argv)
{
	std::cerr
		<< "Usage: " << argv[0] << " <PDB dir> <PDB list>\n\n"
		<< "PDB list is like:\n\n"
		<< "1AX3 A\n"
		<< "2BC1 X\n"
		<< "...\n\n";
}

void init_bin_totals()
{
	for (int aa = 0;aa < Amino::Num;aa++)
	{
		for (int a1 = 0;a1 < TORSION_BINS;a1++)
		{
			for (int a2 = 0;a2 < TORSION_BINS;a2++)
			{
				bin_total[a1][a2][aa] = 0;
			}
		}
	}
}

void add_bin_totals(const Peptide &p)
{
	int phi_bin, psi_bin;

	// (first and last residue do not have both phi & psi angles)

	for (int n = 1;n < p.length() - 1;n++)
	{
		if (Torsion::get_phi_psi_bin(p, n, &phi_bin, &psi_bin))
		{
			int a = p.res(n).amino().num();
			++bin_total[phi_bin][psi_bin][a];
		}
	}
}

void write_bins(const std::string &date_str, const std::string &cmd_str)
{
	std::cout <<  "# dgen data file - "
		<< date_str << "\n# "
		<< cmd_str << "\n";

	for (int aa = 0;aa < Amino::Num;aa++)
	{
		std::cout << "\n" << Amino(aa).abbr() << "\n\n";

		for (int a1 = 0;a1 < TORSION_BINS;a1++)
		{
			for (int a2 = 0;a2 < TORSION_BINS;a2++)
			{
				bool ok = (bin_total[a1][a2][aa] > MIN_RELIABLE_VAL);
				char ch = (ok ? 'x' : '.');

				if (!ok)
				{
					// check surrounding bins
					int a1_1 = (a1 + TORSION_BINS - 1) % TORSION_BINS;
					int a1_2 = (a1 + TORSION_BINS + 1) % TORSION_BINS;
					int a2_1 = (a2 + TORSION_BINS - 1) % TORSION_BINS;
					int a2_2 = (a2 + TORSION_BINS + 1) % TORSION_BINS;

					if (bin_total[a1][a2_1][aa] > MIN_RELIABLE_VAL ||
						bin_total[a1][a2_2][aa] > MIN_RELIABLE_VAL ||
						bin_total[a1_1][a2][aa] > MIN_RELIABLE_VAL ||
						bin_total[a1_1][a2_1][aa] > MIN_RELIABLE_VAL ||
						bin_total[a1_1][a2_2][aa] > MIN_RELIABLE_VAL ||
						bin_total[a1_2][a2][aa] > MIN_RELIABLE_VAL ||
						bin_total[a1_2][a2_1][aa] > MIN_RELIABLE_VAL ||
						bin_total[a1_2][a2_2][aa] > MIN_RELIABLE_VAL)
					{
						ch = '+';
					}
				}

				std::cout << ch;
			}

			std::cout << "\n";
		}
	}
}

// concatenate all of the arguments on the command line
// (with spaces in between)

std::string get_command(int argc, const char *argv[])
{
	char cmd[1000];
	cmd[0] = '\0';

	for (int a = 0;a < argc;a++)
	{
		strcat(cmd, argv[a]);

		if (a < argc - 1)
		{
			int len = strlen(cmd);
			cmd[len] = ' ';
			cmd[len+1] = '\0';
		}
	}

	return std::string(cmd);
}

int main(int argc, const char **argv)
{
	if (argc != 3)
	{
		show_usage(argv);
		exit(0);
	}

	const char *pdb_dir = argv[1];
	const char *pdb_list = argv[2];

	C_File file(pdb_list, "r", "PDB list");
	init_bin_totals();

	const int Max_Len = 1000;
	char buffer[Max_Len];

	while (file.next_line(buffer, Max_Len))
	{
		char pdb_id[Max_Len];
		char chain;

		if (sscanf(buffer, "%s %c", pdb_id, &chain) != 2)
		{
			std::cerr << "Error on line " << file.line_num()
				<< " of " << pdb_list << "\n";
			exit(1);
		}

		char pdb_filename[Max_Len];
		sprintf(pdb_filename, "%s/%s.pdb", pdb_dir, pdb_id);

		Peptide p;
		p.read_pdb(pdb_filename, chain,
			true   // no warnings
			);
		p.conf().calc_torsion_angles();
		std::cerr << "Read " << pdb_filename << "\n";
		add_bin_totals(p);
	}

	write_bins(get_date(), get_command(argc, argv));
	return 0;
}

