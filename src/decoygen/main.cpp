
#include <iostream>
#include <ctime>
#include <cstring>
#include <cstdio>
#include <cctype>
#include <unistd.h>
#include <string>
#include "torsion.h"
#include "amino.h"
#include "peptide.h"
#include "stream_printf.h"
#include "c_file.h"
#include "common.h"
#include "geom.h"
#include "lennard_jones.h"
#include "transform.h"

//#define NUM_DECOYS  20
#define NUM_DECOYS  80
#define MIN_GDT    0.4

const bool show_dots = false;

// Define DO_LOG to create a log file for each chain
//#define DO_LOG

bool valid_bin[TORSION_BINS][TORSION_BINS][Amino::Num];
Lennard_Jones *lj_scorer = NULL;

void call_tmscore(const char *native_filename, const char *decoy_filename,
	double *gdt, double *rmsd)
{
	const char *exec = "./TMscore";
	char command[1000];

	sprintf(command, "%s %s %s", exec, decoy_filename, native_filename);
	FILE *pipe = popen(command, "r");

	if (pipe == NULL)
	{
		std::cerr << "Error: command failed: \"" << command << "\"\n";
		exit(1);
	}

	const int Max_Len = 2000;
	char buffer[Max_Len];

	bool gdt_found = false;
	bool rmsd_found = false;

	while (fgets(buffer, Max_Len, pipe) != NULL && !(gdt_found && rmsd_found))
	{
		const char *gdt_str = "GDT-TS-score=";
		const char *rmsd_str = "RMSD of  the common residues=";
		static const int gdt_str_len = strlen(gdt_str);
		static const int rmsd_str_len = strlen(rmsd_str);

		if (std::string(buffer, gdt_str_len) == gdt_str)
		{
			if (sscanf(buffer + gdt_str_len, "%lf", gdt) != 1)
			{
				std::cerr << "Error: no number after \""
					<< gdt_str << "\" in output of command: \""
					<< command << "\"\n";
				pclose(pipe);
				exit(1);
			}

			gdt_found = true;
		}
		else
		if (std::string(buffer, rmsd_str_len) == rmsd_str)
		{
			if (sscanf(buffer + rmsd_str_len, "%lf", rmsd) != 1)
			{
				std::cerr << "Error: no number after \""
					<< rmsd_str << "\" in output of command: \""
					<< command << "\"\n";
				pclose(pipe);
				exit(1);
			}

			rmsd_found = true;
		}
	}

	if (!gdt_found || !rmsd_found)
	{
		std::cerr
			<< "Error: GDT_TS or RMSD line not found in output of command: \""
			<< command << "\"\n";
		exit(1);
	}

	pclose(pipe);
}

// modify (phi, psi) angles by up to a certain amount
// The resulting angles are guaranteed to be inside a valid bin

void modify_angles(double *phi_val, double *psi_val, int amino_num,
	double amount)
{
	double phi_deg = rad2deg(*phi_val);
	double psi_deg = rad2deg(*psi_val);

	assert(phi_deg >= 0.0 && phi_deg < 360.0);
	assert(psi_deg >= 0.0 && psi_deg < 360.0);

	int phi_bin = (int) (phi_deg / (double) TORSION_SIZE);
	int psi_bin = (int) (psi_deg / (double) TORSION_SIZE);

	int limit = 10000;

	if (!valid_bin[phi_bin][psi_bin][amino_num])
	{
		/*
		std::cerr
			<< "Warning: modify_angle() called with angles not inside a bin:\n"
			<< "phi = " << Printf("%.1f", phi_deg) << "\n"
			<< "psi = " << Printf("%.1f", psi_deg) << "\n"
			<< "amino acid = " << Amino(amino_num).abbr() << "\n";
		*/
		limit = 100;
	}

	// Try random angles until a valid bin is found
	//
	// Should always find a bin, since we know the current
	// bin is valid (and it will be in the current bin if the random
	// numbers are close to 0.0)

	double amount2 = amount * 2.0;

	for (int count = 0; ;count++)
	{
		double h = range_0_360(phi_deg + (drand48() * amount2) - amount);
		double s = range_0_360(psi_deg + (drand48() * amount2) - amount);

		phi_bin = (int) (h / (double) TORSION_SIZE);
		psi_bin = (int) (s / (double) TORSION_SIZE);

		if (valid_bin[phi_bin][psi_bin][amino_num] || count == limit)
		{
			*phi_val = deg2rad(h);
			*psi_val = deg2rad(s);
			break;
		}
	}
}

void create_one_decoy(Peptide &p, const char *outfilename,
	const char *pdb_filename, double target_gdt, double *gdt_ptr,
	double *rmsd_ptr, double deviation_deg, FILE *log_file,
	int *iterations, int *clashes)
{
	Conformation orig_conf = p.conf();
	*gdt_ptr = 999.9;
	*clashes = 0;

#ifdef DO_LOG
	fprintf(log_file, "\n%s\n\n", outfilename);
	fflush(log_file);
#endif // DO_LOG

	int count;
	for (count = 0;*gdt_ptr > target_gdt;count++)
	{
		if (show_dots) { std::cerr << '.'; }

		Conformation conf = p.conf();
		Transform t;
		t.set_to_identity();

		// don't change phi/psi in first/last residue, since we don't have a
		// Ramachandran plot to check for valid angles

		for (int n = 1;n < p.length() - 1;n++)
		{
			double phi = p.conf().phi(n);
			double psi = p.conf().psi(n);
			double new_phi = phi;
			double new_psi = psi;

			assert(phi != TORSION_UNKNOWN);
			assert(psi != TORSION_UNKNOWN);

			modify_angles(&new_phi, &new_psi, p.res(n).amino().num(),
				deviation_deg);

			conf.transform_pos(n, Atom_N, t);
			conf.transform_pos(n, Atom_CA, t);
			conf.transform_pos(n, Atom_CB, t);

			Transform phi_trans(
				p.conf().pos(n, Atom_N),
				p.conf().pos(n, Atom_CA),
				new_phi - phi);

			t = t.times(phi_trans);

			conf.transform_pos(n, Atom_C, t);
			conf.transform_pos(n, Atom_O, t);

			Transform psi_trans(
				p.conf().pos(n, Atom_CA),
				p.conf().pos(n, Atom_C),
				new_psi - psi);

			t = t.times(psi_trans);
		}

		// transform the atoms in the last residue

		int i = p.length() - 1;
		Residue &res = p.res(i);

		for (int a = 0;a < res.num_atoms();a++)
		{
			conf.transform_pos(i, res.atom(a).type().type(), t);
		}

		// give p the new conformation
		p.conf().swap(conf);

		if (lj_scorer->steric_clash(p, 20.0))
		{
			(*clashes)++;
			// revert to previous conformation
			p.conf().swap(conf);
			continue;
		}

		p.write_pdb(outfilename);
		call_tmscore(pdb_filename, outfilename, gdt_ptr, rmsd_ptr);

#ifdef DO_LOG
		fprintf(log_file, "#%d %s: GDT = %f, RMSD = %f",
			count + 1, outfilename, *gdt_ptr, *rmsd_ptr);

		if (*gdt_ptr <= target_gdt)
		{
			fprintf(log_file, ": reached target of %f\n", target_gdt);
		}

		fprintf(log_file, "\n");
		fflush(log_file);
#endif // DO_LOG

	}

#ifdef DO_LOG
	fprintf(log_file, "%s\n\n", get_date().c_str());
#endif // DO_LOG
	
	if (show_dots) { std::cerr << "\n"; }

	// restore the original conformation
	p.conf().swap(orig_conf);
	*iterations = count;
}

void create_decoys(const char *filename, const std::string &pdb_id,
	char chain, FILE *result_file)
{
	Peptide p;

	if (!p.read_pdb(filename, chain, true))
	{
		std::cerr << "Error in PDB file " << filename << "\n";
		exit(1);
	}

	std::cerr << get_date() << "\n";

	chain = p.chain();	// (just in case chain was ' ')
	p.remove_non_backbone_atoms();
	p.conf().calc_torsion_angles();
	double deviation_deg = 50.0 / (double) p.length();

	if (lj_scorer->steric_clash(p))
	{
		std::cerr << "Warning: steric clashes found in " << filename << "\n";
	}

/*
if (self_clash)
{
	std::cout << filename << ": clash!\n";
}
else
{
	std::cout << filename << ": OK!\n";
}
return;
*/

#ifdef DO_LOG
	char log_file_name[1000];
	sprintf(log_file_name, "%s_%c.log", pdb_id.c_str(), chain);
	C_File log_file(log_file_name, "w", "Decoy log file");

	fprintf(log_file, "%s\nCreating decoys for %s\nLength = %d, "
		"deviation = %.3f degrees\n\n",
		get_date().c_str(), filename, p.length(), deviation_deg);
#else
	FILE *log_file = NULL;
#endif // DO_LOG

/*
	if (self_clash)
	{
		fprintf(result_file, "%s*   Steric clash!\n", filename);
		return;
	}
*/
	double gdt, rmsd;

	for (int n = 0;n < NUM_DECOYS;n++)
	{
		double target_gdt = 1.0 -
		 	((double) n / (double) (NUM_DECOYS - 1)) * (1.0 - MIN_GDT);

		/*
		// set the first GDT value to a bit under 1.0
		// (0.997 for 20 decoys with min GDT of 0.4)

		double target_gdt = 1.0 -
			((double) (n + 0.1) /
				(double) (NUM_DECOYS - 0.9)) * (1.0 - MIN_GDT);
		*/

		char outfilename[1000];
		sprintf(outfilename, "%s_%c.Decoy_%02d.pdb", pdb_id.c_str(), chain, n);

		std::cerr << outfilename;
		
		if (show_dots) { std::cerr << ' '; }
		else { std::cerr << '\n'; }

		int count, clash_count;
		create_one_decoy(p, outfilename, filename, target_gdt, &gdt,
			&rmsd, deviation_deg, log_file, &count, &clash_count);

		fprintf(result_file, "%s  GDT %.3f (%.3f)  RMSD %.3f  %d(%d)\n",
			outfilename, gdt, target_gdt, rmsd, count, clash_count);
		fflush(result_file);
	}
}

void read_data(const char *filename)
{
	C_File file(filename, "r", "Data file");

	const int Max_Len = 1000;
	char buffer[Max_Len];

	for (int aa = 0;aa < Amino::Num;aa++)
	{
		if (!file.next_line(buffer, Max_Len) ||
			std::string(buffer, 3) !=  Amino(aa).abbr())
		{
			std::cerr << "Error: expected \""
				<<  Amino(aa).abbr()
				<< "\" on line "
				<< file.line_num()
				<< " of data file " << filename
				<< "\n";
			exit(1);
		}

		for (int a1 = 0;a1 < TORSION_BINS;a1++)
		{
			if (!file.next_line(buffer, Max_Len) ||
				strlen(buffer) != TORSION_BINS + 1)
			{
				std::cerr << "Error: expected row of "
					<< TORSION_BINS
					<< " characters on line "
					<< file.line_num()
					<< " of data file " << filename
					<< "\n";
				exit(1);
			}

			for (int a2 = 0;a2 < TORSION_BINS;a2++)
			{
				char ch = buffer[a2];

				if (ch != '.' && ch != 'x' && ch != '+')
				{
					std::cerr << "Error: illegal character \""
						<< ch
						<< "\" on line "
						<< file.line_num()
						<< " of data file " << filename
						<< "\n";
					exit(1);
				}

				valid_bin[a1][a2][aa] = (ch != '.');
			}
		}
	}
}

std::string extract_pdb_id(const char *pdb_filename)
{
	const char *p = strrchr(pdb_filename, '/');

	if (p == NULL)
	{
		p = pdb_filename;
	}
	else
	{
		p++;
	}

	return std::string(p, 4);
}

void read_list_file(const char *list_filename, const char *pdb_dir,
	FILE *result_file)
{
	C_File list_file(list_filename, "r", "PDB list file");

	const int Max_Len = 1000;
	char buffer[Max_Len];

	while (list_file.next_line(buffer, Max_Len))
	{
		char pdb_filename[1000];
		char chain;

		if (sscanf(buffer, "%s %c", pdb_filename, &chain) != 2)
		{
			std::cerr << "Error on line " << list_file.line_num()
				<< " of PDB list file " << list_filename
				<< "\n";
			exit(1);
		}

		char full_pdb_filename[1000];
		sprintf(full_pdb_filename, "%s/%s", pdb_dir, pdb_filename);

		create_decoys(full_pdb_filename, extract_pdb_id(pdb_filename),
			chain, result_file);
	}

	std::cerr << "Finished at " << get_date() << "\n";
}

void init()
{
	long seed = time(0);
	srand((unsigned) seed);
	srand48(seed);
}

void show_usage(const char **argv)
{
	std::cerr
		<< "Usage: " << argv[0] << " <data file> <PDB dir> <PDB list>\n"
		<< "            or <data file> <PDB file> [<chain>]\n\n"
		<< "Data file is like:\n\n"
		<< "ALA\n\n"
		<< "++xxxxxxxx++.........++xxx+.........\n"
		<< "+xxxxxxxxxx++........+xxx++.........\n"
		<< "...\n\n(36 chars x 36 lines per amino acid)\n\n"
		<< "PDB list is like:\n\n"
		<< "1AX3 A\n"
		<< "2BC1 X\n"
		<< "...\n\n";
}

int main(int argc, const char **argv)
{
	if (argc != 3 && argc != 4)
	{
		show_usage(argv);
		exit(0);
	}

	const char *data_filename = argv[1];

	init();
	read_data(data_filename);

	char result_filename[100];
	sprintf(result_filename, "result_%ld", (long) getpid());
	C_File result_file(result_filename, "w", "Results file");

	lj_scorer = new Lennard_Jones;

	if (argc == 3 || argc == 4 && strlen(argv[3]) == 1)
	{
		const char *pdb_filename = argv[2];
		char chain = toupper(argc == 3 ? ' ' : argv[3][0]);
		create_decoys(pdb_filename, extract_pdb_id(pdb_filename),
			chain, result_file);
	}
	else
	{
		const char *pdb_dir = argv[2];
		const char *pdb_list = argv[3];
		read_list_file(pdb_list, pdb_dir, result_file);
	}

	delete lj_scorer;
	return 0;
}

