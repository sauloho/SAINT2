//#include "string.h" // PG added this
//#include <stdlib.h> // PG added this
#include <cstdlib> // PG added this
#include <cstring> // PG added this
#include <iostream>
#include <string>
#include "geom.h"
#include "c_file.h"
#include "point.h"
#include "geom.h"
#include "conformation.h"
#include "stream_printf.h"

#define MAX_FRAG_SIZE 1000
enum BB_Atom { BB_Atom_N, BB_Atom_CA, BB_Atom_C, Num_BB, BB_Atom_Undef };

const char *to_string(BB_Atom a)
{
	switch (a)
	{
		case BB_Atom_N:	return "N";
		case BB_Atom_CA:	return "CA";
		case BB_Atom_C:	return "C";
		default:		return "???";
	}
}

// (len is the length of the fragment, but the previous and next
// residue is included in pos[] as well, so the actual fragment is
// from 1 to len, with an extra residue at 0 and len + 1)
// (if shorter is true, then there are 10 residues rather than 11 --
// implying that this is at the very start or end of the chain)

void print_fragment(int fragment_num, int tpos, bool shorter, int len,
	const Point pos[MAX_FRAG_SIZE][Num_BB], double score, const char *info)
{
	std::cout
		<< "F " << fragment_num
		<< " P " << tpos
		<< " L " << len
		<< " S " << score
		<< " = " << info
		<< "\n";

	bool missing_first = (shorter && tpos == 0);
	bool missing_last = (shorter && tpos != 0);

	int n1 = (missing_first ? 0 : 1);
	int n2 = n1 + len - 1;

	for (int n = n1;n <= n2;n++)
	{
		double phi, n_angle, c_angle;

		if (n == 0)
		{
			phi = TORSION_UNKNOWN;
			n_angle = BOND_ANGLE_C_N_CA;
		}
		else
		{
			phi = torsion_angle(
				pos[n-1][BB_Atom_C],
				pos[n][BB_Atom_N],
				pos[n][BB_Atom_CA],
				pos[n][BB_Atom_C]);

			n_angle = angle_formed(pos[n-1][BB_Atom_C],
				pos[n][BB_Atom_N], pos[n][BB_Atom_CA]);
		}

		double psi, omega;

		if (n == n2 && missing_last)
		{
			psi = omega = TORSION_UNKNOWN;
			c_angle = BOND_ANGLE_CA_C_N;
		}
		else
		{
			psi = torsion_angle(
				pos[n][BB_Atom_N],
				pos[n][BB_Atom_CA],
				pos[n][BB_Atom_C],
				pos[n+1][BB_Atom_N]);

			omega = torsion_angle(
				pos[n][BB_Atom_CA],
				pos[n][BB_Atom_C],
				pos[n+1][BB_Atom_N],
				pos[n+1][BB_Atom_CA]);

			c_angle = angle_formed(pos[n][BB_Atom_CA],
				pos[n][BB_Atom_C], pos[n+1][BB_Atom_N]);

		}

		double ca_angle = angle_formed(pos[n][BB_Atom_N],
			pos[n][BB_Atom_CA], pos[n][BB_Atom_C]);

		double phi_deg = (phi == TORSION_UNKNOWN ? phi :
			 range_m180_180(rad2deg(phi)));

		double psi_deg = (psi == TORSION_UNKNOWN ? psi :
			 range_m180_180(rad2deg(psi)));

		double omega_deg = (omega == TORSION_UNKNOWN ? omega :
			 range_m180_180(rad2deg(omega)));

		std::cout << n - n1 << ' '
			<< Printf("%6.1f ", phi_deg)
			<< Printf("%6.1f ", psi_deg)
			<< Printf("%6.1f ", omega_deg)
			<< Printf("%6.1f ", rad2deg(n_angle))
			<< Printf("%6.1f ", rad2deg(ca_angle))
			<< Printf("%6.1f ", rad2deg(c_angle))
			<< '\n';
	}
}

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " fragment_library_file\n";
		exit(1);
	}

	static const int Max_Len = 1000;
	char buffer[Max_Len];
	char str[Max_Len];
	int first_tpos = -1;
	bool first_tpos_known = false;
	int prev_tpos = -1;
	int fragment_num = 0;

	const char *filename = argv[1];
	C_File file(filename, "r", "Fragment library file");
	fgets(buffer, Max_Len, file);
	int line_num = 1;
	bool any_shorter_beyond_0 = false;
	int any_shorter_beyond_0_last_line_num = -1;
	int last_frag_line_num = -1;

	while (!feof(file))
	{
		int tpos;
		if (sscanf(buffer, "%s %d", str, &tpos) != 2 ||
			std::string(str) != "TPOS:")
		{
			std::cerr << "Error on line " << line_num
				<< " of " << filename
				<< ": expected \"TPOS:\" followed by a number\n";
			exit(1);
		}

		std::cerr << buffer;

		if (!first_tpos_known)
		{
			first_tpos_known = true;
			first_tpos = tpos;
		}
		else
		{
			if (tpos != prev_tpos + 1)
			{
				std::cerr << "Error on line " << line_num
					<< " of " << filename
					<< ": TPOS = " << tpos
					<< ", expected " << prev_tpos + 1
					<< "\n";
				exit(1);
			}

			if (any_shorter_beyond_0)
			{
				std::cerr << "Error: fragment of length 10 found that is "
					"not at the start or end of the target; line "
					<< any_shorter_beyond_0_last_line_num
					<< " of "
					<< filename
					<< "\n";
				exit(1);
			}
		}

		prev_tpos = tpos;
		fgets(buffer, Max_Len, file);
		line_num++;

		while (!feof(file) && std::string(buffer, 3) == "PDB")
		{
			if (fgets(buffer, Max_Len, file) == NULL)
			{
				std::cerr << "Unexpected end of file #1 in "
					<< filename << "\n";
				exit(1);
			}

			line_num++;
			last_frag_line_num = line_num;

			int len;
			double score;
			char info[Max_Len];
			char pdb_id[100];
			char res[100];
			char rmsd_str[100];
			char chain;

			std::cerr << buffer;
			if (sscanf(buffer, "%s %c %s %d %s %lf",
				pdb_id, &chain, res, &len, rmsd_str, &score) != 6)
			{
				std::cerr << "Error on line " << line_num
					<< " of " << filename
					<< ": expected <pdbid> <chain> <res_seq> <len> "
						"<rmsd> <score>\n";
				exit(1);
			}

			sprintf(info, "%s %c %s R %s", pdb_id, chain, res, rmsd_str);

			if (len + 2 >= MAX_FRAG_SIZE)
			{
				std::cerr << "Error on line " << line_num
					<< " of " << filename
					<< ": fragment size too large ("
					<< len << ", max is "
					<< MAX_FRAG_SIZE - 2
					<< ")\n";
				exit(1);
			}

			Point pos[MAX_FRAG_SIZE][Num_BB];
			bool got_pos[MAX_FRAG_SIZE][Num_BB];
			std::string prev_res_seq;
			int n;

			for (n = 0;n < len + 2;n++)
			{
				for (int a = 0;a < Num_BB;a++)
				{
					got_pos[n][a] = false;
				}
			}

			// not all atoms in the first and last residue are
			// required to calculate the torsion angles, so flag
			// them as already present to avoid warning messages

			/// only need C in first residue
			got_pos[0][BB_Atom_N] = got_pos[0][BB_Atom_CA] = true;

			// only need N and CA in last residue
			got_pos[0][BB_Atom_C] = true;

			n = -1;

			while (fgets(buffer, Max_Len, file) != NULL)
			{
				line_num++;

				if (std::string(buffer, 4) != "ATOM")
				{
					break;
				}

				char res_seq[Max_Len];

				if (strlen(buffer) < 47 ||
					sscanf(buffer + 22, "%s", res_seq) != 1)
				{
					std::cerr << "Missing data on line " << line_num
						<< " of " << filename
						<< "\n";
					exit(1);
				}

				if (std::string(res_seq) != prev_res_seq)
				{
					n++;

					// std::cerr << "Seq " << res_seq << " = " << n << "\n";

					// if (n >= len + 2)
					if (n >= MAX_FRAG_SIZE)
					{
						// error message printed after loop
						break;
					}

					prev_res_seq = res_seq;
				}

				std::string atom = std::string(buffer + 12, 4);
				BB_Atom a;

				if (atom == " CA ") { a = BB_Atom_CA; }
				else if (atom == " C  ") { a = BB_Atom_C; }
				else if (atom == " N  ") { a = BB_Atom_N; }
				else { a = BB_Atom_Undef; }

				if (a != BB_Atom_Undef)
				{
    				if (sscanf(std::string(buffer + 30, 8).c_str(), "%lf",
							&pos[n][a].x) != 1 ||
        				sscanf(std::string(buffer + 38, 8).c_str(), "%lf",
							&pos[n][a].y) != 1 ||
        				sscanf(std::string(buffer + 46, 8).c_str(), "%lf",
							&pos[n][a].z) != 1)
					{
						std::cerr << "Missing coordinates on line "
							<< line_num << " of " << filename
							<< "\n";
						exit(1);
					}

					got_pos[n][a] = true;
				}
			}

			//!!! TO DO: (++n == len)
			//(Sheenal needs to change file format to give actual number
			//of residues, not always 9 or 11)

			bool shorter = false;

			if (++n != len + 2)
			{
				if (n == len + 1)
				{
					shorter = true;
				}
				else
				{
					std::cerr
						<< "Warning: Incorrect number of residues in "
							"fragment before line "
						<< line_num
						<< " of " << filename
						<< " (" << n
						<< " residues, expected "
						<< len << " + 2 = " << len + 2
						<< ")\n";
					continue;
				}
			}

			bool any_missing = false;

			for (int m = 0;m < n && !any_missing;m++)
			{
				for (int a = 0;a < Num_BB;a++)
				{
					if (!got_pos[m][a])
					{
						std::cerr << "Warning: missing "
							<< to_string((BB_Atom) a)
							<< " in residue #"
							<< m + 1
							<< " of fragment before line "
							<< line_num
							<< " of " << filename
							<< "\n";
						any_missing = true;
						break;
					}
				}
			}

			if (any_missing)
			{
				continue;
			}

			print_fragment(fragment_num, tpos - first_tpos, shorter,
				len, pos, score, info);
			fragment_num++;

			if (shorter && tpos != first_tpos)
			{
				any_shorter_beyond_0 = true;
				any_shorter_beyond_0_last_line_num = last_frag_line_num;
			}
		}
	}

	return 0;
}

