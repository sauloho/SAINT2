
#include <iostream>
#include <cmath>
#include <ctime>
#include <cassert>
#include "config.h"
#include "sequence.h"
#include "runner.h"
#include "scorer.h"
#include "reporter.h"
#include "static_init.h"
#include "orientation.h"
//#include "pairwise.h"
#include "torsion.h"
#include "common.h"
#include "c_file.h"
#include "geom.h"

#include "atom_type.h"

//#define MIN_RELIABLE_VAL 15
#define SOLV_RELIABLE_VAL 7
#define TOR_RELIABLE_VAL 7
#define RAPDF_RELIABLE_VAL 7
#define ORIENT_RELIABLE_VAL 7

// definition of "near" = within 10 Angstroms
#define SOLV_DIST 10.0
#define SOLV_BIN_TOP 38
#define FIRST_SOLV_BIN 3

#define RAPDF_TOP_BIN 20
#define FIRST_RAPDF_BIN 2

#define SQUARE(x) ((x)*(x))

long total_rapdf_pairs = 0;

void calc_solvation(const Peptide &p, long total_solv[SOLV_BIN_TOP][Amino::Num],
	const char *filename, char chain)
{
	for (int n = 0;n < p.length();n++)
	{
		const Residue &r = p.res(n);

		// number of CB atoms within SOLV_DIST Angstroms of this residue's CB
		long count = 0;

		Point pos;
		
		if (!p.get_cb_ca_pos(n, &pos))
		{
			continue;
		}

		if (pos.x == 0.0 && pos.y == 0.0 && pos.z == 0.0)
		{
			std::cerr
				<< "Error: Suspicious position (0, 0, 0) for CA/CB of residue "
				<< r.res_seq_str()
				<< " ("
				<< p.res(n).amino().abbr()
				<< ") on line "
				<< r.pdb_line()
				<< " of " << filename
				<< "\n";
			exit(1);
		}

		Point pos2;

		for (int m = 0;m < p.length();m++)
		{
			if (m != n && p.get_cb_ca_pos(m, &pos2))
			{
				if (fabs(pos.x - pos2.x) <= SOLV_DIST &&
					fabs(pos.y - pos2.y) <= SOLV_DIST &&
					fabs(pos.z - pos2.z) <= SOLV_DIST &&
					pos.distance(pos2) <= SOLV_DIST)
				{
					count++;
				}
			}
		}

		if (count == 0)
		{
			std::cout << "Warning: ignoring zero solvation for residue "
				<< r.res_seq_str()
				<< " on line "
				<< r.pdb_line()
				<< " of "
				<< filename
				<< "\n";
			continue;
		}

		if (count >= SOLV_BIN_TOP)
		{
			count = SOLV_BIN_TOP - 1;
		}
		else
		if (count < FIRST_SOLV_BIN)
		{
			count = FIRST_SOLV_BIN;
		}

		++(total_solv[count][r.amino().num()]);
	}
}

void calc_torsion(const Peptide &p,
	long total_torsion[TORSION_BINS][TORSION_BINS][Amino::Num],
	int num_torsion_aminos[TORSION_BINS][TORSION_BINS],
	const char *filename, char chain)
{
	// (a pair of phi/psi angles is required, so the first and last
	// residue (which only have one of these angles) are not included)

	int phi_bin, psi_bin;

	for (int n = 1;n < p.length() - 1;n++)
	{
		if (Torsion::get_phi_psi_bin(p, n, &phi_bin, &psi_bin))
		{
			int a = p.res(n).amino().num();
			
			if (++total_torsion[phi_bin][psi_bin][a] == 7)
			{
				++num_torsion_aminos[phi_bin][psi_bin];
			}
		}
	}
}

void calc_orient(const Peptide &p,
	long total_orient[ORIENT_DISTS][ORIENT_ANGLES][Amino::Num][Amino::Num])
{
	int dist_bin, angle_bin;

	for (int n = 0;n < p.length() - 1;n++)
	{
		int a1 = p.res(n).amino().num();

		for (int m = n + 1;m < p.length();m++)
		{
			if (Orientation::get_bins(p, n, m, &dist_bin, &angle_bin))
			{
				assert(dist_bin >= 0 && dist_bin < ORIENT_DISTS);
				assert(angle_bin >= 0 && angle_bin < ORIENT_ANGLES);

				int a2 = p.res(m).amino().num();

				++(total_orient[dist_bin][angle_bin][a1][a2]);
				++(total_orient[dist_bin][angle_bin][a2][a1]);
			}
		}
	}
}

/*
void calc_pairwise(const Peptide &p,
	long total_pairwise[PAIRWISE_DISTS][PAIRWISE_SEPS][Amino::Num][Amino::Num])
{
	int dist_bin, sep_bin;

	for (int n = 0;n < p.length() - 1;n++)
	{
		int a1 = p.res(n).amino().num();

		for (int m = n + 1;m < p.length();m++)
		{
			if (Pairwise::get_bins(p, n, m, &dist_bin, &sep_bin))
			{
				assert(dist_bin >= 0 && dist_bin < PAIRWISE_DISTS);
				assert(sep_bin >= 0 && sep_bin < PAIRWISE_SEPS);

				int a2 = p.res(m).amino().num();

				++(total_pairwise[dist_bin][sep_bin][a1][a2]);
				++(total_pairwise[dist_bin][sep_bin][a2][a1]);
			}
		}
	}
}
*/

char atom_char1(Atom_Id a)
{
	switch (a)
	{
		case Atom_N:	return 'N';
		case Atom_C:	return 'C';
		case Atom_O:	return 'O';
		case Atom_CA:	return 'A';
		case Atom_CB:	return 'B';
		default:		std::cerr << "Error in atom_char()\n";
						exit(1);
	}
}

void calc_rapdf(const Peptide &p,
	long total[NUM_RAPDF_IDS][NUM_RAPDF_IDS][RAPDF_TOP_BIN])
{
	for (int n1 = 1;n1 < p.length();n1++)
	{
		const Residue &res1 = p.res(n1);

		for (int a1 = 0;a1 < res1.num_atoms();a1++)
		{
			Atom_Id t1 = res1.atom(a1).type().type();
			if (t1 == Atom_Undef) { continue; }

			int id_1 = res1.amino().rapdf_id(t1);

			Point pos1 = p.atom_pos(n1, t1);

			for (int n2 = 0;n2 < n1;n2++)
			{
				const Residue &res2 = p.res(n2);

				// std::cout << n1 << ", " << n2 << std::endl;

				for (int a2 = 0;a2 < res2.num_atoms();a2++)
				{
					Atom_Id t2 = res2.atom(a2).type().type();
					if (t2 == Atom_Undef) { continue; }

					int id_2 = res2.amino().rapdf_id(t2);

					Point pos2 = p.atom_pos(n2, t2);

/*
if (t1 < Num_Backbone && t2 < Num_Backbone)
{
std::cout
<< "R: (" << n2 << ") " << res2.amino().name() << " " << Atom_Type(t2).name() << " id " << id_2
<< " at " << pos2
<< ", (" << n1 << ") " << res1.amino().name() << " " << Atom_Type(t1).name() << " id " << id_1
<< " at " << pos1 << "\n";
}
*/

					if (pos1.closer_than(RAPDF_TOP_BIN, pos2))
					{
						double dist = pos1.distance(pos2);

						/*
						if (dist < 4.0 &&
							t1 < Num_Backbone && t2 < Num_Backbone &&
							abs(n1 - n2) > 1)
						{
							char ch1 = atom_char1(t1);
							char ch2 = atom_char1(t2);

							if (ch1 > ch2)
							{
								char temp = ch1;
								ch1 = ch2;
								ch2 = temp;
							}

							std::cout << "=" << ch1 << ch2 << ' '
								<< dist << "\n";
						}
						*/

						int int_dist = (int) dist;
						assert(int_dist < RAPDF_TOP_BIN);

						if (int_dist < FIRST_RAPDF_BIN)
						{
							int_dist = FIRST_RAPDF_BIN;
						}

/*
if (t1 < Num_Backbone && t2 < Num_Backbone)
{ std::cout << " => dist " << dist << " = bin " << int_dist << "\n"; }
*/

						++total[id_1][id_2][int_dist];
						++total_rapdf_pairs;
					}
				}
			}
		}
	}
}

void fix_missing_data(double val[], bool val_known[], bool val_unreliable[],
	int from, int top, const char *comment[])
{
	// find first known/unreliable values

	int lower, upper;
	for (lower = from;lower < top;lower++)
	{
		if (val_known[lower] || val_unreliable[lower])
		{
			break;
		}
	}

	int n;
	for (n = from;n < top;n++)
	{
		comment[n] = "";
	}

	if (lower == top)
	{
		// all values bad

		for (n = from;n < top;n++)
		{
			val[n] = 0.0;
			comment[n] = "all unknown";
		}

		return;
	}

	// find last known/unreliable value

	for (upper = top - 1;upper > lower;upper--)
	{
		if (val_known[upper] || val_unreliable[upper])
		{
			break;
		}
	}

	// check for unreliable values -- interpolate from nearest values

	int last_known = -1;

	for (n = lower;n <= upper;n++)
	{
		if (val_known[n])
		{
			last_known = n;
		}

		if (val_unreliable[n] || !val_known[n])
		{
			// find the next defined value
			int m;
			for (m = n + 1;m <= upper;m++)
			{
				if (val_known[m])
				{
					break;
				}
			}

			if (last_known == -1)
			{
				if (m > upper)
				{
					// there are no known values, only unreliable ones
					
					for (int x = from;x < top;x++)
					{
						val[x] = 0.0;
						comment[x] = "all sparse data";
					}

					return;
				}

				last_known = m;

				for (int x = n;x < last_known;x++)
				{
					val[x] = val[last_known];
					comment[x] = "sparse data";
				}

				n = last_known;
				continue;
			}

			if (m > upper)
			{
				// there are no more known values

				for (int x = n;x <= upper;x++)
				{
					val[x] = val[last_known];
					comment[x] = "sparse data";
				}

				break;
			}
			else
			{
				// interpolate between val[last_known] and val[m]
				int d1 = m - n;
				int d2 = n - last_known;

				val[n] = ((d1 * val[last_known]) + (d2 * val[m])) / ((double) (m - last_known));
				comment[n] = "interpolated";
			}
		}
	}

	for (n = lower - 1;n >= from;n--)
	{
		val[n] = val[n+1];
		comment[n] = "no data";
	}

	for (n = upper + 1;n < top;n++)
	{
		val[n] = val[n-1];
		comment[n] = "no data";
	}
}

void fix_missing_torsion_data(double val[TORSION_BINS][TORSION_BINS],
	bool val_known[TORSION_BINS][TORSION_BINS])
{
	int phi_bin, psi_bin;

	// set unknown values

	for (phi_bin = 0;phi_bin < TORSION_BINS;phi_bin++)
	{
		for (psi_bin = 0;psi_bin < TORSION_BINS;psi_bin++)
		{
			if (!val_known[phi_bin][psi_bin])
			{
				int best_d_squared = -1;
				int num = 0;
				double t = 0;

				for (int h = 0;h < TORSION_BINS;h++)
				{
					for (int s = 0;s < TORSION_BINS;s++)
					{
						if (!val_known[h][s])
						{
							continue;
						}

						int dh = abs(phi_bin - h);

						if (dh > TORSION_BINS / 2)
						{
							dh = TORSION_BINS - dh;
						}

						int ds = abs(psi_bin - s);

						if (ds > TORSION_BINS / 2)
						{
							ds = TORSION_BINS - ds;
						}

						int d = ((dh * dh) + (ds * ds));

						if (best_d_squared == -1 || d < best_d_squared)
						{
							best_d_squared = d;
							t = val[h][s];
							num = 1;
						}
						else
						if (d == best_d_squared)
						{
							t += val[h][s];
							num++;
						}
					}
				}

				if (num == 0)
				{
					val[phi_bin][psi_bin] = 0.0;
				}
				else
				{
					// t/num is the closest (or average of equal closest)
					// value(s)
					val[phi_bin][psi_bin] = (double) t / (double) num;
				}
			}
		}
	}
}

void generate_solvation(long total_aa[Amino::Num],
	long total_solv[SOLV_BIN_TOP][Amino::Num],
	long total_res, const char *cmd)
{
	static const char *filename = "solv.data";
	C_File out(filename, "w", "Solvation data file");
	std::cerr << "Creating " << filename << "\n";

	fprintf(out, "# Solvation data file - %s - min reliable %d\n# %s\n",
		get_date().c_str(), SOLV_RELIABLE_VAL, cmd);

	fprintf(out, "\n# Min count, max count, radius\n\n%d %d %.1f\n",
		FIRST_SOLV_BIN, SOLV_BIN_TOP - 1, SOLV_DIST);

	long num_in_bin[SOLV_BIN_TOP];
	int bin, a;

	for (bin = FIRST_SOLV_BIN;bin < SOLV_BIN_TOP;bin++)
	{
		long total = 0;

		for (a = 0;a < Amino::Num;a++)
		{
			total += total_solv[bin][a];
		}

		num_in_bin[bin] = total;
	}

	assert(total_res != 0);

	for (a = 0;a < Amino::Num;a++)
	{
		assert(total_aa[a] != 0);

		double val[SOLV_BIN_TOP];
		bool val_known[SOLV_BIN_TOP];
		bool val_unreliable[SOLV_BIN_TOP];

		for (bin = FIRST_SOLV_BIN;bin < SOLV_BIN_TOP;bin++)
		{
			val_known[bin] = false;
			val_unreliable[bin] = false;
		}

		for (bin = FIRST_SOLV_BIN;bin < SOLV_BIN_TOP;bin++)
		{
			if (total_solv[bin][a] >= SOLV_RELIABLE_VAL)
			{
				double ratio1 = total_solv[bin][a] / (double) total_aa[a];
				double ratio2 = num_in_bin[bin] / (double) total_res;
				val[bin] = -log(ratio1 / ratio2);
				val_known[bin] = true;

				/*
				std::cout
					<< Amino(a).abbr() << " " << bin
					<< ": r1 = "
					<< total_solv[bin][a] << " / " << total_aa[a]
					<< " = " << ratio1
					<< ", r2 = "
					<< num_in_bin[bin] << " / " << total_res
					<< " = " << ratio2
					<< ", value = "
					<< val[bin] << "\n";
				*/
			}
			else
			{
				if (total_solv[bin][a] != 0)
				{
					val_unreliable[bin] = true;
				}

				/*
				double ratio1 = total_solv[bin][a] / (double) total_aa[a];
				double ratio2 = num_in_bin[bin] / (double) total_res;

				std::cout
					<< Amino(a).abbr() << " " << bin
					<< ": r1 = "
					<< total_solv[bin][a] << " / " << total_aa[a]
					<< " = " << ratio1
					<< ", r2 = "
					<< num_in_bin[bin] << " / " << total_res
					<< " = " << ratio2
					<< ", value = " << (total_solv[bin][a] == 0 ? "undefined" : "unreliable")
					<< "\n";
				*/
			}
		}

		const char *comment[SOLV_BIN_TOP];
		fix_missing_data(val, val_known, val_unreliable, FIRST_SOLV_BIN,
			SOLV_BIN_TOP, comment);

		fprintf(out, "\n%s\n\n", Amino(a).abbr());

		for (bin = FIRST_SOLV_BIN;bin < SOLV_BIN_TOP;bin++)
		{
			/*
			std::cout << "= " << Amino(a).abbr() << " " << bin << " "
				<< val[bin];
			if (!val_known[bin]) { std::cout << " !"; }
			std::cout << "\n";
			*/

			fprintf(out, "%-2d %.3f", bin, val[bin]);

			if (*comment[bin])
			{
				fprintf(out, "\t# %s", comment[bin]);
			}

			fprintf(out, "\n");
		}
	}
}

void generate_torsion(long total_aa[Amino::Num],
	long total_torsion[TORSION_BINS][TORSION_BINS][Amino::Num],
	int num_torsion_aminos[TORSION_BINS][TORSION_BINS],
	long total_res, const char *cmd)
{
	static const char *filename = "torsion.data";
	C_File out(filename, "w", "Torsion data file");
	std::cerr << "Creating " << filename << "\n";

	fprintf(out, "# Torsion data file - %s\n# %s\n#\n"
		"# Using average closest value for undefined data\n",
		get_date().c_str(), cmd);
	
//fprintf(out, "# Using Kullback-Leibler idea (power of denominator)\n");
//fprintf(out, "# Weighted by number of amino acid types in each bin\n");
	fprintf(out, "# %d x %d total bins (%d x %d degrees), min reliable = %d\n",
		TORSION_BINS, TORSION_BINS,
		360 / TORSION_BINS, 360 / TORSION_BINS,
		TOR_RELIABLE_VAL);

	long num_in_bin[TORSION_BINS][TORSION_BINS];
	int phi_bin, psi_bin, a;

	for (phi_bin = 0;phi_bin < TORSION_BINS;phi_bin++)
	{
		for (psi_bin = 0;psi_bin < TORSION_BINS;psi_bin++)
		{
			long total = 0;

			for (a = 0;a < Amino::Num;a++)
			{
				total += total_torsion[phi_bin][psi_bin][a];
			}

			num_in_bin[phi_bin][psi_bin] = total;
		}
	}

	for (a = 0;a < Amino::Num;a++)
	{
		assert(total_aa[a] != 0);

		double val[TORSION_BINS][TORSION_BINS];
		bool val_known[TORSION_BINS][TORSION_BINS];

		for (phi_bin = 0;phi_bin < TORSION_BINS;phi_bin++)
		{
			for (psi_bin = 0;psi_bin < TORSION_BINS;psi_bin++)
			{
				val_known[phi_bin][psi_bin] = false;
			}
		}

		for (phi_bin = 0;phi_bin < TORSION_BINS;phi_bin++)
		{
			for (psi_bin = 0;psi_bin < TORSION_BINS;psi_bin++)
			{
				double t = total_torsion[phi_bin][psi_bin][a];

				if (t >= TOR_RELIABLE_VAL)
				{
					double ratio1 = t / (double) total_aa[a];
					double ratio2 = num_in_bin[phi_bin][psi_bin] /
						(double) total_res;

				/*
					// Kullback-Liebler distance (Graham's idea):
					// val[phi_bin][psi_bin] = -ratio2 * log(ratio1 / ratio2);

					double val1 = -log(ratio1 / ratio2);
					double val2 = -log(exp(ratio1 * 1.4 - 0.7));

					assert(num_torsion_aminos[phi_bin][psi_bin] > 0);
					double w1 = (num_torsion_aminos[phi_bin][psi_bin] - 1) /
							(double) Amino::Num;
					// #1->0, #2->0.33 ... #5->0.95, #6..20->1.0
					w1 = log(1.0 + w1 * 8.0);
					if (w1 > 1.0) { w1 = 1.0; }
					double w2 = 1.0 - w1;

					val[phi_bin][psi_bin] = (val1 * w1) + (val2 * w2);
				*/
					val[phi_bin][psi_bin] = -log(ratio1 / ratio2);
					val_known[phi_bin][psi_bin] = true;
				}
			}
		}

		fix_missing_torsion_data(val, val_known);

		fprintf(out, "\n%s\n", Amino(a).abbr());

		// print Ramachandran plot

		for (psi_bin = 0;psi_bin < TORSION_BINS;psi_bin++)
		{
			//int s = (psi_bin + TORSION_BINS / 2) % TORSION_BINS;
			//s = TORSION_BINS - 1 - s;
			int s = psi_bin;

			fprintf(out, "\n# ");

			for (phi_bin = 0;phi_bin < TORSION_BINS;phi_bin++)
			{
				//int h = (phi_bin + TORSION_BINS / 2) % TORSION_BINS;
				int h = phi_bin;

				char ch;

				if (val[h][s] == 0.0)
				{
					ch = '.';
				}
				else
				if (val[h][s] >= 0.0)
				{
					ch = ((int) (val[h][s] * 4)) + 'A';
					if (ch > 'Z') { ch = 'Z'; }
				}
				else
				{
					ch = ((int) (-val[h][s] * 4)) + '0';
					if (ch > '9') { ch = '9'; }
				}

				putc(val_known[h][s] ?
					(val[h][s] < 0.0 ? '-' : '+') :
					' ', out);
				putc(ch, out);
			}
		}

		fprintf(out, "\n\n");

		for (phi_bin = 0;phi_bin < TORSION_BINS;phi_bin++)
		{
			for (psi_bin = 0;psi_bin < TORSION_BINS;psi_bin++)
			{
				fprintf(out,
					"%-2d %-2d %.3f   # (%ld / %ld) / (%ld / %ld)%s [%d]\n",
					phi_bin, psi_bin,
					val[phi_bin][psi_bin],
					total_torsion[phi_bin][psi_bin][a],
					total_aa[a],
					num_in_bin[phi_bin][psi_bin],
					total_res,
					(val_known[phi_bin][psi_bin] ? "" : " ?"),
					num_torsion_aminos[phi_bin][psi_bin]
					);
			}
			
			fprintf(out, "\n");
		}
	}
}

void generate_orient(
	long total_orient[ORIENT_DISTS][ORIENT_ANGLES][Amino::Num][Amino::Num],
	const char *cmd)
{
	static const char *filename = "orient.data";
	C_File out(filename, "w", "Orientation data file");
	std::cerr << "Creating " << filename << "\n";

	fprintf(out, "# Orientation data file - %s - min reliable %d\n# %s\n",
		get_date().c_str(), ORIENT_RELIABLE_VAL, cmd);

	fprintf(out, "# %d distances, %d angles\n", ORIENT_DISTS, ORIENT_ANGLES);

	long num_in_bin[ORIENT_DISTS][ORIENT_ANGLES];
	long num_in_all_bins = 0;
	int a1, a2, a, d;

	for (d = 0;d < ORIENT_DISTS;d++)
	{
		for (a = 0;a < ORIENT_ANGLES;a++)
		{
			long t = 0;

			for (a1 = 0;a1 < Amino::Num;a1++)
			{
				for (a2 = 0;a2 < Amino::Num;a2++)
				{
					t += total_orient[d][a][a1][a2];
				}
			}

			num_in_bin[d][a] = t;
			num_in_all_bins += t;
		}
	}

	for (a1 = 0;a1 < Amino::Num;a1++)
	{
		for (a2 = a1;a2 < Amino::Num;a2++)
		{
			fprintf(out, "\n%s %s\n", Amino(a1).abbr(), Amino(a2).abbr());

			long total = 0;

			for (d = 0;d < ORIENT_DISTS;d++)
			{
				for (a = 0;a < ORIENT_ANGLES;a++)
				{
					total += total_orient[d][a][a1][a2];
				}
			}

			for (a = 0;a < ORIENT_ANGLES;a++)
			{
				const char *desc = "?";

				switch (a)
				{
					case 0: desc = "a1 < 30, a2 < 30"; break;
					case 1: desc = "a1 < 30, a2 30-90"; break;
					case 2: desc = "a1 < 30, a2 90+"; break;
					case 3: desc = "a1 30-90, a2 < 30"; break;
					case 4: desc = "a1 > 90, a2 < 30"; break;
					case 5: desc = "a1 > 90, a2 > 90"; break;
					case 6: desc = "a1 > 90, a2 30-90"; break;
					case 7: desc = "a1 30-90, a2 > 90"; break;
					case 8: desc = "a1 30-90, a2 30-90, torsion < 90"; break;
					case 9: desc = "a1 30-90, a2 30-90, torsion > 90"; break;
				}

				fprintf(out, "\n# %s\n\n", desc);

				for (d = 0;d < ORIENT_DISTS;d++)
				{
					double r1 = 0.0;
					
					if (total != 0 &&
						total_orient[d][a][a1][a2] >= ORIENT_RELIABLE_VAL)
					{
						r1 = total_orient[d][a][a1][a2] / (double) total;
					}

					double r2 = 0.0;
					
					if (num_in_bin[d][a] >= ORIENT_RELIABLE_VAL)
					{
						r2 = num_in_bin[d][a] / (double) num_in_all_bins;
					}

					double val = 0.0;
					bool kn = true;

					if (r1 != 0.0 && r2 != 0.0)
					{
						val = -log(r1 / r2);
					}
					else
					{
						val = 0.0;
						kn = false;
					}

					fprintf(out, "%-2d %6.3f \t# ((%d / %d) / (%d / %d))%s\n",
						d + 3, val,
						total_orient[d][a][a1][a2], total,
						num_in_bin[d][a], num_in_all_bins,
						(kn ? "" : " ?"));

					/*
					std::cout
						<< Amino(a1).abbr() << " " 
						<< Amino(a2).abbr()
						<< " a " << a
						<< " d " << d
						<< ": r1 " << total_orient[d][a][a1][a2]
						<< " / " << total
						<< " r2 "
						<< num_in_bin[d][a]
						<< " / "
						<< num_in_all_bins
						<< " -log(" << r1 << " / " << r2 << ")"
						<< " = "
						<< val
						<< "\n";
					*/
				}
			}
		}
	}
}

/*
void generate_pairwise(
	long total_pairwise[PAIRWISE_DISTS][PAIRWISE_SEPS][Amino::Num][Amino::Num],
	const char *cmd)
{
	static const char *filename = "pairwise.data";
	C_File out(filename, "w", "Pairwise data file");
	std::cerr << "Creating " << filename << "\n";

	fprintf(out, "# Pairwise data file - %s\n# %s\n",
		get_date().c_str(), cmd);

	long num_in_bin[PAIRWISE_DISTS][PAIRWISE_SEPS];
	long num_in_all_bins = 0;
	int a1, a2, s, d;

	for (d = 0;d < PAIRWISE_DISTS;d++)
	{
		for (s = 0;s < PAIRWISE_SEPS;s++)
		{
			long t = 0;

			for (a1 = 0;a1 < Amino::Num;a1++)
			{
				for (a2 = 0;a2 < Amino::Num;a2++)
				{
					t += total_pairwise[d][s][a1][a2];
				}
			}

			num_in_bin[d][s] = t;
			num_in_all_bins += t;
		}
	}

	for (a1 = 0;a1 < Amino::Num;a1++)
	{
		for (a2 = a1;a2 < Amino::Num;a2++)
		{
			fprintf(out, "\n%s %s\n\n", Amino(a1).abbr(), Amino(a2).abbr());

			long total = 0;

			for (d = 0;d < PAIRWISE_DISTS;d++)
			{
				for (s = 0;s < PAIRWISE_SEPS;s++)
				{
					total += total_pairwise[d][s][a1][a2];
				}
			}

			for (s = 1;s < PAIRWISE_SEPS;s++)
			{
				for (d = 0;d < PAIRWISE_DISTS;d++)
				{
					double r1 = 0.0;
					
					if (total != 0 &&
						total_pairwise[d][s][a1][a2] >= ORIENT_RELIABLE_VAL)
					{
						r1 = total_pairwise[d][s][a1][a2] / (double) total;
					}

					double r2 = 0.0;
					
					if (num_in_bin[d][s] >= ORIENT_RELIABLE_VAL)
					{
						r2 = num_in_bin[d][s] / (double) num_in_all_bins;
					}

					double val = 0.0;
					bool kn = true;

					if (r1 != 0.0 && r2 != 0.0)
					{
						val = -log(r1 / r2);
					}
					else
					{
						val = 0.0;
						kn = false;
					}

					fprintf(out,
						"%-2d %-2d %6.3f \t# ((%d / %d) / (%d / %d))%s\n",
						s, d + 3, val,
						total_pairwise[d][s][a1][a2], total,
						num_in_bin[d][s], num_in_all_bins,
						(kn ? "" : " ?"));

					std::cout
						<< Amino(a1).abbr() << " " 
						<< Amino(a2).abbr()
						<< " s " << s
						<< " d " << d
						<< ": r1 " << total_pairwise[d][s][a1][a2]
						<< " / " << total
						<< " r2 "
						<< num_in_bin[d][s]
						<< " / "
						<< num_in_all_bins
						<< " -log(" << r1 << " / " << r2 << ")"
						<< " = "
						<< val
						<< "\n";
				}
				
				if (s != PAIRWISE_SEPS - 1)
				{
					putc('\n', out);
				}
			}
		}
	}
}
*/

namespace
{
	double rapdf_val[NUM_RAPDF_IDS][NUM_RAPDF_IDS][RAPDF_TOP_BIN];
	bool rapdf_val_known[NUM_RAPDF_IDS][NUM_RAPDF_IDS][RAPDF_TOP_BIN];
	long rapdf_num_in_bin[NUM_RAPDF_IDS][NUM_RAPDF_IDS];
}

void generate_rapdf(long total[NUM_RAPDF_IDS][NUM_RAPDF_IDS][RAPDF_TOP_BIN],
	const char *cmd)
{
	static const char *filename = "rapdf.data";
	C_File out(filename, "w", "RAPDF data file");
	std::cerr << "Creating " << filename << "\n";

	fprintf(out, "# RAPDF data file - %s - min reliable %d\n# %s\n",
		get_date().c_str(), RAPDF_RELIABLE_VAL, cmd);

	// fprintf(out, "# total atom pairs within 20 = %ld\n", total_rapdf_pairs);
	fprintf(out, "\n# First bin, top bin\n\n%d %d\n", FIRST_RAPDF_BIN,
		RAPDF_TOP_BIN);

	long num_at_dist[RAPDF_TOP_BIN];
	double total_score[RAPDF_TOP_BIN];
	int total_nonzero[RAPDF_TOP_BIN];
	long num_in_all_bins = 0;

	int b1, b2, d;

	for (b1 = 0;b1 < NUM_RAPDF_IDS;b1++)
	{
		for (b2 = 0;b2 < NUM_RAPDF_IDS;b2++)
		{
			rapdf_num_in_bin[b1][b2] = 0;
		}
	}

	for (d = FIRST_RAPDF_BIN;d < RAPDF_TOP_BIN;d++)
	{
		num_at_dist[d] = 0;
		total_score[d] = 0.0;
		total_nonzero[d] = 0;
	}

// std::cout << "======================\n";

	for (d = FIRST_RAPDF_BIN;d < RAPDF_TOP_BIN;d++)
	{
		for (b1 = 0;b1 < NUM_RAPDF_IDS;b1++)
		{
			for (b2 = 0;b2 < NUM_RAPDF_IDS;b2++)
			{
				long num = total[b1][b2][d];
				rapdf_num_in_bin[b1][b2] += num;
				num_at_dist[d] += num;
				num_in_all_bins += num;
				rapdf_val_known[b1][b2][d] = false;
// std::cout << b1 << ' ' << b2 << ' ' << d << " = " << num << "\n";
			}
		}
	}
// std::cout << ":::::::::::::::::::::::::::\n";

	for (b1 = 0;b1 < NUM_RAPDF_IDS;b1++)
	{
		for (b2 = 0;b2 < NUM_RAPDF_IDS;b2++)
		{
			for (d = FIRST_RAPDF_BIN;d < RAPDF_TOP_BIN;d++)
			{
				double t = total[b1][b2][d];

				if (t != 0)
				{
					/*
					// ratio1 = proportion of [b1][b2] at distance d
					// ratio2 = proportion of all at distance d
					double ratio1 = t / (double) rapdf_num_in_bin[b1][b2];
					double ratio2 = num_at_dist[d] / (double) num_in_all_bins;
					bad */

					/* good:
					*/
					double ratio1 = t / (double) num_at_dist[d];
					double ratio2 =
						rapdf_num_in_bin[b1][b2] / (double) num_in_all_bins;

					assert(ratio2 != 0.0);
					double s = -log(ratio1 / ratio2);

					if (t >= RAPDF_RELIABLE_VAL)
					{
						rapdf_val[b1][b2][d] = s;
						rapdf_val_known[b1][b2][d] = true;
					}

					total_score[d] += s;
					total_nonzero[d]++;
				}
			}
		}
	}

	for (b1 = 0;b1 < NUM_RAPDF_IDS;b1++)
	{
		for (b2 = 0;b2 <= b1;b2++)
		{
			fprintf(out, "\n%s %s  %s %s  # (%ld)\n\n",
				Amino::rapdf_amino(b1).abbr(),
				Amino::rapdf_atom(b1).name(),
				Amino::rapdf_amino(b2).abbr(),
				Amino::rapdf_atom(b2).name(),
				rapdf_num_in_bin[b1][b2]);

			for (d = FIRST_RAPDF_BIN;d < RAPDF_TOP_BIN;d++)
			{
				if (!rapdf_val_known[b1][b2][d])
				{
					// !!!
					if (false && total_nonzero[d] != 0)
					{
						// use the average score, including scores from
						// unreliable frequencies
						rapdf_val[b1][b2][d] =
							total_score[d] / (double) total_nonzero[d];
					}
					else
					{
						rapdf_val[b1][b2][d] = 0.0;
					}

					// find the closest value and use the maximum
					// of that value and the value calculated above

					int dd;
					bool found = false;

					if (d < (RAPDF_TOP_BIN + FIRST_RAPDF_BIN) / 2)
					{
						for (dd = d + 1;dd < RAPDF_TOP_BIN;dd++)
						{
							if (rapdf_val_known[b1][b2][dd])
							{
								found = true;
								break;
							}
						}
					}
					else
					{
						for (dd = d - 1;dd >= FIRST_RAPDF_BIN;dd--)
						{
							if (rapdf_val_known[b1][b2][dd])
							{
								found = true;
								break;
							}
						}
					}

					if (found &&
						rapdf_val[b1][b2][dd] > rapdf_val_known[b1][b2][d])
					{
						rapdf_val_known[b1][b2][d] = rapdf_val[b1][b2][dd];
					}
				}

				/*
				fprintf(out, "%d %.3f%s\n",
					d, rapdf_val[b1][b2][d],
					(rapdf_val_known[b1][b2][d] ? "" :
						(total[b1][b2][d] == 0 ? " \t# ?" : " \t# low")));
				*/

				fprintf(out, "%d %.3f \t# ((%ld / %ld) / (%ld / %ld))%s\n",
					d, rapdf_val[b1][b2][d],
					
					total[b1][b2][d],
					num_at_dist[d],
					rapdf_num_in_bin[b1][b2],
					num_in_all_bins,

					(rapdf_val_known[b1][b2][d] ? "" : " ?"));
			}
		}
	}
}

void count_amino_acids(const Peptide &p, long total_aa[])
{
	for (int n = 0;n < p.length();n++)
	{
		int a = p.res(n).amino().num();
		assert(a >= 0 && a < 20);
		total_aa[a]++;
	}
}

void get_command(int argc, const char *argv[], char cmd[1000])
{
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
}

namespace
{

// some of these arrrays are very large; caused a stack overrun when they were
// local variables inside main(), so they have been made static instead

long total_solv[SOLV_BIN_TOP][Amino::Num];
long total_torsion[TORSION_BINS][TORSION_BINS][Amino::Num];
long total_orient[ORIENT_DISTS][ORIENT_ANGLES][Amino::Num][Amino::Num];
//long total_pairwise[PAIRWISE_DISTS][PAIRWISE_SEPS][Amino::Num][Amino::Num];
long total_rapdf[NUM_RAPDF_IDS][NUM_RAPDF_IDS][RAPDF_TOP_BIN];

int num_torsion_aminos[TORSION_BINS][TORSION_BINS];
}

int main(int argc, const char *argv[])
{
	if (argc != 3)
	{
		std::cerr << "Usage: " << argv[0] << " <dir> <PDB chain list file>\n"
			<< "\nChain list file contains lines of the form:\n\n"
			<< "<PDB id> <chain>\n\n";
		exit(1);
	}

	const char *dir = argv[1];
	C_File file(argv[2], "r", "Chain list file");
	char s[80], pdb_id[80];
	char pdb_filename[1000];
	char chain_id;
	int line_num = 0;
	long total_aa[Amino::Num];
	long total_res = 0;

	for (int a = 0;a < Amino::Num;a++)
	{
		total_aa[a] = 0;

		for (int solv_bin = 0;solv_bin < SOLV_BIN_TOP;solv_bin++)
		{
			total_solv[solv_bin][a] = 0;
		}

		for (int phi_bin = 0;phi_bin < TORSION_BINS;phi_bin++)
		{
			for (int psi_bin = 0;psi_bin < TORSION_BINS;psi_bin++)
			{
				total_torsion[phi_bin][psi_bin][a] = 0;
			}
		}

		for (int b = 0;b < Amino::Num;b++)
		{
			for (int od = 0;od < ORIENT_DISTS;od++)
			{
				for (int oa = 0;oa < ORIENT_ANGLES;oa++)
				{
					total_orient[od][oa][a][b] = 0;
				}
			}

			/*
			for (int pd = 0;pd < PAIRWISE_DISTS;pd++)
			{
				for (int pn = 0;pn < PAIRWISE_SEPS;pn++)
				{
					total_pairwise[pd][pn][a][b] = 0;
				}
			}
			*/
		}
	}

	for (int phi_bin = 0;phi_bin < TORSION_BINS;phi_bin++)
	{
		for (int psi_bin = 0;psi_bin < TORSION_BINS;psi_bin++)
		{
			num_torsion_aminos[phi_bin][psi_bin] = 0;
		}
	}

	for (int b1 = 0;b1 < NUM_RAPDF_IDS;b1++)
	{
		for (int b2 = 0;b2 < NUM_RAPDF_IDS;b2++)
		{
			for (int d = FIRST_RAPDF_BIN;d < RAPDF_TOP_BIN;d++)
			{
				total_rapdf[b1][b2][d] = 0;
			}
		}
	}

	std::cerr << "Started " << get_date() << "\n";

	while (fgets(s, 80, file) != NULL)
	{
		line_num++;

		chain_id = ' ';
		if (sscanf(s, "%s %c", pdb_id, &chain_id) < 1)
		{
			std::cerr << "Error on line " << line_num
				<< " of " << file.name()
				<< "\n";
			exit(1);
		}

		sprintf(pdb_filename, "%s/%s_%c.pdb", dir, pdb_id, chain_id);

		if (!file_exists(pdb_filename))
		{
			sprintf(pdb_filename, "%s/%s.pdb", dir, pdb_id);
		}

		Peptide p;
		
		if (!p.read_pdb(pdb_filename, chain_id,
			true	// no warnings
			))
		{
			std::cerr << "Errors found in PDB file " << pdb_filename
				<< "\n";
			exit(1);
		}

		std::cerr << "Read " << pdb_id << " chain " << chain_id
			<< ", length = " << p.length() << "\n";
		std::cout.flush();

		total_res += p.length();
		count_amino_acids(p, total_aa);
		p.conf().calc_torsion_angles();

		calc_torsion(p, total_torsion, num_torsion_aminos,
			pdb_filename, chain_id);
		calc_solvation(p, total_solv, pdb_filename, chain_id);
		calc_orient(p, total_orient);
		calc_rapdf(p, total_rapdf);

		//////// calc_pairwise(p, total_pairwise);
	}

	char cmd[1000];
	get_command(argc, argv, cmd);

	generate_torsion(total_aa, total_torsion, num_torsion_aminos,
		total_res, cmd);
	generate_solvation(total_aa, total_solv, total_res, cmd);
	generate_orient(total_orient, cmd);
	generate_rapdf(total_rapdf, cmd);

	//////// generate_pairwise(total_pairwise, cmd);

	std::cerr << "Finished " << get_date() << "\n";

	return 0;
}

