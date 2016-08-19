
#include <iostream>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>

// scoring terms

enum Term
{
	Solv, Orient, LJ, RAPDF, HBond, Torsion,
	Num_Terms
};

typedef std::vector<double> DVec;
DVec base_score[Num_Terms];
DVec gdt_score;
DVec w_score;
double gdt_total = 0.0;
double gdt_mean = 0.0;
double gdt_msq_sum = 0.0;

inline double square(double x)
{
	return x * x;
}

double weighted_correlation(double s, double o, double lj, double r,
	double h, double t, int num)
{
	double total_score = 0.0;

	int n;
	for (n = 0;n < num;n++)
	{
		total_score += (w_score[n] =
			s * base_score[Solv][n] +
			o * base_score[Orient][n] +
			lj * base_score[LJ][n] +
			r * base_score[RAPDF][n] +
			h * base_score[HBond][n] +
			t * base_score[Torsion][n]
			);
	}

	double mean_score = total_score / (double) num;
	double both_sum = 0.0;
	double score_msq_sum = 0.0;

	for (n = 0;n < num;n++)
	{
		double w = w_score[n] - mean_score;
		both_sum += w * (gdt_score[n] - gdt_mean);
		score_msq_sum += w * w;
	}

	return (both_sum / sqrt(score_msq_sum * gdt_msq_sum));
}

std::string weights_to_str(double best_weight[Num_Terms])
{
	char buffer[1000];

	sprintf(buffer, "%f,%f,%f,%f,%f,%f",
		best_weight[Solv],
		best_weight[Orient],
		best_weight[LJ],
		best_weight[RAPDF],
		best_weight[HBond],
		best_weight[Torsion]);
	return std::string(buffer);
}

typedef std::map<std::string,double> SDMap;
SDMap known_corr;

double downhill_search(double best_weight[Num_Terms], double best_c)
{
	double dist = 0.1;
	int dist_level = 0;
	double weight[Num_Terms];

	known_corr[weights_to_str(best_weight)] = best_c;
	int num = (int) gdt_score.size();
	int count = 0;

	// try all combinations of weights +/-dist or unchanged

	for ( ; ; )
	{
		int best_ds, best_dor, best_dlj, best_dr, best_dp, best_dh, best_dt;
		best_ds = best_dor = best_dlj = best_dr = best_dp = best_dh =
			best_dt = 0;
		double new_best_c = best_c;

		int dds = -1;
		if (best_weight[Solv] + dds * dist < 0.0) { dds = 0; }

		for (int ds = dds;ds <= 1;ds++)
		{
			weight[Solv] = best_weight[Solv] + ds * dist;

			int ddor = -1;
			if (best_weight[Orient] + ddor * dist < 0.0) { ddor = 0; }

			for (int dor = ddor;dor <= 1;dor++)
			{
				weight[Orient] = best_weight[Orient] + dor * dist;

				int ddlj = -1;
				if (best_weight[LJ] + ddlj * dist < 0.0) { ddlj = 0; }

				for (int dlj = ddlj;dlj <= 0;dlj++)
				//for (int dlj = 0;dlj <= 0;dlj++)
				{
					weight[LJ] = best_weight[LJ] + dlj * dist;

					int ddr = -1;
					if (best_weight[RAPDF] + ddr * dist < 0.0) { ddr = 0; }

					for (int dr = ddr;dr <= 1;dr++)
					{
						weight[RAPDF] = best_weight[RAPDF] + dr * dist;

						int ddh = -1;
						if (best_weight[HBond] + ddh * dist < 0.0) { ddh = 0; }

						for (int dh = ddh;dh <= 1;dh++)
						{
							weight[HBond] = best_weight[HBond] + dh * dist;

							int ddt = -1;
							if (best_weight[Torsion] + ddt * dist < 0.0) { ddt = 0; }
							//for (int dt = ddt;dt <= 1;dt++)
							for (int dt = 0;dt <= 0;dt++)
							{
							weight[Torsion] = best_weight[Torsion] + dt * dist;
//							if (weight[Torsion] < 0.0) { weight[Torsion] = 0.0; }

							std::string s = weights_to_str(weight);
							SDMap::iterator i = known_corr.find(s);

							if (i != known_corr.end())
							{
//std::cout << dist << " Old: " << s << " = " << i->second << "\n";
								// already seen this one
								continue;
							}

							double c = weighted_correlation(
								weight[Solv], weight[Orient], 
								weight[LJ], weight[RAPDF],
								weight[HBond], weight[Torsion],
								num);

							known_corr[s] = c;

//std::cout << dist << " New: " << s << " = " << c << "\n";

							if (fabs(c) > fabs(new_best_c))
							{
								std::cerr << c << "\n";

								new_best_c = c;
								best_ds = ds;
								best_dlj = dlj;
								best_dr = dr;
								best_dor = dor;
								best_dh = dh;
								best_dt = dt;
							}
							}
						}
					}
				}
			}
		}
		
		if (best_ds == 0 && best_dor == 0 && best_dlj == 0 &&
			best_dr == 0 && best_dp == 0 && best_dh == 0 && best_dt == 0)
		{
			if (++dist_level == 3)
			{
				return best_c;
			}

			dist *= 0.1;
			count = 0;
		}
		else
		{
			best_c = new_best_c;
			best_weight[Solv]    += best_ds * dist;
			best_weight[Orient]  += best_dor * dist;
			best_weight[LJ] += best_dlj * dist;
			best_weight[RAPDF]   += best_dr * dist;
			best_weight[HBond]   += best_dh * dist;
			best_weight[Torsion] += best_dt * dist;

			/*
			if (best_weight[Solv] < 0.0) { best_weight[Solv] = 0.0; }
			if (best_weight[Orient] < 0.0) { best_weight[Orient] = 0.0; }
			if (best_weight[LJ] < 0.0) { best_weight[LJ] = 0.0; }
			if (best_weight[RAPDF] < 0.0) { best_weight[RAPDF] = 0.0; }
			if (best_weight[HBond] < 0.0) { best_weight[HBond] = 0.0; }
			if (best_weight[Torsion] < 0.0) { best_weight[Torsion] = 0.0; }
			*/

			count++;
			std::cerr << dist << " #" << count
				<< " (" << best_c << ")"
				<< " s " << best_weight[Solv]
				<< " o " << best_weight[Orient]
				<< " l " << best_weight[LJ]
				<< " r " << best_weight[RAPDF]
				<< " h " << best_weight[HBond]
				<< " t " << best_weight[Torsion]
				<< "\n";

		}
	}
}

void normalise_weights(double w[])
{
	double t = 0.0;

	int n;
	for (n = 0;n < Num_Terms;n++)
	{
		t += w[n];
	}

	for (n = 0;n < Num_Terms;n++)
	{
		w[n] /= t;
	}
}

void find_weights()
{
	//static const double inc = 0.1;
	//static const double maxval = 1.20001;
	//static const int num_inc = (int) (maxval / inc);
	//static const int num_inc_cubed = num_inc * num_inc * num_inc;
	double best_weight[Num_Terms];
	double best_c = 0.0;

	int num = (int) gdt_score.size();
	w_score.resize(num);

	int g;
	for (g = 0;g < num;g++)
	{
		gdt_total += gdt_score[g];
	}

	gdt_mean = gdt_total / (double) num;

	for (g = 0;g < num;g++)
	{
		gdt_msq_sum += square(gdt_score[g] - gdt_mean);
	}

/*
	int count = 0;

	for (double s = inc;s < maxval;s += inc)
	{
		for (double o = inc;o < maxval;o += inc)
		{
			for (double t = inc;t < maxval;t += inc)
			{
				count++;
				std::cerr << ((double) count / (double) (num_inc_cubed)) * 100
					<< " %\n";

				for (double r = inc;r < maxval;r += inc)
				{
						for (double h = inc;h < maxval;h += inc)
						{
							double c = weighted_correlation(s, o, t, r, p, h,
								num);

							if (fabs(c) > fabs(best_c))
							{
								best_c = c;
								best_weight[Solv]    = s;
								best_weight[Orient]  = o;
								best_weight[LJ] = t;
								best_weight[RAPDF]   = r;
								best_weight[HBond]   = h;

								std::cerr << best_c << "\n";
							}
						}
				}
			}
		}
	}
*/

	best_weight[Solv]    = 25.0;
	best_weight[LJ]      = 25.0;
	best_weight[RAPDF]   = 25.0;
	best_weight[Orient]  = 5.0;
	best_weight[HBond]   = 5.0;
	//best_weight[Torsion] = 15.0;
	best_weight[Torsion] = 0.0;

	/*
	best_weight[Solv]    = 20.0;
	best_weight[Orient]  = 5.0;
	best_weight[LJ]		 = 45.0;
	best_weight[RAPDF]   = 15.0;
	best_weight[HBond]   = 0.0;
	best_weight[Torsion] = 10.0;
	*/

	/*
	best_weight[Solv]    = drand48() * 0.8 + 1.0;
	best_weight[Orient]  = drand48() * 0.8 + 1.0;
	best_weight[LJ]		 = drand48() * 0.8 + 1.0;
	best_weight[RAPDF]   = drand48() * 0.8 + 1.0;
	best_weight[HBond]   = drand48() * 0.8 + 1.0;
	best_weight[Torsion] = drand48() * 0.8 + 1.0;
	*/

	best_c = weighted_correlation(
		best_weight[Solv], best_weight[Orient], 
		best_weight[LJ], best_weight[RAPDF],
		best_weight[HBond], best_weight[Torsion],
		num);

std::cout << "INITIAL:\n";
std::cout << "\n"
	"weight_rapdf       = " << best_weight[RAPDF] << "\n"
	"weight_solvation   = " << best_weight[Solv] << "\n"
	"weight_lj          = " << best_weight[LJ] << "\n"
	"weight_hbond       = " << best_weight[HBond] << "\n"
	"weight_orientation = " << best_weight[Orient] << "\n"
	"weight_torsion     = " << best_weight[Torsion] << "\n"
	"# correlation      = " << best_c << "\n";
std::cout << "\n\n";

	best_c = downhill_search(best_weight, best_c);
	normalise_weights(best_weight);

	std::cout << "\n"
		"weight_rapdf       = " << best_weight[RAPDF] << "\n"
		"weight_solvation   = " << best_weight[Solv] << "\n"
		"weight_lj          = " << best_weight[LJ] << "\n"
		"weight_hbond       = " << best_weight[HBond] << "\n"
		"weight_orientation = " << best_weight[Orient] << "\n"
		"weight_torsion     = " << best_weight[Torsion] << "\n"
		"# correlation      = " << best_c << "\n";
}

void read_data(const char *filename)
{
	static const int Max_Length = 1000;
	char buffer[Max_Length];
	FILE *file;

	if ((file = fopen(filename, "r")) == NULL)
	{
		std::cerr << "Cannot open " << filename << "\n";
		exit(1);
	}

	for (int n = 0;fgets(buffer, Max_Length, file) != NULL;n++)
	{
		for (int i = 0;i < Num_Terms;i++)
		{
			base_score[i].resize(n + 1);
		}

		gdt_score.resize(n + 1);

		int c = sscanf(buffer, "%*s %lf %*d %lf %lf %lf %lf %lf %lf",
			&gdt_score[n],
			&base_score[Solv][n],
			&base_score[Orient][n],
			&base_score[LJ][n],
			&base_score[RAPDF][n],
			&base_score[HBond][n],
			&base_score[Torsion][n]);

		if (c != 7)
		{
			std::cerr << "Error on line " << n + 1
				<< " of " << filename
				<< ": scanf returned " << c
				<< "\n";
			fclose(file);
			exit(1);
		}
	}

	fclose(file);
}

void show_usage(const char *argv0)
{
	std::cout << "Usage: " << argv0
		<< " <modsubname_gdt_len_s_o_l_r_h_t>\n\n"
		"Fields are:\n\n"
		"1) Decoy name (ignored)\n"
		"2) GDT_TS of decoy and native\n"
		"3) Protein length (ignored)\n"
		"4) Solvation score (decoy - native)\n"
		"5) Orientation score (decoy - native)\n"
		"6) LJ score (decoy - native)\n"
		"7) RAPDF score (decoy - native)\n"
		"8) Hydrogen bond score (decoy - native)\n"
		"9) Torsion score (decoy - native)\n";
}

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		show_usage(argv[0]);
		exit(1);
	}

	srand(time(0));
	read_data(argv[1]);
	find_weights();

	return 0;
}

