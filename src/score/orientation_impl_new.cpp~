
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>

#include "peptide.h"
#include "residue.h"
#include "atom.h"
#include "amino.h"
#include "geom.h"
#include "c_file.h"
#include "scorer_combined.h"
#include "orientation.h"
#include "orientation_impl.h"

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678119
#endif

Orientation_impl::Orientation_impl() :
	m_data_loaded(false)
{
}

Orientation_impl::~Orientation_impl()
{
}

void Orientation_impl::set_data_file(const std::string &filename)
{
	m_filename = filename;
}

void Orientation_impl::load_data()
{
	if (m_data_loaded)
	{
		return;
	}

	if (m_filename.empty())
	{
		std::cerr << "Error: Orientation data file undefined\n";
		exit(1);
	}

//std::cout << "LOADING ORIENTATION DATA\n";

	// exits with error message if the file does not exist
	C_File file(m_filename, "r", "Orientation data file");

	static const int Max_Len = 1000;
	char buffer[Max_Len];
	double val;

	for (int a1 = 0;a1 < Amino::Num;a1++)
	{
		for (int a2 = a1;a2 < Amino::Num;a2++)
		{
			if (!file.next_line(buffer, Max_Len) ||
				std::string(buffer, 3) != Amino(a1).abbr() ||
				std::string(buffer + 4, 3) != Amino(a2).abbr())
			{
				std::cerr << "Error: expected \""
					<< Amino(a1).abbr() << ' ' << Amino(a2).abbr()
					<< "\" on line "
					<< file.line_num()
					<< " of orientation data file "
					<< m_filename
					<< "\n";
				exit(1);
			}

			for (int angle = 0;angle < ORIENT_ANGLES;angle++)
			{
				for (int dist = 0;dist < ORIENT_DISTS;dist++)
				{
					int d;

					if (!file.next_line(buffer, Max_Len) ||
						sscanf(buffer, "%d %lf", &d, &val) != 2 ||
						d != dist + 3)
					{
						std::cerr << "Error: expected \"" << dist + 3
							<< "\" followed by value on line "
							<< file.line_num()
							<< " of orientation data file "
							<< m_filename
							<< "\n";
						exit(1);
					}

/*
std::cout << ": "
<< dist << " "
<< angle << " "
<< a1 << " " << a2 << " "
<< " = " << val
<< "\n";
*/
					m_data[dist][angle][a1][a2] = val;
					m_data[dist][angle][a2][a1] = val;
				}
			}
		}
	}

	if (file.next_line(buffer, Max_Len))
	{
		std::cerr << "Warning: ignoring extra data on line "
			<< file.line_num()
			<< " of orientation data file "
			<< m_filename
			<< "\n";
	}

	m_data_loaded = true;
}

double Orientation_impl::score(const Peptide& p, bool verbose)
{
	// (does nothing if already loaded)
	load_data();
	bool a1_lt_90, a1_lt_45, a2_lt_90, a2_lt_45;
	int dist_bin, angle_bin, num,a1,a2;
	double total = 0.0, t,d,t2;
	static const double NinetyDeg = deg2rad(90.0);
	Point ca1,ca2,s1,s2;
	
	for (int n = p.start() + 2;n <= p.end();n++)
	{
		a1 = p.res(n).amino().num();
		t = 0.0;		num = 0;

		for (int m = p.start();m < n - 1;m++)
		{
			a2 = p.res(m).amino().num();
			
			if (!(p.atom_exists(m, Atom_CA) && p.atom_exists(n, Atom_CA))) continue;

			ca1 = p.atom_pos2(m, Atom_CA);
			ca2 = p.atom_pos2(n, Atom_CA);
			d=sqrt(pow((ca1.x-ca2.x),2)+pow((ca1.y-ca2.y),2)+pow((ca1.z-ca2.z),2));
	
			dist_bin = (d < 3.0 ? 0 : (int) (d - 2.0));


			if (dist_bin >= ORIENT_DISTS) dist_bin = ORIENT_DISTS - 1;
	
			if (!(p.get_side_chain_pos(m, &s1) && p.get_side_chain_pos(n, &s2) &&
				  p.atom_exists(m, Atom_C) && p.atom_exists(n, Atom_C)))
				continue;

			if (d < 0.1)
			{
				a1_lt_90 = a1_lt_45 = a2_lt_90 = a2_lt_45 = true;
			}
			else
			{
				Point ca1_s1 = (s1.minus(ca1)).normalised();
				Point ca2_s2 = (s2.minus(ca2)).normalised();
				Point ca1_ca2 = (ca2.minus(ca1)).normalised();
				Point ca2_ca1 = ca1_ca2.negated();

				double dp1 = ca1_ca2.dot_product(ca1_s1);
				double dp2 = ca2_ca1.dot_product(ca2_s2);

				a1_lt_90 = (dp1 > 0);
				a1_lt_45 = (dp1 > M_SQRT1_2);
				a2_lt_90 = (dp2 > 0);
				a2_lt_45 = (dp2 > M_SQRT1_2);
			}		


			// summary of angle values:
			//
			// val		a1		a2		torsion
			//
			//   0		<45		<45
			//   1		<45		45-90
			//   2		<45		90+
			//   3		45-90	<45
			//   4		90+		<45
			//   5		90+		90+
			//   6		90+		45-90
			//   7		45-90	90+	
			//   8		45-90	45-90	<90
			//   9		45-90	45-90	90+
		
			if (a1_lt_45)
			{
				if (a2_lt_45)	   { angle_bin = 0; }
				else if (a2_lt_90) { angle_bin = 1; }
				else			   { angle_bin = 2; }
			}
			else
			{ 
				if (a2_lt_45)
				{
					if (a1_lt_90)	   { angle_bin = 3; }
					else			   { angle_bin = 4; }
				}
				else
				{
					if (!a1_lt_90)
					{
						if (!a2_lt_90) angle_bin = 5;
						else 		   angle_bin = 6;
					}
					else
					{
						if (!a2_lt_90) angle_bin = 7;
						else
						{
							if (d < 0.1) t2 = 0.0;
							else t2 = fabs(torsion_angle(s1, ca1, ca2, s2));

							if (t2 < NinetyDeg)  angle_bin = 8;
							else 				 angle_bin = 9;
						}
					}
				}
			}

			assert(angle_bin >= 0 && angle_bin < ORIENT_ANGLES);
			t += m_data[dist_bin][angle_bin][a1][a2];
			num++;
		}
		total += t;
	}

#ifndef RAW_SCORE
	int len = p.length();
	total /= (double) (len + 20);
	// normalise so that all score types have approximately the same range
	total *= 42.0;
#endif // RAW_SCORE

	return total;
}

void Orientation_impl::dump(std::ostream &out /*=std::cout*/)
{
}

