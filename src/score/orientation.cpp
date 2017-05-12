
#include "orientation.h"
#include "orientation_impl.h"
#include "peptide.h"
#include "geom.h"
#include "scorer_combined.h"

Orientation::Orientation()
{
	m_short = new Orientation_impl;
	m_long = new Orientation_impl;
}

Orientation::~Orientation()
{
	delete m_short;
	delete m_long;
}

void Orientation::set_short_data_file(const std::string &filename)
{
	m_short->set_data_file(filename);
}

void Orientation::set_long_data_file(const std::string &filename)
{
	m_long->set_data_file(filename);
}

double Orientation::score(const Peptide &p, bool verbose, bool continuous)
{
	if (p.length() <= SHORT_PEPTIDE || continuous)
	{
		return m_short->score(p, verbose, continuous);
	}
	else
	{
		return m_long->score(p, verbose);
	}
}

bool Orientation::get_bins(const Peptide &p, int n, int m,
	int *dist_bin, int *angle_bin)
{
	//static const double MaxOrientDist = ORIENT_DISTS + 2.0;

	// static const double FortyFiveDeg = deg2rad(45.0);
	static const double NinetyDeg = deg2rad(90.0);

	int r1, r2;

	// set (r1 and r2) to (n and m) such that
	// r1's amino num < r2's amino num

	if (p.res(n).amino().num() <
		p.res(m).amino().num())
	{
		r1 = n;
		r2 = m;
	}
	else
	{
		r1 = m;
		r2 = n;
	}

	if (!(p.atom_exists(r1, Atom_CA) && p.atom_exists(r2, Atom_CA)))
	{
		return false;
	}

	Point ca1 = p.atom_pos(r1, Atom_CA);
	Point ca2 = p.atom_pos(r2, Atom_CA);

	// if (ca1.further_than(10.0 - 0.001, ca2))

	/*
	if (ca1.further_than(MaxOrientDist, ca2))
	{
		return false;
	}
	*/

	double d = ca1.distance(ca2);

	/*
	// check for edge case (exactly at MaxOrientDist)
	if (d >= MaxOrientDist)
	{
		return false;
	}
	*/

	// distance 0-3 = bin 0, distance 3-4 = bin 1 ...
	// distance (MaxOrientDist - 0.001) = bin (ORIENT_DISTS - 1)

	*dist_bin = (d < 3.0 ? 0 : (int) (d - 2.0));
	//assert(*dist_bin < ORIENT_DISTS);

	if (*dist_bin >= ORIENT_DISTS)
	{
		*dist_bin = ORIENT_DISTS - 1;
	}

	Point s1, s2;

	if (!(p.get_side_chain_pos(r1, &s1) &&
		  p.get_side_chain_pos(r2, &s2) &&
		  p.atom_exists(r1, Atom_C) &&
		  p.atom_exists(r2, Atom_C)))
	{
		return false;
	}

	bool a1_lt_90, a1_lt_45, a2_lt_90, a2_lt_45;

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

	// a1 = angle_formed(s1, ca1, ca2);
	// a2 = angle_formed(s2, ca2, ca1);

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
		if (a2_lt_45)	   { *angle_bin = 0; }
		else if (a2_lt_90) { *angle_bin = 1; }
		else			   { *angle_bin = 2; }
	}
	else
	if (a2_lt_45)
	{
		if (a1_lt_90)	   { *angle_bin = 3; }
		else			   { *angle_bin = 4; }
	}
	else
	if (!a1_lt_90)
	{
		if (!a2_lt_90)
		{
			*angle_bin = 5;
		}
		else
		{
			*angle_bin = 6;
		}
	}
	else
	if (!a2_lt_90)
	{
		*angle_bin = 7;
	}
	else
	{
		double t;

		if (d < 0.1)
		{
			t = 0.0;
		}
		else
		{
			t = fabs(torsion_angle(s1, ca1, ca2, s2));

			// std::cout << "! t " << rad2deg(t) << "\n";
		}

		if (t < NinetyDeg)
		{
			*angle_bin = 8;
		}
		else
		{
			*angle_bin = 9;
		}
	}

	assert(*angle_bin >= 0 && *angle_bin < ORIENT_ANGLES);
	return true;
}

