#include <cstdlib> // PG added this
#include <iostream>
#include <cmath>
#include "matrix.h"
#include "point.h"
#include "geom.h"

// (M_PI and M_TWOPI should be defined in math.h, but just in case)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_TWOPI
#define M_TWOPI (M_PI * 2.0)
#endif

double rad2deg(double radians)
{
	return range_0_360(radians * (180.0 / M_PI));
}

double deg2rad(double degrees)
{
	return range_0_2pi(degrees * (M_PI / 180.0));
}

double range_0_360(double val)
{
	// (avoid calls to fmod() if possible)

	if (val < 0.0)
	{
		if (val < 360.0)
		{
			val = 360.0 - fmod(-val, 360.0);
			return (val == 360.0 ? 0.0 : val);
		}
		else
		{
			return val + 360.0;
		}
	}
	else
	if (val >= 360.0)
	{
		if (val >= 720.0)
		{
			return fmod(val, 360.0);
		}
		else
		{
			return val - 360.0;
		}
	}
	else
	{
		return val;
	}
}

double range_m180_180(double val)
{
	val = range_0_360(val);
	return (val >= 180.0 ? val - 360.0 : val);
}

double range_0_2pi(double val)
{
	// (avoid calls to fmod() if possible)

	if (val < 0.0)
	{
		if (val < M_TWOPI)
		{
			val = M_TWOPI - fmod(-val, M_TWOPI);
			return (val == M_TWOPI ? 0.0 : val);
		}
		else
		{
			return val + M_TWOPI;
		}
	}
	else
	if (val >= M_TWOPI)
	{
		if (val >= M_TWOPI * 2.0)
		{
			return fmod(val, M_TWOPI);
		}
		else
		{
			return val - M_TWOPI;
		}
	}
	else
	{
		return val;
	}
}

double range_mpi_pi(double val)
{
	val = range_0_2pi(val);
	return (val >= M_PI ? val - M_PI * 2.0 : val);
}

double torsion_angle(const Point& p1, const Point& p2,
	const Point& p3, const Point& p4)
{
   double xij = p1.x - p2.x;
   double yij = p1.y - p2.y;
   double zij = p1.z - p2.z;

   double xkj = p3.x - p2.x;
   double ykj = p3.y - p2.y;
   double zkj = p3.z - p2.z;

   double xkl = p3.x - p4.x;
   double ykl = p3.y - p4.y;
   double zkl = p3.z - p4.z;

   // normal to plane 1

   double dxi = yij * zkj - zij * ykj;
   double dyi = zij * xkj - xij * zkj;
   double dzi = xij * ykj - yij * xkj;

   // normal to plane 2

   double gxi = zkj * ykl - ykj * zkl;
   double gyi = xkj * zkl - zkj * xkl;
   double gzi = ykj * xkl - xkj * ykl;

   // lengths of normals

   double bi = sqrt(dxi * dxi + dyi * dyi + dzi * dzi);
   double bk = sqrt(gxi * gxi + gyi * gyi + gzi * gzi);
   
   double ct = (dxi * gxi) + (dyi * gyi) + (dzi * gzi);
   
   double z1 = 1.0 / bi;
   double z2 = 1.0 / bk;

   ct = ct * z1 * z2;

   if (ct > 1.0)
   {
      ct = 1.0;
   }

   if (ct < (-1.0))
   {
      ct = -1.0;
   }
  
   double ap = acos(ct);
   
   double s = (xkj * (dzi * gyi - dyi * gzi) +
               ykj * (dxi * gzi - dzi * gxi) +
               zkj * (dyi * gxi - dxi * gyi));

   if (s < 0.0)
   {
      ap = -ap;
   }
   
   ap = (ap > 0.0) ? M_PI - ap : - (M_PI + ap);
   
   return ap;
}

double angle_formed(const Point& p1, const Point& p2, const Point& p3)
{
	static const double TinyDist = 0.001;

	Point vec_1 = p1.minus(p2);
	Point vec_2 = p3.minus(p2);

	double len_1 = vec_1.length();
	double len_2 = vec_2.length();

	if (len_1 < TinyDist || len_2 < TinyDist)
	{
		return 0.0;
	}

	double d = vec_1.dot_product(vec_2);
	double x = (d / (len_1 * len_2));

	// std::cout << " { dot = " << d << ", lengths = " << len_1
	//	<< " & " << len_2 << ", x = " << x << " } ";

	// value should never be < -1 or > 1, but check for floating
	// point inaccuracies

	if (x < -1.0) { x = -1.0; std::cout << "CORRECTING angle_formed: < -1\n"; }
	else if (x > 1.0) { x = 1.0; std::cout << "CORRECTING angle_formed: > 1\n"; }

	return acos(x);
}

// This function depends on the fact that the four atoms bonded to a CA atom
// (C, N, CB and H atoms) form a tetrahedral shape (ie. all bond angles are
// 109.5 degrees).

Point estimate_CB_pos(const Point &ca_pos, const Point &n_pos,
	const Point &c_pos)
{
	// (see code below for meaning of these constants)
	static const double PlaneOffset = 0.28488;
	static const double DistMidToCB = 2.1496;

	static const double LenFactor = BOND_LENGTH_C_C / BOND_LENGTH_N_CA;

	// vector from CA to C
	Point ca_to_c = c_pos.minus(ca_pos);

	// vector from CA to N, adjusted to average Carbon-Carbon bond length
	Point ca_to_n = n_pos.minus(ca_pos);
	ca_to_n.scale(LenFactor);

	Point adjusted_n = ca_pos.plus(ca_to_n);

	if (ca_to_n.is_tiny() || ca_to_c.is_tiny())
	{
		std::cerr << "Error: illegal atom positions in estimate_CB_pos(): "
			<< "CA = " << ca_pos
			<< ", N = " << n_pos
			<< ", C = " << c_pos
			<< "\n";
		exit(1);
	}

	// midpoint of CA and adjusted N position

	Point mid(
		(c_pos.x + adjusted_n.x) * 0.5,
		(c_pos.y + adjusted_n.y) * 0.5,
		(c_pos.z + adjusted_n.z) * 0.5);

	// find point p in the plane containing { C, adjusted N, predicted CB }
	// such that the vector from mid to p points towards CB

	Point p = ca_to_n.cross_product(ca_to_c);
	p.scale(PlaneOffset);
	p.add(ca_pos);

	// move p to the predicted CB position
	// (DistMidToDB is the expected distance from mid to CB)

	p.subtract(mid);
	p.normalise();
	p.scale(DistMidToCB);
	p.add(mid);

	// Could also make sure the distance from CA to p is exactly 1.52,
	// but it makes almost no difference - testing with real protein CB
	// positions vs the value calcuated by this function showed that the
	// improvement was at most 0.1 Angstroms; 99% of the time less than 0.04
	//
	// The code would be:
	//
	// p.subtract(ca_pos);
	// p.normalise();
	// p.scale(BOND_LENGTH_C_C);
	// p.add(ca_pos);

	return p;
}

Point estimate_O_pos(const Point &ca_pos, const Point &c_pos, const Point &next_n_pos)
{
	// the C, CA, N and O atoms are all in a plane (more or less),
	// so just use vector arithmetic to calculate the position
	//
	// Note that the N-C-O bond angle is slightly larger than the CA-C-O bond angle
	// (about 125 vs 121 degrees), but this doesn't make a significant difference
	// to the estimated position (less than 0.05 Angstroms)

	Point v1 = c_pos.minus(next_n_pos);
	Point v2 = c_pos.minus(ca_pos);

	v1.normalise();
	v2.normalise();

	Point p = v1.plus(v2);

	p.scale(0.953 * BOND_LENGTH_C_O);

	// the line above does the same thing as:
	//
	// p.normalise();
	// p.scale(BOND_LENGTH_C_O);

	return c_pos.plus(p);
}

Point torsion_to_coord(const Point &p1, const Point &p2, const Point &p3,
	double bond_length, double bond_angle, double torsion_angle,
	double prev_bond_length)
{
	// (This function is based on "Practical Conversion from Torsion Space to
	// Cartesian Space for In Silico Protein Synthesis", Parsons et al, 2005).
	//
	// p1, p2 and p3 correspond to points A, B and C in the paper.
	//
	// Note that the x coordinate of D2 in the top right of p. 1066 should
	// be R cos(-theta), not R cos(theta).

	double cos_bond = cos(bond_angle);
	double sin_bond = sin(bond_angle);
	double cos_tor = cos(torsion_angle);
	double sin_tor = sin(torsion_angle);

	Point p(
		bond_length * -cos_bond,
		bond_length * sin_bond * cos_tor,
		bond_length * sin_bond * sin_tor
	);

	// find the columns of the rotation matrix (mx, my and mz)

	Point vec_1_2 = p2.minus(p1);
	Point mx = p3.minus(p2);
	Point mz = vec_1_2.cross_product(mx);
	mx.scale(1.0 / prev_bond_length);		// normalise to length 1.0
	mz.normalise();
	Point my = mz.cross_product(mx);

	Matrix_3_3 m(mx, my, mz);
	return m.times(p).plus(p3);
}

Point torsion_to_ideal(const Point &p1, const Point &p2, const Point &p3,
	double torsion_angle, bool is_ca, bool is_n)
{
	static double x_CA = BOND_LENGTH_N_CA * -cos(BOND_ANGLE_C_N_CA);
	static double yz_CA = BOND_LENGTH_N_CA * sin(BOND_ANGLE_C_N_CA);
	static double mx_scale_CA = 1.0 / BOND_LENGTH_C_N;

	static double x_C = BOND_LENGTH_C_C * -cos(BOND_ANGLE_N_CA_C);
	static double yz_C = BOND_LENGTH_C_C * sin(BOND_ANGLE_N_CA_C);
	static double mx_scale_C = 1.0 / BOND_LENGTH_N_CA;

	static double x_N = BOND_LENGTH_C_N * -cos(BOND_ANGLE_CA_C_N);
	static double yz_N = BOND_LENGTH_C_N * sin(BOND_ANGLE_CA_C_N);
	static double mx_scale_N = 1.0 / BOND_LENGTH_C_C;

	double cos_tor = cos(torsion_angle);
	double sin_tor = sin(torsion_angle);

	double x, yz;
	double mx_scale;

	if (is_ca)
	{
		x = x_CA;
		yz = yz_CA;
		mx_scale = mx_scale_CA;
	}
	else
	if (is_n)
	{
		x = x_N;
		yz = yz_N;
		mx_scale = mx_scale_N;
	}
	else   // C atom
	{
		x = x_C;
		yz = yz_C;
		mx_scale = mx_scale_C;
	}

	Point p(x, yz * cos_tor, yz * sin_tor);

	// find the columns of the rotation matrix (mx, my and mz)

	Point vec_1_2 = p2.minus(p1);
	Point mx = p3.minus(p2);
	Point mz = vec_1_2.cross_product(mx);
	mx.scale(mx_scale);		// normalise to length 1.0
	mz.normalise();
	Point my = mz.cross_product(mx);

	Matrix_3_3 m(mx, my, mz);
	return m.times(p).plus(p3);
}

void get_initial_ideal(Point *n, Point *ca, Point *c,
	double angle /*=BOND_ANGLE_N_CA_C*/)
{
	*n = Point(-BOND_LENGTH_N_CA, 0, 0);
	*ca = Point(0, 0, 0);
	*c = Point(-cos(angle), sin(angle), 0);
	c->normalise();
	c->scale(BOND_LENGTH_C_C);
	
	// assert(fabs(angle_formed(*n, *ca, *c) - BOND_ANGLE_N_CA_C) < 0.0001);
}

#ifdef TEST

#include <iostream>
#include <cstdlib>
#include <ctime>

void test_ideal()
{
	static const int Num = 20;
	Point n_pos[Num], ca_pos[Num], c_pos[Num];
	double phi[Num], psi[Num];
	double omega = deg2rad(180.0);

	get_initial_ideal(&n_pos[0], &ca_pos[0], &c_pos[0]);

	double phi_deg = -35.0;
	double psi_deg = 150.0;
	double prev_psi_deg = 160.0;
	double prev_omega = omega;

	phi[0] = psi[0] = phi[Num-1] = psi[Num-1] = 0.0;

	int n;
	for (n = 1;n < Num;n++)
	{
		n_pos[n] = torsion_to_ideal(n_pos[n-1], ca_pos[n-1], c_pos[n-1],
			deg2rad(prev_psi_deg),
			false, true);

		psi[n-1] = prev_psi_deg;
		
		ca_pos[n] = torsion_to_ideal(ca_pos[n-1], c_pos[n-1], n_pos[n],
			deg2rad(prev_omega),
			true, false);

		c_pos[n] = torsion_to_ideal(c_pos[n-1], n_pos[n], ca_pos[n],
			deg2rad(phi_deg),
			false, false);

		phi[n] = phi_deg;

		prev_psi_deg = psi_deg;
		prev_omega = omega;
		phi_deg += 10.0;
		psi_deg -= 15.0;
	}

	for (n = 0;n < Num;n++)
	{
		std::cout
			<< n << ") "
			<< phi[n] << "\t" << psi[n] << "\t"
			<< "N " << n_pos[n]
			<< "\tCA " << ca_pos[n]
			<< "\tC " << c_pos[n]
			<< "\n";

		if (n > 0)
		{
			double val1 = range_0_360(phi[n] + 0.01);
			double val2 = range_0_360(rad2deg(torsion_angle(c_pos[n-1], n_pos[n], ca_pos[n], c_pos[n])) + 0.01);

			if (fabs(val1 - val2) >= 0.0001)
			{
				std::cout << "val1 " << val1 << " val2 " << val2 << "\n";
			}

			assert(fabs(val1 - val2) < 0.0001);
			assert(fabs(omega - rad2deg(torsion_angle(ca_pos[n-1], c_pos[n-1], n_pos[n], ca_pos[n]))) < 0.0001);
			assert(fabs(n_pos[n].distance(c_pos[n-1]) - BOND_LENGTH_C_N) < 0.0001);
			assert(fabs(angle_formed(c_pos[n-1], n_pos[n], ca_pos[n]) - BOND_ANGLE_C_N_CA) < 0.0001);
		}

		assert(fabs(ca_pos[n].distance(n_pos[n]) - BOND_LENGTH_N_CA) < 0.0001);
		assert(fabs(c_pos[n].distance(ca_pos[n]) - BOND_LENGTH_C_C) < 0.0001);
		assert(fabs(angle_formed(n_pos[n], ca_pos[n], c_pos[n]) - BOND_ANGLE_N_CA_C) < 0.0001);

		if (n < Num - 1)
		{
			double val1 = range_0_360(psi[n] + 0.01);
			double val2 = range_0_360(rad2deg(torsion_angle(n_pos[n], ca_pos[n], c_pos[n], n_pos[n+1])) + 0.01);
			assert(fabs(val1 - val2) < 0.0001);

			assert(fabs(angle_formed(ca_pos[n], c_pos[n], n_pos[n+1]) - BOND_ANGLE_CA_C_N) < 0.0001);
		}
	}

	std::cout << "test_ideal(): OK\n";
}

void test_torsion()
{
	srand48(1);

	for (int n = 0;n < 10000;n++)
	{
		Point p1(-4, 4, 1);
		Point p2(-2, 2, 1);
		Point p3(3, 2, 1);
		Point wanted(7, 5, 4);

		// create a random rotation matrix
		Point mx(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);
		mx.normalise();

		Point r(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);
		Point mz = mx.cross_product(r);
		mz.normalise();

		Point my = mz.cross_product(mx);

		Matrix_3_3 m(mx, my, mz);
		m.verify_orthonormal();

		Point d(drand48() * 10.0 - 5.0, drand48() * 10.0 - 5.0, drand48() * 10.0 - 5.0);

		p1 = m.times(p1).plus(d);
		p2 = m.times(p2).plus(d);
		p3 = m.times(p3).plus(d);
		wanted = m.times(wanted).plus(d);

		Point pred = torsion_to_coord(p1, p2, p3,
			p3.distance(wanted),			// bond length
			angle_formed(p2, p3, wanted),	// bond angle
			deg2rad(45),					// torsion angle
			5.0);							// prev. bond length

		if (pred.distance(wanted) > 0.0001)
		{
			std::cout << "Error in test_torsion (n = " << n << ")\n";
			exit(1);
		}
	}

	std::cout << "test_torsion(): OK\n";
}

void test_torsion2()
{
	srand48(1);

	for (int n = 0;n < 10000;n++)
	{
		Point p1(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);
		Point p2(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);
		Point p3(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);
		Point p4(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);

		double t = torsion_angle(p1, p2, p3, p4);

		Point pred = torsion_to_coord(p1, p2, p3,
			p3.distance(p4),				// bond length
			angle_formed(p2, p3, p4),		// bond angle
			t,								// torsion angle
			p2.distance(p3));				// prev. bond length

		double dist = pred.distance(p4);

		if (dist > 0.0001)
		{
			std::cout << "Error in test_torsion2 (n = " << n << ")\n";
			exit(1);
		}
	}

	std::cout << "test_torsion2(): OK\n";
}

void test_1hkx_torsion()
{
	// points from protein 1HKX

	Point c_335(3.101, -8.304, 41.379);

	Point n_336(3.433, -7.021, 41.473);
	Point ca_336(3.177, -6.071, 40.395);
	Point c_336(4.444, -5.636, 39.652);

	Point n_337(5.491, -6.457, 39.711);
	Point ca_337(6.746, -6.140, 39.031);
	Point c_337(6.494, -5.895, 37.540);

	Point n_338(6.114, -6.939, 36.810);

	double phi = torsion_angle(c_335, n_336, ca_336, c_336);
	double psi = torsion_angle(n_336, ca_336, c_336, n_337);
	double omega = torsion_angle(ca_336, c_336, n_337, ca_337);

	std::cout << "336 phi = " <<
		range_m180_180(rad2deg(phi)) << " (expected -107.589)\n";
	std::cout << "336 psi = " <<
		range_m180_180(rad2deg(psi)) << " (expected 23.0184)\n";
	std::cout << "336 omega = " <<
		range_m180_180(rad2deg(omega)) << "\n";

	std::cout << "337 phi = " <<
		range_m180_180(rad2deg(torsion_angle(c_336, n_337, ca_337, c_337)))
		<< " (expected -55.719)\n";
	std::cout << "337 psi = " <<
		range_m180_180(rad2deg(torsion_angle(n_337, ca_337, c_337, n_338)))
		<< " (expected -67.612)\n\n";

	//////////////////////////////

	// predict position of c_336

	double ba1 = angle_formed(n_336, ca_336, c_336);

	std::cout << "Bond angle = " << rad2deg(ba1)
		<< ", torsion angle = " << range_m180_180(rad2deg(phi))
		<< "\n\n";

	Point pred = torsion_to_coord(c_335, n_336, ca_336,
		ca_336.distance(c_336), ba1, phi,
		n_336.distance(ca_336));

	std::cout << "Actual = " << c_336 << "\n"
			  << "Pred   = " << pred << "\n";
}

void test_cb()
{
    Point c(-0.471405, 0.816497, -0.333333);
    Point cb(0.942809, 0, -0.333333);
    Point n(-0.471405, -0.816497, -0.333333);
    c.scale(1.52);
    cb.scale(1.52);
    n.scale(1.45);

	Point ca(10, -20, 30);
	c.add(ca);
	cb.add(ca);
	n.add(ca);

	Point p = estimate_CB_pos(ca, n, c);
	
	std::cout << "CB " << cb
		<< " p " << p
		<< " dist " << p.distance(cb)
		<< "\n";
}

void test_angles()
{
	Point v1(0,  0, 1);
	Point v2(0,  1, 1);
	Point v3(-1, 0, 1);
	Point v4(1,  0, 1);

	v2.normalise();
	v3.normalise();
	v4.normalise();

	std::cout << "Dot product = " << v1.dot_product(v2) << "\n";
	std::cout << "Dot product = " << v1.dot_product(v3) << "\n";
	std::cout << "Dot product = " << v1.dot_product(v4) << "\n";
}

int main(int argc, char **argv)
{
	// test_cb();
	// test_torsion();
	// test_1hkx_torsion();
	// test_ideal();
	test_angles();
}

#endif // TEST
