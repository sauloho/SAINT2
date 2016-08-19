
#ifndef GEOM_H_INCLUDED
#define GEOM_H_INCLUDED

// Geometrical functions.

#include <cmath>
#include "point.h"

class Matrix_3_3;
class Peptide;

// convert radians to degrees (0.0 <= return value < 360.0)
double rad2deg(double radians);

// convert degrees to radians (0 <= return value < 2 pi)
double deg2rad(double degrees);

// force a degree value into the range 0 <= val < 360
double range_0_360(double val);

// force a degree value into the range -180 <= val < 180
double range_m180_180(double val);

// force a radians value into the range 0 <= val < 2 pi
double range_0_2pi(double val);

// force a radians value into the range -pi <= val < pi
double range_mpi_pi(double val);

// torsion (dihedral) angle (in radians, from -PI to PI)
//
// Equals the angle formed between p2p1 and p3p4 while looking down p2p3.
//
// phi is torsion_angle(prev C, N, CA, C)
// psi is torsion_angle(N, CA, C, next N)
// omega is torsion_angle(CA, C, next N, next CA)

double torsion_angle(const Point& p1, const Point& p2,
	const Point& p3, const Point& p4);

// find the angle p1-p2-p3 (in radians)
// ie. the angle between vectors (p2 to p1) and (p2 to p3).
// Returns 0 if p1 or p3 is very close to p2.

double angle_formed(const Point& p1, const Point& p2, const Point& p3);

// estimate the CB position from backbone CA, N and C positions

Point estimate_CB_pos(const Point &ca_pos, const Point &n_pos,
	const Point &c_pos);

// estimate the oxygen position from backbone C, CA and N positions
// (C and CA are in the same residue as the O, N is in the next residue)
Point estimate_O_pos(const Point &ca_pos, const Point &c_pos,
	const Point &next_n_pos);

// calculate the position of the next backbone atom given the last three
// backbone positions, the bond length & angle, the torsion angle
// and the previous bond length (ie. distance p2-p3)
// (angles are in radians).

Point torsion_to_coord(const Point &p1, const Point &p2, const Point &p3,
	double bond_length, double bond_angle, double torsion_angle,
	double prev_bond_length);

// same as torsion_to_coord(), but assumes ideal bond lengths and angles
// (torsion angles is in radians).

Point torsion_to_ideal(const Point &p1, const Point &p2, const Point &p3,
	double torsion_angle,
	bool is_ca,		// whether this is a CA atom
	bool is_n		// whether this is an N atom
);

// constants

namespace
{
	// ideal bond lengths and angles (taken from
	// "Proteins: Structures and Molecular Properties", T. E. Creighton)

	static const double BOND_LENGTH_C_C = 1.52;		// any two C atoms
	static const double BOND_LENGTH_C_N = 1.33;		// backbone C and N
	static const double BOND_LENGTH_N_CA = 1.45;	// N and CA
	static const double BOND_LENGTH_C_O = 1.23;		// backbone C and O

	static const double BOND_ANGLE_N_CA_C = deg2rad(109.5);
	static const double BOND_ANGLE_CA_C_N = deg2rad(115.6);
	static const double BOND_ANGLE_C_N_CA = deg2rad(121.9);
}

// assign ideal positions to the atoms in the first residue
//
// The N atom is placed on the (negative) X axis
// The CA atom is placed at (0, 0, 0)
// The C atom is placed on the XY plane

void get_initial_ideal(Point *n, Point *ca, Point *c,
	double angle = BOND_ANGLE_N_CA_C);

#endif // GEOM_H_INCLUDED

