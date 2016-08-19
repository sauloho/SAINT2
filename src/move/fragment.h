
#ifndef FRAGMENT_H_INCLUDED
#define FRAGMENT_H_INCLUDED

#include <vector>
#include <assert.h> // PG added this
// Torsion and bond angles for a single residue

struct Residue_Angles
{
	double phi, psi, omega;				// torsion angles
	double n_angle, ca_angle, c_angle;	// bond angles

	Residue_Angles() :
		phi(0.0), psi(0.0), omega(0.0),
		n_angle(0.0), ca_angle(0.0), c_angle(0.0)
	{ }

	Residue_Angles(double phi_val, double psi_val, double omega_val,
		double n_val, double ca_val, double c_val) :
		phi(phi_val), psi(psi_val), omega(omega_val),
		n_angle(n_val), ca_angle(ca_val), c_angle(c_val)
	{ }
};

class Fragment
{
public:
	static const double Min_Score;

	// add a residue to the fragment (angles are in radians)
	void add(double phi_val, double psi_val, double omega_val,
		double n_val, double ca_val, double c_val)
	{
		m_angle.push_back(Residue_Angles(
			phi_val, psi_val, omega_val, n_val, ca_val, c_val));
	}

	void set_score(double val)
	{
		assert(val >= Min_Score);
		m_score = val;
	}

	// get the fragment's size
	int length() const
	{ return m_angle.size(); }
	
	// get the fragment's score
	double score() const
	{ return m_score; }

	// get the phi angle for the nth residue (first index is 0)
	// Angle is in radians
	double phi(int n) const
	{
		assert(n >= 0 && n < length());
		return m_angle[n].phi;
	}

	// get the psi angle for the nth residue (first index is 0)
	// Angle is in radians
	double psi(int n) const
	{
		assert(n >= 0 && n < length());
		return m_angle[n].psi;
	}

	// get the omega angle for the nth residue (first index is 0)
	// Angle is in radians
	double omega(int n) const
	{
		assert(n >= 0 && n < length());
		return m_angle[n].omega;
	}

	double n_angle(int n) const
	{
		assert(n >= 0 && n < length());
		return m_angle[n].n_angle;
	}

	double ca_angle(int n) const
	{
		assert(n >= 0 && n < length());
		return m_angle[n].ca_angle;
	}

	double c_angle(int n) const
	{
		assert(n >= 0 && n < length());
		return m_angle[n].c_angle;
	}

private:
	double m_score;
	std::vector<Residue_Angles> m_angle;
};

typedef std::vector<Fragment> Fragment_Vec;

#endif // FRAGMENT_H_INCLUDED

