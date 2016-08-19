#ifndef CONFORMATION_H_INCLUDED
#define CONFORMATION_H_INCLUDED

// This class holds all of the position data (and related information
// such as torsion angles) for a Peptide.
//
// Space is always allocated for backbone atoms;
// if set_pos() is called with a non-backbone atom, space is
// also allocated for side chain atoms.

#include <cassert>
#include <vector>
#include "atom_id.h"
#include "point.h"

#define TORSION_UNKNOWN 999.0

class Peptide;
class Transform;

class Conformation
{
	friend class CORE_impl;
public:
	Conformation();
	~Conformation();

	// associate the conformation with a particular Peptide
	// (calls set_num_res() with the Peptide's current length)
	Conformation(const Peptide *p);

	Conformation(const Conformation &other);
	Conformation& operator = (const Conformation &other);

	// swap all values with another conformation
	void swap(Conformation &other);

	// set the number of residues in the conformation
	void set_num_res(int n);

	void clear();

	// get the position of an atom in residue n
	Point &pos(int n, Atom_Id atom);

	// const version of above function
	const Point &pos(int n, Atom_Id atom) const
	{
		// call non-const version
		return (const_cast<Conformation*>(this))->pos(n, atom);
	}
	
	// set the position of an atom in residue n
	void set_pos(int n, Atom_Id atom, const Point &p);

	// transform an atom position
	void transform_pos(int n, Atom_Id a, const Transform &t);

	// forget positions of all backbone atoms
	void remove_non_backbone_atoms();

	// calculate torsion angles for all residues
	void calc_torsion_angles();

	// make sure torsion angles are correct
	void verify_torsion_angles() const;

	// torsion angles in radians (returns TORSION_UNKNOWN if value has not been
	// calculated, including due to missing backbone atoms)
	double phi(int n) const;
	double psi(int n) const;
	double omega(int n) const;

	// set torsion angles (in radians)
	void set_phi(int n, double angle);
	void set_psi(int n, double angle);
	void set_omega(int n, double angle);

	// check if the angles for residue n, n - 1 and n + 1
	// are consistent with an alpha helix
	bool in_helix(int n) const;

	const Peptide *peptide() const
	{ return m_peptide; }

	void set_peptide(const Peptide *p)
	{ m_peptide = p; }

private:

	struct ResData
	{
		// torsion angles in radians (between -pi and pi, or TORSION_UNKNOWN)
		double phi, psi, omega;

		ResData() :
			phi(TORSION_UNKNOWN), psi(TORSION_UNKNOWN), omega(TORSION_UNKNOWN)
		{
		}
	};

	typedef std::vector<ResData> ResData_Vec;

	// peptide associated with this conformation
	const Peptide *m_peptide;

	// For the sake of allocation speed, use a single vector to
	// store all backbone atom positions.
	//
	// The index for residue n, Atom_Id a is (n * Num_Backbone + a)
	Point_Vec m_backbone;

	// Single array for all non-backbone atom positions.
	// Only allocated if needed (ie. size 0 if there are no
	// non-backbone atoms).
	//
	// The index for residue n, Atom_Id a is
	// (n * Num_Non_Backbone + a - Num_Backbone)

	Point_Vec m_non_backbone;

	// Other data for each residue (torsion angles etc.)
	ResData_Vec m_res_data;
};

typedef std::vector<Conformation> Conf_Vec;

#endif // CONFORMATION_H_INCLUDED

