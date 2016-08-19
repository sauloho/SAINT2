#ifndef ATOM_ID_H_INCLUDED
#define ATOM_ID_H_INCLUDED

// Types of atom (wrapped by class Atom_Type).
//
// (The order of these atom types must be kept in sync with the
// values in class Atom_Type_Helper in atom_type.cpp. The first
// five (the backbone atoms) should never be changed, since they
// are used as indexes in the Residue class).

enum Atom_Id
{
	// backbone atoms (including beta carbon)

	Atom_N, Atom_CA, Atom_C, Atom_O, Atom_CB,
	Num_Backbone,

	// side chain carbon

	Atom_CG = Num_Backbone, Atom_CG1, Atom_CG2,
	Atom_CD, Atom_CD1, Atom_CD2,
	Atom_CE, Atom_CE1, Atom_CE2, Atom_CE3,
	Atom_CZ, Atom_CZ2, Atom_CZ3,
	Atom_CH2,

	// side chain nitrogen

	Atom_ND1, Atom_ND2,
	Atom_NE, Atom_NE1, Atom_NE2,
	Atom_NZ,
	Atom_NH1, Atom_NH2,

	// side chain oxygen

	Atom_OG, Atom_OG1,
	Atom_OD1, Atom_OD2,
	Atom_OE1, Atom_OE2,
	Atom_OH,

	// side chain sulphur

	Atom_SG,
	Atom_SD,

	// unknown atom type

	Atom_Undef,
	Num_Atom_Ids,

	Num_NonBackbone = Atom_Undef - Num_Backbone
};

#endif // ATOM_ID_H_INCLUDED

