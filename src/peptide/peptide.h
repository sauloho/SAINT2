
#ifndef PEPTIDE_H_INCLUDED
#define PEPTIDE_H_INCLUDED

// Class representing a protein peptide (ie. a conformation).
//
// The peptide may be initialised by calling create_from_sequence()
// with an amino acid sequence or by reading the whole peptide from
// a PDB file using read_pdb().
//
// The peptide may be extended just be changing its length
// (using set_length() or add_length()); do not use the append_residue()
// function for extending (which is intended for creating things like
// fragments one residue at a time).

#include <vector>
#include <iostream>
#include <string>
#include "residue.h"
#include "transform.h"
#include "conformation.h"
#include "common.h"

// forward declarations
class Sequence;

class Peptide
{
	friend class CORE_impl;
public:
	// constructor
	Peptide();

	// destructor
	~Peptide();

	// copy constructor
	Peptide(const Peptide& other);

	// assignment operator
	Peptide& operator = (const Peptide& other);

	// delete all residues
	void clear();

	// initialise all residues in the peptide to the amino acids in
	// the sequence specified (and the codons, if known).
	// The current length (#residues extruded) is set to zero
	void create_from_sequence(const Sequence &seq);

	// read a chain from a PDB file (if the chain is ' ', the first
	// chain found is used).
	// Exits with an error message if the file cannot be opened.
	// Returns false and prints a message if there is an error;
	// in this case, the peptide will not be valid
	bool read_pdb(const char *filename, char chain = ' ',
		bool no_warnings = false);

	// get the filename used in the last call to read_pdb() or write_pdb()
	const char *get_filename() const;

	// write to a PDB file
	// (if filename is an empty string, writes to stdout)
	void write_pdb(const char *filename, bool backbone_only = false,
		bool atom_records_only = false) const;
	
	void forget_is_from_pdb()
	{ m_from_pdb = false; }

	// add a residue to the C terminal end and return a pointer to it;
	// none of the residue values are set.
	//
	// This function is intended to be used for things like fragments
	// rather than the main peptide object (which should use set_length()
	// or add_length() to perform extending).
	//
	// NOTE: sets current length to full length.
	Residue *append_residue();

	// add side chain atoms using SCWRL
	//void call_scwrl();

	// check if the protein has been fully synthesised
	bool full_grown() const;

	// get the length of the full grown peptide
	// (ie. the length of the amino acid sequence)
	int full_length() const
	{ return m_full_length; }

	// get the current length (number of residues extruded)
	int length() const
	{ return m_length; }

	// change the current length (cannot be beyond full_length())
	// (newly extruded residues will need coordinates to be set, etc.)
	void set_length(int len);
	void add_length(int amount);

	// lowest residue index (of the extruded residues)
	int start() const
	{ return (reverseSaint ? m_full_length - m_length : 0); }

	// highest residue index (of the extruded residues)
	int end() const
	{ return (reverseSaint ? m_full_length - 1 : m_length - 1); }

	// get the current conformation of the peptide (all of the atom positions)
	Conformation &conf() { return m_conf; }
	const Conformation &conf() const { return m_conf; }

	// convenience functions for getting and setting atom positions

	bool atom_exists(int n, Atom_Id atom) const
	{ return m_res[n].atom_exists(atom); }

	Point atom_pos(int n, Atom_Id atom)
	{
		assert(m_res[n].atom_exists(atom));
		return m_conf.pos(n, atom);
	}

	Point atom_pos2(int n, Atom_Id atom) const
	{
		return m_conf.pos(n, atom);
	}

	Point atom_pos(int n, Atom_Id atom) const
	{
		assert(m_res[n].atom_exists(atom));
		return m_conf.pos(n, atom);
	}

	void set_atom_pos(int n, Atom_Id atom, Point p)
	{ m_conf.set_pos(n, atom, p); }

	void transform_pos(int n, Atom_Id atom, const Transform &t)
	{ m_conf.transform_pos(n, atom, t); }

	double atom_dist(int n1, Atom_Id atom1, int n2, Atom_Id atom2) const
	{ return atom_pos(n1, atom1).distance(atom_pos(n2, atom2)); }

	// delete all atoms apart from CA, CB and backbone C, N, O
	void remove_non_backbone_atoms();

	// check if there are a large number of non-backbone atoms present
	// (if more than half the non-glycine residues have any non-backbone atoms)
	bool many_non_backbone() const;

	// find the proportion of residues that are in contact with a previous
	// residue (ie. further towards the N terminus).
	// Assumes non-backbone atoms are present and torsion angles
	// have been calculated.
	//
	// Note that it tries to detect if the protein starts with a helix;
	// if so, the region up to the end of the helix is counted as all
	// 100% contacts. This makes the algorithm very sensitive to what
	// happens at the beginning of the protein (a bend in the helix
	// just beyond the "helix" angle limit will have a big effect;
	// so will starting a helix a bit further from the N terminus
	// than the limit (approx. the 10th residue instead of the 9th)).
	double previous_contact_proportion() const;

	// find the proportion of residues that are in contact with and
	// between a pair of previous residues (further towards the N terminus).
	// Assumes non-backbone atoms are present.
	double interleaving_proportion() const;

	// deletes the nth residue
	// (the only normal reason to call this function is when the
	// peptide is being checked and an invalid residue is found)
	// void remove_residue(int n);

	// Remove backbone atoms, find torsion angles, and
	// set all bond lengths and angles to ideal values.
	void idealise();

	// Set all bond lengths to ideal values, but keep the current
	// torsion and bond angles
	void idealise_bond_lengths();

	// get the radius of gyration of all the CB atoms in the peptide
	// (CA atom for glycine) (= square root of the average squared
	// distance from the mean)
	// If up_to_res is specified, RG is only calculated from the
	// N terminus up to that residue (inclusive)
	double radius_of_gyr(int up_to_res = -1) const;

	// get the largest distance between any pair of atoms in the peptide
	double diameter() const;

	// get the index of the most recently extruded residue
	int latest() const
	{ return m_length - 1; }

	// get the nth residue (N terminus has index 0)
	const Residue &res(int n) const;
	Residue &res(int n);

	// convenience function to check if residue n is glycine
	bool is_glycine(int n) const
	{ return res(n).amino().is_glycine(); }

	void set_chain(char ch)
	{ m_chain = ch; }

	char chain() const
	{ return m_chain; }

	// get the number of gaps between residues
	int num_missing_res() const;

	// print a list of missing residues (gaps in residue sequence
	// numbers or a long distance between residues)
	void print_missing_res(const char *filename) const;

	// get the number of missing backbone atoms
	int num_missing_backbone() const;

	// debugging function
	void dump(std::ostream &out = std::cout) const;

	//static void set_scwrl_exec(const std::string &path)
	//{ m_scwrl_exec = path; }

	// get the CB position of residue n, or CA for glycine
	bool get_cb_ca_pos(int n, Point *pos) const;

	// get the CB position of residue n, or where the CB would be for glycine
	bool get_side_chain_pos(int n, Point *pos) const;

	// verify that all bond lengths and angles are ideal
	void verify_ideal() const;

	// verify that all of the bond lengths, bond angles and torsion
	// angles are reasonable
	void verify() const;

	// find the RMSD to another peptide
	double calc_rmsd(const Peptide &p) const;

	void dump_backbone(std::ostream &out = std::cout) const;

    void alloc_satisfied_con(int length);

    void set_satisfied_con(int pos, int value);

    // void print_satisfied_con();

    // This variable is used to store which contacts are satisfied in a given conformation
    

private:
	// helper function for copy constructor and assignment operator
	void copy(const Peptide &other);

	// check bond lengths and angles of backbone N, C and CA atoms
	// (atoms with bad positions are removed).
	//
	// If the N atom in one residue is too far from the C atom in
	// the previous residue, set_missing_after(true) is called for
	// the previous residue.
	void check_N_CA_C_positions(const char *filename, bool no_warnings);

	// Check existing CB positions and add missing ones based on the
	// positions of the N, CA and C atoms in each residue
	// (should be called after check_N_CA_positions, otherwise it
	// might base the estimated CB position on incorrect atom positions)
	void check_CB_positions(const char *filename, bool no_warnings);

	// Check existing backbone O positions and add missing ones based on
	// the position of the CA, C and N atoms
	// (should be called after check_N_CA_positions, otherwise it
	// might base the estimated O position on incorrect atom positions)
	void check_O_positions(const char *filename, bool no_warnings);

	// check for missing side chain atoms
	void check_side_chain_atoms(const char *filename, bool no_warnings);

	// make sure m_scwrl_exec is valid
	//static bool verify_scwrl_exec(bool exit_on_err = true);

private:
	//static std::string m_scwrl_exec;	// path of SCWRL executable

	int m_length;				// number of residues extruded so far
	int m_full_length;			// maximum ("full grown") length
    int *m_satisfied_con;
	Residue_Vec m_res;			// data for each residue (indexes correspond
								// to indexes in m_seq)
	Conformation m_conf;		// atom positions and related information
	char m_chain;				// chain code

	bool m_from_pdb;			// whether peptide was read from a PDB file
	mutable std::string m_filename;		// filename in last call to read_pdb()
										// or write_pdb
};

// a vector of peptides
typedef std::vector<Peptide> Peptide_Vec;

#endif // PEPTIDE_H_INCLUDED

