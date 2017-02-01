
#ifndef MOVER_FRAGMENT_REV_H_INCLUDED
#define MOVER_FRAGMENT_REV_H_INCLUDED

// A Mover that is based on fragment replacement
// with extrusions in reverse (C to N instead of N to C)

#include <string>
#include "mover.h"
#include "mover_fragment.h"
#include "point.h"
#include "fragment.h"
#include "distribution.h"
#include "conformation.h"

// forward declarations
class Peptide;

class Mover_Fragment_Rev : public Mover_Fragment
{
public:
	// constructor
	Mover_Fragment_Rev();

	// destructor
	virtual ~Mover_Fragment_Rev();

	// set a peptide to its initial state (sequential)
	// Note: since this class uses fragment replacement, it wouldn't
	// make sense to have an initial length of less than any fragment's
	// size, so initial_length is treated as a minimum value
	virtual void init_sequential(Peptide &s, int initial_length,
		Run_Observer *observer);

	// start a sequential run from a pdb file which has been read in
	virtual void init_sequential_from_segment(Peptide &s, int initial_length,
		Run_Observer *observer);

	// create a set of structures from a peptide; each one is a random
	// move away from the original structure
	virtual void do_random_move(Peptide &p, int num,
		bool exhaustive_for_pos, Conf_Vec &result,
		Run_Observer *observer, int curr_length_moves);

	// replace a random fragment in the peptide (may be still growing
	// or fully grown)
	void do_random_move(Peptide &p, Run_Observer *observer, int curr_length_moves);

	// extend the peptide by the requested number of residues
	virtual void extend(Peptide &s, int num_res,
		bool ribosome_wall, Run_Observer *observer);

	// dump internal state (debugging function)
	virtual void dump();

private:
	// disable copy and assignment by making them private
	Mover_Fragment_Rev(const Mover_Fragment_Rev&);
	Mover_Fragment_Rev &operator = (const Mover_Fragment_Rev&);

	// initialise vector m_end_fragment[]
	void init_end_fragments();

	// initialise probability distributions
	void init_distributions();

	// pick a random start (C terminal) fragment with a minimum length
	Fragment *get_starting_fragment(int min_length);

	// pick a random fragment starting at the specified position
	Fragment *random_fragment_starting_at(int start_pos);

	// pick a random fragment (evenly distributed over all positions
	// >= min_start_pos; for all possible positions, use 0)
	Fragment *random_fragment(int min_start_pos, int *start_pos);

	// set the torsion and bond angles in the peptide to the angles in the
	// fragment (ending at p_end_index)
	void change_angles(Peptide &p, int p_end_index, const Fragment *f,
			int f_start_index, int f_end_index, int verbose);

	// add a new fragment
	virtual Fragment *add_fragment(int start_pos, int length);

	// called at end of load_fragments()
	virtual void after_fragments_loaded(int c_terminus);

private:
	// list of fragments starting at each starting position
	Fragment_Vec_Vec m_fragment;

	// probability distribution for fragments at each starting position
	Distribution_Vec m_frag_distrib;

	// list of fragments ending at the C terminus (pointers into m_fragment)
	Fragment_Ptr_Vec m_end_fragment;

	// probability distribution for m_end_fragment[]
	Distribution m_end_frag_distrib;

	// last available start position in m_fragment (ie. last index
	// in m_fragment for which m_fragment[index].size() != 0)
	int m_last_start_pos;

	// index of C terminus
	int m_c_terminus;
};

#endif // MOVER_FRAGMENT_REV_H_INCLUDED

