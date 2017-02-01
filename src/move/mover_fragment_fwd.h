
#ifndef MOVER_FRAGMENT_FWD_H_INCLUDED
#define MOVER_FRAGMENT_FWD_H_INCLUDED

// A Mover that is based on fragment replacement.

#include <string>
#include "mover.h"
#include "mover_fragment.h"
#include "point.h"
#include "fragment.h"
#include "distribution.h"
#include "conformation.h"

// forward declarations
class Peptide;
class Matrix_3_3;

class Mover_Fragment_Fwd : public Mover_Fragment
{
public:
	// constructor
	Mover_Fragment_Fwd();

	// destructor
	virtual ~Mover_Fragment_Fwd();

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
	Mover_Fragment_Fwd(const Mover_Fragment_Fwd&);
	Mover_Fragment_Fwd &operator = (const Mover_Fragment_Fwd&);

	// initialise vector m_start_fragment[]
	void init_start_fragments();

	// initialise probability distributions
	void init_distributions();

	// pick a random start (N terminal) fragment with a minimum length
	Fragment *get_starting_fragment(int min_length);

	// pick a random fragment ending at the specified position
	Fragment *random_fragment_ending_at(int end_pos);

	// pick a random fragment (evenly distributed over all positions
	// <= max_end_pos; for all possible positions, use Peptide length - 1)
	Fragment *random_fragment(int max_end_pos, int *end_pos);

	// set the torsion and bond angles in the peptide to the angles in the
	// fragment (starting from p_start_index)
	void change_angles(Peptide &p, int p_start_index, const Fragment *f,
			int f_start_index, int f_end_index);

	// add a new fragment
	virtual Fragment *add_fragment(int start_pos, int length);

	// called at end of load_fragments()
	virtual void after_fragments_loaded(int c_terminus);

private:
	// list of fragments ending at each ending position
	Fragment_Vec_Vec m_fragment;

	// probability distribution for fragments at each ending position
	Distribution_Vec m_frag_distrib;

	// list of fragments starting at position 0 (pointers into m_fragment)
	Fragment_Ptr_Vec m_start_fragment;

	// probability distribution for fragments starting at position 0
	Distribution m_start_frag_distrib;

	// first available ending position in m_fragment (ie. first index
	// in m_fragment for which m_fragment[index].size() != 0)
	int m_build_from_pos;
};

#endif // MOVER_FRAGMENT_FWD_H_INCLUDED

