
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <iostream>
#include <sstream>		// for std::ostringstream
#include "mover.h"
#include "mover_fragment_rev.h"
#include "run_observer.h"
#include "c_file.h"
#include "config.h"
#include "peptide.h"
#include "conformation.h"
#include "runner.h"
#include "point.h"
#include "random.h"
#include "geom.h"
#include "transform.h"
#include "stream_printf.h"

Mover_Fragment_Rev::Mover_Fragment_Rev()
	: m_c_terminus(-1)
{
}

Mover_Fragment_Rev::~Mover_Fragment_Rev()
{
}

Fragment *Mover_Fragment_Rev::random_fragment(int min_start_pos, int *start_pos)
{
	assert(min_start_pos <= *p_m_build_from_pos);
	int num;
	int count = 0;

	if (min_start_pos < 0)
	{
		min_start_pos = 0;
	}

	int num_possible = *p_m_build_from_pos - min_start_pos + 1;

	// find a position with at least one fragment in it

	do
	{
		// avoid infinite loop
		if (++count > 100000)
		{
			std::cerr << "Infinite loop: no fragments found in "
				"Mover_Fragment_Rev::random_fragment()\n";
			exit(1);
		}

		*start_pos = rand() % num_possible + min_start_pos;
		num = (int) m_fragment[*start_pos].size();
	} while (num == 0);

	int n = m_frag_distrib[*start_pos].select();
	return &(m_fragment[*start_pos][n]);
}

// DONE
Fragment *Mover_Fragment_Rev::random_fragment_starting_at(int start_pos)
{
	assert(start_pos >= 0 && start_pos < (int) m_fragment.size());
	int num = (int) m_fragment[start_pos].size();
	
	if (num == 0)
	{
		std::cerr << "Error: no fragments starting at position "
			<< start_pos << "\n";
		exit(1);
	}

	int n = m_frag_distrib[start_pos].select();
	return &m_fragment[start_pos][n];
}

// DONE
Fragment *Mover_Fragment_Rev::get_starting_fragment(int min_length)
{
	if (m_end_fragment.size() == 0)
	{
		std::cerr << "Error: no fragments at C terminus\n";
		exit(1);
	}

	Fragment *frag = NULL;
	
	int count;
	for (count = 0;count < 100;count++)
	{
		int n = m_end_frag_distrib.select();
		Fragment *f = m_end_fragment[n];

		if (frag == NULL || f->length() < frag->length())
		{
			frag = f;

			if (frag->length() <= min_length)
			{
				break;
			}
		}
	}

	// no fragment found with length <= min_length, so just return
	// shortest one found
	return frag;
}

void Mover_Fragment_Rev::init_sequential(Peptide &p, int initial_length,
	Run_Observer *observer)
{
	std::cout << "In init_sequential -- m_build_from_pos: " << *p_m_build_from_pos << "\n";
	assert(initial_length > 0);
	assert(initial_length <= p.full_length());
	load_fragments(observer);

	// m_c_terminus should have been set by load_fragments()
	assert(m_c_terminus != -1);

	if (m_c_terminus != p.full_length() - 1)
	{
		std::cerr
			<< "Error: fragment library implies C terminus is at position "
			<< m_c_terminus + 1 << ", but sequence length is "
			<< p.full_length() << "\n";
		exit(1);
	}

	Fragment *f = get_starting_fragment(initial_length);
	p.set_length(f->length());
	//change_angles(p, p.end(), f);
	change_angles(p, p.end() - initial_length + 1, f, 0, f->length() - 1, 0);

	if (p.length() < initial_length)
	{
		extend(p, initial_length - p.length(),
			true,	// ribosome wall -- no harm done if it isn't actually true
			observer);
	}
}

void Mover_Fragment_Rev::init_sequential_from_segment(Peptide &p, int initial_length,
	Run_Observer *observer)
{
	// don't limit fragment library use unless we have a long start segment
	// m_build_from_pos is -1 by default
	if (initial_length > 9)
	{
		// we want the last fragment that we are using to have its end at the end
		// of the segment. Minus 1 allows segment to be used also.
		*p_m_build_from_pos = p.full_length() - initial_length + 1;
	}

	std::cout << "In init_sequential_from_segment -- m_build_from_pos: " << *p_m_build_from_pos << "\n";
	
	assert(initial_length > 0);
	assert(initial_length <= p.full_length());
	assert(initial_length == p.length());
	load_fragments(observer);

	// m_c_terminus should have been set by load_fragments()
	assert(m_c_terminus != -1);

	if (m_c_terminus != p.full_length() - 1)
	{
		std::cerr
			<< "Error: fragment library implies C terminus is at position "
			<< m_c_terminus + 1 << ", but sequence length is "
			<< p.full_length() << "\n";
		exit(1);
	}

	/*Fragment *f = get_starting_fragment(initial_length);
	p.set_length(f->length());
	//change_angles(p, p.end(), f);
	change_angles(p, p.end() - initial_length + 1, f);

	if (p.length() < initial_length)
	{
		extend(p, initial_length - p.length(),
			true,	// ribosome wall -- no harm done if it isn't actually true
			observer);
	}*/
}

// DONE
void Mover_Fragment_Rev::do_random_move(Peptide &p, int num,
	bool /*exhaustive_for_pos*/, Conf_Vec &result, Run_Observer *observer, int curr_length_moves)
{
	result.clear();

	for (int n = 0;n < num;n++)
	{
		result.push_back(p.conf());

		p.conf().swap(result.back());
		do_random_move(p, observer, curr_length_moves);
		p.conf().swap(result.back());
	}
}

// DONE
void Mover_Fragment_Rev::do_random_move(Peptide &p,
	Run_Observer *observer, int curr_length_moves)
{
	int verbose = 0;
	assert(m_fragments_loaded);

	int start, end;
	Fragment *f = random_fragment(p.start(), &start);
	end = start + f->length() - 1;
	//change_angles(p, end, f);
	//change_angles(p, start, f);
	//std::cout << "p.length(): " << p.length() << "\n";
	if (p.length() == 221 && curr_length_moves == 23)
	{
		if (!verbose)
		{
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~\nReached length 113\n~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		}
		verbose = 1;
	}
	
	if (verbose)
	{
		std::cout << "starting do_random_move change angles\n";
	}
	// Is the fragment overlapping with the segment?
	if (end <= *p_m_build_from_pos)
	{
		if (verbose)
		{
			std::cout << "frag doesn't overlap\n";
		}
		// if not, we're fine, use whole fragment
		change_angles(p, start, f, 0, f->length() - 1, verbose);
	}
	else
	{
		if (verbose)
		{
			std::cout << "frag overlaps\n";
			std::cout << "frag overlap is: " << (end - *p_m_build_from_pos) << "\n";
			std::cout << "start: " << start << "  f->length(): " << f->length() << "  end: " << end << "\n";
		}
		// overlapping, so truncate fragment
		// overlap is (end - *p_m_build_from_pos)
		change_angles(p, start, f, 0, f->length() - 1 - (end - *p_m_build_from_pos), verbose);
		//change_angles(p, start, f, 0, f->length() - 1);
	}
	//std::cout << "finished do_random_move change angles\n";

	if (m_double_replacement_prob != 0.0)
	{
		if (m_double_replacement_prob == 1.0 ||
			Random::rnd(1.0) < m_double_replacement_prob)
		{
			f = random_fragment(p.start(), &start);
			end = start + f->length() - 1;
			//change_angles(p, end, f);
			//change_angles(p, start, f);
			// Is the fragment overlapping with the segment?
			if (end <= *p_m_build_from_pos)
			{
				// if not, we're fine, use whole fragment
				change_angles(p, start, f, 0, f->length() - 1, verbose);
			}
			else
			{
				// overlapping, so truncate fragment
				// overlap is (end - *p_m_build_from_pos)
				change_angles(p, start, f, 0, f->length() - 1 - (end - *p_m_build_from_pos), verbose);
				//change_angles(p, start, f, 0, f->length() - 1);
			}
		}
	}
}

// DONE
void Mover_Fragment_Rev::extend(Peptide &p, int num_res,
	bool ribosome_wall, Run_Observer *observer)
{
	assert(m_fragments_loaded);
	assert(num_res >= 0);

	if (p.length() + num_res > p.full_length())
	{
		num_res = p.full_length() - p.length();
		assert(p.start() - num_res == 0);
	}

	if (num_res == 0)
	{
		return;
	}

//std::cout << "Mover_Fragment_Rev::extend() by " << num_res << "\n";
	int old_start = p.start();
	int new_start = p.start() - num_res;
	int longest_found = 0;
	int count = 0;
	int p_end_index;
	Fragment *f = NULL;

	bool fragment_ok = false;

	while (!fragment_ok)
	{
		if (++count > 10000)
		{
			std::cerr
				<< "Infinite loop: could not find a fragment to extend "
					"from position "
				<< p.start()
				<< " to "
				<< new_start;

			if (longest_found != 0)
			{
				std::cerr
					<< " (no single fragment was long enough; longest "
						"found was "
					<< longest_found
					<< ")";
			}

			std::cerr << "\n";
			exit(1);
		}

		f = random_fragment_starting_at(new_start);

		// find the index of the residue in p corresponding
		// to the fragment's last residue
		p_end_index = new_start + f->length() - 1;

		// make sure the fragment isn't longer than the number of
		// residues extruded (ie. won't reach beyond the C terminus)

		if (p_end_index >= p.full_length())
		{
			continue;
		}

		// need at least one overlapping residue between the fragment
		// and the current peptide

		if (p_end_index < old_start)
		{
			if (f->length() > longest_found)
			{
				longest_found = f->length();
			}

			continue;
		}

		fragment_ok = true;
	}

	assert(p_end_index >= old_start);
	assert(p_end_index - f->length() + 1 == new_start);

	p.add_length(num_res);
	//change_angles(p, p_end_index, f);
	//change_angles(p, new_start, f);
	
	//std::cout << "about to do extend change angles\n";
	//std::cout << "p_end_index: " << p_end_index << "\n";
	//std::cout << "*p_m_build_from_pos: " << *p_m_build_from_pos << "\n";
	//std::cout << "new_start: " << new_start << "\n";
	// Is the fragment overlapping with the segment?
	if (p_end_index <= *p_m_build_from_pos)
	{
		// if not, we're fine, use whole fragment
		change_angles(p, new_start, f, 0, f->length() - 1, 0);
	}
	else
	{
		// overlapping, so truncate fragment
		// overlap is (p_end_index - *p_m_build_from_pos)
		// I really don't know why there's a minus 1 on this line in the fwd version:
		//std::cout << "p_end_index - *p_m_build_from_pos: " << p_end_index - *p_m_build_from_pos << "\n";
		change_angles(p, new_start, f, 0, f->length() - 1 - (p_end_index - *p_m_build_from_pos), 0);
		//change_angles(p, new_start, f, 0, f->length() - 1);
	}
	//std::cout << "finished extend change angles\n";

	if (ribosome_wall && !p.full_grown())
	{
		reorient_for_ribosome(p);
	}

	//p.conf().verify_torsion_angles();
}


void Mover_Fragment_Rev::change_angles(Peptide &p, int p_start_index,
	const Fragment *f, int f_start_index, int f_end_index, int verbose)
{
	//std::cout << "starting assertions\n";
	assert(p_start_index >= 0);
	//std::cout << "passed assertion0\n";
	int p_end_index = p_start_index + f_end_index - f_start_index;
	assert(p_end_index < p.full_length());
	//std::cout << "passed assertion1\n";
	//int f_part_length = f_end_index - f_start_index + 1;
	//assert(p_end_index - f_part_length + 1 >= *p_m_build_from_pos);
	//std::cout << "passed assertion2\n";

	bool realign_before = (p_start_index > p.start());
	Point old_C, old_N, old_CA;

	if (realign_before)
	{
		old_C = p.atom_pos(p_start_index - 1, Atom_C);
		old_N = p.atom_pos(p_start_index, Atom_N);
		old_CA = p.atom_pos(p_start_index, Atom_CA);
//std::cout << "checkpoint 1\n";
//fflush(stdout);
	}

	int n = f_end_index;			// position in fragment
	int i = p_start_index + n;		// position in p

	if (verbose)
	{
		Point n_pos;
		Point ca_pos;
		Point c_pos;
		Point next_n_pos;

		Point c_to_n;
		Point n_to_ca;
		Point ca_to_c;

		std::cout << "############ whole thing -------------------------------\n";
		std::cout << "p.full_length()" << p.full_length() << "\n";
		for (int counter = 0; counter < p.full_length() - 1; counter++)
		{
			n_pos = p.atom_pos(counter, Atom_N);
			ca_pos = p.atom_pos(counter, Atom_CA);
			c_pos = p.atom_pos(counter, Atom_C);
			next_n_pos = p.atom_pos(counter+1, Atom_N);

			n_to_ca = n_pos.minus(ca_pos);
			std::cout << n_to_ca.length() << "\t";

			ca_to_c = ca_pos.minus(c_pos);
			std::cout << ca_to_c.length() << "\t";
			
			if (counter != p.full_length() - 1)
			{
				c_to_n = c_pos.minus(next_n_pos);
				std::cout << c_to_n.length() << "\t";
			}

			std::cout << "\n";
		}
	}

	if (i == p.end() - 1)
	{
std::cout << "checkpoint 2\n";
fflush(stdout);
		Point n_pos, ca_pos, c_pos, pos;
		get_initial_ideal(&n_pos, &ca_pos, &c_pos, f->ca_angle(f_end_index));
		p.set_atom_pos(i, Atom_N, n_pos),
		p.set_atom_pos(i, Atom_CA, ca_pos),
		p.set_atom_pos(i, Atom_C, c_pos);

		if (!p.is_glycine(i))
		{
			p.set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
		}

		// calculate where next N would be to get position for O
//std::cout << "checkpoint 3\n";
//fflush(stdout);

		Point est_n_pos = torsion_to_coord(n_pos, ca_pos, c_pos,
			BOND_LENGTH_C_N,
			BOND_ANGLE_CA_C_N,  // (or f->c_angle(n), but is just an estimate)
			f->psi(n), BOND_LENGTH_C_C);
		p.set_atom_pos(i, Atom_O, estimate_O_pos(ca_pos, c_pos, est_n_pos));

		// also need to find C position for previous residue

		assert(i > 0);
		p.set_atom_pos(i - 1, Atom_C,
			torsion_to_coord(p.atom_pos(i, Atom_C), p.atom_pos(i, Atom_CA),
				p.atom_pos(i, Atom_N), BOND_LENGTH_C_N,
				f->n_angle(n), f->phi(n), BOND_LENGTH_N_CA));
		p.conf().set_phi(i, f->phi(n));

		n--;
		i--;
//std::cout << "checkpoint 4\n";
//fflush(stdout);
	}

	Point ca_pos = p.atom_pos(i + 1, Atom_CA);
	Point n_pos = p.atom_pos(i + 1, Atom_N);
	Point c_pos = p.atom_pos(i, Atom_C);
	Point c_to_n;
	Point n_to_ca;
	Point ca_to_c;

	if (verbose) {
		std::cout << "-------------------------------\n";
	}
	for ( ;n >= f_start_index;n--, i--)
	{
//std::cout << "checkpoint 5\n";
//fflush(stdout);
		if (verbose) {
			n_to_ca = n_pos.minus(ca_pos);
			std::cout << n_to_ca.length() << "\t";
		}
		if (verbose) {
			char pdb_out[150];
std::cout << "checkpoint 5.1.1\n";
fflush(stdout);
			sprintf(pdb_out,"output/1_part%d_%dpre",p.length(),i);
std::cout << "checkpoint 5.1.2\n";
fflush(stdout);
			p.write_pdb(pdb_out);
std::cout << "checkpoint 5.1.3\n";
fflush(stdout);
		}
		if (false)
		{
			ca_pos = torsion_to_coord(ca_pos, n_pos, c_pos, BOND_LENGTH_C_C,
				0.1, f->omega(n), BOND_LENGTH_C_N);
std::cout << "LOOK AT ME: " << f->c_angle(n) << "\n";
fflush(stdout);

		}
		else
		{
			ca_pos = torsion_to_coord(ca_pos, n_pos, c_pos, BOND_LENGTH_C_C,
				f->c_angle(n), f->omega(n), BOND_LENGTH_C_N);
		}
		p.set_atom_pos(i, Atom_CA, ca_pos);
		p.conf().set_omega(i, f->omega(n));
std::cout << "checkpoint 5.1\n";
fflush(stdout);

		if (verbose) {
			char pdb_out[150];
std::cout << "checkpoint 5.1.1\n";
fflush(stdout);
			sprintf(pdb_out,"output/1_part%d_%dbefore",p.length(),i);
std::cout << "checkpoint 5.1.2\n";
fflush(stdout);
			p.write_pdb(pdb_out);
std::cout << "checkpoint 5.1.3\n";
fflush(stdout);
		}
std::cout << "checkpoint 5.2\n";
fflush(stdout);

		p.set_atom_pos(i, Atom_O, estimate_O_pos(ca_pos, c_pos, n_pos));

std::cout << "checkpoint 5.3\n";
fflush(stdout);
		if (verbose) {
			char pdb_out[150];
			sprintf(pdb_out,"output/1_part%d_%dafter",p.length(),i);
			p.write_pdb(pdb_out);
		}
std::cout << "checkpoint 5.4\n";
fflush(stdout);

		if (verbose) {
			c_to_n = c_pos.minus(n_pos);
			std::cout << c_to_n.length() << "\t";
		}
std::cout << "checkpoint 5.5\n";
fflush(stdout);
		n_pos = torsion_to_coord(n_pos, c_pos, ca_pos, BOND_LENGTH_N_CA,
			f->ca_angle(n), f->psi(n), BOND_LENGTH_C_C);
		p.set_atom_pos(i, Atom_N, n_pos);
		p.conf().set_psi(i, f->psi(n));

		if (!p.is_glycine(i))
		{
			p.set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
		}
//std::cout << "checkpoint 6\n";
//fflush(stdout);

		if (i > 0)
		{
//std::cout << "checkpoint 7\n";
//fflush(stdout);
			if (verbose) {
				ca_to_c = ca_pos.minus(c_pos);
				std::cout << ca_to_c.length();
			}
			c_pos = torsion_to_coord(c_pos, ca_pos, n_pos, BOND_LENGTH_C_N,
				f->n_angle(n), f->phi(n), BOND_LENGTH_N_CA);
			p.set_atom_pos(i - 1, Atom_C, c_pos);
			p.conf().set_phi(i, f->phi(n));
		}
		if (verbose) {
			std::cout << "\n";
		}
	}

	if (realign_before)
	{
//std::cout << "checkpoint 8\n";
//fflush(stdout);
		Point new_C = p.atom_pos(p_start_index - 1, Atom_C);
		Point new_N = p.atom_pos(p_start_index, Atom_N);
		Point new_CA = p.atom_pos(p_start_index, Atom_CA);
//std::cout << "checkpoint 8.0\n";
//fflush(stdout);

		// Because the bond lengths are ideal, t.find_alignment() will
		// align the old C and N atoms exactly onto the new positions.
		// The CA atom will be placed in a position which causes
		// the old omega torsion angle (for residue p_start_index - 1)
		// to be retained.
		//
		// The overall effect is:
		//
		// The torsion and bond angles for residue (p_start_index - 1)
		// are the same as the old values.
		// The torsion and bond angles for residue (p_start_index)
		// come from the fragment.
if (verbose) {
std::cout << old_C.x << "\t" << old_C.y << "\t" << old_C.z << "\n";
std::cout << old_N.x << "\t" << old_N.y << "\t" << old_N.z << "\n";
std::cout << old_CA.x << "\t" << old_CA.y << "\t" << old_CA.z << "\n";
std::cout << new_C.x << "\t" << new_C.y << "\t" << new_C.z << "\n";
std::cout << new_N.x << "\t" << new_N.y << "\t" << new_N.z << "\n";
std::cout << new_CA.x << "\t" << new_CA.y << "\t" << new_CA.z << "\n";
fflush(stdout);
}
		Transform t;
		t.find_alignment(old_C, old_N, old_CA, new_C, new_N, new_CA);
//std::cout << "checkpoint 8.01\n";

		for (i = 0;i < p_start_index - 1;i++)
		{
//std::cout << "checkpoint 8.1: " << i << "\n";
//fflush(stdout);
			for (int a = 0;a < Num_Backbone;a++)
			{
				if (!(p.is_glycine(i) && (Atom_Id) a == Atom_CB))
				{
					p.transform_pos(i, (Atom_Id) a, t);
				}
			}
		}

		// i equals (p_start_index - 1)
//std::cout << "checkpoint 8.2\n";
//fflush(stdout);

		p.transform_pos(i, Atom_N, t);
		p.transform_pos(i, Atom_CA, t);
		// Atom_C has already been changed
//std::cout << "checkpoint 8.3\n";
//fflush(stdout);

		if (!p.is_glycine(i))
		{
			p.set_atom_pos(i, Atom_CB,
				estimate_CB_pos(
					p.atom_pos(i, Atom_CA),
					p.atom_pos(i, Atom_N),
					p.atom_pos(i, Atom_C)));
		}
//std::cout << "checkpoint 8.4\n";
//fflush(stdout);

		p.set_atom_pos(i, Atom_O,
			estimate_O_pos(
				p.atom_pos(i, Atom_CA),
				p.atom_pos(i, Atom_C),
				p.atom_pos(i + 1, Atom_N)));
//std::cout << "checkpoint 9\n";
//fflush(stdout);
	}
}

// DONE
Fragment *Mover_Fragment_Rev::add_fragment(int start_pos, int length)
{
	if (m_fragment.size() <= (unsigned) start_pos)
	{
		m_fragment.resize(start_pos + 1);
	}

	int num = (int) m_fragment[start_pos].size();
	m_fragment[start_pos].resize(num + 1);

	return &m_fragment[start_pos][num];
}

// DONE
void Mover_Fragment_Rev::after_fragments_loaded(int c_terminus)
{
	assert(c_terminus != -1);
	m_c_terminus = c_terminus;
	init_end_fragments();
	init_distributions();
}

// DONE
void Mover_Fragment_Rev::init_end_fragments()
{
	m_end_fragment.clear();
	// m_build_from_pos = -1;

	// find all fragments that end at the C terminus
	// (ie. end at position "pos" and have length "pos + 1")

	for (int pos = (int) m_fragment.size() - 1;pos >= 0;pos--)
	{
		if (*p_m_build_from_pos == -1 && m_fragment[pos].size() != 0)
		{
			*p_m_build_from_pos = pos;
		}

		for (int n = 0;n < (int) m_fragment[pos].size();n++)
		{
			if (pos + m_fragment[pos][n].length() - 1 == m_c_terminus)
			{
				m_end_fragment.push_back(&m_fragment[pos][n]);
			}
		}
	}

	if (m_end_fragment.size() == 0)
	{
		std::cerr << "Error: no fragments found ending at C terminus\n";
		exit(1);
	}
}

// DONE
void Mover_Fragment_Rev::init_distributions()
{
	unsigned int n, pos;

	m_frag_distrib.resize(m_fragment.size());

	for (pos = 0;pos <= (unsigned) *p_m_build_from_pos;pos++)
	{
		for (n = 0;n < m_fragment[pos].size();n++)
		{
			double prob = m_fragment[pos][n].score();
			assert(prob >= Fragment::Min_Score);
			m_frag_distrib[pos].add(prob, n);
		}
	}

	for (n = 0;n < m_end_fragment.size();n++)
	{
		double prob = m_end_fragment[n]->score();
		assert(prob >= Fragment::Min_Score);
		m_end_frag_distrib.add(prob, n);
	}
}

// DONE
void Mover_Fragment_Rev::dump()
{
	std::cout << "Dumping fragments:\n\n";

	for (unsigned int pos = 0;pos < m_fragment.size();pos++)
	{
		std::cout << m_fragment[pos].size()
			<< " fragments starting at " << pos << "\n\n";

		assert((int) m_fragment[pos].size() == m_frag_distrib[pos].num());

		for (unsigned int n = 0;n < m_fragment[pos].size();n++)
		{
			const Fragment &f = m_fragment[pos][n];
			std::cout << pos << " to " << pos + f.length() - 1
				<< " (" << f.score() << ")\n";
		}
	}

	std::cout << "\n" << m_end_fragment.size()
		<< " fragments ending at C terminus\n\n";

	assert((int) m_end_fragment.size() == m_end_frag_distrib.num());

	for (unsigned int n = 0;n < m_end_fragment.size();n++)
	{
		const Fragment &f = *(m_end_fragment[n]);
		std::cout << m_c_terminus - f.length() + 1 << " to "
			<< m_c_terminus << " (" << f.score() << ")\n";
	}
}

