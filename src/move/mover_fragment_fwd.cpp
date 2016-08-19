
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <iostream>
#include <sstream>		// for std::ostringstream
#include "mover.h"
#include "mover_fragment_fwd.h"
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

Mover_Fragment_Fwd::Mover_Fragment_Fwd()
{
}

Mover_Fragment_Fwd::~Mover_Fragment_Fwd()
{
}

Fragment *Mover_Fragment_Fwd::random_fragment(int max_end_pos, int *end_pos)
{
	assert(max_end_pos > 0);
	int num;
	int count = 0;
	int max_frag = (int) m_fragment.size() - 1;

	if (max_end_pos > max_frag)
	{
		max_end_pos = max_frag;
	}

	int num_possible = max_end_pos - m_first_end_pos + 1;

	// find a position with at least one fragment in it

	do
	{
		// avoid infinite loop
		if (++count > 100000)
		{
			std::cerr << "Infinite loop: no fragments found in "
				"Mover_Fragment_Fwd::random_fragment()\n";
			exit(1);
		}

		*end_pos = rand() % num_possible + m_first_end_pos;
		num = (int) m_fragment[*end_pos].size();
	} while (num == 0);

	int n = m_frag_distrib[*end_pos].select();
	return &(m_fragment[*end_pos][n]);
}

Fragment *Mover_Fragment_Fwd::random_fragment_ending_at(int end_pos)
{
	assert(end_pos >= 0 && end_pos < (int) m_fragment.size());
	int num = (int) m_fragment[end_pos].size();
	
	if (num == 0)
	{
		std::cerr << "Error: no fragments ending at position "
			<< end_pos << "\n";
		exit(1);
	}

	// int n = rand() % num;
	int n = m_frag_distrib[end_pos].select();
	return &m_fragment[end_pos][n];
}

Fragment *Mover_Fragment_Fwd::get_starting_fragment(int min_length)
{
	if (m_start_fragment.size() == 0)
	{
		std::cerr << "Error: no fragments at N terminus\n";
		exit(1);
	}

	Fragment *frag = NULL;
	
	int count;
	for (count = 0;count < 100;count++)
	{
		int n = m_start_frag_distrib.select();
		Fragment *f = m_start_fragment[n];

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

void Mover_Fragment_Fwd::init_sequential(Peptide &p, int initial_length,
	Run_Observer *observer)
{
	assert(initial_length > 0);
	assert(initial_length <= p.full_length());
	load_fragments(observer);

	Fragment *f = get_starting_fragment(initial_length);
	p.set_length(f->length());
	change_angles(p, 0, f);

	if (p.length() < initial_length)
	{
		extend(p, initial_length - p.length(),
			true,	// ribosome wall -- no harm done if it isn't actually true
			observer);
	}
}

void Mover_Fragment_Fwd::do_random_move(Peptide &p, int num,
	bool exhaustive_for_pos, Conf_Vec &result, Run_Observer *observer)
{
	result.clear();

	for (int n = 0;n < num;n++)
	{
		result.push_back(p.conf());

		p.conf().swap(result.back());
		do_random_move(p, observer);
		p.conf().swap(result.back());
	}
}

void Mover_Fragment_Fwd::do_random_move(Peptide &p,
	Run_Observer *observer)
{
	assert(m_fragments_loaded);

	int start, end;
	Fragment *f = random_fragment(p.length() - 1, &end);
	start = end - f->length() + 1;
	change_angles(p, start, f);

	if (m_double_replacement_prob != 0.0)
	{
		if (m_double_replacement_prob == 1.0 ||
			Random::rnd(1.0) < m_double_replacement_prob)
		{
			f = random_fragment(p.length() - 1, &end);
			start = end - f->length() + 1;
			change_angles(p, start, f);
		}
	}
}

void Mover_Fragment_Fwd::extend(Peptide &p, int num_res,
	bool ribosome_wall, Run_Observer *observer)
{
	assert(m_fragments_loaded);
	assert(num_res >= 0);

	if (p.length() + num_res > p.full_length())
	{
		num_res = p.full_length() - p.length();
	}

	if (num_res == 0)
	{
		return;
	}

	int old_end = p.length() - 1;
	int new_end = p.length() + num_res - 1;
	int longest_found = 0;
	int count = 0;
	int p_start_index;
	Fragment *f = NULL;

	bool fragment_ok = false;

	while (!fragment_ok)
	{
		if (++count > 10000)
		{
			std::cerr
				<< "Infinite loop: could not find a fragment to extend "
					"from position "
				<< p.length() - 1
				<< " to "
				<< new_end;

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

		f = random_fragment_ending_at(new_end);

		// find the index of the residue in p corresponding
		// to the fragment's residue 0
		p_start_index = new_end - f->length() + 1;

		// make sure the fragment isn't longer than the number of
		// residues extruded (ie. won't reach beyond the N terminus)

		if (p_start_index < 0)
		{
			continue;
		}

		// need at least one overlapping residue between the fragment
		// and the current peptide

		if (p_start_index > old_end)
		{
			if (f->length() > longest_found)
			{
				longest_found = f->length();
			}

			continue;
		}

		fragment_ok = true;
		// fragment_ok = wont_hit_ribosome(f, p,
		//    p_start_index, frag_rot, frag_d);
	}

	assert(p_start_index <= old_end);
	assert(p_start_index + f->length() - 1 == new_end);

	p.add_length(num_res);
	change_angles(p, p_start_index, f);

	if (ribosome_wall && !p.full_grown())
	{
		reorient_for_ribosome(p);
	}
}

void Mover_Fragment_Fwd::change_angles(Peptide &p, int p_start_index,
	const Fragment *f)
{
	assert(p_start_index >= 0);

	//int p_end_index = p_start_index + f->length() - 1;
	//assert(p_end_index < p.length());

	bool realign_before = (p_start_index > 0);
	Point old_C, old_N, old_CA;

	if (realign_before)
	{
		old_C = p.atom_pos(p_start_index - 1, Atom_C);
		old_N = p.atom_pos(p_start_index, Atom_N);
		old_CA = p.atom_pos(p_start_index, Atom_CA);
	}

	int n = f->length() - 1;		// position in fragment
	int i = p_start_index + n;		// position in p

	if (i == p.length() - 1)
	{
		Point n_pos, ca_pos, c_pos, pos;
		get_initial_ideal(&n_pos, &ca_pos, &c_pos, f->ca_angle(n));
		p.set_atom_pos(i, Atom_N, n_pos),
		p.set_atom_pos(i, Atom_CA, ca_pos),
		p.set_atom_pos(i, Atom_C, c_pos);

		if (!p.is_glycine(i))
		{
			p.set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
		}

		// calculate where next N would be to get position for O

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
	}

	Point ca_pos = p.atom_pos(i + 1, Atom_CA);
	Point n_pos = p.atom_pos(i + 1, Atom_N);
	Point c_pos = p.atom_pos(i, Atom_C);

	for ( ;n >= 0;n--, i--)
	{
		ca_pos = torsion_to_coord(ca_pos, n_pos, c_pos, BOND_LENGTH_C_C,
			f->c_angle(n), f->omega(n), BOND_LENGTH_C_N);
		p.set_atom_pos(i, Atom_CA, ca_pos);
		p.conf().set_omega(i, f->omega(n));

		p.set_atom_pos(i, Atom_O, estimate_O_pos(ca_pos, c_pos, n_pos));

		n_pos = torsion_to_coord(n_pos, c_pos, ca_pos, BOND_LENGTH_N_CA,
			f->ca_angle(n), f->psi(n), BOND_LENGTH_C_C);
		p.set_atom_pos(i, Atom_N, n_pos);
		p.conf().set_psi(i, f->psi(n));

		if (!p.is_glycine(i))
		{
			p.set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
		}

		if (i > 0)
		{
			c_pos = torsion_to_coord(c_pos, ca_pos, n_pos, BOND_LENGTH_C_N,
				f->n_angle(n), f->phi(n), BOND_LENGTH_N_CA);
			p.set_atom_pos(i - 1, Atom_C, c_pos);
			p.conf().set_phi(i, f->phi(n));
		}
	}

	if (realign_before)
	{
		Point new_C = p.atom_pos(p_start_index - 1, Atom_C);
		Point new_N = p.atom_pos(p_start_index, Atom_N);
		Point new_CA = p.atom_pos(p_start_index, Atom_CA);

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

		Transform t;
		t.find_alignment(old_C, old_N, old_CA, new_C, new_N, new_CA);

		for (i = 0;i < p_start_index - 1;i++)
		{
			for (int a = 0;a < Num_Backbone;a++)
			{
				if (!(p.is_glycine(i) && (Atom_Id) a == Atom_CB))
				{
					p.transform_pos(i, (Atom_Id) a, t);
				}
			}
		}

		// i equals (p_start_index - 1)

		p.transform_pos(i, Atom_N, t);
		p.transform_pos(i, Atom_CA, t);
		// Atom_C has already been changed

		if (!p.is_glycine(i))
		{
			p.set_atom_pos(i, Atom_CB,
				estimate_CB_pos(
					p.atom_pos(i, Atom_CA),
					p.atom_pos(i, Atom_N),
					p.atom_pos(i, Atom_C)));
		}

		p.set_atom_pos(i, Atom_O,
			estimate_O_pos(
				p.atom_pos(i, Atom_CA),
				p.atom_pos(i, Atom_C),
				p.atom_pos(i + 1, Atom_N)));
	}
}

Fragment *Mover_Fragment_Fwd::add_fragment(int start_pos, int length)
{
	int end_pos = start_pos + length - 1;

	if (m_fragment.size() <= (unsigned) end_pos)
	{
		m_fragment.resize(end_pos + 100);
	}

	int num = (int) m_fragment[end_pos].size();
	m_fragment[end_pos].resize(num + 1);

	return &m_fragment[end_pos][num];
}

void Mover_Fragment_Fwd::after_fragments_loaded(int /*c_terminus*/)
{
	init_start_fragments();
	init_distributions();
}

void Mover_Fragment_Fwd::init_start_fragments()
{
	m_start_fragment.clear();
	m_first_end_pos = -1;

	// find all fragments that start at position 0
	// (ie. end at position "pos" and have length "pos + 1")

	for (int pos = 0;pos < (int) m_fragment.size();pos++)
	{
		if (m_first_end_pos == -1 && m_fragment[pos].size() != 0)
		{
			m_first_end_pos = pos;
		}

		for (int n = 0;n < (int) m_fragment[pos].size();n++)
		{
			if (m_fragment[pos][n].length() == pos + 1)
			{
				m_start_fragment.push_back(&m_fragment[pos][n]);
			}
		}
	}

	if (m_start_fragment.size() == 0)
	{
		std::cerr << "Error: no fragments found at first position\n";
		exit(1);
	}

	/*
	std::cerr << "(" << m_start_fragment.size()
		<< " start fragments)\n";
	*/
}

void Mover_Fragment_Fwd::init_distributions()
{
	unsigned int n, pos;

	m_frag_distrib.resize(m_fragment.size());

	for (pos = m_first_end_pos;pos < m_fragment.size();pos++)
	{
		for (n = 0;n < m_fragment[pos].size();n++)
		{
			double prob = m_fragment[pos][n].score();
			assert(prob >= Fragment::Min_Score);
			m_frag_distrib[pos].add(prob, n);
		}
	}

	for (n = 0;n < m_start_fragment.size();n++)
	{
		double prob = m_start_fragment[n]->score();
		assert(prob >= Fragment::Min_Score);
		m_start_frag_distrib.add(prob, n);
	}
}

void Mover_Fragment_Fwd::dump()
{
	std::cout << "Dumping fragments:\n\n";

	for (unsigned int pos = 0;pos < m_fragment.size();pos++)
	{
		std::cout << m_fragment[pos].size()
			<< " fragments ending at " << pos << "\n\n";

		assert((int) m_fragment[pos].size() == m_frag_distrib[pos].num());

		for (unsigned int n = 0;n < m_fragment[pos].size();n++)
		{
			const Fragment &f = m_fragment[pos][n];
			std::cout << (int) pos - f.length() + 1
				<< " to " << pos
				<< " (" << f.score() << ")\n";
		}
	}

	std::cout << "\n" << m_start_fragment.size()
		<< " fragments starting at 0\n\n";

	assert((int) m_start_fragment.size() == m_start_frag_distrib.num());

	for (unsigned int n = 0;n < m_start_fragment.size();n++)
	{
		const Fragment &f = *(m_start_fragment[n]);
		std::cout << "0 to "
			<< f.length() - 1
			<< " (" << f.score() << ")\n";
	}
}

