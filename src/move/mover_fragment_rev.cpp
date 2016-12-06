
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

// DONE
Fragment *Mover_Fragment_Rev::random_fragment(int min_start_pos, int *start_pos)
{
	assert(min_start_pos <= m_last_start_pos);
	int num;
	int count = 0;

	if (min_start_pos < 0)
	{
		min_start_pos = 0;
	}

	int num_possible = m_last_start_pos - min_start_pos + 1;

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

// DONE
void Mover_Fragment_Rev::init_sequential(Peptide &p, int initial_length,
	Run_Observer *observer)
{
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
	change_angles(p, p.end(), f);

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
	change_angles(p, p.end(), f);

	if (p.length() < initial_length)
	{
		extend(p, initial_length - p.length(),
			true,	// ribosome wall -- no harm done if it isn't actually true
			observer);
	}
}

// DONE
void Mover_Fragment_Rev::do_random_move(Peptide &p, int num,
	bool /*exhaustive_for_pos*/, Conf_Vec &result, Run_Observer *observer)
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

// DONE
void Mover_Fragment_Rev::do_random_move(Peptide &p,
	Run_Observer *observer)
{
	assert(m_fragments_loaded);

	int start, end;
	Fragment *f = random_fragment(p.start(), &start);
	end = start + f->length() - 1;
	change_angles(p, end, f);

	if (m_double_replacement_prob != 0.0)
	{
		if (m_double_replacement_prob == 1.0 ||
			Random::rnd(1.0) < m_double_replacement_prob)
		{
			f = random_fragment(p.start(), &start);
			end = start + f->length() - 1;
			change_angles(p, end, f);
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
	change_angles(p, p_end_index, f);

	if (ribosome_wall && !p.full_grown())
	{
		reorient_for_ribosome(p);
	}

	//p.conf().verify_torsion_angles();
}

void Mover_Fragment_Rev::change_angles(Peptide &p, int p_end_index,
	const Fragment *f)
{
	assert(p_end_index < p.full_length());
	assert(p_end_index - f->length() + 1 >= p.start());

	bool realign_after = (p_end_index < p.end());
	Point old_C, old_N, old_CA;
	double old_n_angle = 0.0;

	if (realign_after)
	{
		old_CA = p.atom_pos(p_end_index + 1, Atom_CA);
		old_N = p.atom_pos(p_end_index + 1, Atom_N);
		old_C = p.atom_pos(p_end_index, Atom_C);
		old_n_angle = angle_formed(old_CA, old_N, old_C);

/*
double curr_phi = torsion_angle(
	p.atom_pos(p_end_index, Atom_C),
	p.atom_pos(p_end_index + 1, Atom_N),
	p.atom_pos(p_end_index + 1, Atom_CA),
	p.atom_pos(p_end_index + 1, Atom_C));
std::cout << "ORIG phi: " << curr_phi << "  " << p.conf().phi(p_end_index + 1) << "\n";
assert(true && approx_equal_angle(curr_phi, p.conf().phi(p_end_index + 1)));
*/

	}

	int n = 0;								// position in fragment
	int i = p_end_index - f->length() + 1;	// position in p

//std::cout << "CHANGE_ANGLES: i = " << i << ", p.start() = " << p.start() << ", p.end() = " << p.end() << ", realign_after = " << (realign_after ? "true" : "false") << "\n";

	if (i == p.start())
	{
		Point n_pos, ca_pos, c_pos, pos;
		get_initial_ideal(&n_pos, &ca_pos, &c_pos, f->ca_angle(0));
		p.set_atom_pos(i, Atom_N, n_pos),
		p.set_atom_pos(i, Atom_CA, ca_pos),
		p.set_atom_pos(i, Atom_C, c_pos);

//std::cout << "INITIAL IDEAL:" << " N = " << n_pos
//<< " CA= " << ca_pos << " C = " << c_pos << "\n";

		if (!p.is_glycine(i))
		{
			p.set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
		}

		Point next_n_pos = torsion_to_coord(n_pos, ca_pos, c_pos,
			BOND_LENGTH_C_N, f->c_angle(0), f->psi(0), BOND_LENGTH_C_C);
		p.set_atom_pos(i, Atom_O, estimate_O_pos(ca_pos, c_pos, next_n_pos));

		p.set_atom_pos(i + 1, Atom_N, next_n_pos);
		p.conf().set_phi(i, f->phi(0));
		p.conf().set_psi(i, f->psi(0));

		n++;
		i++;
	}

	Point ca_pos = p.atom_pos(i - 1, Atom_CA);
	Point c_pos = p.atom_pos(i - 1, Atom_C);
	Point n_pos = p.atom_pos(i, Atom_N);

	for ( ;n < f->length();n++, i++)
	{
//std::cout << "@ n = " << n << " i = " << i << "\n";
		double omega_val;

		if (n > 0)
		{
			omega_val = f->omega(n - 1);
			p.conf().set_omega(i - 1, f->omega(n - 1));
		}
		else
		{
			omega_val = p.conf().omega(i - 1);
		}

		ca_pos = torsion_to_coord(ca_pos, c_pos, n_pos, BOND_LENGTH_N_CA,
			f->n_angle(n), omega_val, BOND_LENGTH_C_N);
		p.set_atom_pos(i, Atom_CA, ca_pos);

		c_pos = torsion_to_coord(c_pos, n_pos, ca_pos, BOND_LENGTH_C_C,
			f->ca_angle(n), f->phi(n), BOND_LENGTH_N_CA);
		p.set_atom_pos(i, Atom_C, c_pos);
		p.conf().set_phi(i, f->phi(n));

		if (!p.is_glycine(i))
		{
			p.set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
		}

		if (i < p.end())
		{
			n_pos = torsion_to_coord(n_pos, ca_pos, c_pos, BOND_LENGTH_C_N,
				f->c_angle(n), f->psi(n), BOND_LENGTH_C_C);
			p.set_atom_pos(i + 1, Atom_N, n_pos);
			p.conf().set_psi(i, f->psi(n));

			p.set_atom_pos(i, Atom_O, estimate_O_pos(ca_pos, c_pos, n_pos));
		}
		else
		{
			Point est_n_pos = torsion_to_coord(n_pos, ca_pos, c_pos,
				BOND_LENGTH_C_N, BOND_ANGLE_CA_C_N,
				deg2rad(180.0), BOND_LENGTH_C_C);

			 p.set_atom_pos(i, Atom_O,
			 	estimate_O_pos(ca_pos, c_pos, est_n_pos));
		}
	}

	if (realign_after)
	{
//std::cout << "## realign_after: p_end_index = " << p_end_index
//<< ", end = " << p.end() << "\n";
		p.conf().set_omega(p_end_index, f->omega(f->length() - 1));

		ca_pos = torsion_to_coord(ca_pos, c_pos, n_pos, BOND_LENGTH_N_CA,
			old_n_angle, f->omega(f->length() - 1), BOND_LENGTH_C_N);
		p.set_atom_pos(p_end_index + 1, Atom_CA, ca_pos);

		Point new_CA = ca_pos;
		Point new_N = p.atom_pos(p_end_index + 1, Atom_N);
		Point new_C = p.atom_pos(p_end_index, Atom_C);

		// Because the bond lengths are ideal, t.find_alignment() will
		// align the old CA and N atoms exactly onto the new positions.
		// The C atom will be placed in a position which causes
		// the old phi torsion angle (for residue p_end_index + 1)
		// to be retained.
		//
		// The overall effect is:
		//
		// The torsion and bond angles for residue (p_end_index + 1)
		// are the same as the old values.
		// The torsion and bond angles for residue (p_end_index)
		// come from the fragment.

		Transform t;
		t.find_alignment(old_CA, old_N, old_C, new_CA, new_N, new_C);
		//t.find_alignment(old_C, old_N, old_CA, new_C, new_N, new_CA);

		for (i = p_end_index + 2;i <= p.end();i++)
		{
			for (int a = 0;a < Num_Backbone;a++)
			{
				if (!(p.is_glycine(i) && (Atom_Id) a == Atom_CB))
				{
					p.transform_pos(i, (Atom_Id) a, t);
				}
			}
		}

		i = p_end_index + 1;
		p.transform_pos(i, Atom_C, t);

/*
double curr_phi = torsion_angle(
	p.atom_pos(p_end_index, Atom_C),
	p.atom_pos(p_end_index + 1, Atom_N),
	p.atom_pos(p_end_index + 1, Atom_CA),
	p.atom_pos(p_end_index + 1, Atom_C));
std::cout << "phi: " << curr_phi << "  " << p.conf().phi(p_end_index + 1) << "\n";
assert(approx_equal_angle(curr_phi, p.conf().phi(p_end_index + 1)));
*/

		// Atom_N and Atom_CA have already been changed

		if (!p.is_glycine(i))
		{
			p.set_atom_pos(i, Atom_CB,
				estimate_CB_pos(
					p.atom_pos(i, Atom_CA),
					p.atom_pos(i, Atom_N),
					p.atom_pos(i, Atom_C)));
		}

		Point next_n_pos;

		if (i < p.end())
		{
			next_n_pos = p.atom_pos(i + 1, Atom_N);
		}
		else
		{
			next_n_pos = torsion_to_coord(
				p.atom_pos(i, Atom_N),
				p.atom_pos(i, Atom_CA),
				p.atom_pos(i, Atom_C),
				BOND_LENGTH_C_N, BOND_ANGLE_CA_C_N,
				deg2rad(180.0), BOND_LENGTH_C_C);
		}

		p.set_atom_pos(i, Atom_O,
			estimate_O_pos(
				p.atom_pos(i, Atom_CA),
				p.atom_pos(i, Atom_C),
				next_n_pos));
	}
	
/*
	std::cout << "Verify after change angles: start = "
		<< p.start() << " end = " << p.end()
		<< " frag = " << p_end_index - f->length() + 1 
		<< " to " << p_end_index << std::endl;
	//p.conf().verify_torsion_angles();
	//p.verify();
*/
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
	m_last_start_pos = -1;

	// find all fragments that end at the C terminus
	// (ie. end at position "pos" and have length "pos + 1")

	for (int pos = (int) m_fragment.size() - 1;pos >= 0;pos--)
	{
		if (m_last_start_pos == -1 && m_fragment[pos].size() != 0)
		{
			m_last_start_pos = pos;
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

	for (pos = 0;pos <= (unsigned) m_last_start_pos;pos++)
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

