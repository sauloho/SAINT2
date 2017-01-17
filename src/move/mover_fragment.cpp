
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <iostream>
#include <sstream>		// for std::ostringstream
#include "mover.h"
#include "mover_fragment.h"
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

// static data members
const char *Mover_Fragment::c_type = "fragment";
const char *Mover_Fragment::c_param_lib = "lib";
const char *Mover_Fragment::c_param_double_replacement_prob =
	"double_replacement_prob";

Mover_Fragment::Mover_Fragment()
	: m_double_replacement_prob(0.0), m_fragments_loaded(false),
	  m_build_from_pos(-1)
{
}

Mover_Fragment::~Mover_Fragment()
{
}

void Mover_Fragment::set_library(const std::string &lib)
{
	if (lib != m_lib)
	{
		m_fragments_loaded = false;
	}

	m_lib = lib;
}

bool Mover_Fragment::parse_parameter(const std::string &name,
	const std::string &value)
{
	if (name == c_param_lib)
	{
		set_library(value);
		return true;
	}

	if (name == c_param_double_replacement_prob)
	{
		std::string full_name = Mover_Fragment::config_section();
		full_name += " ";
		full_name += c_type;
		full_name += " ";
		full_name += name;

		m_double_replacement_prob = parse_double(value, full_name, 0.0);
		return true;
	}

	return false;
}

void Mover_Fragment::verify_parameters()
{
    if (m_lib.empty())
    {
        Config::missing_parameter(Mover::config_section(),
			c_type, c_param_lib);
    }
}

void Mover_Fragment::init_non_sequential(Peptide &p, bool random_coil,
	Run_Observer *observer)
{
	load_fragments(observer);

	if (random_coil)
	{
		std::cerr << "Mover_Fragment::init_non_sequential(): "
			"random coil not implemented\n";
		exit(1);
	}

	p.set_length(p.full_length());
	int i = p.length() - 1;

	const double phi_val = deg2rad(180.0);
	const double psi_val = phi_val;
	const double omega_val = phi_val;

	Point n_pos, ca_pos, c_pos, pos;
	get_initial_ideal(&n_pos, &ca_pos, &c_pos);
	p.set_atom_pos(i, Atom_N, n_pos);
	p.set_atom_pos(i, Atom_CA, ca_pos);
	p.set_atom_pos(i, Atom_C, c_pos);

	if (!p.is_glycine(i))
	{
		p.set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
	}

	// calculate where next N would be to get position for O

	Point est_n_pos = torsion_to_coord(n_pos, ca_pos, c_pos,
		BOND_LENGTH_C_N, BOND_ANGLE_CA_C_N, psi_val, BOND_LENGTH_C_C);
	p.set_atom_pos(i, Atom_O, estimate_O_pos(ca_pos, c_pos, est_n_pos));
	p.conf().set_phi(i, phi_val);

	// create fully extended using ideal bond lengths and angles 

	for (i--;i >= 0;i--)
	{
		c_pos = torsion_to_coord(c_pos, ca_pos, n_pos,
			BOND_LENGTH_C_N, BOND_ANGLE_C_N_CA, phi_val, BOND_LENGTH_N_CA);
		p.set_atom_pos(i, Atom_C, c_pos);

		ca_pos = torsion_to_coord(ca_pos, n_pos, c_pos,
			BOND_LENGTH_C_C, BOND_ANGLE_CA_C_N, omega_val, BOND_LENGTH_C_N);
		p.set_atom_pos(i, Atom_CA, ca_pos);

		p.set_atom_pos(i, Atom_O, estimate_O_pos(ca_pos, c_pos, n_pos));

		n_pos = torsion_to_coord(n_pos, c_pos, ca_pos,
			BOND_LENGTH_N_CA, BOND_ANGLE_N_CA_C, psi_val, BOND_LENGTH_C_C);
		p.set_atom_pos(i, Atom_N, n_pos);

		if (!p.is_glycine(i))
		{
			p.set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
		}

		p.conf().set_psi(i, psi_val);
		p.conf().set_omega(i, omega_val);

		if (i > 0)
		{
			p.conf().set_phi(i, phi_val);
		}
	}

	p.verify_ideal();
}

void Mover_Fragment::init_from_peptide(Peptide &p, Run_Observer *observer)
{
	load_fragments(observer);
	//p.verify_ideal_bond_lengths();
}

void Mover_Fragment::reorient_for_ribosome(Peptide &p)
{
	Point zero(0.0, 0.0, 0.0);

	if (reverseSaint)
	{
		assert(p.atom_pos(p.start(), Atom_CA).close_to(zero));
	}
	else
	{
		assert(p.atom_pos(p.length() - 1, Atom_CA).close_to(zero));
	}

	// find average C alpha position

	Point a(0.0, 0.0, 0.0);

	int n;
	for (n = p.start();n <= p.end();n++)
	{
		a.add(p.atom_pos(n, Atom_CA));
	}

	a.scale(1.0 / (p.end() - p.start() + 1.0));

	if (a.distance(zero) < 1.0)
	{
		// no reason to perform any particular rotation, the centre of
		// gravity is very close to the most recently extruded residue
		return;
	}

	// any old point not colinear with "zero" and "a"
	Point b(a.x, a.y * 2, a.z * 3);

	// find transformation which keeps the most recently extruded residue
	// at (0, 0, 0) and puts "a" on the -x axis

	Transform t;
	t.find_alignment(
		zero, a, b,
		zero, Point(1.0, 0.0, 0.0), Point(0.0, 0.0, 1.0));

	if (!t.orthonormal())
	{
		// error, but keep going
		std::cout
			<< "@Error in Mover_Fragment::reorient_for_ribosome():"
				" transformation = " << t << "\n";
		exit(1);
		return;
	}

	for (n = p.start();n <= p.end();n++)
	{
		Residue &r = p.res(n);

		for (int i = 0;i < r.num_atoms();i++)
		{
			Atom_Id atom_id = r.atom(i).type().type();

			if (atom_id != Atom_Undef)
			{
				p.transform_pos(n, atom_id, t);
			}
		}
	}
}

void Mover_Fragment::print_template(std::ostream &out,
	bool commented /*= true */)
{
	const char *c = (commented ? "#" : "");

    out << c << "type = " << c_type << "\n"
        << c << c_param_lib << " = ...\t\t\t# fragment library location\n"
        << "\n";
}

void Mover_Fragment::load_fragments(Run_Observer *observer)
{
	if (m_fragments_loaded)
	{
		return;
	}

	// (checked in verify_parameters())
	assert(!m_lib.empty());

	observer->msg(std::string("Loading fragments from ") + m_lib);

	// file format is like:
	//
	// F <fragment#> P <pos> L <length> S <score> # <info>
	// 0 <phi> <psi> <omega> <N angle> <CA angle> <C angle>
	// 1 <phi> <psi> <omega> <N angle> <CA angle> <C angle>
	// ...
	// F ...
	
	C_File file(m_lib.c_str(), "r", "Fragment library file");

	static const int Max_Len = 1000;
	char buffer[Max_Len];
	int frag_num = -1;
	int c_terminus = -1;

	for (int f = 0;file.next_line(buffer, Max_Len);f++)
	{
		int frag_pos, frag_len;
		double frag_score;

		frag_num++;

		if (sscanf(buffer, "F %d P %d L %d S %lf",
			&frag_num, &frag_pos, &frag_len, &frag_score) != 4 ||
			frag_num != f)
		{
			std::cerr << "Error on line "
				<< file.line_num()
				<< " of fragment library file "
				<< m_lib
				<< ": expected \"F " << f
				<< " P <pos> L <len> S <score>\"\n";
			exit(1);
		}

		Fragment *frag = add_fragment(frag_pos, frag_len);

		if (c_terminus == -1 || (frag_pos + frag_len - 1 > c_terminus))
		{
			c_terminus = frag_pos + frag_len - 1;
		}

		if (frag_score < Fragment::Min_Score)
		{
			std::cerr << "Error on line "
				<< file.line_num()
				<< " of fragment library file "
				<< m_lib
				<< ": illegal score value "
				<< frag_score
				<< " (minimum is "
				<< Fragment::Min_Score
				<< ")\n";
			exit(1);
		}

		frag->set_score(frag_score);
		
		for (int n = 0;n < frag_len;n++)
		{
			if (!file.next_line(buffer, Max_Len))
			{
				std::cerr << "Error: unexpected end of file in fragment "
					"library file "
					<< m_lib << "\n";
				exit(1);
			}
			
			int index;
			double phi, psi, omega, n_angle, ca_angle, c_angle;

			if (sscanf(buffer, "%d %lf %lf %lf %lf %lf %lf", &index,
				&phi, &psi, &omega, &n_angle, &ca_angle, &c_angle) != 7 ||
				index != n)
			{
				std::cerr << "Error on line "
					<< file.line_num()
					<< " of fragment library file "
					<< m_lib
					<< ": expected \"" << n
					<< " <phi> <psi> <omega> <N angle> <CA angle> <C angle>\n";
				exit(1);
			}

			// (some values may be 999 ie. undefined, but no harm in
			// calling deg2rad on them, since they will never be used)

			frag->add(deg2rad(phi), deg2rad(psi), deg2rad(omega),
				deg2rad(n_angle), deg2rad(ca_angle), deg2rad(c_angle));
		}
	}

	m_fragments_loaded = true;
	after_fragments_loaded(c_terminus);

	std::ostringstream msg;
	msg << frag_num << " fragments read";

	observer->msg(msg.str());
}

