
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <sys/stat.h>	// (for stat())
#include <sys/types.h>	// (for stat())
#include "peptide.h"
#include "sequence.h"
#include "c_file.h"
#include "temp_file.h"
#include "stream_printf.h"
#include "pdb_atom_rec.h"
#include "atom.h"
#include "parse.h"
#include "rmsd.h"
#include "geom.h"

// static member variables
//std::string Peptide::m_scwrl_exec;

Peptide::Peptide() :
	m_length(0), m_full_length(0), m_conf(this), m_chain('A'),
	m_from_pdb(false)
{
}

Peptide::~Peptide()
{
	// do not delete m_seq (not owner)
}

Peptide::Peptide(const Peptide &other)
{
	// assert(!"called Peptide copy constructor!");

	copy(other);
}

Peptide& Peptide::operator = (const Peptide &other)
{
	// assert(!"Peptide assignment!");

	copy(other);
	return *this;
}

void Peptide::copy(const Peptide &other)
{
	// assert(!"called Peptide::copy()!");

	m_length = other.m_length;
	m_full_length = other.m_full_length;
	m_res = other.m_res;	// copies whole vector
	m_conf = other.m_conf;
	m_conf.set_peptide(this);
	m_chain = other.m_chain;
	m_from_pdb = other.m_from_pdb;
	m_filename = other.m_filename;
}

void Peptide::clear()
{
	m_length = m_full_length = 0;
	m_conf.clear();
	m_chain = 'A';
	m_from_pdb = false;
	
	// (using a temporary variable like this reduces the vector's
	// capacity to 0; just using erase() would not deallocate any memory)
	Residue_Vec empty_vec;
	m_res.swap(empty_vec);
}

void Peptide::alloc_satisfied_con(int length)
{
    m_satisfied_con = (int *)malloc(sizeof(int)*(length+1));
}

void Peptide::set_satisfied_con(int pos,int value)
{
    m_satisfied_con[pos]=value;
}
/***
void Peptide::remove_residue(int n)
{
	assert(n >= 0 && n < (int) m_res.size());
	assert(m_full_length = m_res.size());

	if (n > 0)
	{
		m_res[n-1].set_missing_after(true);
	}

	for ( ;n < m_full_length - 1;n++)
	{
		m_res[n] = m_res[n+1];
	}

	m_full_length--;
	m_res.resize((unsigned) m_full_length);
}
***/

bool Peptide::full_grown() const
{
	// (don't return true if the peptide has not been initialised yet
	// (ie. if m_full_length is 0))

	return ((m_length == m_full_length) && m_full_length != 0);
}

void Peptide::set_length(int len)
{
	assert(len >= 0 && len <= m_full_length);
	m_length = len;
}

void Peptide::add_length(int amount)
{
	set_length(m_length + amount);
}

void Peptide::create_from_sequence(const Sequence &seq)
{
	m_full_length = seq.length();
	m_length = 0;
	m_res.resize(m_full_length);
	m_from_pdb = false;

	int n;
	for (n = 0;n < m_full_length;n++)
	{
		m_res[n].set_amino(seq.amino(n));
		m_res[n].allocate_backbone_atoms();
	}

	if (seq.codons_known())
	{
		for (n = 0;n < m_full_length;n++)
		{
			m_res[n].set_codon(seq.codon(n));
		}
	}

	m_conf.set_num_res(m_full_length);
}

bool Peptide::read_pdb(const char *filename, char chain /*=' '*/,
	bool no_warnings /*=false*/)
{
	// reset all data in the Peptide
	clear();

	m_from_pdb = true;

	bool use_first_chain = (chain == ' ');
	bool found_chain = false;
	bool removed_first_res = false;

	// (exits with a message if the file cannot be opened)
	const char *desc = "PDB file";
	C_File file(filename, "r", desc);

	m_filename = filename;

	// residue currently being read
	Residue *res = NULL;

	// index of current residue
	int res_index = -1;

	// sequence number and i-code of the current residue
	int curr_res_seq = -99999;
	char curr_i_code = ' ';

	const int Max_Len = 10000;
	char buffer[Max_Len];
	int line_num = 0;

	while (fgets(buffer, Max_Len, file))
	{
		line_num++;
		std::string record(buffer, 6);

		if (record == "ENDMDL")
		{
			// found end of model
			break;
		}

		// check for terminal residue

		if (record == "TER   " && found_chain)
		{
			// treat the line the same way as an "ATOM" record;
			// some things will be missing (eg. xyz positions),
			// but it can still parse the line partially
			// (to get the chain, res_seq and i_code)

			PDB_Atom_Rec rec;
			std::string err_msg;

			rec.parse(buffer, &err_msg);

			if (rec.chain_id == chain &&
				rec.res_seq == curr_res_seq &&
				rec.i_code == curr_i_code)
			{
				// found terminal residue
				break;
			}
		}

		if (record == "ATOM  ")
		{
			PDB_Atom_Rec rec;
			std::string err_msg;
			
			if (!rec.parse(buffer, &err_msg))
			{
				if (!no_warnings)
				{
					std::cerr << "Warning: "
						<< err_msg
						<< " on line " << line_num
						<< " of " << filename
						<< "\n";
				}

				continue;
			}

			// ignore terminal oxygen
			if (rec.name == " OXT")
			{
				continue;
			}

			// if chain was not specified, set it to the first one found

			if (chain == ' ')
			{
				chain = rec.chain_id;
			}

			if (rec.chain_id == chain)
			{
				found_chain = true;

				// extract the amino acid type
				Amino aa(rec.res_name.c_str());

				if (!aa.standard())
				{
					if (!no_warnings)
					{
						std::cerr
							<< "Warning: ignoring unknown amino acid type "
							<< rec.res_name
							<< " on line "
							<< line_num
							<< " of "
							<< filename
							<< "\n";
					}

					continue;
				}

				// check if a new residue has been started

				if (!(rec.res_seq == curr_res_seq &&
					  rec.i_code == curr_i_code))
				{
// std::cout << "NEW RES " << rec.res_seq << "\n";
					// check for missing residues (gaps in sequence numbering)
					//
					// (Note that two residues in a row can have the same
					// sequence number but different i-codes)

					bool id_gap = (curr_res_seq != -99999 &&
						rec.res_seq != curr_res_seq &&
						rec.res_seq != curr_res_seq + 1 &&
						res != NULL);

					bool remove_init_res = false;
					bool at_init_res =
						(m_res.size() == 1 && !removed_first_res);
					const char *remove_msg;

					if (id_gap)
					{
						if (at_init_res)
						{
							remove_init_res = true;
							remove_msg = " with immediate gap";
						}
						else
						{
							res->set_missing_after(true);
						}
					}
					else
					if (at_init_res &&
						!(res->atom_exists(Atom_N) &&
			  			  res->atom_exists(Atom_CA) &&
			  			  res->atom_exists(Atom_C)))
					{
						remove_init_res = true;
						remove_msg = " with missing atoms";
					}

					if (remove_init_res)
					{
						if (!no_warnings)
						{
							std::cout
								<< "Warning: ignoring initial residue "
								<< m_res[0].res_seq_str()
								<< remove_msg
								<< " on line "
								<< m_res[0].pdb_line()
								<< " of "
								<< filename
								<< "\n";
						}

						m_res.erase(m_res.begin());
						removed_first_res = true;
					}

					res_index = (int) m_res.size();
					m_res.resize(res_index + 1);
					m_conf.set_num_res(res_index + 1);

					curr_res_seq = rec.res_seq;
					curr_i_code = rec.i_code;

					res = &m_res[res_index];
					res->clear();
					res->set_amino(aa);
					res->set_res_seq(curr_res_seq, curr_i_code);
					res->set_pdb_line(line_num);
				}
				else   // same residue as last ATOM record
				{
					if (aa != res->amino())
					{
						/***
						if (!no_warnings)
						{
							std::cerr
								<< "Warning: ignoring multiple occupancy atom "
									"on line "
								<< line_num
								<< " of " << filename
								<< "\n";
						}
						***/

						continue;
					}
				}

				// add the new atom to the residue
				Atom_Type t(rec.name.c_str());
			
				if (t.undefined())
				{
					if (!no_warnings)
					{
						// probably hydrogen if contains an "H";
						// don't bother reporting
						// ("OXT" is terminal oxygen)

						if (strchr(rec.name.c_str(), 'H') == NULL &&
							rec.name != " OXT")
						{
							std::cerr << "Warning: ignoring unknown atom type "
								<< rec.name
								<< " on line " << line_num
								<< " of " << filename
								<< "\n";
						}
					}

					continue;
				}

				if (aa.rapdf_id(t) == -1)
				{
					if (!no_warnings)
					{
						std::cerr << "Warning: ignoring "
							<< rec.name
							<< " found in "
							<< aa.name()
							<< " on line " << line_num
							<< " of " << filename
							<< "\n";
					}
				}
				else
				if (!res->atom_exists(t))
				{
// std::cout << "ADD ATOM " << rec.name << " TO RES " << rec.res_seq << "\n";
					Atom *a = res->add_atom(t);
					a->set_pdb_line(line_num);
					set_atom_pos(res_index, t.type(), rec.pos);
				}
			}
		}
	}

	if (m_res.size() == 0)
	{
		if (use_first_chain)
		{
			std::cerr << "Error: no valid residues found in first chain in "
				<< filename << "\n";
		}
		else
		{
			std::cerr << "Error: chain '" << chain
				<< "' not found in " << filename
				<< "\n";
		}

		return false;
	}

	if (m_res.size() > 1)
	{
		assert(res != NULL);

		bool remove_final = false;
		const char *desc;

		if (m_res[m_res.size()-2].missing_after())
		{
			remove_final = true;
			desc = " after id gap";
		}
		else
		if (!(res->atom_exists(Atom_N) &&
			  res->atom_exists(Atom_CA) &&
			  res->atom_exists(Atom_C)))
		{
			remove_final = true;
			desc = " with missing atoms";
		}

		if (remove_final)
		{
			if (!no_warnings)
			{
				std::cout
					<< "Warning: ignoring final residue "
					<< res->res_seq_str()
					<< desc
					<< " on line "
					<< res->pdb_line()
					<< " of "
					<< filename
					<< "\n";
			}

			m_res.pop_back();
		}
	}

	if (m_res.size() > 1)
	{
		m_res[m_res.size()-1].set_missing_after(false);
	}

	m_chain = chain;
	m_length = m_full_length = (int) m_res.size();

	check_N_CA_C_positions(filename, no_warnings);
	check_CB_positions(filename, no_warnings);
	check_O_positions(filename, no_warnings);
	// check_side_chain_atoms(filename, no_warnings);
	if (!no_warnings) { print_missing_res(filename); }

	return true;
}

const char *Peptide::get_filename() const
{
	return m_filename.c_str();
}

inline bool outside_length_tolerance(double dist, double ideal_dist)
{
	return (dist < ideal_dist * 0.8 || dist > ideal_dist * 1.2);
}

inline bool outside_angle_tolerance(double angle, double ideal_angle)
{
	static const double AngleTolerance = deg2rad(20.0);

	return (fabs(angle - ideal_angle) > AngleTolerance);
}

struct IndexSet
{
	int n1, n2, n3;	// (n3 is -1 if unused)

	IndexSet(int v1, int v2, int v3)
		: n1(v1), n2(v2), n3(v3)
	{
	}

	bool operator == (const IndexSet &rhs)
	{
		return (n1 == rhs.n1 && n2 == rhs.n2 && n3 == rhs.n3);
	}
};

typedef std::list<IndexSet> IndexSetList;
typedef std::vector<IndexSetList> IndexSetListVec;

void remove_index_set(IndexSetList &the_list, const IndexSet &to_remove)
{
	for (IndexSetList::iterator i = the_list.begin();
		 i != the_list.end();++i)
	{
		if (*i == to_remove)
		{
			the_list.erase(i);
			return;
		}
	}

	assert(!"Couldn't find IndexSet in remove_index_set");
}

void Peptide::check_N_CA_C_positions(const char *filename, bool no_warnings)
{
	static const double bond_length_from[3] =
	{
		BOND_LENGTH_N_CA,
		BOND_LENGTH_C_C,	// from CA to C
		BOND_LENGTH_C_N
	};

	static const double bond_angle_from[3] =
	{
		BOND_ANGLE_N_CA_C,
		BOND_ANGLE_CA_C_N,
		BOND_ANGLE_C_N_CA
	};

	// maximum distance between two residues (C in 1st & N in 2nd)
	// before a "gap" is assumed and set_missing_after(true) is called
	static const double MaxResGap = BOND_LENGTH_C_N * 2.0;

	// this function depends on the AtomId values Atom_N, Atom_CA and Atom_C
	// being in consecutive order

	assert(Atom_CA == (Atom_Id) (Atom_N + 1));
	assert(Atom_C == (Atom_Id) (Atom_N + 2));

	int num = (int) m_res.size() * 3;

	// transfer all atom positions to a single array
	// (at the same time, check if any pair of residues is a long way
	// apart - if so, call set_missing_after(true) for the first one)

	bool *exists = new bool[num+2];		// +2 to avoid bounds checking later
	bool *gap_after = new bool[num+2];
	Point *pos = new Point[num];
	bool *bad = new bool[num];

	int n, a;
	for (n = a = 0;n < (int) m_res.size();n++)
	{
		Residue &r = m_res[n];

		for (Atom_Id t = Atom_N;t <= Atom_C;t = (Atom_Id) (t + 1), a++)
		{
			bool atom_exists = exists[a] = r.atom_exists(t);

			if (atom_exists)
			{
				pos[a] = atom_pos(n, t);
			}
			else
			{
				if (!no_warnings)
				{
					std::cout << "Warning: missing "
						<< Atom_Type(t).name()
						<< " in residue "
						<< r.res_seq_str()
						<< " on line "
						<< r.pdb_line()
						<< " of "
						<< filename
						<< "\n";
				}
			}

			if (t == Atom_C)
			{
				gap_after[a] = r.missing_after();

				// check if there is a large gap between the C
				// in this residue and the N in the next one

				if (!gap_after[a] &&
					atom_exists &&
					n < (int) m_res.size() - 1  &&
					m_res[n+1].atom_exists(Atom_N) &&
					pos[a].distance(atom_pos(n + 1, Atom_N)) > MaxResGap)
				{
					gap_after[a] = true;
					r.set_missing_after(true);

					if (!no_warnings)
					{
						std::cout << "Warning: large distance between residue "
							<< r.res_seq_str()
							<< " (" << r.amino().abbr()
							<< ") and "
							<< m_res[n+1].res_seq_str()
							<< " (" << m_res[n+1].amino().abbr()
							<< ") on line "
							<< m_res[n+1].pdb_line()
							<< " of " << filename
							<< "\n";
					}
				}
			}
			else
			{
				gap_after[a] = false;
			}
		}
	}

	assert(a == num);
	exists[num] = exists[num+1] = false;

	IndexSetListVec v;
	v.resize(num);

	for (a = 0;a < num;a++)
	{
		int i = a % 3;

		if (exists[a] && exists[a+1] && !gap_after[a])
		{
			double length = pos[a].distance(pos[a+1]);

			if (outside_length_tolerance(length, bond_length_from[i]))
			{
				IndexSet x(a, a + 1, -1);
				v[a].push_back(x);
				v[a+1].push_back(x);

				if (!no_warnings)
				{
					int nn = a / 3;
					std::cout << "Warning: wrong bond length "
						<< Atom_Type((Atom_Id) (i + Atom_N)).name()
						<< "-"
						<< Atom_Type((Atom_Id) ((i + 1) % 3 + Atom_N)).name()
						<< " = "
						<< Printf("%.2f", length)
						<< " (vs ideal "
						<< Printf("%.2f", bond_length_from[i])
						<< ") in residue "
						<< m_res[nn].res_seq_str()
						<< " (" << m_res[nn].amino().abbr()
						<< ") on line "
						<< m_res[nn].atom((Atom_Id) (i + Atom_N)).pdb_line()
						<< " of " << filename
						<< "\n";
				}
			}
		}

		if (exists[a] && exists[a+1] && exists[a+2] &&
			!gap_after[a] && !gap_after[a+1])
		{
			double angle = angle_formed(pos[a], pos[a+1], pos[a+2]);

			if (outside_angle_tolerance(angle, bond_angle_from[i]))
			{
				IndexSet x(a, a + 1, a + 2);
				v[a].push_back(x);
				v[a+1].push_back(x);
				v[a+2].push_back(x);

				if (!no_warnings)
				{
					int nn = a / 3;
					std::cout << "Warning: wrong bond angle "
						<< Atom_Type((Atom_Id) (i + Atom_N)).name()
						<< "-"
						<< Atom_Type((Atom_Id) ((i + 1) % 3 + Atom_N)).name()
						<< "-"
						<< Atom_Type((Atom_Id) ((i + 2) % 3 + Atom_N)).name()
						<< " = "
						<< Printf("%.1f", rad2deg(angle))
						<< " (vs ideal "
						<< Printf("%.1f", rad2deg(bond_angle_from[i]))
						<< ") in residue "
						<< m_res[nn].res_seq_str()
						<< " (" << m_res[nn].amino().abbr()
						<< ") near line "
						<< m_res[nn].pdb_line()
						<< " of " << filename
						<< "\n";
				}
			}
		}

		bad[a] = false;
	}

	// inefficient, but this function should rarely be called

	for (unsigned int k = 5;k >= 2;k--)
	{
		for (a = 0;a < num;a++)
		{
			if (v[a].size() >= k)
			{
				bad[a] = true;
			}
		}

		for (a = 0;a < num;a++)
		{
			if (bad[a] && v[a].size() != 0)
			{
				for (IndexSetList::iterator i = v[a].begin();
					 i != v[a].end();++i)
				{
					if (i->n1 != a)
					{
						remove_index_set(v[i->n1], *i);
					}

					if (i->n2 != a)
					{
						remove_index_set(v[i->n2], *i);
					}

					if (i->n3 != a && i->n3 != -1)
					{
						remove_index_set(v[i->n3], *i);
					}
				}

				v[a].erase(v[a].begin(), v[a].end());
			}
		}
	}

	for (a = 0;a < num;a++)
	{
		if (v[a].size() != 0)
		{
			IndexSetList::iterator i = v[a].begin();

			if (i->n3 == -1)
			{
				// can't tell if 1st or 2nd atom is bad
				bad[i->n1] = bad[i->n2] = true;
			}
			else
			{
				// assume the centre atom it the bad one
				bad[i->n2] = true;
			}
		}
	}

	for (a = 0;a < num;a++)
	{
		if (bad[a])
		{
			int n = a / 3;
			Residue &r = m_res[n];
			Atom_Id t = (Atom_Id) (Atom_N + (a % 3));

			if (!no_warnings)
			{
				std::cout << "=> Ignoring bad "
					<< Atom_Type(t).name()
					<< " in residue "
					<< r.res_seq_str()
					<< " ("
					<< r.amino().abbr()
					<< ") on line "
					<< r.atom(t).pdb_line()
					<< " of " << filename
					<< "\n";
			}

			r.remove_atom(t);
		}
	}

	delete [] exists;
	delete [] gap_after;
	delete [] pos;
	delete [] bad;
}

void Peptide::check_CB_positions(const char *filename, bool no_warnings)
{
	for (int n = 0;n < (int) m_res.size();n++)
	{
		Residue &r = m_res[n];

		if (!r.amino().is_glycine())
		{
			if (r.atom_exists(Atom_CA) &&
				r.atom_exists(Atom_C) &&
				r.atom_exists(Atom_N))
			{
				Point ca_pos = atom_pos(n, Atom_CA);
				Point n_pos = atom_pos(n, Atom_N);
				Point c_pos = atom_pos(n, Atom_C);

				// (should never happen, since bond lengths should have
				// already been checked by check_N_CA_C_positions())

				if (n_pos.close_to(ca_pos) ||
					c_pos.close_to(ca_pos))
				{
					if (!no_warnings)
					{
						std::cout
							<< "Warning: N or C position same as CA position"
							   " on line "
							<< r.atom(Atom_CA).pdb_line()
							<< " of " << filename
							<< "\n";
					}

					continue;
				}

				Point pred_cb = estimate_CB_pos(ca_pos, n_pos, c_pos);
			
				if (r.atom_exists(Atom_CB))
				{
					static const double DistTolerance = 0.5;
					double dist = pred_cb.distance(atom_pos(n, Atom_CB));

					if (dist > DistTolerance)
					{
						set_atom_pos(n, Atom_CB, pred_cb);

						if (!no_warnings)
						{
							std::cout
								<< "Warning: fixing incorrect CB position for "
									"residue "
								<< r.res_seq_str()
								<< " (" << r.amino().abbr()
								<< ") on line "
								<< r.atom(Atom_CB).pdb_line()
								<< " of "
								<< filename
								<< " (" << Printf("%.2f", dist)
								<< " Angstroms from expected position)\n";
						}
					}
				}
				else
				{
					Atom *a = r.add_atom(Atom_CB);
					a->set_pdb_line(r.pdb_line());
					set_atom_pos(n, Atom_CB, pred_cb);

					if (!no_warnings)
					{
						std::cout << "Warning: estimating missing CB position "
							"for residue "
							<< r.res_seq_str()
							<< " (" << r.amino().abbr()
							<< ") on line "
							<< r.pdb_line()
							<< " of "
							<< filename
							<< "\n";
					}
				}
			}
			else
			if (!r.atom_exists(Atom_CB))
			{
				if (!no_warnings)
				{
					std::cout << "Warning: missing CB in residue "
						<< r.res_seq_str()
						<< " on line "
						<< r.pdb_line()
						<< " of "
						<< filename
						<< "\n";
				}
			}
		}
	}
}

void Peptide::check_O_positions(const char *filename, bool no_warnings)
{
	// (don't let n = last residue, since this code references m_res[n+1])

	for (int n = start();n < end();n++)
	{
		Residue &r = m_res[n];

		if (r.atom_exists(Atom_CA) &&
			r.atom_exists(Atom_C))
		{
			Point ca_pos = atom_pos(n, Atom_CA);
			Point c_pos = atom_pos(n, Atom_C);

			// need the the N position from the next residue

			if (r.missing_after() || !m_res[n+1].atom_exists(Atom_N))
			{
				if (!r.atom_exists(Atom_O))
				{
					if (!no_warnings)
					{
						std::cout << "Warning: cannot fix missing O atom in "
							"residue "
							<< r.res_seq_str()
							<< " (" << r.amino().abbr()
							<< ") on line "
							<< r.pdb_line()
							<< " of " << filename
							<< "\n";
					}
				}
				else
				{
					// all we can do is check is the bond length

					Point o_pos = atom_pos(n, Atom_O);
					double dist = c_pos.distance(o_pos);

					if (outside_length_tolerance(dist, BOND_LENGTH_C_O))
					{
						if (!no_warnings)
						{
							std::cout
								<< "Warning: cannot fix O with incorrect bond "
								<< "length " << Printf("%.2f", dist)
								<< " (vs ideal "
								<< Printf("%.2f", BOND_LENGTH_C_O)
								<< ") in residue " << r.res_seq_str()
								<< " (" << r.amino().abbr()
								<< ") on line "
								<< r.atom(Atom_O).pdb_line()
								<< " of " << filename
								<< "; ignoring atom\n";
						}

						r.remove_atom(Atom_O);
					}
				}
			}
			else
			{
				Point n_pos = atom_pos(n + 1, Atom_N);
				Point pred = estimate_O_pos(ca_pos, c_pos, n_pos);

				if (r.atom_exists(Atom_O))
				{
					static const double DistTolerance = 0.5;
					Point o_pos = atom_pos(n, Atom_O);
					double dist = o_pos.distance(pred);

					if (dist > DistTolerance)
					{
						// replace with expected position
						set_atom_pos(n, Atom_O, pred);

						if (!no_warnings)
						{
							std::cout
								<< "Warning: fixing incorrect O position for "
									"residue "
								<< r.res_seq_str()
								<< " (" << r.amino().abbr()
								<< ") on line "
								<< r.atom(Atom_O).pdb_line()
								<< " of "
								<< filename
								<< " (" << Printf("%.2f", dist)
								<< " Angstroms from expected position)\n";
						}
					}
				}
				else
				{
					Atom *a = r.add_atom(Atom_O);
					a->set_pdb_line(r.pdb_line());
					set_atom_pos(n, Atom_O, pred);

					if (!no_warnings)
					{
						std::cout << "Warning: estimating missing O position "
							"for residue "
							<< r.res_seq_str()
							<< " (" << r.amino().abbr()
							<< ") on line "
							<< r.pdb_line()
							<< " of "
							<< filename
							<< "\n";
					}
				}
			}
		}
		else
		if (!r.atom_exists(Atom_O))
		{
			if (!no_warnings)
			{
				std::cout << "Warning: missing O in residue "
					<< r.res_seq_str()
					<< " on line "
					<< r.pdb_line()
					<< " of "
					<< filename
					<< "\n";
			}
		}
	}
}

void Peptide::check_side_chain_atoms(const char *filename, bool no_warnings)
{
	if (no_warnings)
	{
		return;
	}

	for (int n = 0;n < (int) m_res.size();n++)
	{
		Residue &r = m_res[n];

		for (int i = 0;i < r.amino().expected_atoms();i++)
		{
			Atom_Id j = r.amino().expected_atom(i);

			if (j < Num_Backbone)
			{
				continue;
			}

			if (!r.atom_exists(j))
			{
				std::cout
					<< "Warning: missing side chain atom "
					<< Atom_Type(j).name()
					<< " in residue "
					<< r.res_seq_str()
					<< " (" << r.amino().abbr()
					<< ") on line "
					<< r.pdb_line()
					<< " of "
					<< filename
					<< "\n";
			}
		}
	}
}

Residue *Peptide::append_residue()
{
	m_length = ++m_full_length;

	m_res.resize(m_res.size() + 1);
	return &(m_res.back());
}

void Peptide::write_pdb(const char *filename,
	bool backbone_only /*= false*/,
	bool atom_records_only /*= false*/) const
{
	C_File file(filename, "w", "PDB output file");
	m_filename = filename;

	char chain = m_chain;

	if (chain == ' ')
	{
		chain = 'A';
	}

	if (!atom_records_only)
	{
		// header

		fprintf(file, "HEADER                                            DD-MMM-YY   ????              \n");
		fprintf(file, "TITLE     SAINT 2 OUTPUT                                                        \n");
		// sequence

		// Note: PDB files can have residues in the SEQRES list that
		// do not appear in ATOM records; if this peptide was originally
		// read from a PDB file, the SEQRES lines output by this function
		// may not be identical to the original PDB file.

		int sr_line = 1;

		for (int n = start();n <= end();n++)
		{
			if (n % 13 == 0)
			{
				if (n != 0)
				{
					fprintf(file, "\n");
				}
		
				fprintf(file, "SEQRES %3d %c %4d ", sr_line, chain, m_length);
				sr_line++;
			}

			fprintf(file, " %s", m_res[n].amino().abbr());
		}

		fprintf(file, "\n");
	}

	// atoms

	int serial = 0;

	for (int n = start();n <= end();n++)
	{
		const Residue &r = m_res[n];

		for (int i = 0;i < r.num_atoms();i++)
		{
			const Atom &a = r.atom(i);

			if (a.undefined())
			{
				continue;
			}

			serial++;

			if (backbone_only && !a.type().is_backbone())
			{
				continue;
			}

			int res_seq;
			char i_code;

			if (m_from_pdb)
			{
				res_seq = r.res_seq();
				i_code = r.i_code();
			}
			else
			{
				res_seq = n + 1;
				i_code = ' ';
			}

			PDB_Atom_Rec rec;
			rec.set_values(
				serial,			// atom number
				a.type(),
				res_seq,
				i_code,
				r.amino(),
				atom_pos(n, a.id()),
				chain
			);

			char buffer[100];
			rec.format(buffer);

			fputs(buffer, file);
			putc('\n', file);
		}
	}

	// "ATOM

	fprintf(file, "END                                                                             \n");
}

/******************

bool Peptide::check_bond_values(Residue &r, Residue_Vec &r_vec, const char *filename,
	int line_num, bool no_warnings)
{
	// bond angle tolerance
	static const double AngleTolerance = deg2rad(20.0);

	assert(r.atom_exists(Atom_N));
	assert(r.atom_exists(Atom_CA));
	assert(r.atom_exists(Atom_C));

	Point n_pos = r.atom(Atom_N).pos();
	Point ca_pos = r.atom(Atom_CA).pos();
	Point c_pos = r.atom(Atom_C).pos();

	double dist_n_ca = n_pos.distance(ca_pos);

	if (dist_n_ca < BOND_LENGTH_N_CA * 0.5 || dist_n_ca > BOND_LENGTH_N_CA * 1.5)
	{
		if (!no_warnings)
		{
			std::cout << "Warning: Illegal N-CA bond length "
				<< Printf("%.2f", dist_n_ca)
				<< " (expected " << Printf("%.2f", BOND_LENGTH_N_CA)
				<< ") in " << r.amino().abbr() << " near line " << line_num
				<< " of " << filename << "\n";
		}

		return false;
	}

	double dist_ca_c = ca_pos.distance(c_pos);

	if (dist_ca_c < BOND_LENGTH_C_C * 0.5 || dist_ca_c > BOND_LENGTH_C_C * 1.5)
	{
		if (!no_warnings)
		{
			std::cout << "Warning: Illegal CA-C bond length "
				<< Printf("%.2f", dist_ca_c)
				<< " (expected " << Printf("%.2f", BOND_LENGTH_C_C)
				<< ") in " << r.amino().abbr() << " near line " << line_num
				<< " of " << filename << "\n";
		}

		return false;
	}

	double angle_n_ca_c = angle_formed(n_pos, ca_pos, c_pos);

	if (fabs(angle_n_ca_c - BOND_ANGLE_N_CA_C) > AngleTolerance)
	{
		if (!no_warnings)
		{
			std::cout << "Warning: Illegal N-CA-C bond angle "
				<< Printf("%.1f", rad2deg(angle_n_ca_c))
				<< " (expected " << Printf("%.1f", rad2deg(BOND_ANGLE_N_CA_C))
				<< ") in " << r.amino().abbr() << " near line " << line_num
				<< " of " << filename << "\n";
		}

		return false;
	}

	if (r_vec.size() != 0)
	{
		Residue &prev = r_vec[r_vec.size()-1];

		if (!prev.missing_after())
		{
			assert(prev.atom_exists(Atom_C));
			assert(prev.atom_exists(Atom_CA));

			Point prev_c_pos = prev.atom(Atom_C).pos();
			Point prev_ca_pos = prev.atom(Atom_CA).pos();
			double dist_c_n = prev_c_pos.distance(n_pos);

			if (dist_c_n < BOND_LENGTH_N_C * 0.5 || dist_c_n > BOND_LENGTH_N_C * 1.5)
			{
				if (!no_warnings)
				{
					std::cout << "Warning: Illegal C-N bond length "
						<< Printf("%.2f", dist_c_n)
						<< " (expected " << Printf("%.2f", BOND_LENGTH_N_C)
						<< ") in " << r.amino().abbr() << " near line " << line_num
						<< " of " << filename << "\n";
				}

				return false;
			}

			double angle_ca_c_n = angle_formed(prev_ca_pos, prev_c_pos, n_pos);

			if (fabs(angle_ca_c_n - BOND_ANGLE_CA_C_N) > AngleTolerance)
			{
				std::cout << "Warning: Illegal CA-C-N bond angle "
					<< Printf("%.1f", rad2deg(angle_ca_c_n))
					<< " (expected " << Printf("%.1f", rad2deg(BOND_ANGLE_CA_C_N))
					<< ") in " << r.amino().abbr() << " near line " << line_num
					<< " of " << filename << "\n";
				return false;
			}

			double angle_c_n_ca = angle_formed(prev_c_pos, n_pos, ca_pos);

			if (fabs(angle_c_n_ca - BOND_ANGLE_C_N_CA) > AngleTolerance)
			{
				std::cout << "Warning: Illegal C-N-CA bond angle "
					<< Printf("%.1f", rad2deg(angle_c_n_ca))
					<< " (expected " << Printf("%.1f", rad2deg(BOND_ANGLE_C_N_CA))
					<< ") in " << r.amino().abbr() << " near line " << line_num
					<< " of " << filename << "\n";
				return false;
			}
		}
	}

	return true;
}
******************/

void Peptide::remove_non_backbone_atoms()
{
	for (int n = 0;n < m_full_length;n++)
	{
		m_res[n].remove_non_backbone_atoms();
	}

	m_conf.remove_non_backbone_atoms();
}

bool Peptide::many_non_backbone() const
{
	int non_gly = 0;
	int total = 0;

	for (int n = 0;n < m_full_length;n++)
	{
		if (!is_glycine(n))
		{
			non_gly++;

			if (m_res[n].num_atoms() > Num_Backbone)
			{
				total++;
			}
		}
	}

	// return true if more than half of the non-glycine residues
	// have non-backbone atoms present

	return (total > non_gly / 2);
}

//#define DEBUG_CONTACT

double Peptide::previous_contact_proportion() const
{
	assert(!reverseSaint || (m_length == m_full_length));

	static const double Contact_Dist = 5.0;	// 4.0, 6.0
	static const int min_offset = 4;
	int total = 0;

	std::vector<bool> contact;
	std::vector<bool> helix;
	contact.resize(m_length);
	helix.resize(m_length);

	for (int i = 0;i < m_length;i++)
	{
		helix[i] = m_conf.in_helix(i);
		contact[i] = false;
	}

	// ignore helix near the beginning

	int start = min_offset;

	if (m_length > 5)
	{
		// find the longest group of helices near the beginning

		bool in_helix = false;
		int helix_start = -1;
		int longest_helix = 0;
		int longest_helix_end = -1;

		for (int x = min_offset;x < m_length;x++)
		{
			if (helix[x])
			{
				if (!in_helix)
				{
					in_helix = true;
					helix_start = x;
				}
			}
			else   // not in helix
			{
				if (in_helix)
				{
					int len = x - helix_start;

					if (len > longest_helix)
					{
						longest_helix = len;
						longest_helix_end = x - 1;
					}

					in_helix = false;
				}

				if (x > min_offset + 4)
				{
					break;
				}
			}
		}

		if (longest_helix > 1)
		{
			start = longest_helix_end + 1;
			total += (start - min_offset);

#ifdef DEBUG_CONTACT
			std::cout << "Included " << total
				<< " helix values up to start ("
				<< min_offset
				<< " to "
				<< start - 1<< ")\n";
#endif // DEBUG_CONTACT
		}
	}

	int n;
	for (n = start;n < m_length;n++)
	{
		int offset = min_offset;

		if (helix[n])
		{
#ifdef DEBUG_CONTACT
			std::cout << "Helix " << n << "\n";
#endif // DEBUG_CONTACT
			offset = 8;
		}

		// check if any previous residues are close to this one

		for (int i = 0;i < m_res[n].num_atoms();i++)
		{
			Atom n_atom = m_res[n].atom(i);
			if (n_atom.undefined()) { continue; }

			Point n_pos = m_conf.pos(n, n_atom.type().type());

			for (int m = 0;m <= n - offset;m++)
			{
				for (int j = 0;j < m_res[m].num_atoms();j++)
				{
					Atom m_atom = m_res[m].atom(j);
					if (m_atom.undefined()) { continue; }

					Point m_pos = m_conf.pos(m, m_atom.type().type());

					if (n_pos.closer_than(Contact_Dist, m_pos))
					{
#ifdef DEBUG_CONTACT
						std::cout << "Contact: " << n << " " << m
						 	<< "  " << n_pos.distance(m_pos) << "\n";
#endif // DEBUG_CONTACT
						contact[n] = true;
						total++;
						goto Next_n;
					}
				}
			}
		}

		Next_n: ;
	}

	for (n = start;n < m_length;n++)
	{
		if (helix[n] && !contact[n])
		{
			if (contact[n - 1] || (n < m_length - 1 && contact[n + 1]))
			{
				// count as a contact, since the helix as a whole
				// touches a previous residue
				total++;

#ifdef DEBUG_CONTACT
				std::cout << "Extra helix contact " << n << "\n";
#endif // DEBUG_CONTACT
			}
			else
			{
#ifdef DEBUG_CONTACT
				std::cout << "Missing helix contact " << n << "\n";
#endif // DEBUG_CONTACT
			}
		}
#ifdef DEBUG_CONTACT
		else
		if (!contact[n])
		{
			std::cout << "Missing contact " << n << "\n";
		}
#endif // DEBUG_CONTACT
	}

	return (double) total / (double) (m_length - min_offset);
}

//#define DEBUG_INTERLEAVE

double Peptide::interleaving_proportion() const
{
	assert(!reverseSaint || (m_length == m_full_length));

	static const Atom_Id atomtype[3] = { Atom_N, Atom_CA, Atom_C };

	// less than 15 (three parallel beta strands * 1.5)
	static const double Max_Dist = 14.0; // 13.0
	static const double Close_To_Line_Dist = 2.0;
	static const int min_offset = 6;  // 4
	static const int offset2 = 6; // 4

	std::vector<bool> helix;
	helix.resize(m_length);
	for (int i = 0;i < m_length;i++)
	{
		helix[i] = m_conf.in_helix(i);
	}

	int total = 0;
	int num_in_a_row = 0;

	for (int n = min_offset;n < m_length;n++)
	{
		bool success = false;

		for (int m1 = offset2;m1 <= n - min_offset;m1++)
		{
			// check if there are helices all the way from m1 to n
			bool all_helix = true;

			int a;
			for (a = m1 + 1;a < n;a++)
			{
				if (!helix[a])
				{
					all_helix = false;
					break;
				}
			}

			if (all_helix) { continue; }

			for (int m2 = 0;m2 <= m1 - offset2;m2++)
			{
				// check if there are helices all the way from m2 to m1
				all_helix = true;

				for (a = m2 + 1;a < m1;a++)
				{
					if (!helix[a])
					{
						all_helix = false;
						break;
					}
				}

				if (all_helix) { continue; }

				for (int s = 0;s < 3;s++)
				{
					Point n_pos = m_conf.pos(n, atomtype[s]);
	
					for (int s1 = 0;s1 < 3;s1++)
					{
						Point m1_pos = m_conf.pos(m1, atomtype[s1]);

						for (int s2 = 0;s2 < 3;s2++)
						{
							Point m2_pos = m_conf.pos(m2, atomtype[s2]);

							if (m1_pos.closer_than(Max_Dist, m2_pos))
							{
								Point closest;
								double param = n_pos.closest_point_param(
									m1_pos, m2_pos, &closest, true);

								if (param >= 0.0 && param <= 1.0)
								{
									double dist = n_pos.distance(closest);

									if (dist < Close_To_Line_Dist)
									{
#ifdef DEBUG_INTERLEAVE
										std::cout << "BETWEEN: "
											<< n << " ("
											<< m1 << " " << m2 << ")\n"
											<< " param = " << param
											<< ", dist to line = "
											<< dist
											<< "\n " << n << ") " << n_pos
											<< "\n " << m1 << ") " << m1_pos
											<< "\n " << m2 << ") " << m2_pos
											<< "\n dist apart "
												<< m1_pos.distance(m2_pos)
											<< "\n";
#endif // DEBUG_INTERLEAVE
										success = true;
										goto Next_N;
									}
								}
							}
						}
					}
				}
			}
		}

	Next_N:

		if (success)
		{
			num_in_a_row++;
		}
		else
		{
			if (num_in_a_row > 1)
			{
				total += num_in_a_row;
			}
#ifdef DEBUG_INTERLEAVE
			else
			if (num_in_a_row == 1)
			{
				std::cout << "-- ignoring 1 in a row\n";
			}
#endif // DEBUG_INTERLEAVE

			num_in_a_row = 0;
		}
	}

	if (num_in_a_row > 1)
	{
		total += num_in_a_row;
	}
#ifdef DEBUG_INTERLEAVE
	else
	if (num_in_a_row == 1)
	{
		std::cout << "-- ignoring 1 in a row\n";
	}
#endif // DEBUG_INTERLEAVE

	return (double) total / (double) (m_length - min_offset);
}

/*
double Peptide::interleaving_proportion() const
{
	static const double Contact_Dist = 5.0;
	static const int min_offset = 4;
	static const int offset2 = 3;

	int total = 0;

	for (int n = min_offset;n < m_length;n++)
	{
		for (int i = 0;i < m_res[n].num_atoms();i++)
		{
			Atom n_atom = m_res[n].atom(i);
			if (n_atom.undefined()) { continue; }

			Point n_pos = m_conf.pos(n, n_atom.type().type());

			for (int m1 = offset2;m1 <= n - min_offset;m1++)
			{
				for (int j1 = 0;j1 < m_res[m1].num_atoms();j1++)
				{
					Atom m1_atom = m_res[m1].atom(j1);
					if (m1_atom.undefined()) { continue; }

					Point m1_pos = m_conf.pos(m1, m1_atom.type().type());

					if (n_pos.closer_than(Contact_Dist, m1_pos))
					{
						for (int m2 = 0;m2 <= m1 - offset2;m2++)
						{
							for (int j2 = 0;j2 < m_res[m2].num_atoms();j2++)
							{
								Atom m2_atom = m_res[m2].atom(j2);
			                    if (m2_atom.undefined()) { continue; }

								Point m2_pos = m_conf.pos(m2,
									m2_atom.type().type());

								if (n_pos.closer_than(Contact_Dist, m2_pos))
								{
									// check if n_pos is between m1_pos
									// and m2_pos

									Point n_to_m1 = m1_pos.minus(n_pos);
									Point n_to_m2 = m2_pos.minus(n_pos);
									n_to_m1.normalise();
									n_to_m2.normalise();

									if (n_to_m1.dot_product(n_to_m2) < -0.8)
									{
										std::cout << "BETWEEN: "
											<< n << " ("
											<< m1 << " " << m2 << ")\n";
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return (double) total / (double) (m_length - min_offset);
}
*/

/*
double Peptide::previous_contact_proportion() const
{
	static const double Contact_Dist = 10.0;
	static const int min_offset = 6;
	int total = 0;

	for (int n = min_offset;n < m_length;n++)
	{
		// check if any previous residues are close to this one

		double closest_dist = 100.0;
		int closest = -1;

		for (int i = 0;i < m_res[n].num_atoms();i++)
		{
			Atom n_atom = m_res[n].atom(i);
			if (n_atom.undefined()) { continue; }

			Point n_pos = m_conf.pos(n, n_atom.type().type());

			for (int m = 0;m <= n - min_offset;m++)
			{
				for (int j = 0;j < m_res[m].num_atoms();j++)
				{
					Atom m_atom = m_res[m].atom(j);
					if (m_atom.undefined()) { continue; }

					Point m_pos = m_conf.pos(m, m_atom.type().type());

					if (n_pos.closer_than(Contact_Dist, m_pos))
					{
						double d = n_pos.distance(m_pos);

						if (d < closest_dist)
						{
							closest_dist = d;
							closest = m;
						}
					}
				}
			}
		}

		std::cout << "Contact " << n << " & " << closest << ": "
			<< closest_dist << "\n";

		// std::cout << "RG up to " << n << " = " << radius_of_gyr(n) << "\n";
	}

	return (double) total / (double) (m_length - min_offset);
}
*/

/*
void Peptide::call_scwrl()
{
	if (!verify_scwrl_exec(false))
	{
		std::cerr << "Warning: SCWRL not found; outputting backbone only\n";
		return;
	}

	std::cerr << "Calling SCWRL ...\n\n";

	Temp_File tmp_in("scwrl_in_");
	Temp_File tmp_out("scwrl_out_");

	// save and restore m_filename (modified by write_pdb() and read_pdb())
	std::string old_filename = m_filename;

	write_pdb(tmp_in.name(),
		true,	// backbone only
		true	// atom records only
	);

	char command[1000];
	sprintf(command, "%s -h -i %s -o %s >/dev/null",
		m_scwrl_exec.c_str(), tmp_in.name(), tmp_out.name());
//std::cerr << "SCWRL command:\n" << command << "\n";
	system(command);
//std::cerr << "FINISHED SCWRL\n";

	if (!read_pdb(tmp_out.name(), ' ',
		true	// don't print warning messages
	))
	{
		std::cerr << "Errors in PDB file read back from SCWRL\n";
		exit(1);
	}

	// restore old filename
	m_filename = old_filename;

	// temp files are deleted automatically
}
*/

/*
bool Peptide::verify_scwrl_exec(
	bool exit_on_err //=true
)
{
	if (m_scwrl_exec.empty())
	{
		if (!exit_on_err)
		{
			return false;
		}

		std::cerr << "Error: SCWRL executable is undefined "
			"(should be defined in config file)\n";
		exit(1);
	}
	else
	{
		struct stat s;
		
		if (stat(m_scwrl_exec.c_str(), &s) != 0 ||
			((s.st_mode & (0x40|0x08|0x01)) == 0) || // "executable" mode bits
			S_ISDIR(s.st_mode))
		{
			if (!exit_on_err)
			{
				return false;
			}

			std::cerr << "Error: SCWRL executable is incorrect: \""
				<< m_scwrl_exec
				<< "\"\n";
			exit(1);
		}
	}

	return true;
}
*/

const Residue &Peptide::res(int n) const
{
	assert(n >= 0 && n < m_full_length);
	return m_res[(int)n];
}

Residue &Peptide::res(int n)
{
	assert(n >= 0 && n < m_full_length);
	return m_res[(int)n];
}

bool Peptide::get_cb_ca_pos(int n, Point *pos) const
{
	bool gly = res(n).amino().is_glycine();
	Atom_Id a = (gly ? Atom_CA : Atom_CB);

	if (!atom_exists(n, a))
	{
		return false;
	}

	*pos = atom_pos(n, a);
	return true;
}

bool Peptide::get_side_chain_pos(int n, Point *pos) const
{
	if (res(n).amino().is_glycine())
	{
		if (!(atom_exists(n, Atom_CA) &&
			  atom_exists(n, Atom_N) &&
			  atom_exists(n, Atom_C)))
		{
			return false;
		}

		*pos = estimate_CB_pos(
			atom_pos(n, Atom_CA),
			atom_pos(n, Atom_N),
			atom_pos(n, Atom_C));
	}
	else
	{
		if (!atom_exists(n, Atom_CB))
		{
			return false;
		}

		*pos = atom_pos(n, Atom_CB);
	}

	return true;
}

void Peptide::idealise()
{
	assert(m_length == m_full_length);

	if (num_missing_res() != 0 ||
		num_missing_backbone() != 0)
	{
		std::cerr << "Error: cannot idealise backbone, since there are "
			"missing backbone atoms or residues\n";
		exit(1);
	}

	int i = m_length - 1;

	Point n_pos, ca_pos, c_pos, pos;
	get_initial_ideal(&n_pos, &ca_pos, &c_pos);
	set_atom_pos(i, Atom_N, n_pos),
	set_atom_pos(i, Atom_CA, ca_pos);
	set_atom_pos(i, Atom_C, c_pos);

	if (!is_glycine(i))
	{
		set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
	}

	// calculate where next N would be to get position for O

	Point est_n_pos = torsion_to_coord(n_pos, ca_pos, c_pos,
		BOND_LENGTH_C_N, BOND_ANGLE_CA_C_N, deg2rad(180.0), BOND_LENGTH_C_C);
	set_atom_pos(i, Atom_O, estimate_O_pos(ca_pos, c_pos, est_n_pos));

	for (i--;i >= 0;i--)
	{
		double phi_val = m_conf.phi(i + 1);
		c_pos = torsion_to_coord(c_pos, ca_pos, n_pos,
			BOND_LENGTH_C_N, BOND_ANGLE_C_N_CA, phi_val, BOND_LENGTH_N_CA);
		set_atom_pos(i, Atom_C, c_pos);

		double omega_val = m_conf.omega(i);
		ca_pos = torsion_to_coord(ca_pos, n_pos, c_pos,
			BOND_LENGTH_C_C, BOND_ANGLE_CA_C_N, omega_val, BOND_LENGTH_C_N);
		set_atom_pos(i, Atom_CA, ca_pos);

		set_atom_pos(i, Atom_O, estimate_O_pos(ca_pos, c_pos, n_pos));

		double psi_val = m_conf.psi(i);
		n_pos = torsion_to_coord(n_pos, c_pos, ca_pos,
			BOND_LENGTH_N_CA, BOND_ANGLE_N_CA_C, psi_val, BOND_LENGTH_C_C);
		set_atom_pos(i, Atom_N, n_pos);

		if (!is_glycine(i))
		{
			set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
		}
	}

	verify_ideal();
	m_conf.verify_torsion_angles();
}

void Peptide::idealise_bond_lengths()
{
	assert(m_length == m_full_length);

	if (num_missing_res() != 0 ||
		num_missing_backbone() != 0)
	{
		std::cerr << "Error: cannot idealise cond lengths, since there are "
			"missing backbone atoms or residues\n";
		exit(1);
	}

	// get the current bond angles for all main chain atoms

	std::vector<double> angle_n;
	std::vector<double> angle_ca;
	std::vector<double> angle_c;

	angle_n.resize(m_length);
	angle_ca.resize(m_length);
	angle_c.resize(m_length);

	Point n_pos, ca_pos, c_pos;

	for (int n = 0;n < m_length;n++)
	{
		n_pos = atom_pos(n, Atom_N);
		ca_pos = atom_pos(n, Atom_CA);
		c_pos = atom_pos(n, Atom_C);

		if (n == 0)
		{
			angle_n[n] = BOND_ANGLE_C_N_CA;
		}
		else
		{
			angle_n[n] = angle_formed(atom_pos(n - 1, Atom_C), n_pos, ca_pos);
		}

		angle_ca[n] = angle_formed(n_pos, ca_pos, c_pos);

		if (n == m_length - 1)
		{
			angle_c[n] = BOND_ANGLE_CA_C_N;
		}
		else
		{
			angle_c[n] = angle_formed(ca_pos, c_pos, atom_pos(n + 1, Atom_N));
		}

		/*
		std::cout << "BA " << n << " "
			<< rad2deg(angle_n[n]) << " "
			<< rad2deg(angle_ca[n]) << " "
			<< rad2deg(angle_c[n]) << "\n";
		*/
	}

	int i = m_length - 1;

	// Handle residue i (the C terminus)
	// Leave CA where it is, and rescale bond lengths to N and C atoms

	ca_pos = atom_pos(i, Atom_CA);
	Point ca_to_n = atom_pos(i, Atom_N).minus(ca_pos);
	Point ca_to_c = atom_pos(i, Atom_C).minus(ca_pos);
	ca_to_n.normalise();
	ca_to_c.normalise();
	ca_to_n.scale(BOND_LENGTH_N_CA);
	ca_to_c.scale(BOND_LENGTH_C_C);
	n_pos = ca_pos.plus(ca_to_n);
	c_pos = ca_pos.plus(ca_to_c);

	// c_pos, n_pos and ca_pos are now all correct

	set_atom_pos(i, Atom_N, n_pos);
	set_atom_pos(i, Atom_C, c_pos);

	// calculate where next N would be to get position for O

	double pretend_psi = deg2rad(180.0);
	Point est_n_pos = torsion_to_coord(n_pos, ca_pos, c_pos,
        BOND_LENGTH_C_N, BOND_ANGLE_CA_C_N, pretend_psi, BOND_LENGTH_C_C);
    set_atom_pos(i, Atom_O, estimate_O_pos(ca_pos, c_pos, est_n_pos));

	if (!is_glycine(i))
	{
		set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
	}

	// calculate main chain atom positions for the rest of the residues

	for (i--;i >= 0;i--)
	{
		c_pos = torsion_to_coord(c_pos, ca_pos, n_pos,
			BOND_LENGTH_C_N, angle_n[i + 1], m_conf.phi(i + 1),
			BOND_LENGTH_N_CA);
		set_atom_pos(i, Atom_C, c_pos);

		ca_pos = torsion_to_coord(ca_pos, n_pos, c_pos,
			BOND_LENGTH_C_C, angle_c[i], m_conf.omega(i), BOND_LENGTH_C_N);
		set_atom_pos(i, Atom_CA, ca_pos);

		set_atom_pos(i, Atom_O, estimate_O_pos(ca_pos, c_pos, n_pos));

		n_pos = torsion_to_coord(n_pos, c_pos, ca_pos,
			BOND_LENGTH_N_CA, angle_ca[i], m_conf.psi(i), BOND_LENGTH_C_C);
		set_atom_pos(i, Atom_N, n_pos);

		if (!is_glycine(i))
		{
			set_atom_pos(i, Atom_CB, estimate_CB_pos(ca_pos, n_pos, c_pos));
		}
	}
}

double Peptide::radius_of_gyr(int up_to_res /*= -1*/) const
{
	if (up_to_res == -1)
	{
		up_to_res = end();
	}

	double total = 0.0;
	int num = 0;

	int n;
	for (n = start() + 1;n <= up_to_res;n++)
	{
		if (m_res[n].atom_exists(Atom_CA))
		{
			Point p1 = atom_pos(n, Atom_CA);

			for (int m = start();m < n;m++)
			{
				if (m_res[m].atom_exists(Atom_CA))
				{
					Point p2 = atom_pos(m, Atom_CA);
					total +=
						square(p1.x - p2.x) +
						square(p1.y - p2.y) +
						square(p1.z - p2.z);
					num++;
				}
			}
		}
	}

	return sqrt(total / (double) num);
}

double Peptide::diameter() const
{
	assert(!reverseSaint || (m_length == m_full_length));

	double max_dist = 0.0;

	for (int n = start() + 1;n <= end();n++)
	{
		const Residue &res_n = m_res[n];

		for (int nn = 0;nn < res_n.num_atoms();nn++)
		{
			if (res_n.atom(nn).type().undefined()) { continue; }

			Point p1 = atom_pos(n, res_n.atom(nn).type().type());

			for (int m = start();m < n;m++)
			{
				const Residue &res_m = m_res[m];

				for (int mm = 0;mm < res_m.num_atoms();mm++)
				{
					if (res_m.atom(mm).type().undefined()) { continue; }

					Point p2 = atom_pos(m, res_m.atom(mm).type().type());

					if (fabs(p1.x - p2.x) > max_dist ||
						fabs(p1.y - p2.y) > max_dist ||
						fabs(p1.z - p2.z) > max_dist)
					{
						double d = p1.distance(p2);

						if (d > max_dist)
						{
							max_dist = d;
						}
					}
				}
			}
		}
	}

	return max_dist;
}

int Peptide::num_missing_res() const
{
	int total = 0;

	for (int n = 0;n < (int) m_res.size() - 1;n++)
	{
		if (m_res[n].missing_after())
		{
			total++;
		}
	}

	return total;
}

void Peptide::print_missing_res(const char *filename) const
{
	for (int n = 0;n < (int) m_res.size() - 1;n++)
	{
		if (m_res[n].missing_after())
		{
			std::cout << "=> gap between residues "
				<< m_res[n].res_seq_str()
				<< " and "
				<< m_res[n+1].res_seq_str()
				<< " (# "
				<< n << " / " << m_res.size()
				<< " ) on line "
				<< m_res[n+1].pdb_line()
				<< " of "
				<< filename
				<< "\n";
		}
	}
}

int Peptide::num_missing_backbone() const
{
	int total = 0;

	for (int n = 0;n < (int) m_res.size();n++)
	{
		total += m_res[n].num_missing_backbone();
	}

	return total;
}

void Peptide::verify_ideal() const
{
	assert(!reverseSaint || (m_length == m_full_length));

	// check if all atoms exist

	int i;
	for (i = 0;i < m_length;i++)
	{
		const Residue &res = m_res[i];
		bool is_gly = res.amino().is_glycine();

		for (int a = 0;a < Num_Backbone;a++)
		{
			if (!(is_gly && (Atom_Id) a == Atom_CB))
			{
				assert(res.atom_exists((Atom_Id) a));
			}
		}

		assert(res.num_atoms() == (is_gly ? Num_Backbone - 1 : Num_Backbone));
	}

	if (m_length != 0)
	{
		// check all bond lengths, bond angles and torsion angles

		Point n_pos = atom_pos(0, Atom_N);
		Point ca_pos = atom_pos(0, Atom_CA);
		Point c_pos = atom_pos(0, Atom_C);

		assert(fabs(ca_pos.distance(n_pos) - BOND_LENGTH_N_CA) < 0.001);
		assert(fabs(c_pos.distance(ca_pos) - BOND_LENGTH_C_C) < 0.001);
		assert(fabs(angle_formed(n_pos, ca_pos, c_pos) - BOND_ANGLE_N_CA_C)
			< 0.001);

		for (i = 1;i < m_length;i++)
		{
			Point n_pred = torsion_to_coord(n_pos, ca_pos, c_pos,
				BOND_LENGTH_C_N, BOND_ANGLE_CA_C_N,
				m_conf.psi(i - 1), BOND_LENGTH_C_C);
			n_pos = atom_pos(i, Atom_N);

			assert(n_pred.close_to(n_pos, 0.1));
			assert(fabs(n_pos.distance(c_pos) - BOND_LENGTH_C_N) < 0.01);
			assert(fabs(angle_formed(ca_pos, c_pos, n_pos) - BOND_ANGLE_CA_C_N)
				< 0.01);

			Point ca_pred = torsion_to_coord(ca_pos, c_pos, n_pos,
				BOND_LENGTH_N_CA, BOND_ANGLE_C_N_CA,
				m_conf.omega(i - 1), BOND_LENGTH_C_N);
			ca_pos = atom_pos(i, Atom_CA);

			assert(ca_pred.close_to(ca_pos, 0.1));
			assert(fabs(ca_pos.distance(n_pos) - BOND_LENGTH_N_CA) < 0.01);
			assert(fabs(angle_formed(c_pos, n_pos, ca_pos) - BOND_ANGLE_C_N_CA)
				< 0.01);

			Point c_pred = torsion_to_coord(c_pos, n_pos, ca_pos,
				BOND_LENGTH_C_C, BOND_ANGLE_N_CA_C,
				m_conf.phi(i), BOND_LENGTH_N_CA);
			c_pos = atom_pos(i, Atom_C);

			assert(c_pred.close_to(c_pos, 0.1));
			assert(fabs(c_pos.distance(ca_pos) - BOND_LENGTH_C_C) < 0.01);
			assert(fabs(angle_formed(n_pos, ca_pos, c_pos) - BOND_ANGLE_N_CA_C)
				< 0.01);
		}

		// check CB and O positions

		for (i = 0;i < m_length;i++)
		{
			n_pos = atom_pos(i, Atom_N);
			ca_pos = atom_pos(i, Atom_CA);
			c_pos = atom_pos(i, Atom_C);

			if (!is_glycine(i))
			{
				Point cb_pos = atom_pos(i, Atom_CB);
				Point cb_pred = estimate_CB_pos(ca_pos, n_pos, c_pos);
				assert(cb_pred.close_to(cb_pos, 0.1));
			}

			Point o_pos = atom_pos(i, Atom_O);

			if (i < m_length - 1)
			{
				Point next_n_pos = atom_pos(i + 1, Atom_N);
				Point o_pred = estimate_O_pos(ca_pos, c_pos, next_n_pos);
				assert(o_pred.close_to(o_pos, 0.1));
			}
			else
			{
				// don't know where O should be, so just check bond length
				assert(fabs(c_pos.distance(o_pos) - BOND_LENGTH_C_O) < 0.3);
			}
		}
	}
}

void Peptide::verify() const
{
	// check if all atoms exist

	int i;
	for (i = start();i <= end();i++)
	{
		const Residue &res = m_res[i];
		bool is_gly = res.amino().is_glycine();

		for (int a = 0;a < Num_Backbone;a++)
		{
			if (!(is_gly && (Atom_Id) a == Atom_CB))
			{
				assert(res.atom_exists((Atom_Id) a));
			}
		}

		assert(res.num_atoms() == (is_gly ? Num_Backbone - 1 : Num_Backbone));
	}

	if (m_length != 0)
	{
		static const double thirtyDeg = deg2rad(30.0);
		// check all bond lengths, bond angles and torsion angles

		Point n_pos = atom_pos(start(), Atom_N);
		Point ca_pos = atom_pos(start(), Atom_CA);
		Point c_pos = atom_pos(start(), Atom_C);

		assert(fabs(ca_pos.distance(n_pos) - BOND_LENGTH_N_CA) < 0.5);
		assert(fabs(c_pos.distance(ca_pos) - BOND_LENGTH_C_C) < 0.5);
		assert(fabs(angle_formed(n_pos, ca_pos, c_pos) - BOND_ANGLE_N_CA_C)
			< thirtyDeg);

		for (i = start() + 1;i <= end();i++)
		{
			Point n_pred = torsion_to_coord(n_pos, ca_pos, c_pos,
				BOND_LENGTH_C_N, BOND_ANGLE_CA_C_N,
				m_conf.psi(i - 1), BOND_LENGTH_C_C);
			//Point old_n = n_pos;
			n_pos = atom_pos(i, Atom_N);

			assert(n_pred.close_to(n_pos, 1.0));
			assert(fabs(n_pos.distance(c_pos) - BOND_LENGTH_C_N) < 0.5);
			assert(fabs(angle_formed(ca_pos, c_pos, n_pos) - BOND_ANGLE_CA_C_N)
				< thirtyDeg);

			Point ca_pred = torsion_to_coord(ca_pos, c_pos, n_pos,
				BOND_LENGTH_N_CA, BOND_ANGLE_C_N_CA,
				m_conf.omega(i - 1), BOND_LENGTH_C_N);
			ca_pos = atom_pos(i, Atom_CA);

			assert(ca_pred.close_to(ca_pos, 1.0));
			assert(fabs(ca_pos.distance(n_pos) - BOND_LENGTH_N_CA) < 0.5);
			assert(fabs(angle_formed(c_pos, n_pos, ca_pos) - BOND_ANGLE_C_N_CA)
				< thirtyDeg);

			Point c_pred = torsion_to_coord(c_pos, n_pos, ca_pos,
				BOND_LENGTH_C_C, BOND_ANGLE_N_CA_C,
				m_conf.phi(i), BOND_LENGTH_N_CA);
			c_pos = atom_pos(i, Atom_C);

			assert(c_pred.close_to(c_pos, 1.0));
			assert(fabs(c_pos.distance(ca_pos) - BOND_LENGTH_C_C) < 0.5);
			assert(fabs(angle_formed(n_pos, ca_pos, c_pos) - BOND_ANGLE_N_CA_C)
				< thirtyDeg);
		}

		// check CB and O positions

		for (i = start();i <= end();i++)
		{
			n_pos = atom_pos(i, Atom_N);
			ca_pos = atom_pos(i, Atom_CA);
			c_pos = atom_pos(i, Atom_C);

			if (!is_glycine(i))
			{
				Point cb_pos = atom_pos(i, Atom_CB);
				Point cb_pred = estimate_CB_pos(ca_pos, n_pos, c_pos);
				assert(cb_pred.close_to(cb_pos, 1.0));
			}

			Point o_pos = atom_pos(i, Atom_O);

			if (i < end())
			{
				Point next_n_pos = atom_pos(i + 1, Atom_N);
				Point o_pred = estimate_O_pos(ca_pos, c_pos, next_n_pos);
				assert(o_pred.close_to(o_pos, 1.0));
			}
			else
			{
				// don't know where O should be, so just check bond length
				assert(fabs(c_pos.distance(o_pos) - BOND_LENGTH_C_O) < 0.5);
			}
		}
	}
}

/***/
namespace
{
	// variables used in calc_rmsd()

	const int RMSD_Max_Atoms = 10000;
	double RMSD_pos1[RMSD_Max_Atoms][3];
	double RMSD_pos2[RMSD_Max_Atoms][3];
}

double Peptide::calc_rmsd(const Peptide &p) const
{
	if (m_length != p.m_length)
	{
		std::cerr << "Error in Peptide::calc_rmsd(): peptides have different "
			"lengths (" << m_length << " & " << p.m_length << ")\n";
		exit(1);
	}

	int num = 0;

	for (int n = start();n <= end();n++)
	{
		const Residue &r = m_res[n];
		const Residue &r2 = p.m_res[n];

		for (int i = 0;i < r.num_atoms();i++)
		{
			const Atom &a = r.atom(i);

			if (r2.atom_exists(a.type()))
			{
				Point p1 = atom_pos(n, a.type().type());
				Point p2 = p.atom_pos(n, a.type().type());

				if (num >= RMSD_Max_Atoms)
				{
					std::cerr << "Error: Too many atoms to compare in "
						"Peptide::calc_rmsd()\n";
					exit(1);
				}

				RMSD_pos1[num][0] = p1.x;
				RMSD_pos1[num][1] = p1.y;
				RMSD_pos1[num][2] = p1.z;

				RMSD_pos2[num][0] = p2.x;
				RMSD_pos2[num][1] = p2.y;
				RMSD_pos2[num][2] = p2.z;

				num++;
			}
		}
	}

	return get_rmsd(RMSD_pos1, RMSD_pos2, num);
}
/***/

void Peptide::dump(std::ostream &out /*= std::cout*/) const
{
	out << "Dumping peptide: current length = " << m_length
		<< ", full length = " << m_full_length
		<< ", " << m_res.size() << " residues in vector\n\n";

	for (int n = start();n <= end();n++)
	{
		const Residue &r = m_res[n];

		out << n << ". " << r.amino().code()
			<< " (" << r.amino().abbr() << " / " << r.amino().name()
			<< "), codon " << r.codon().str()
			<< ", "
			<< r.num_atoms() << " atoms\n\n";

		for (int i = 0;i < r.num_atoms();i++)
		{
			const Atom &a = r.atom(i);

			out << "   " << i << ". id = "
				<< a.rapdf_id() << " (";

			a.type().dump(out);

			out << ")";

			if (a.undefined())
			{
				out << " - missing";
			}

			Point p = atom_pos(n, a.id());

			out << ", pos = ("
				<< p.x << ", "
				<< p.y << ", "
				<< p.z << ")";

			out << "\n";
		}

		out << "\n";
	}
}

void Peptide::dump_backbone(std::ostream &out /*= std::cout*/) const
{
	double total_n_ca_c = 0.0;
	double total_ca_c_n = 0.0;
	double total_c_n_ca = 0.0;

	for (int n = start();n <= end();n++)
	{
		// bond lengths

		Point n_pos = atom_pos(n, Atom_N);
		Point ca_pos = atom_pos(n, Atom_CA);
		Point c_pos = atom_pos(n, Atom_C);

		double n_ca = n_pos.distance(ca_pos);
		double ca_c = ca_pos.distance(c_pos);

		out << n << ") N-CA "
			<< Printf("%.3f ", n_ca)
			<< Printf("%.4f\n", n_ca - BOND_LENGTH_N_CA);

		out << n << ") CA-C "
			<< Printf("%.3f ", ca_c)
			<< Printf("%.4f\n", ca_c - BOND_LENGTH_C_C);

		if (n != end())
		{
			Point next_n_pos = atom_pos(n + 1, Atom_N);
			double c_next_n = c_pos.distance(next_n_pos);

			out << n << ") C-N  "
				<< Printf("%.3f ", c_next_n)
				<< Printf("%.4f\n", c_next_n - BOND_LENGTH_C_N);
		}

		// bond angles
		double diff;

		if (n != start())
		{
			Point prev_c_pos = atom_pos(n - 1, Atom_C);
			double c_n_ca = rad2deg(angle_formed(prev_c_pos, n_pos, ca_pos));
			diff = c_n_ca - 121.9;
			if (diff > 180.0) { diff -= 360.0; }
			if (diff < -180.0) { diff += 360.0; }

			total_c_n_ca += fabs(diff);

			out << n << ") angle C-N-CA "
				<< Printf("%.2f ", c_n_ca)
				<< Printf("%.3f\n", diff);
		}

		double n_ca_c = rad2deg(angle_formed(n_pos, ca_pos, c_pos));

		diff = n_ca_c - 109.5;
		if (diff > 180.0) { diff -= 360.0; }
		if (diff < -180.0) { diff += 360.0; }

		total_n_ca_c += fabs(diff);

		out << n << ") angle N-CA-C "
			<< Printf("%.2f ", n_ca_c)
			<< Printf("%.3f\n", diff);

		if (n != end())
		{
			Point next_n_pos = atom_pos(n + 1, Atom_N);
			double ca_c_n = rad2deg(angle_formed(ca_pos, c_pos, next_n_pos));

			diff = ca_c_n - 115.6;
			if (diff > 180.0) { diff -= 360.0; }
			if (diff < -180.0) { diff += 360.0; }

			total_ca_c_n += fabs(diff);

			out << n << ") angle CA-C-N "
				<< Printf("%.2f ", ca_c_n)
				<< Printf("%.3f\n", diff);
		}
	}

	out << "\nAverage abs difference:\n\n"
		<< "C-N-CA " << total_c_n_ca / (double) m_length << "\n"
		<< "N-CA-C " << total_n_ca_c / (double) m_length << "\n"
		<< "CA-C-N " << total_ca_c_n / (double) m_length << "\n";
}

