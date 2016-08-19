
#include <cassert>
#include <cmath>
#include "peptide.h"
#include "geom.h"
#include "common.h"
#include "transform.h"
#include "conformation.h"

Conformation::Conformation()
	: m_peptide(NULL)
{
}

Conformation::Conformation(const Peptide *p)
	: m_peptide(p)
{
	set_num_res(p->full_length());
}

Conformation::~Conformation()
{
}

Conformation::Conformation(const Conformation &other)
	: m_peptide(other.m_peptide),
	  m_backbone(other.m_backbone),
	  m_non_backbone(other.m_non_backbone),
	  m_res_data(other.m_res_data)
{
}

Conformation& Conformation::operator = (const Conformation &other)
{
	m_peptide = other.m_peptide;
	m_backbone = other.m_backbone;
	m_non_backbone = other.m_non_backbone;
	m_res_data = other.m_res_data;

	return *this;
}

void Conformation::swap(Conformation &other)
{
	m_backbone.swap(other.m_backbone);
	m_non_backbone.swap(other.m_non_backbone);
	m_res_data.swap(other.m_res_data);

	const Peptide *p = m_peptide;
	m_peptide = other.m_peptide;
	other.m_peptide = p;
}

void Conformation::set_num_res(int n)
{
	m_backbone.resize(n * Num_Backbone);

	if (m_non_backbone.size() != 0)
	{
		m_non_backbone.resize(n * Num_NonBackbone);
	}

	m_res_data.resize(n);
}

void Conformation::clear()
{
    // (using a temporary variable like this reduces the vector's
    // capacity to 0; just using erase() would not deallocate any memory)

   	Point_Vec empty_vec;
   	m_backbone.swap(empty_vec);
	remove_non_backbone_atoms();

	// TO DO: clear m_res_data as well
}

void Conformation::remove_non_backbone_atoms()
{
	if (m_non_backbone.size() != 0)
	{
    	// (using a temporary variable like this reduces the vector's
    	// capacity to 0; just using erase() would not deallocate any memory)

    	Point_Vec empty_vec;
    	m_non_backbone.swap(empty_vec);
	}
}

Point& Conformation::pos(int n, Atom_Id atom)
{
	assert(n * Num_Backbone < (int) m_backbone.size());

	if (atom < Num_Backbone)
	{
		return m_backbone[n * Num_Backbone + atom];
	}
	else
	{
		assert(m_non_backbone.size() != 0);
		return m_non_backbone[(n * Num_NonBackbone) + atom - Num_Backbone];
	}
}

void Conformation::set_pos(int n, Atom_Id atom, const Point &p)
{
	assert(n * Num_Backbone < (int) m_backbone.size());

	if (atom < Num_Backbone)
	{
		m_backbone[n * Num_Backbone + atom] = p;
	}
	else
	{
		if (m_non_backbone.size() == 0)
		{
			m_non_backbone.resize(
				(m_backbone.size() / Num_Backbone) * Num_NonBackbone);
		}

		int i = (n * Num_NonBackbone) + atom - Num_Backbone;
		m_non_backbone[i] = p;
	}
}

void Conformation::transform_pos(int n, Atom_Id atom, const Transform &t)
{
	assert(n * Num_Backbone < (int) m_backbone.size());

	if (atom < Num_Backbone)
	{
		int i = n * Num_Backbone + atom;
		m_backbone[i] = t.times(m_backbone[i]);
	}
	else
	{
		if (m_non_backbone.size() == 0)
		{
			m_non_backbone.resize(
				(m_backbone.size() / Num_Backbone) * Num_NonBackbone);
		}

		int i = (n * Num_NonBackbone) + atom - Num_Backbone;
		m_non_backbone[i] = t.times(m_non_backbone[i]);
	}
}

void Conformation::calc_torsion_angles()
{
	assert(m_peptide != NULL);
	m_res_data.resize(m_peptide->full_length());

	Point prev_c_pos;
	bool prev_c_missing = true;

	Point next_n_pos, next_ca_pos;
	bool next_n_missing, next_ca_missing;

	bool n_missing = !m_peptide->atom_exists(0, Atom_N);
	Point n_pos = pos(0, Atom_N);

	bool ca_missing = !m_peptide->atom_exists(0, Atom_CA);
	Point ca_pos = pos(0, Atom_CA);

	bool gap_prev = true;
	bool gap_next;

	int num_res = m_peptide->length();

	for (int n = 0;n < m_peptide->length();n++)
	{
		gap_next = m_peptide->res(n).missing_after();

		bool c_missing = !m_peptide->atom_exists(n, Atom_C);
		Point c_pos = pos(n, Atom_C);

		if (prev_c_missing || n_missing || ca_missing || c_missing || gap_prev)
		{
			m_res_data[n].phi = TORSION_UNKNOWN;
		}
		else
		{
			m_res_data[n].phi =
				torsion_angle(prev_c_pos, n_pos, ca_pos, c_pos);
		}

		next_n_missing = ((n == num_res - 1) ||
			!m_peptide->atom_exists(n + 1, Atom_N));

		if (!next_n_missing)
		{
			next_n_pos = pos(n + 1, Atom_N);
		}

		if (n_missing || ca_missing || c_missing || next_n_missing ||
			gap_next)
		{
			m_res_data[n].psi = TORSION_UNKNOWN;
		}
		else
		{
			m_res_data[n].psi =
				torsion_angle(n_pos, ca_pos, c_pos, next_n_pos);
		}

		next_ca_missing = ((n == num_res - 1) ||
			!m_peptide->atom_exists(n + 1, Atom_CA));

		if (!next_ca_missing)
		{
			next_ca_pos = pos(n + 1, Atom_CA);
		}

		if (ca_missing || c_missing || next_n_missing || next_ca_missing ||
			gap_next)
		{
			m_res_data[n].omega = TORSION_UNKNOWN;
		}
		else
		{
			m_res_data[n].omega = 
				torsion_angle(ca_pos, c_pos, next_n_pos, next_ca_pos);
		}

		n_missing = next_n_missing;
		n_pos = next_n_pos;

		ca_missing = next_ca_missing;
		ca_pos = next_ca_pos;

		prev_c_missing = c_missing;
		prev_c_pos = c_pos;

		gap_prev = gap_next;
	}
}

void Conformation::verify_torsion_angles() const
{
	Point prev_c_pos;
	bool prev_c_missing = true;

	Point next_n_pos, next_ca_pos;
	bool next_n_missing, next_ca_missing;

	int i = m_peptide->start();
	bool n_missing = !m_peptide->atom_exists(i, Atom_N);
	Point n_pos = pos(i, Atom_N);

	bool ca_missing = !m_peptide->atom_exists(i, Atom_CA);
	Point ca_pos = pos(i, Atom_CA);

	bool gap_prev = true;
	bool gap_next;

	int num_res = m_peptide->length();

	for (int n = m_peptide->start();n <= m_peptide->end();n++)
	{
		gap_next = m_peptide->res(n).missing_after();

		bool c_missing = !m_peptide->atom_exists(n, Atom_C);
		Point c_pos = pos(n, Atom_C);

		if (prev_c_missing || n_missing || ca_missing || c_missing || gap_prev)
		{
			//assert(m_res_data[n].phi == TORSION_UNKNOWN);
		}
		else
		{
			/*
			std::cout << "VERIFY_TORSION_ANGLES: n = " << n
				<< " phi = " << m_res_data[n].phi
				<< " real = " << torsion_angle(prev_c_pos, n_pos, ca_pos, c_pos)
				<< " pC= " << prev_c_pos
				<< " N = " << n_pos
				<< " CA= " << ca_pos
				<< " C = " << c_pos
				<< "\n";
			*/
			assert(approx_equal_angle(m_res_data[n].phi,
				torsion_angle(prev_c_pos, n_pos, ca_pos, c_pos)));
		}

		next_n_missing = ((n == num_res - 1) ||
			!m_peptide->atom_exists(n + 1, Atom_N));

		if (!next_n_missing)
		{
			next_n_pos = pos(n + 1, Atom_N);
		}

		if (n_missing || ca_missing || c_missing || next_n_missing ||
			gap_next)
		{
			//assert(m_res_data[n].psi == TORSION_UNKNOWN);
		}
		else
		{
			/*
			std::cout << "VERIFY_TORSION_ANGLES: n = " << n
				<< " psi = " << m_res_data[n].psi
				<< " real = " << torsion_angle(n_pos, ca_pos, c_pos, next_n_pos)
				<< " N = " << n_pos
				<< " CA= " << ca_pos
				<< " C = " << c_pos
				<< " N2= " << next_n_pos
				<< "\n";
			*/
			assert(approx_equal_angle(m_res_data[n].psi,
				torsion_angle(n_pos, ca_pos, c_pos, next_n_pos)));
		}

		next_ca_missing = ((n == num_res - 1) ||
			!m_peptide->atom_exists(n + 1, Atom_CA));

		if (!next_ca_missing)
		{
			next_ca_pos = pos(n + 1, Atom_CA);
		}

		if (ca_missing || c_missing || next_n_missing || next_ca_missing ||
			gap_next)
		{
			//assert(m_res_data[n].omega == TORSION_UNKNOWN);
		}
		else
		{
			/*
			std::cout << "VERIFY_TORSION_ANGLES: n = " << n
				<< " omega = " << m_res_data[n].omega
				<< " real = " << torsion_angle(ca_pos, c_pos, next_n_pos, next_ca_pos)
				<< " diff = " << fabs(m_res_data[n].omega - torsion_angle(ca_pos, c_pos, next_n_pos, next_ca_pos))
				<< " CA = " << ca_pos
				<< " C  = " << c_pos
				<< " N2 = " << next_n_pos
				<< " CA2= " << next_ca_pos
				<< "\n";
			*/
			assert(approx_equal_angle(m_res_data[n].omega,
				torsion_angle(ca_pos, c_pos, next_n_pos, next_ca_pos)));
		}

		n_missing = next_n_missing;
		n_pos = next_n_pos;

		ca_missing = next_ca_missing;
		ca_pos = next_ca_pos;

		prev_c_missing = c_missing;
		prev_c_pos = c_pos;

		gap_prev = gap_next;
	}
}

double Conformation::phi(int n) const
{
	if (n >= (int) m_res_data.size())
	{
		return TORSION_UNKNOWN;
	}
	
	return m_res_data[n].phi;
}

double Conformation::psi(int n) const
{
	if (n >= (int) m_res_data.size())
	{
		return TORSION_UNKNOWN;
	}

	return m_res_data[n].psi;
}

double Conformation::omega(int n) const
{
	if (n >= (int) m_res_data.size())
	{
		return TORSION_UNKNOWN;
	}

	return m_res_data[n].omega;
}

void Conformation::set_phi(int n, double angle)
{
	m_res_data[n].phi = range_mpi_pi(angle);
}

void Conformation::set_psi(int n, double angle)
{
	m_res_data[n].psi = range_mpi_pi(angle);
}

void Conformation::set_omega(int n, double angle)
{
	m_res_data[n].omega = range_mpi_pi(angle);
}

bool Conformation::in_helix(int n) const
{
	static const double phi_max = range_mpi_pi(deg2rad(-35.0));
	static const double phi_min = range_mpi_pi(deg2rad(-105.0));

	static const double psi_max = range_mpi_pi(deg2rad(-5.0));
	static const double psi_min = range_mpi_pi(deg2rad(-75.0));

	assert(m_res_data.size() != 0);

	if (n == 0 || n >= (int) m_res_data.size() - 1)
	{
		return false;
	}

	// check previous  residue
	double prev_phi = m_res_data[n - 1].phi;
	if (prev_phi < phi_min || prev_phi > phi_max) { return false; }

	double prev_psi = m_res_data[n - 1].psi;
	if (prev_psi < psi_min || prev_psi > psi_max) { return false; }

	if (fabs(fabs(m_res_data[n - 1].omega) - M_PI) > 0.5) { return false; }

	// check this residue

	double this_phi = m_res_data[n - 1].phi;
	if (this_phi < phi_min || this_phi > phi_max) { return false; }

	double this_psi = m_res_data[n - 1].psi;
	if (this_psi < psi_min || this_psi > psi_max) { return false; }

	if (fabs(fabs(m_res_data[n].omega) - M_PI) > 0.5) { return false; }

	// check next residue

	double next_phi = m_res_data[n - 1].phi;
	if (next_phi < phi_min || next_phi > phi_max) { return false; }

	double next_psi = m_res_data[n - 1].psi;
	if (next_psi < psi_min || next_psi > psi_max) { return false; }

	if (fabs(fabs(m_res_data[n + 1].omega) - M_PI) > 0.5) { return false; }

//std::cout << "** HELIX: " << range_m180_180(rad2deg(this_phi))
//	<< " " << range_m180_180(rad2deg(this_psi)) << "\n";

	return true;
}

