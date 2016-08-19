
#include <cmath>
#include <iostream>
#include "matrix.h"

Matrix_3_3::Matrix_3_3()
{
}

Matrix_3_3::Matrix_3_3(const Point &a, const Point &b, const Point &c)
{
	m_val[0][0] = a.x;
	m_val[1][0] = a.y;
	m_val[2][0] = a.z;

	m_val[0][1] = b.x;
	m_val[1][1] = b.y;
	m_val[2][1] = b.z;

	m_val[0][2] = c.x;
	m_val[1][2] = c.y;
	m_val[2][2] = c.z;
}

Matrix_3_3::Matrix_3_3(const Matrix_3_3 &other)
{
	copy(other);
}

Matrix_3_3& Matrix_3_3::operator = (const Matrix_3_3 &other)
{
	copy(other);
	return *this;
}

void Matrix_3_3::copy(const Matrix_3_3 &other)
{
	for (int row = 0;row < Mat_Size;row++)
	{
		for (int col = 0;col < Mat_Size;col++)
		{
			m_val[row][col] = other.m_val[row][col];
		}
	}
}

void Matrix_3_3::set_to_zero()
{
	for (int row = 0;row < Mat_Size;row++)
	{
		for (int col = 0;col < Mat_Size;col++)
		{
			m_val[row][col] = 0.0;
		}
	}
}

void Matrix_3_3::set_to_identity()
{
	for (int row = 0;row < Mat_Size;row++)
	{
		for (int col = 0;col < Mat_Size;col++)
		{
			m_val[row][col] = (row == col ? 1.0 : 0.0);
		}
	}
}

void Matrix_3_3::find_inverse(Matrix_3_3 *result,
	double determinant_val /*= 0.0*/) const
{
	if (determinant_val == 0.0)
	{
		determinant_val = determinant();
	}

	assert(fabs(determinant_val) > 0.000000001);

	double det1 = 1.0 / determinant_val;

	// formula derived from http://www.dr-lex.be/random/matrix_inv.html

	result->m_val[0][0] =
		((m_val[2][2] * m_val[1][1]) - (m_val[2][1] * m_val[1][2])) * det1;
	result->m_val[1][0] =
		((m_val[2][0] * m_val[1][2]) - (m_val[2][2] * m_val[1][0])) * det1;
	result->m_val[2][0] =
		((m_val[2][1] * m_val[1][0]) - (m_val[2][0] * m_val[1][1])) * det1;

	result->m_val[0][1] =
		((m_val[2][1] * m_val[0][2]) - (m_val[2][2] * m_val[0][1])) * det1;
	result->m_val[1][1] =
		((m_val[2][2] * m_val[0][0]) - (m_val[2][0] * m_val[0][2])) * det1;
	result->m_val[2][1] =
		((m_val[2][0] * m_val[0][1]) - (m_val[2][1] * m_val[0][0])) * det1;

	result->m_val[0][2] =
		((m_val[1][2] * m_val[0][1]) - (m_val[1][1] * m_val[0][2])) * det1;
	result->m_val[1][2] =
		((m_val[1][0] * m_val[0][2]) - (m_val[1][2] * m_val[0][0])) * det1;
	result->m_val[2][2] =
		((m_val[1][1] * m_val[0][0]) - (m_val[1][0] * m_val[0][1])) * det1;
}

double Matrix_3_3::determinant() const
{
	/*
	std::cout << "\ndeterminant():\n"
		<< m_val[0][0] << " * (("
			<< m_val[2][2] << " * " << m_val[1][1] << ") - ("
			<< m_val[2][1] << " * " <<  m_val[1][2] << ")) +\n"
		<< m_val[1][0] << " * (("
			<< m_val[2][1] << " * " << m_val[0][2] << ") - ("
			<< m_val[2][2] << " * " << m_val[0][1] << ")) +\n"
		<< m_val[2][0] << " * (("
			<< m_val[1][2] << " * " << m_val[0][1] << ") - ("
			<< m_val[1][1] << " * " << m_val[0][2] << "))\n";
	*/

	// formula derived from http://www.dr-lex.be/random/matrix_inv.html

	return m_val[0][0] *
			((m_val[2][2] * m_val[1][1]) - (m_val[2][1] * m_val[1][2])) +
		m_val[1][0] *
			((m_val[2][1] * m_val[0][2]) - (m_val[2][2] * m_val[0][1])) +
		m_val[2][0] *
			((m_val[1][2] * m_val[0][1]) - (m_val[1][1] * m_val[0][2]));
}

void Matrix_3_3::multiply(const Matrix_3_3 &m, Matrix_3_3 *result) const
{
	for (int row = 0;row < Mat_Size;row++)
	{
		for (int col = 0;col < Mat_Size;col++)
		{
			result->m_val[row][col] =
				m_val[row][0] * m.m_val[0][col] +
				m_val[row][1] * m.m_val[1][col] +
				m_val[row][2] * m.m_val[2][col];
		}
	}
}

void Matrix_3_3::multiply(const Point &p, Point *result) const
{
	*result = this->times(p);
}

Point Matrix_3_3::times(const Point &p) const
{
	return Point(
		m_val[0][0] * p.x + m_val[0][1] * p.y + m_val[0][2] * p.z,
		m_val[1][0] * p.x + m_val[1][1] * p.y + m_val[1][2] * p.z,
		m_val[2][0] * p.x + m_val[2][1] * p.y + m_val[2][2] * p.z);
}

void Matrix_3_3::verify_orthonormal() const
{
	// a, b and c are the three columns of the matrix
	Point a(m_val[0][0], m_val[1][0], m_val[2][0]);
	Point b(m_val[0][1], m_val[1][1], m_val[2][1]);
	Point c(m_val[0][2], m_val[1][2], m_val[2][2]);

	assert(fabs(a.length() - 1.0) < 0.01);
	assert(fabs(b.length() - 1.0) < 0.01);
	assert(fabs(c.length() - 1.0) < 0.01);

	assert(fabs(a.dot_product(b) < 0.01));
	assert(fabs(a.dot_product(c) < 0.01));
	assert(fabs(b.dot_product(c) < 0.01));
}

void Matrix_3_3::find_rotation_matrix(const Point &vec, double angle)
{
	assert(!vec.is_tiny());

	// formula from http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	// /ArbitraryAxisRotation.html ("Rotations about the origin")

	double u = vec.x;
	double v = vec.y;
	double w = vec.z;
	double u2 = u * u;
	double v2 = v * v;
	double w2 = w * w;
	double denom = u2 + v2 + w2;
	double rs = sqrt(denom) * sin(angle);
	double c = cos(angle);
	double cc = 1.0 - c;

	m_val[0][0] = (u2 + (v2 + w2) * c) / denom;
	m_val[1][0] = ((u * v * cc) + (w * rs)) / denom;
	m_val[2][0] = ((u * w * cc) - (v * rs)) / denom;

	m_val[0][1] = ((u * v * cc) - (w * rs)) / denom;
	m_val[1][1] = (v2 + (u2 + w2) * c) / denom;
	m_val[2][1] = ((v * w * cc) + (u * rs)) / denom;

	m_val[0][2] = ((u * w * cc) + (v * rs)) / denom;
	m_val[1][2] = ((v * w * cc) - (u * rs)) / denom;
	m_val[2][2] = (w2 + (u2 + v2) * c) / denom;
}

std::ostream& operator << (std::ostream &out, const Matrix_3_3 &m)
{
	out << '[';

	for (int row = 0;row < Matrix_3_3::Mat_Size;row++)
	{
		if (row > 0)
		{
			out << " ;";
		}

		for (int col = 0;col < Matrix_3_3::Mat_Size;col++)
		{
			out << ' ' << m.m_val[row][col];
		}
	}

	return out << " ]";
}

