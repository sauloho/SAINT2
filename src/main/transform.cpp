
#include <cmath>
#include <iostream>
#include "matrix.h"
#include "transform.h"

Transform::Transform()
{
}

Transform::Transform(const Transform &other)
{
	copy(other);
}

Transform& Transform::operator = (const Transform &other)
{
	copy(other);
	return *this;
}

void Transform::copy(const Transform &other)
{
	for (int row = 0;row < Mat_Rows;row++)
	{
		for (int col = 0;col < Mat_Cols;col++)
		{
			m_val[row][col] = other.m_val[row][col];
		}
	}
}

Transform::Transform(const Point &p1, const Point &p2, double angle)
{
	find_rotation(p1, p2, angle);
}

void Transform::create_transform(const Matrix_3_3 &rotation,
	Point displacement)
{
	for (int row = 0;row < 3;row++)
	{
		for (int col = 0;col < 3;col++)
		{
			m_val[row][col] = rotation(row, col);
		}
	}

	m_val[0][3] = displacement.x;
	m_val[1][3] = displacement.y;
	m_val[2][3] = displacement.z;
}

void Transform::set_to_identity()
{
	for (int row = 0;row < Mat_Rows;row++)
	{
		for (int col = 0;col < Mat_Cols;col++)
		{
			m_val[row][col] = (row == col ? 1.0 : 0.0);
		}
	}
}

void Transform::multiply(const Transform &m, Transform *result) const
{
	int row;

	for (row = 0;row < Mat_Rows;row++)
	{
		for (int col = 0;col < Mat_Cols;col++)
		{
			result->m_val[row][col] =
				m_val[row][0] * m.m_val[0][col] +
				m_val[row][1] * m.m_val[1][col] +
				m_val[row][2] * m.m_val[2][col];
		}
	}

	for (row = 0;row < Mat_Rows;row++)
	{
		result->m_val[row][3] += m_val[row][3];
	}
}

void Transform::multiply(const Point &p, Point *result) const
{
	*result = this->times(p);
}

Transform Transform::times(const Transform &m)
{
	Transform result;
	multiply(m, &result);
	return result;
}

Point Transform::times(const Point &p) const
{
	return Point(
	 (m_val[0][0] * p.x + m_val[0][1] * p.y + m_val[0][2] * p.z + m_val[0][3]),
	 (m_val[1][0] * p.x + m_val[1][1] * p.y + m_val[1][2] * p.z + m_val[1][3]),
	 (m_val[2][0] * p.x + m_val[2][1] * p.y + m_val[2][2] * p.z + m_val[2][3]));
}

void Transform::find_rotation(const Point &p1, const Point &p2, double angle)
{
	assert(!p1.close_to(p2));

	// formula from http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	// /ArbitraryAxisRotation.html ("Rotations about an arbitrary line")

	double a = p1.x;
	double b = p1.y;
	double c = p1.z;
	double u = p2.x - p1.x;
	double v = p2.y - p1.y;
	double w = p2.z - p1.z;
	double u2 = u * u;
	double v2 = v * v;
	double w2 = w * w;
	double au = a * u;
	double bv = b * v;
	double cw = c * w;
	double denom = u2 + v2 + w2;
	double rs = sqrt(denom) * sin(angle);
	double ca = cos(angle);
	double ca2 = 1.0 - ca;

	m_val[0][0] = (u2 + (v2 + w2) * ca) / denom;
	m_val[1][0] = ((u * v * ca2) + (w * rs)) / denom;
	m_val[2][0] = ((u * w * ca2) - (v * rs)) / denom;

	m_val[0][1] = ((u * v * ca2) - (w * rs)) / denom;
	m_val[1][1] = (v2 + (u2 + w2) * ca) / denom;
	m_val[2][1] = ((v * w * ca2) + (u * rs)) / denom;

	m_val[0][2] = ((u * w * ca2) + (v * rs)) / denom;
	m_val[1][2] = ((v * w * ca2) - (u * rs)) / denom;
	m_val[2][2] = (w2 + (u2 + v2) * ca) / denom;

	m_val[0][3] = (a * (v2 + w2) - u * (bv + cw) +
		(u * (bv + cw) - a * (v2 + w2)) * ca +
		(b * w - c * v) * rs) / denom;

	m_val[1][3] = (b * (u2 + w2) - v * (au + cw) +
		(v * (au + cw) - b * (u2 + w2)) * ca +
		(c * u - a * w) * rs) / denom;

	m_val[2][3] = (c * (u2 + v2) - w * (au + bv) +
		(w * (au + bv) - c * (u2 + v2)) * ca +
		(a * v - b * u) * rs) / denom;
}

void Transform::find_alignment(Point p1, Point p2, Point p3,
	Point q1, Point q2, Point q3)
{
	// vectors t, u, v

	Point t = p2.minus(p1);
	t.normalise();

	Point v = t.cross_product(p3.minus(p1));
	v.normalise();

	Point u = v.cross_product(t);

	assert(fabs(t.dot_product(u)) < 0.001);
	assert(fabs(t.dot_product(v)) < 0.001);
	assert(fabs(u.dot_product(v)) < 0.001);

	// vectors x, y, z

	Point x = q2.minus(q1);
	x.normalise();

	Point z = x.cross_product(q3.minus(q1));
	z.normalise();

	Point y = z.cross_product(x);

	assert(fabs(x.dot_product(y)) < 0.001);
	assert(fabs(x.dot_product(z)) < 0.001);
	assert(fabs(y.dot_product(z)) < 0.001);

	// put vectors into a matrix

	Matrix_3_3 m1(t, u, v);
	Matrix_3_3 m2(x, y, z);

	Matrix_3_3 m1_inv, rotation;
	m1.find_inverse(&m1_inv);

	m2.multiply(m1_inv, &rotation);
	Point disp = q1.minus(rotation.times(p1));

	create_transform(rotation, disp);

	/*
	Point t_to = rotation->times(t);
	Point u_to = rotation->times(u);
	Point v_to = rotation->times(v);
	assert(t_to.close_to(x));
	assert(u_to.close_to(y));
	assert(v_to.close_to(z));

	verify_transformation(*rotation, *displacement, p1, p2, p3, q1, q2, q3);
	*/
}

inline bool length_one(double a, double b, double c)
{
	return (fabs((a*a + b*b + c*c) - 1.0) < 0.001);
}

bool Transform::orthonormal() const
{
	return (
		length_one(m_val[0][0], m_val[0][1], m_val[0][2]) &&
		length_one(m_val[1][0], m_val[1][1], m_val[1][2]) &&
		length_one(m_val[2][0], m_val[2][1], m_val[2][2]) &&
		length_one(m_val[0][0], m_val[1][0], m_val[2][0]) &&
		length_one(m_val[0][1], m_val[1][1], m_val[2][1]) &&
		length_one(m_val[0][2], m_val[1][2], m_val[2][2]));
}

std::ostream& operator << (std::ostream &out, const Transform &m)
{
	out << '[';

	for (int row = 0;row < Transform::Mat_Rows;row++)
	{
		if (row > 0)
		{
			out << " ;";
		}

		for (int col = 0;col < Transform::Mat_Cols;col++)
		{
			out << ' ' << m.m_val[row][col];
		}
	}

	return out << " ]";
}

#ifdef TEST

#include <iostream>
#include <string>
#include "geom.h"
#include "config.h"
#include "stream_printf.h"

#ifndef M_TWOPI
#define M_TWOPI (M_PI * 2.0)
#endif

// (these variables and functions are here just to stop having to link
// everything)

std::string Config::m_cmd;

double rad2deg(double radians)
{
	return range_0_360(radians * (180.0 / M_PI));
}

double deg2rad(double degrees)
{
	return range_0_2pi(degrees * (M_PI / 180.0));
}

double range_0_360(double val)
{
	if (val < 0.0)
	{
		val = 360.0 - fmod(-val, 360.0);
		return (val == 360.0 ? 0.0 : val);
	}
	else
	if (val >= 360.0)
	{
		return fmod(val, 360.0);
	}
	else
	{
		return val;
	}
}

double range_0_2pi(double val)
{
	if (val < 0.0)
	{
		val = M_TWOPI - fmod(-val, M_TWOPI);
		return (val == M_TWOPI ? 0.0 : val);
	}
	else
	if (val >= M_TWOPI)
	{
		return fmod(val, M_TWOPI);
	}
	else
	{
		return val;
	}
}

void test_transform(const Point &p1, const Point &p2, const Point &p)
{
	std::cout << "\nTransforming " << p << " around vector "
		<< p2.minus(p1) << " from " << p1 << "\n\n";

	for (double a = 0.0;a < 360.0;a += 30.0)
	{
		Transform t(p1, p2, deg2rad(a));

		std::cout << Printf("%.1f", a)
			<< "  "
			<< t.times(p)
			<< "\n";
	}
}

int main(int argc, char **argv)
{
	Point p(3, 4, 5);

	test_transform(p, p.plus(Point(10, 0, 0)), p.plus(Point(0, 1, 0)));
	test_transform(p, p.plus(Point(0, 10, 0)), p.plus(Point(0, 0, 1)));
	test_transform(p, p.plus(Point(0, 0, 10)), p.plus(Point(1, 0, 0)));
	return 0;
}

/* Output:

Transforming (3, 5, 5) around vector (10, 0, 0) from (3, 4, 5)

0.0  (3, 5, 5)
30.0  (3, 4.86603, 5.5)
60.0  (3, 4.5, 5.86603)
90.0  (3, 4, 6)
120.0  (3, 3.5, 5.86603)
150.0  (3, 3.13397, 5.5)
180.0  (3, 3, 5)
210.0  (3, 3.13397, 4.5)
240.0  (3, 3.5, 4.13397)
270.0  (3, 4, 4)
300.0  (3, 4.5, 4.13397)
330.0  (3, 4.86603, 4.5)

Transforming (3, 4, 6) around vector (0, 10, 0) from (3, 4, 5)

0.0  (3, 4, 6)
30.0  (3.5, 4, 5.86603)
60.0  (3.86603, 4, 5.5)
90.0  (4, 4, 5)
120.0  (3.86603, 4, 4.5)
150.0  (3.5, 4, 4.13397)
180.0  (3, 4, 4)
210.0  (2.5, 4, 4.13397)
240.0  (2.13397, 4, 4.5)
270.0  (2, 4, 5)
300.0  (2.13397, 4, 5.5)
330.0  (2.5, 4, 5.86603)

Transforming (4, 4, 5) around vector (0, 0, 10) from (3, 4, 5)

0.0  (4, 4, 5)
30.0  (3.86603, 4.5, 5)
60.0  (3.5, 4.86603, 5)
90.0  (3, 5, 5)
120.0  (2.5, 4.86603, 5)
150.0  (2.13397, 4.5, 5)
180.0  (2, 4, 5)
210.0  (2.13397, 3.5, 5)
240.0  (2.5, 3.13397, 5)
270.0  (3, 3, 5)
300.0  (3.5, 3.13397, 5)
330.0  (3.86603, 3.5, 5)
*/

#endif

