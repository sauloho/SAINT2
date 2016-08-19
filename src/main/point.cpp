
#include <cmath>
#include <cassert>
#include "point.h"

Point::Point()
	: x(0.0), y(0.0), z(0.0)
{
}

Point::Point(double x_val, double y_val, double z_val)
	: x(x_val), y(y_val), z(z_val)
{
}

void Point::normalise()
{
	double len = length();

	if (!(fabs(len) > 0.000000000001))
	{
		int z;
		z++;
	}

	assert(fabs(len) > 0.000000000001);

	double one_over_len = 1.0 / len;

	x *= one_over_len;
	y *= one_over_len;
	z *= one_over_len;
}

Point Point::normalised() const
{
	Point p = *this;
	p.normalise();
	return p;
}

double Point::length() const
{
	return sqrt(x * x + y * y + z * z);
}

double Point::dot_product(const Point &p2) const
{
	return ((x * p2.x) + (y * p2.y) + (z * p2.z));
}

Point Point::cross_product(const Point &p2) const
{
	return Point(
		(y * p2.z) - (z * p2.y),
		(z * p2.x) - (x * p2.z),
		(x * p2.y) - (y * p2.x));
}

void Point::add(const Point &p2)
{
	x += p2.x;
	y += p2.y;
	z += p2.z;
}

Point Point::plus(const Point &p2) const
{
	return Point(x + p2.x, y + p2.y, z + p2.z);
}

void Point::subtract(const Point &p2)
{
	x -= p2.x;
	y -= p2.y;
	z -= p2.z;
}

Point Point::minus(const Point &p2) const
{
	return Point(x - p2.x, y - p2.y, z - p2.z);
}

void Point::negate()
{
	x = -x;
	y = -y;
	z = -z;
}

Point Point::negated() const
{
	return Point(-x, -y, -z);
}

void Point::scale(double factor)
{
	x *= factor;
	y *= factor;
	z *= factor;
}

bool Point::is_tiny(double tiny_val /*= 0.001*/) const
{
	return (fabs(x) < tiny_val && fabs(y) < tiny_val && fabs(z) < tiny_val);
}

bool Point::close_to(const Point &p2, double dist /*= 0.001*/) const
{
	return (fabs(x - p2.x) < dist &&
			fabs(y - p2.y) < dist &&
			fabs(z - p2.z) < dist);
}

inline double square(double x)
{
	return (x * x);
}

double Point::distance(const Point &p2) const
{
	return sqrt(square(x - p2.x) + square(y - p2.y) + square(z - p2.z));
}

class IntDist
{
public:
	static const int Max_Dist_Squared = 10000;
	int m_sqrt[Max_Dist_Squared];

	IntDist()
	{
		int s = 0;
		for (int n = 0;n < Max_Dist_Squared;n++)
		{
			if (square(s + 1) == n) { s++; }

			m_sqrt[n] = s;
		}
	}

	int distance(const Point &p1, const Point &p2)
	{
		int d_squared = (int) (
			square(p1.x - p2.x) +
			square(p1.y - p2.y) +
			square(p1.z - p2.z));

		if (d_squared < Max_Dist_Squared)
		{
			return m_sqrt[d_squared];
		}

		return (int) sqrt((double) d_squared);
	}
};

int Point::int_distance(const Point &p2) const
{
	static IntDist intDist;
	return intDist.distance(*this, p2);
}

bool Point::further_than(double dist, const Point &p2) const
{
	double dx, dy, dz;
	
	if ((dx = fabs(x - p2.x)) > dist ||
		(dy = fabs(y - p2.y)) > dist ||
		(dz = fabs(z - p2.z)) > dist)
	{ return true; }

	return (square(dx) + square(dy) + square(dz) > square(dist));
}

bool Point::closer_than(double dist, const Point &p2) const
{
	double dx, dy, dz;
	
	if ((dx = fabs(x - p2.x)) > dist ||
		(dy = fabs(y - p2.y)) > dist ||
		(dz = fabs(z - p2.z)) > dist)
	{ return false; }

	return (square(dx) + square(dy) + square(dz) < square(dist));
}

double Point::closest_point_param(const Point &p1, const Point &p2,
	Point *closest /*= NULL*/, bool only_if_in_segment /*= false*/) const
{
	Point p1_to_this = this->minus(p1);
	Point p1_to_p2 = p2.minus(p1);
	double dist_p1_p2_squared =
		square(p1.x - p2.x) + square(p1.y - p2.y) + square(p1.z - p2.z);

	assert(dist_p1_p2_squared != 0.0);
	double param = p1_to_this.dot_product(p1_to_p2) / dist_p1_p2_squared;

	if (closest != NULL &&
		!(only_if_in_segment && (param < 0.0 || param > 1.0)))
	{
		closest->x = p2.x * param + p1.x * (1.0 - param);
		closest->y = p2.y * param + p1.y * (1.0 - param);
		closest->z = p2.z * param + p1.z * (1.0 - param);
	}

	return param;
}

double Point::dist_to_line(const Point &p1, const Point &p2) const
{
	Point closest;
	(void) closest_point_param(p1, p2, &closest);
	return this->distance(closest);
}

std::ostream& operator << (std::ostream &out, const Point &p)
{
	return out << '(' << p.x << ", " << p.y << ", " << p.z << ')';
}

