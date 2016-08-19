#ifndef POINT_H_INCLUDED
#define POINT_H_INCLUDED

#include <iostream>
#include <vector>

/// @brief A 3D point or vector.

struct Point
{
	double x;
	double y;
	double z;

	/// @brief Default constructor (initialised to (0, 0, 0)).
	Point();

	/// @brief Constructor.
	Point(double x_val, double y_val, double z_val);

	// don't need to define copy and assignment,
	// since implicit versions are OK (memberwise copy).

	/// @brief Set to a length of one, keeping the same direction.
	void normalise();

	/// @brief Return a vector with length one, keeping the same direction.
	Point normalised() const;

	/// @brief Find the length of the vector.
	double length() const;

	/// @brief Find the dot product with another vector
	/// (zero if the vectors are perpendicular).
	double dot_product(const Point &p2) const;

	// @brief Find the cross product with another vector
	/// (perpendicular to both of them).
	Point cross_product(const Point &p2) const;

	/// @brief Add a vector to this point.
	void add(const Point &p2);

	/// @brief Get the result of adding a vector to this point.
	Point plus(const Point &p2) const;

	/// @brief Subtract a vector from this point.
	void subtract(const Point &p2);

	/// @brief Get the result of subtracting a vector from this point.
	Point minus(const Point &p2) const;

	/// @brief Negate all three coordinates.
	void negate();

	/// @brief Same point with all three coordinates negated.
	Point negated() const;

	/// @brief Multiply all coordinates by the specified value.
	void scale(double factor);

	bool operator == (const Point &p) const
	{ return (x == p.x && y == p.y && z == p.z); }

	bool operator != (const Point &p) const
	{ return (x != p.x || y != p.y || z != p.z); }

	/// @brief Check if all three coordinates are close to zero
	/// @param tiny_val Distance tolerance.
	bool is_tiny(double tiny_val = 0.001) const;

	/// @brief Check if all three coordinates of this point and p2 are
	/// within a certain distance of each other (eg. as a quick check
	/// to see if the points are approximately equal)
	/// @param dist Distance tolerance.
	bool close_to(const Point &p2, double dist = 0.001) const;

	// get the distance to another point
	double distance(const Point &p2) const;

	// get the distance to another point, rounded down to the nearest
	// integer (faster than distance())
	int int_distance(const Point &p2) const;

	// check if a point is more than a certain distance away
	// (faster than using distance())
	bool further_than(double dist, const Point &p2) const;

	// check if a point is closer than a certain distance away
	// (faster than using distance())
	bool closer_than(double dist, const Point &p2) const;

	// get the parameter value representing the position of the
	// closest point to the line through p1 and p2.
	//
	// A value of 0 means at p1
	// A value of 1 means at p2
	// A value between 0 and 1 means somewhere along the line
	//
	// If closest is not NULL, it is set to the closest point
	// (unless only_if_in_segment is true and the value is not between 0 and 1)
	double closest_point_param(const Point &p1, const Point &p2,
		Point *closest = NULL, bool only_if_in_segment = false) const;

	// find the distance to the closest point on the line
	// through p1 and p2
	double dist_to_line(const Point &p1, const Point &p2) const;
};

inline Point midpoint(const Point &p1, const Point &p2)
{
	return Point((p1.x + p2.x) * 0.5,
				 (p1.y + p2.y) * 0.5,
				 (p1.z + p2.z) * 0.5);
}

/// @brief A std::vector of Points.
typedef std::vector<Point> Point_Vec;

/// @brief Print a point to a stream (format is like "(1.5, 2.3, -0.9")).
/// Example:
/// <pre>
/// Point p(1.5, 2.3, -0.9);
/// std::cout << p << "\n";
/// </pre>

std::ostream& operator << (std::ostream &out, const Point &p);

#endif // POINT_H_INCLUDED

