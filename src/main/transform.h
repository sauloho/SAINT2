#ifndef TRANSFORM_H_INCLUDED
#define TRANSFORM_H_INCLUDED

#include <iostream>
#include <cassert>
#include "matrix.h"
#include "point.h"

/// @brief A transformation matrix (used for rotating a point in 3D).
///
/// The matrix looks like:
///
/// A B C X
/// D E F Y
/// G H J Z
//  0 0 0 1
///
/// but the last row is not stored, since it is constant.
///
/// The top left 3x3 matrix is a rotation, and (X Y Z) is a translation.

class Transform
{
	friend std::ostream& operator << (std::ostream &out, const Transform &m);

public:
	/// @brief Matrix size.
	static const int Mat_Cols = 4;
	static const int Mat_Rows = 3;

	/// @brief Constructor (does not initialise values).
	Transform();

	/// @brief Copy constructor.
	Transform(const Transform &other);

	/// @brief Assignment operator.
	Transform& operator = (const Transform &other);

	/// @brief Calls find_rotation() to create the transformation.
	Transform(const Point &p1, const Point &p2, double angle);

	/// @brief Create transformation with the specified rotation and
	/// displacement.
	void create_transform(const Matrix_3_3 &rotation, Point displacement);

	/// @brief Set to the identity matrix.
	void set_to_identity();

	/// @brief Multiply another matrix by this one.
	///
	/// @param m Other matrix.
	/// @param result [out] This matrix times \a m.
	void multiply(const Transform &m, Transform *result) const;

	/// @brief Same as multiply(), but returns the matrix.
	Transform times(const Transform &m);

	/// @brief Multiply a point (or vector) by this matrix.
	///
	/// @param p Point/vector.
	/// @param result [out] This matrix times \a p.
	void multiply(const Point &p, Point *result) const;

	/// @brief Same as multiply(), but returns the point/vector.
	///
	/// @param p Point/vector.
	/// @return This matrix times \a p.
	Point times(const Point &p) const;

	/// @brief Access an element of the matrix, eg. mat(0, 2).
	/// The result is an lvalue, so can be used in an assignment statement.
	///
	/// @param row Row (0-2).
	/// @param col Column (0-2).
	/// @return Value at [row][col].
	double &operator() (int row, int col)
	{ 
		assert(col >= 0 && col < Mat_Cols && row >= 0 && row < Mat_Rows);
		return m_val[row][col];
	}

	/// @brief As above (const version).
	const double &operator() (int row, int col) const
	{ 
		assert(col >= 0 && col < Mat_Cols && row >= 0 && row < Mat_Rows);
		return m_val[row][col];
	}

	/// @breief Find the rotation matrix about a line through two points.
	/// The angle (in radians) is clockwise relative to the direction the
	/// vector from p1 to p2 is pointing (eg. if the direction is (0, 0, -1),
	/// it will be clockwise looking down on the xy plane).
	void find_rotation(const Point &p1, const Point &p2, double angle);

	// find the transformation which moves point p1 to q1, changes
	// vector directions (p1-p2) to (q1-q2), and ensures that
	// the plane of (p1-p2-p3) corresponds to the plane of (q1-q2-q3).

	void find_alignment(Point p1, Point p2, Point p3,
		Point q1, Point q2, Point q3);

	// returns true if the rotational component of the transformation
	// is orthonormal (which should always be true), ie. each row
	// and column treated as a vector has a length of exactly one
	bool orthonormal() const;

private:
	/// @brief Copy another matrix to this one.
	void copy(const Transform &other);

private:
	/// Matrix values.
	double m_val[Mat_Rows][Mat_Cols];
};

/// @brief Print the matrix on a single line;
/// the format is like "[ 1 2 3 4 ; 5 6 7 8 ; 9 10 11 12 ]".
std::ostream& operator << (std::ostream &out, const Transform &m);

#endif // TRANSFORM_H_INCLUDED

