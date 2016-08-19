#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <iostream>
#include <cassert>
#include "point.h"

/// @brief A 3x3 matrix (used for rotating a point in 3D).

class Matrix_3_3
{
	friend std::ostream& operator << (std::ostream &out, const Matrix_3_3 &m);

public:
	/// @brief Matrix size.
	static const int Mat_Size = 3;

	/// @brief Constructor (does not initialise values).
	Matrix_3_3();

	/// @brief Constructor from three vectors. The x, y and z values of
	/// each vector (Point) correspond to row 0, 1 and 2 of the matrix
	/// respectively.
	///
	/// @param a First column.
	/// @param b Second column.
	/// @param c Third column.
	Matrix_3_3(const Point &a, const Point &b, const Point &c);

	/// @brief Copy constructor.
	Matrix_3_3(const Matrix_3_3 &other);

	/// @brief Assignment operator.
	Matrix_3_3& operator = (const Matrix_3_3 &other);

	/// @brief Clear the matrix to all zeroes.
	void set_to_zero();

	/// @brief Set to the identity matrix.
	void set_to_identity();

	/// @brief Calculate the determinant of the matrix.
	///
	/// @return The determinant.
	double determinant() const;

	/// @brief Find the inverse of this matrix.
	///
	/// @param result [out] Inverse.
	/// @param determinant_val The determinant (if zero, the value is
	/// calculated)
	void find_inverse(Matrix_3_3 *result,
		double determinant_val = 0.0) const;

	/// @brief Multiply another matrix by this one.
	///
	/// @param m Other matrix.
	/// @param result [out] This matrix times \a m.
	void multiply(const Matrix_3_3 &m, Matrix_3_3 *result) const;

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
		assert(col >= 0 && col < Mat_Size && row >= 0 && row < Mat_Size);
		return m_val[row][col];
	}

	/// @brief As above (const version).
	const double &operator() (int row, int col) const
	{ 
		assert(col >= 0 && col < Mat_Size && row >= 0 && row < Mat_Size);
		return m_val[row][col];
	}

	/// @brief Call a set of assert() statements to verify that
	/// all three columns are orthogonal to each other.
	void verify_orthonormal() const;

	// find the rotation matrix about a vector passing through
	// the origin.
	//
	// To rotate a point p about an arbitrary vector (p1 to p2):
	//
	// Matrix_3_3 m;
	// m.find_rotation_matrix(p2.minus(p1), angle);
	// new_p = m.times(p.minus(p1)).plus(p1);
	void find_rotation_matrix(const Point &vec, double angle);

private:
	/// @brief Copy another matrix to this one.
	void copy(const Matrix_3_3 &other);

private:
	/// Matrix values.
	double m_val[Mat_Size][Mat_Size];
};

/// @brief Print the matrix on a single line;
/// the format is like "[ 1 2 3 ; 4 5 6 ; 7 8 9 ]".
///
/// Example: 
/// <pre>
/// Matrix m;
/// m.set_to_identity();
/// std::cout << m << '\n';
/// </pre>
///
/// @param out Output stream.
/// @param m Matrix to print.
std::ostream& operator << (std::ostream &out, const Matrix_3_3 &m);

#endif // MATRIX_H_INCLUDED

