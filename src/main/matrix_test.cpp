
#include <iostream>
#include "point.h"
#include "matrix.h"

void test_cross(const Point &p1, const Point &p2)
{
	std::cout << p1 << " x " << p2 << " = "
		<< p1.find_cross_product(p2) << "\n";
}

int main(int argc, char **argv)
{
	Point p(5.0, 7.0, -5.0);

	std::cout << p << ", length " << p.length() << "\n";
	p.normalise();
	std::cout << "normalised = " << p << ", length " << p.length() << "\n\n";

	test_cross(Point(5, 3, 0), Point(2, 7, 0));
	test_cross(Point(2, 7, 0), Point(5, 3, 0));
	test_cross(Point(0, 0, 1), Point(0, 1, 0));
	test_cross(Point(0, 0, 1), Point(1, 0, 0));
	test_cross(Point(0, 1, 0), Point(0, 0, 1));
	test_cross(Point(0, 1, 0), Point(1, 0, 0));
	test_cross(Point(1, 0, 0), Point(0, 0, 1));
	test_cross(Point(1, 0, 0), Point(0, 1, 0));
	test_cross(Point(-1, 0, 0), Point(0, 1, 0));

	std::cout << "\n";
	Matrix_3_3 m;

	m.set_to_zero();
	std::cout << "Zero = " << m << "\n";
	m.set_to_identity();
	std::cout << "Id   = " << m << "\n";

	m(0, 0) = 5.0;  m(0, 1) = 8.0; m(0, 2) = 9.0;
	m(1, 0) = 7.0;  m(1, 1) = 2.0; m(1, 2) = 3.0;
	m(2, 0) = 11.0; m(2, 1) = 1.0; m(2, 2) = 6.0;

	std::cout << "\nDet " << m << " = " << m.determinant()
		<< " (should be -162)\n";

	Matrix_3_3 r, r2;
	m.find_inverse(&r, 1.0);	// pretend the determinant is 1
	std::cout << "Inverse = " << r << " / " << m.determinant() << "\n"
			  << "Expected: [ 9 -39 6 ; -9 -69 48 ; -15 83 -46 ] / -162\n\n";

	m(0, 0) = 2.0;  m(0, 1) = -3.0;  m(0, 2) = -5.0;
	m(1, 0) = 7.0;  m(1, 1) = 11.0; m(1, 2) = 13.0;
	m(2, 0) = 17.0; m(2, 1) = 19.0; m(2, 2) = 23.0;

	std::cout << "Det " << m << " = " << m.determinant()
		<< " (should be 102)\n";
	m.find_inverse(&r, 1.0);	// pretend the determinant is 1
	std::cout << "Inverse = " << r << " / " << m.determinant() << "\n"
			  << "Expected: [ 6 -26 16 ; 60 131 -61 ; -54 -89 43 ] / 102\n\n";

	// use proper determinant
	m.find_inverse(&r);
	std::cout << "M = " << m << "\n";
	std::cout << "R = " << r << "\n";
	m.multiply(r, &r2);
	std::cout << "m * inverse(m) = " << r2 << "\n\n";

	Point p2;
	p = Point(1.0, 2.0, 3.0);
	std::cout << "Transform p " << p;
	m.multiply(p, &p2);
	std::cout << " => " << p2 << "\n";

	return 0;
}

