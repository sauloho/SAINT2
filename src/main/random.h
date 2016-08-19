#ifndef RANDOM_H_INCLUDED
#define RANDOM_H_INCLUDED

/// @brief A wrapper for the standard C random number functions.
///
/// Use this class in preference to calling rand(), etc. directly.
///
/// Example:
/// <pre>
/// Random::set_seed(time(0));
/// int n = Random::rnd(10);	// a random number from 0 to 9
/// </pre>

class Random
{
public:
	/// @brief Set the random number seed. Call this before using the
	/// other random number functions.
	static void set_seed(long seed);

	/// @brief Get the random number seed.
	static long get_seed();

	/// @brief Generate a random double (0.0 <= n < \a max_val).
	static double rnd(double max_val);

	/// @brief Generate a random integer (0 <= n < \a max_val).
	static int rnd(int max_val);

private:
	/// The random number seed.
	static long m_seed;
};

#endif // RANDOM_H_INCLUDED

