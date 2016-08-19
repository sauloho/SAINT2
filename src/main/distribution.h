#ifndef DISTRIBUTION_H_INCLUDED
#define DISTRIBUTION_H_INCLUDED

#include <vector>
#include <iostream>

/// @brief For selecting values using a probability distribution.
///
/// Example:
/// <pre>
/// Distribution d;
/// d.add(0.5, 0);
/// d.add(1.0, 1);
/// d.add(0.001, 2);
/// int x = d.select();
/// </pre>
///
/// This example will assign x a value of 0 just under 1/3 of the time
/// (0.5 / (0.5 + 1.0 + 0.001)), 1 just under 2/3 of the time,
/// and 2 a very small amount of the time.

class Distribution
{
public:
	/// Constructor.
	Distribution();

	/// @brief Make the distribution empty.
	void clear();

	/// @brief Add an element to the distribution.
	/// @param probability Proportional probability of selecting this element.
	/// @param value Value returned by select() when this element is chosen.
	void add(double probability, int value);

	/// @brief Get the number of elements in the distribution.
	int num() const
	{ return (int) m_val.size(); }

	/// @brief Select a random element from the distribution.
	/// @return The value of the element (as specified in add()).
	int select();

	/// @brief Print the distribution (for debugging).
	/// Note that sum_to_here values will be undefined until the
	/// probability values are sorted (after select() is called).
	void dump(std::ostream &out = std::cout);

	static void self_test();

private:
	/// @brief Sort m_val by ascending probability, and calculate
	/// m_prob_total and the intermediate prob_sum values.
	void sort_values();

	/// @brief Called by select() to do the actual work.
	int binary_search(double val) const;

private:
	/// A single element in the distribution.
	struct Dist_Element
	{
		/// relative probability
		double prob;

		/// sum of probabilities in m_val up to & including this one
		double sum_to_here;

		/// value returned by select() for this element
		int val;

		Dist_Element()
		{ }

		Dist_Element(double p, int v)
			: prob(p), val(v)
		{ }

		// (required for sorting)
		bool operator < (const Dist_Element &other) const
		{ return prob < other.prob; }
	};

	typedef std::vector<Dist_Element> Element_Vec;

	/// Values and their probabilities. After select() is called,
	/// the values are sorted.
	Element_Vec m_val;

	/// Total of all probabilities in m_val.
	double m_prob_total;

	/// Whether m_val is currently sorted.
	bool m_sorted;
};

typedef std::vector<Distribution> Distribution_Vec;

#endif // DISTRIBUTION_H_INCLUDED
