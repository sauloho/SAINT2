#ifndef LENNARD_JONES_H_INCLUDED
#define LENNARD_JONES_H_INCLUDED

class Peptide;

/**
 * @brief Lennard-Jones (short range) potential.
 */
class Lennard_Jones
{
	friend class LJ_Static_Init;
public:
	Lennard_Jones();
	~Lennard_Jones();

	// score a peptide (low scores are better)
	double score(const Peptide& peptide, bool verbose = false);

	// check if any single LJ score is beyond a certain threshold
	bool steric_clash(const Peptide &peptide,
		double max_total = -1.0	// maximum total score (-1 means no limit)
	);

protected:
	struct LJ_Params
	{
		double c12;
		double c6;
		double max_dist;

		LJ_Params();
		LJ_Params(double c12_val, double c6_val);
	};

	static LJ_Params m_lj[Num_Backbone][Num_Backbone];
};

#endif // LENNARD_JONES_H_INCLUDED
