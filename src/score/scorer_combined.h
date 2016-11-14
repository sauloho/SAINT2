
#ifndef SCORER_COMBINED_H_INCLUDED
#define SCORER_COMBINED_H_INCLUDED

// The default type of Scorer.

#include "scorer.h"
#include <iostream>

// a peptide is "short" if it is less than or equal to this length
#define SHORT_PEPTIDE 150

// forward declarations
class Config;
class Peptide;
class Atom;
class Solvation;
class Orientation;
class Lennard_Jones;
class RAPDF;
class HBond;
class Saulo;
class Eleanor;
class CORE;
class PredSS;
class Rgyr;
class Contact;
class Crowding;
class Randomscr;
class Torsion;
class PredTor;
class Ribosome;

enum Score_Term
{
	SC_SOLV, SC_ORIENT, SC_LJ, SC_RAPDF, SC_HBOND, SC_SAULO, SC_ELEANOR, SC_CORE, SC_PREDSS, SC_RGYR, SC_CONTACT, SC_CROWD, SC_RANDSCR, SC_TOR, SC_PREDTOR,
	SC_RIBO, SC_NUM
};

class Scorer_Combined : public Scorer
{
public:
	// constructor
	Scorer_Combined();

	// destructor
	virtual ~Scorer_Combined();

	/// Parse a config file parameter. Returns false if the parameter
	/// is not recognised.
	virtual bool parse_parameter(const std::string &name,
		const std::string &value);

	// verify that the parameters are consistent and complete
	// (otherwise exits with an error message)
	virtual void verify_parameters();

	// score a peptide (low scores ate better)
	virtual double score(const Peptide &p, double progress = 1.0,
		double *progress1_score = NULL);

	// print a brief description of the type of scoring
	virtual void print_desc(std::ostream &out);

    // print sample config file parameters
    static void print_template(std::ostream &out, bool commented = true);

    // "type" value in config file
    static const char *type()
    { return c_type; }

	// set the "long peptide" weight for one of the scoring terms
	void set_long_weight(Score_Term s, double val);

	// set the "short peptide" weight for one of the scoring terms
	void set_short_weight(Score_Term s, double val);

	// returns true if the ribosome weight is not zero
	bool ribosome_wall() const
	{ return (m_weight[SC_RIBO] != 0.0); }

	// dump the values used for scoring
	void dump(std::ostream &out = std::cout);

private:
    // disable copy and assignment by making them private
	Scorer_Combined(const Scorer_Combined&);
	Scorer_Combined &operator = (const Scorer_Combined&);

private:
	// "type" parameter name
	static const char *c_type;

	static const char *c_score_name[SC_NUM];

	// parameter names

	//static const char *c_param_rapdf_type;
	//static const char *c_rapdf_default_type;

	static const char *c_param_filename[SC_NUM];
	static const char *c_param_short_filename[SC_NUM];

	static const char *c_param_weight[SC_NUM];
	static const char *c_param_short_weight[SC_NUM];

	RAPDF *m_rapdf;
	Solvation *m_solvation;
	Torsion *m_torsion;
	PredTor *m_predtor;
	HBond *m_hbond;
	Saulo *m_saulo;
	Eleanor *m_eleanor;
	CORE *m_core;
	PredSS *m_predss;
	Rgyr *m_rgyr;
	Contact *m_contact;
	Crowding *m_crowding;
	Randomscr *m_randomscr;	
	Orientation *m_orientation;
	Lennard_Jones *m_lj;
	Ribosome *m_ribosome;
	

	double m_weight[SC_NUM];
	double m_short_weight[SC_NUM];
};

#endif // SCORER_COMBINED_INCLUDED
