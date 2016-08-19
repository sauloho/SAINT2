#ifndef RIBOSOME_H_INCLUDED
#define RIBOSOME_H_INCLUDED

// forward declarations
class Peptide;

class Ribosome
{
public:
	// constructor
	Ribosome();

	// destructor
	~Ribosome();

	// score a peptide (low scores ate better)
	double score(const Peptide &p, bool verbose = false);

private:
    // disable copy and assignment by making them private
	Ribosome(const Ribosome&);
	Ribosome &operator = (const Ribosome&);
};

#endif // RIBOSOME_INCLUDED
