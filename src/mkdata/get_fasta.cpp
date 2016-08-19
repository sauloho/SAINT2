#include <cstdlib> // PG added this
#include <iostream>
#include <cstring>
#include <string>
#include "sequence.h"
#include "peptide.h"

using namespace std;

int main(int argc, const char *argv[])
{
	if (argc != 2 && argc != 3)
	{
		cerr << "Usage: " << argv[0] << " pdbfile [chain]\n";
		exit(1);
	}

	const char *filename = argv[1];
	char chain = ' ';

	if (argc == 3)
	{
		if (strlen(argv[2]) != 1)
		{
			cerr << "Chain must be a single character\n";
			exit(1);
		}

		chain = argv[2][0];
	}

	Peptide p;

	if (!p.read_pdb(filename, chain))
	{
		cerr << "Failed to read " << filename << "\n";
		exit(1);
	}

	const char *pdb = strrchr(filename, '/');

	if (pdb == NULL) { pdb = filename; }
	else { pdb++; }

	std::cout << "> seq_" << p.chain() << ' '
		<< std::string(pdb, 4) << '\n';
	Sequence seq;
	seq.create_from_peptide(p);
	seq.print(std::cout, 60);
	std::cout << '\n';
	
	return 0;
}

