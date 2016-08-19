#include <cstdlib> // PG added this
#include <iostream>
#include <cstring>
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

	p.forget_is_from_pdb();
	p.write_pdb("");

	return 0;	// success
}

