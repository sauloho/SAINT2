
// reads all of the PDB files in a list (eg. a list of PISCES chains)
// and outputs a list of the ones with no errors and length >= 20

#include <iostream>
#include <cmath>
#include <cstring>
#include <cassert>
#include "sequence.h"
#include "static_init.h"
#include "peptide.h"
#include "c_file.h"
#include "point.h"
#include "geom.h"

void process_file(const char *filename, const char *pdb_id, char chain_id)
{
	Peptide p;

	// if (!p.read_pdb(filename, chain_id, true))
	if (!p.read_pdb(filename, chain_id, false))
	{
		std::cerr << "Ignoring: " << pdb_id << " " << chain_id << "\n";
	}
	else
	if (p.length() < 20)
	{
		std::cerr << "Short: " << pdb_id << " " << p.chain() << "\n";
	}
	else
	{
		int missing_res = p.num_missing_res();
		int missing_atoms = p.num_missing_backbone();

		std::cout
			<< "PDB: "
			<< pdb_id << " " << p.chain()
			<< " " << p.length()
			<< " " << missing_res
			<< " " << missing_atoms;

		if (missing_res == 0 && missing_atoms == 0)
		{
			std::cout << " ok";
		}

		std::cout << "\n";
		std::cout.flush();
	}
}

int main(int argc, const char *argv[])
{
	if (argc != 2 && argc != 3)
	{
		std::cerr << "Usage: " << argv[0]
			<< " <pdb file> [<chain>] | <dir> <PDB chain list file>\n"
			<< "\nChain list file contains lines of the form:\n\n"
			<< "<PDB id> <chain>\n\n";
		exit(1);
	}

	bool from_file;
	char chain_id;

	if (argc == 2)
	{
		from_file = false;
		chain_id = ' ';
	}
	else
	{
		if (strlen(argv[2]) == 1)
		{
			from_file = false;
			chain_id = argv[2][0];
		}
		else
		{
			from_file = true;
		}
	}

	if (from_file)
	{
		const char *dir = argv[1];
		C_File file(argv[2], "r", "Chain list file");
		char s[80], pdb_id[80];
		char pdb_filename[1000];
		int line_num = 0;

		while (fgets(s, 80, file) != NULL)
		{
			line_num++;

			chain_id = ' ';
			if (sscanf(s, "%s %c", pdb_id, &chain_id) < 1)
			{
				std::cerr << "Error on line " << line_num
					<< " of " << file.name()
					<< "\n";
				exit(1);
			}

			sprintf(pdb_filename, "%s/%s.pdb", dir, pdb_id);
			process_file(pdb_filename, pdb_id, chain_id);
		}
	}
	else
	{
		char pdb_id[5];

		// find PDB id in filename
		const char *p = strrchr(argv[1], '/');

		if (p == NULL)
		{
			p = argv[1];
		}
		else
		{
			p++;
		}

		strncpy(pdb_id, p, 4);
		pdb_id[4] = '\0';
		process_file(argv[1], pdb_id, chain_id);
	}

	return 0;
}

