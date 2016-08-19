
// input format:
// 1f86,A,10,9,cplmvkvld,2ilk,A,152,9,aymtmkirn,89.38,0,89.47,10,4.53129

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <cstring>
#include <map>
#include "c_file.h"

using namespace std;

#define PDB_HOME "/home/jellis/data/pdb/data/structures/divided/pdb/"

struct PDB_key
{
	std::string pdb;
	char chain;

	PDB_key() : chain(' ')
	{
	}

	PDB_key(const char *pdb_val, char chain_val)
		: pdb(pdb_val), chain(chain_val)
	{
	}
		
	bool operator == (const PDB_key &other) const
	{
		return (pdb == other.pdb && chain == other.chain);
	}

	bool operator < (const PDB_key &other) const
	{
		if (pdb < other.pdb) { return true; }
		if (pdb > other.pdb) { return false; }
		return (chain < other.chain);
	}
};

struct PDB_data
{
	int start_res;
	int end_res;

	PDB_data()
	{
	}

	PDB_data(int start, int end)
		: start_res(start), end_res(end)
	{
	}
};

typedef map<PDB_key,PDB_data> PDB_Map;
PDB_Map pdb_map;

void dump_pdb_map()
{
	cerr << pdb_map.size() << " PDBs used:";

	int count = 0;

	for (PDB_Map::const_iterator i = pdb_map.begin();i != pdb_map.end();++i)
	{
		cerr << ' ' << i->first.pdb << i->first.chain;
		
		if (++count > 50) { cerr << " ... etc ..."; break; }
	}

	cerr << '\n';
}

void get_pdb_values(const char *pdb, char chain, int *start_res,
	int *end_res, const char *override_filename = NULL)
{
	string filename;
	
	if (override_filename == NULL)
	{
		filename = PDB_HOME;
		filename += string(pdb + 1, 2);
		filename += "/pdb";
		filename += pdb;
		filename += ".ent";
	}
	else
	{
		filename = override_filename;
	}

	C_File file(filename.c_str(), "r", "PDB file");

	// cout << filename << "\n";

	char buffer[1000];
	bool first_res_found = false;

	while (file.next_line(buffer, 1000))
	{
		if (buffer[21] == chain &&
			buffer[0] == 'A' &&
			buffer[1] == 'T' &&
			buffer[2] == 'O' &&
			buffer[3] == 'M' &&
			buffer[4] == ' ')
		{
			int res;

			if (sscanf(buffer + 22, "%d", &res) == 1)
			{
				if (!first_res_found)
				{
					*start_res = res;
					first_res_found = true;
				}

				*end_res = res;
			}
		}
	}

	// cout << "@@ " << pdb << " " << chain << " = "
	//	<< *start_res << " to " << *end_res << "\n";
}

bool fragment_ok(const char *pdb, char chain, int pos, int len,
	bool at_start_of_target, bool at_end_of_target)
{
	PDB_key key(pdb, chain);
	PDB_Map::iterator i = pdb_map.find(key);

	if (i == pdb_map.end())
	{
		int start_res, end_res;

		get_pdb_values(pdb, chain, &start_res, &end_res);
		pdb_map[key] = PDB_data(start_res, end_res);
		i = pdb_map.find(key);
	}

	bool at_start = (pos <= i->second.start_res);
	bool at_end = (pos + len - 1 >= i->second.end_res);

	if ((at_start && !at_start_of_target) ||
		(at_end && !at_end_of_target))
	{
		return false;
	}

	/*
	cerr << "OK: " << pdb << " " << chain << " "
		<< i->second.start_res << " < pos " << pos
		<< " & pos + len " << pos + len - 1 
		<< " < " << i->second.end_res
		<< "\n";
	*/

	return true;

	// return (pos > i->second.start_res && pos + len - 1 < i->second.end_res);
}

void replace(char *s, char from, char to)
{
	for ( ;*s;s++)
	{
		if (*s == from) { *s = to; }
	}
}

int main(int argc, char **argv)
{
	if ((argc != 4 && argc != 5 && argc != 6) || strlen(argv[3]) != 1)
	{
		cerr << "Usage: " << argv[0]
			<< " CSV_filename target_pdb target_chain [target_filename|"
			<< "(start_res_id end_res_id)]\n";
		cerr <<
			"\nThis program filters out fragments in the input file that are"
			"\nat the very start or end of the chain (unless the target"
			"\nposition is also the start or end of the chain).\n"
			"\nThe format of CSV_filename is like:\n"
			"\n1f86,A,10,9,cplmvkvld,1gte,A,175,9,...\n"
			"\nBy default, the location of the file with PDB id WXYZ is assumed to be:\n"
			"\n" << PDB_HOME << "XY/pdbWXYZ.ent\n"
			"\nThe values start_res_id and end_res_id are the id of the first"
			"\nand last residue in the target respectively. If target_filename"
			"\nis specified, these values are obtained from the file; if"
			"\nneither is specified, the target file is located in the same"
			"\nway as the the other PDB files.\n\n";
		exit(1);
	}

	C_File file(argv[1], "r", "CSV file");

	string target_pdb = argv[2];
	char target_chain = argv[3][0];
	int target_start, target_end;

	if (argc == 6)
	{
		if (sscanf(argv[4], "%d", &target_start) != 1 ||
			sscanf(argv[5], "%d", &target_end) != 1)
		{
			cerr << argv[0] << ": Illegal start / end residue id: "
				<< argv[4] << " / " << argv[5] << "\n";
			exit(1);
		}
	}
	else
	{
		const char *target_filename = (argc == 5 ? argv[4] : NULL);

		get_pdb_values(target_pdb.c_str(), target_chain,
			&target_start, &target_end, target_filename);
	}

	char buffer[1000];
	char pdb1[100], pdb2[100], chain1, chain2;
	int pos1, len1, pos2, len2;
	int ignore_count = 0;
	int keep_count = 0;

	while (file.next_line(buffer, 1000))
	{
		replace(buffer, ',', ' ');

		// cout << "[" << buffer << "]\n";

		int num = sscanf(buffer, "%s %c %d %d %*s %s %c %d %d",
			pdb1, &chain1, &pos1, &len1,
			pdb2, &chain2, &pos2, &len2);

		// cout << "Num = " << num << ": [" << pdb1 << "] [" << chain1
		//	<< "] [" << pos1 << "] [" << len1 << "], ["
		//	<< pdb2 << "] [" << chain2 << "] [" << pos2
		//	<< "] [" << len2 << "]\n";

		replace(buffer, ' ', ',');

		if (num != 8)
		{
			cerr << "Error on line " << file.line_num()
				<< " of " << file.name() << "\n";
			exit(1);
		}

		char *the_pdb;
		char the_chain;
		int the_pos, the_len;
		bool at_start = false, at_end = false;

		if (target_pdb == pdb1 && target_chain == chain1)
		{
			at_start = (pos1 == target_start);
			at_end = (pos1 + len1 - 1 == target_end);

			the_pdb = pdb2;
			the_chain = chain2;
			the_pos = pos2;
			the_len = len2;
		}
		else
		if (target_pdb == pdb2 && target_chain == chain2)
		{
			at_start = (pos2 == target_start);
			at_end = (pos2 + len2 - 1 == target_end);

			the_pdb = pdb1;
			the_chain = chain1;
			the_pos = pos1;
			the_len = len1;
		}
		else
		{
			cerr << "Error: target pdb \"" << target_pdb
				<< " " << target_chain << "\" not found on line "
				<< file.line_num() << " of " << file.name() << "\n"
				<< buffer << "\n";
			exit(1);
		}

		if (fragment_ok(the_pdb, the_chain, the_pos, the_len,
			at_start, at_end))
		{
			cout << buffer;
			keep_count++;
		}
		else
		{
			// cout << "# " << buffer;
			ignore_count++;
		}
	}

	cerr << "Total " << keep_count << " lines kept, "
		<< ignore_count << " ignored\n";

	dump_pdb_map();
	return 0;
}

