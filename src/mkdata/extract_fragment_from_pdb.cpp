
// extract the PDB "ATOM" records for a single fragment from
// a PDB file
//
// Includes the previous and next reidue (ie. outputs 11 residues).

#include <cstdio>
#include <cstring>
#include <string>
#include <cctype>
#include <iostream>

using namespace std;

const int ChainPos = 21;
const int ResPos = 22;
const int ResLen = 5;

void print_previous_residue(char prev_line[50][1000], int n)
{
	int m = (n - 1 + 50) % 50;

	for ( ; ; )
	{
		/*
		cout << "Prev. comparing (" << n << ") '"
			<< string(prev_line[n] + ResPos, ResLen)
			<< "' & (" << m << ") '"
			<< string(prev_line[m] + ResPos, ResLen)
			<< "'\n";
		*/

		if (strncmp(
				prev_line[n] + ResPos,
				prev_line[m] + ResPos, ResLen) != 0)
		{
			break;
		}

		m = (m - 1 + 50) % 50;

		if (m == n) { break; }
	}

	m = (m + 1) % 50;

	for ( ; ; )
	{
		cout << prev_line[m];

		if (m == n) { break; }

		m = (m + 1) % 50;
	}
}

int main(int argc, char **argv)
{
	if (argc != 4)
	{
		cerr << "Usage: " << argv[0] << " <file> <chain> <start residue>\n";
		exit(1);
	}

	const char *filename = argv[1];
	char chain = argv[2][0];

	// residue can have a letter after it (eg. "237A")
	char start_res[100];

	if (isalpha(argv[3][strlen(argv[3]) - 1]))
	{
		sprintf(start_res, "%5s", argv[3]);
	}
	else
	{
		sprintf(start_res, "%4s ", argv[3]);
	}

	FILE *file;

	file = fopen(filename, "r");

	if (file == NULL)
	{
		cerr << argv[0] << ": cannot open " << filename << "\n";
		exit(1);
	}

	char buffer[1000];
	char prev_line[50][1000];
	bool found = false;
	int res_count = 0;
	int p = 0;
	int prev_p = 49;
	bool any_prev = false;
	bool didnt_run_out = false;

	for (int n = 0;n < 50;n++)
	{
		prev_line[n][0] = '\0';
		// fill up the rest with spaces, so the comparison in
		// print_previous_residue (with the middle of the string) is safe

		for (int m = 1;m < 50;m++) { prev_line[n][m] = ' '; }
	}

	while (fgets(buffer, 1000, file))
	{
		if (strncmp(buffer, "ATOM ", 5) == 0 &&
			(chain == '-' || buffer[ChainPos] == chain))
		{
			if (!found)
			{
				if (strncmp(buffer + ResPos, start_res, ResLen) == 0)
				{
					found = true;
					res_count = 0;

					if (!any_prev)
					{
						// cerr << "Error: no previous residue before " << start_res << "\n";
						// exit(1);
					}
					else
					{
						print_previous_residue(prev_line, prev_p);
					}
				}
				else
				{
// cout << "Adding (" << p << ") " << buffer;
					strcpy(prev_line[p], buffer);
					prev_p = p;
					p = (p + 1) % 50;
					any_prev = true;
					continue;
				}
			}

			if (found)
			{
				/*
				cout << "Comparing '"
					<< string(buffer + ResPos, ResLen)
					<< "' & '"
					<< string(prev_line[prev_p] + ResPos, ResLen)
					<< "'\n";
				*/

				if (strncmp(prev_line[prev_p] + ResPos,
					buffer + ResPos, ResLen) != 0)
				{
					res_count++;
					// cout << "!! res count = " << res_count << "\n";

					if (res_count > 10)
					{
						didnt_run_out = true;
						break;
					}
				}

				cout << buffer;

				strcpy(prev_line[p], buffer);
				prev_p = p;
				p = (p + 1) % 50;
			}
		}
		else
		{
			//if (found) { break; }
		}
	}

	fclose(file);

	if (!found)
	{
		cerr << "Residue " << start_res
			<< "not found in " << filename
			<< " chain " << chain
			<< "\n";
		exit(1);
	}

	if (!didnt_run_out) { res_count++; }

	if (res_count != 11 && res_count != 10)
	{
		cerr << "Error: not enough residues after position " << start_res
			<< "in" << filename
			<< " chain " << chain << "\n";
		exit(1);
	}

	return 0;
}

