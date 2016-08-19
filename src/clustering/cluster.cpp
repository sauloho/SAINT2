
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include "c_file.h"
#include "cluster_template.h"

using namespace std;

typedef vector<string> StringVec;
typedef vector<double> DoubleVec;
typedef vector<vector <double> > DoubleVecVec;

int num_obj = 0;
StringVec obj_name;		// name of each object
DoubleVecVec obj_dist;	// distances between all pairs of objects
DoubleVec temp_vec;		// for get_centroid() (global to avoid reallocation)

const int Max_Objects = 10000;

void usage_and_exit(const char *argv0)
{
	cerr << "Usage: " << argv0
		<< " <filename> <cutoff distance> (-h for help)\n";
	exit(1);
}

void verbose_usage_and_exit(const char *argv0)
{
	cout
		<< "Usage: " << argv0 << " filename\n\n"
		<< "Perform hierarchical clustering with average linkage on a list of\n"
		"objects with known pairwise distances (eg. RMSD between proteins).\n"
		"\n"
		"The \"centroid\" of each cluster is the object with the minimum\n"
		"total squared distance to all of the other objects in the cluster.\n"
		"\n"
		"Input file format:\n"
		"\n"
		"<number of objects>\n"
		"\n"
		"(optional part)\n"
		"\n"
		"0=name of object zero\n"
		"1=name of object one\n"
		"...\n"
		"\n"
		"(mandatory part)\n"
		"\n"
		"0 1 <distance>\n"
		"0 2 <distance>\n"
		"0 3 <distance>\n"
		"...\n"
		"1 2 <distance>\n"
		"1 3 <distance>\n"
		"...\n\n";
	exit(0);
}

void read_data(const string &filename)
{
	C_File file(filename.c_str(), "r", "input file");

	const int Buff_Len = 1000;
	char buffer[Buff_Len];

	if (!file.next_line(buffer, Buff_Len))
	{
		cerr << "Unexpected end of file in " << filename << "\n";
		exit(1);
	}

	if (sscanf(buffer, "%d", &num_obj) != 1)
	{
		cerr << "Error in " << filename
			<< ": expected first line to be the number of objects\n";
		exit(1);
	}

	if (num_obj < 1 || num_obj > Max_Objects)
	{
		cerr << "Error in " << filename
			<< ": number of objects must from 1 to " << Max_Objects << "\n";
		exit(1);
	}

	obj_dist.resize(num_obj);
	obj_name.resize(num_obj);
	temp_vec.resize(num_obj);

	for (int k = 0;k < num_obj;k++)
	{
		obj_dist[k].resize(num_obj);
		obj_dist[k][k] = 0.0;
	}

	// next pair of indexes expected is (i j)
	int i = 0;
	int j = 1;

	while (file.next_line(buffer, Buff_Len))
	{
		char *p = strchr(buffer, '=');

		if (p != NULL)
		{
			int n;
			if (sscanf(buffer, "%d", &n) != 1)
			{
				cerr << "Error on line " << file.line_num()
					<< " of " << filename
					<< ": expected integer index at start of line\n";
				exit(1);
			}

			if (n < 0 || n >= num_obj)
			{
				cerr << "Error on line " << file.line_num()
					<< " of " << filename
					<< ": index must be from 0 to #objects minus one ("
					<< num_obj - 1 << ")\n";
				exit(1);
			}

			// remove trailing newline
			int len = strlen(p);
			p[len - 1] = '\0';

			// find first non-space character after "="
			for (p++;isspace(*p);p++)
			{ }

			obj_name[n] = p;
			// cout << "name[" << n << "] = \"" << p << "\"\n";
		}
		else  // a distance value
		{
			if (i >= num_obj - 1)
			{
				cerr << "Extra data on line " << file.line_num()
					<< " of " << filename << "\n";
				exit(1);
			}

			int n1, n2;
			double d;

			if (sscanf(buffer, "%d %d %lf", &n1, &n2, &d) != 3 ||
				n1 != i || n2 != j)
			{
				cerr << "Error on line " << file.line_num()
					<< " of " << filename
					<< ": expected \"" << i << " " << j << " <distance>\"\n";
				exit(1);
			}

			obj_dist[i][j] = obj_dist[j][i] = d;

			if (++j >= num_obj)
			{
				++i;
				j = i + 1;
			}
		}
	}

	if (i != num_obj - 1)
	{
		cerr << "Not enough data in " << filename
			<< " (last index should be \""
			<< num_obj - 2 << " " << num_obj - 1 << "\")\n";
		exit(1);
	}

	// check for unnamed objects

	for (int z = 0;z < num_obj;z++)
	{
		if (obj_name[z].empty())
		{
			char name[100];

			sprintf(name, "#%d", z);
			obj_name[z] = name;
		}
	}
}

// find the "central representative" of a set of objects
// (the one with the smallest sum of squared distances to all
// the other objects in the set)

int get_centroid(const StringVec &objname, const vector<int> &index)
{
	unsigned i, j;

	// set temp_vec[] to zero for all values in index[]

	for (i = 0;i < index.size();i++)
	{
		temp_vec[i] = 0.0;
	}

	// find the total squared distance between all values in index[]

	for (i = 1;i < index.size();i++)
	{
		int n1 = index[i];

		for (j = 0;j < i;j++)
		{
			int n2 = index[j];
			double d = obj_dist[n1][n2];
			double d_squared = d * d;

//cout << ":" << n1 << " & " << n2 << " += " << d << " squared\n";
			temp_vec[i] += d_squared;
			temp_vec[j] += d_squared;
		}
	}

	// find the index with the lowest total

	int best = 0;

//cout << "CENTROID of:\n";
//cout << index[0] << " (" << temp_vec[0] << ")\n";

	for (i = 1;i < index.size();i++)
	{
//cout << index[i] << " (" << temp_vec[i] << ")\n";
		if (temp_vec[i] < temp_vec[best])
		{
			best = i;
		}
	}
//cout << "=> " << index[best] << "\n";

	return index[best];
}

void print_central_cluster(vector<vector<int> > &result, vector<int> &centroid)
{
	vector<double> t;
	t.resize(centroid.size());

	int i;
	for (i = 0;i < (int) t.size();i++)
	{
		t[i] = 0.0;
	}

	for (int n1 = 1;n1 < (int) result.size();n1++)
	{
		int c1 = centroid[n1];

		for (int n2 = 0;n2 < n1;n2++)
		{
			int c2 = centroid[n2];

			double d = obj_dist[c1][c2];
			double d_squared = d * d;

			t[n1] += d_squared;
			t[n2] += d_squared;
		}
	}

	int best = 0;

	for (i = 1;i < (int) t.size();i++)
	{
		if (t[i] < t[best])
		{
			best = i;
		}
	}

	cout << "central cluster " << best << ' '
		<< obj_name[centroid[best]] << "\n\n";
}

void do_clustering(double cutoff)
{
	vector<vector<int> > result;
	vector<int> centroid;
	Cluster<string> clust;
	clust.cluster(obj_name, obj_dist, get_centroid, cutoff, result, centroid);

	cout << result.size() << " clusters, largest = " << result[0].size()
		<< " (cut off = " << cutoff << ")\n\n";

	print_central_cluster(result, centroid);

	unsigned n;
	for (n = 0;n < result.size();n++)
	{
		int c = centroid[n];

		cout << "cluster " << n
			<< " size " << result[n].size()
			<< " centre " << obj_name[c] << '\n';
	}

	cout << "\n";

	for (n = 0;n < result.size();n++)
	{
		for (unsigned m = 0;m < result[n].size();m++)
		{
			cout << "c" << n << ' ' << obj_name[result[n][m]] << '\n';
		}
	}
}

int main(int argc, char **argv)
{
	if (argc == 2)
	{
		if (argv[1][0] == '-' &&
			(argv[1][1] == 'h' || argv[1][1] == '?'))
		{
			verbose_usage_and_exit(argv[0]);
		}

		usage_and_exit(argv[0]);
	}
	else
	if (argc != 3)
	{
		usage_and_exit(argv[0]);
	}

	string filename = argv[1];
	double cutoff;

	if (sscanf(argv[2], "%lf", &cutoff) != 1)
	{
		cerr << "Illegal cutoff value \"" << argv[2] << "\"\n";
		exit(1);
	}

	read_data(filename);
	do_clustering(cutoff);
	return 0;
}

