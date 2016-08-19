
#ifndef CLUSTER_H_INCLUDED
#define CLUSTER_H_INCLUDED

// Author: David Fisher   21/10/2010
//
// A template class that performs hierarchical clustering using average linkage.
//
// This version of the class represents the centroid of a cluster by the most
// "central" object in the cluster (rather than creating a new centroid object).
//
// The following function needs to be defined to use this class
// (it can have any name):
//
// int getCentroid(const vector<T> &objects, const vector<int> &index)
//
// - this returns one of the integers in vector "index", which
//   are indexes in the "objects" vector. The index of the most central
//   object should be returned.

/////////////////////////////////////////////////////////////////////
// Example using protein fragments (class Fragment) for type T:

// int getCentroid(const vector<Fragment> &frag, const vector<int> &index)
// {
//    ...
//    return index[...];
// }
//
// ...
//    vector<Fragment> fragments;		// the fragments to cluster
//    vector<vector<double> > distance;	// RMSDs between fragments
// ...
//    Cluster<Fragment> clust;
//    vector<vector<int> > result;
//    const double threshold = 2.5;
//
//    clust.cluster(fragments, distance, getCentroid, threshold, result);
//
//    cout << result.size() << " clusters found\n";

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <cassert>

using namespace std;

//#define DEBUG
//#define CAN_PRINT_OBJECT

template <class T>
class Cluster
{
public:
	typedef int (*GetCentroidFn)(const vector<T> &obj,
		const vector<int> &index);

	Cluster()
	{ }

	~Cluster()
	{ }

	// Perform hierarchical clustering, given a vector of objects,
	// a function to see if an object is closer than a certain distance
	// from another one, and a distance cut off value (threshold).
	//
	// A list of object indexes for each cluster is put in a separate
	// vector in "result"; the indexes of the centroid objects are put
	// in result_centres
	void cluster(const vector<T> &object,
		const vector<vector<double> > &distance,
		GetCentroidFn getCentroidFn, double threshold,
		vector<vector<int> > &result, vector<int> &result_centroids);

private:
	struct Node
	{
		int m_index;		// index in vector<object>, or -1 for created nodes
		int m_val;			// index in vector<object> (centroid or leaf node)
		Node *m_closest;
		double m_closest_dist;
		Node *m_child_1;
		Node *m_child_2;
		int m_num;			// number of objects represented by this node

		Node(int index, int val) : m_index(index), m_val(val),
			m_closest(NULL), m_closest_dist(9e99),
			m_child_1(NULL), m_child_2(NULL), m_num(1)
		{ }


		// comparison operator for sorting
		bool operator < (const Node &other) const
		{
			// order by distance to closest node
			return (m_closest_dist < other.m_closest_dist);
		}
	};

	struct NodePtr
	{
		Node *node_ptr;

		NodePtr()
		{ }

		NodePtr(Node *ptr) : node_ptr(ptr)
		{ }

		// comparison operator for sorting
		bool operator < (const NodePtr &other) const
		{
			// reverse order by number of children
			return (node_ptr->m_num > other.node_ptr->m_num);
		}
	};

private:
	typedef list<Node> NodeList;
	NodeList m_node;
	NodeList m_child_node;

private:
	void find_child_indexes(Node *node, vector<int> &child_vec);
	void dump_node_list(const NodeList &nodes);
};

// since this is a template class, the implementation needs to be
// in the header file

template <class T>
void Cluster<T>::cluster(const vector<T> &object,
	const vector<vector<double> > &distance,
	GetCentroidFn getCentroidFn, double threshold,
	vector<vector<int> > &result, vector<int> &result_centroids)
{
	m_node.erase(m_node.begin(), m_node.end());
	m_child_node.erase(m_child_node.begin(), m_child_node.end());

	// create m_node and find the closest node to each node

	int n;
	for (n = 0;n < (int) object.size();n++)
	{
		m_node.push_back(Node(n, n));

		typename NodeList::iterator i;
		typename NodeList::iterator curr = m_node.end();
		--curr;

		curr->m_closest_dist = 9e99;

		for (i = m_node.begin();i != curr;++i)
		{
			double d = distance[i->m_val][curr->m_val];

			if (d < curr->m_closest_dist)
			{
				curr->m_closest_dist = d;
				curr->m_closest = &(*i);
			}

			if (d < i->m_closest_dist)
			{
				i->m_closest_dist = d;
				i->m_closest = &(*curr);
			}
		}
	}

	m_node.sort();

#ifdef DEBUG
	dump_node_list(m_node);
#endif // DEBUG

	while (m_node.size() > 1 &&
		   m_node.front().m_closest_dist <= threshold)
	{
		// merge the front node with its closest node, put the result
		// in sorted order in m_node (unless its closest node is beyond
		// the threshold cut off value, in which case don't bother sorting,
		// just add it to the end), and move the two nodes into the
		// m_child_node list

		Node *node = &m_node.front();
		int new_val;

		vector<int> obj_indexes;

		find_child_indexes(node, obj_indexes);
		find_child_indexes(node->m_closest, obj_indexes);

		new_val = getCentroidFn(object, obj_indexes);
		assert(new_val >= 0 && new_val < (int) object.size());

		// move the closest node from m_node to m_child_node

		typename NodeList::iterator i = m_node.begin();

		for (++i;i != m_node.end() && &(*i) != node->m_closest;++i)
		{ }

		if (i == m_node.end())
		{
			cerr << "Error #1 during clustering\n";
			exit(1);
		}

		m_child_node.splice(m_child_node.end(), m_node, i);
		assert(&m_child_node.back() == node->m_closest);
		Node *child_2 = node->m_closest;

		m_child_node.splice(m_child_node.end(), m_node, m_node.begin());
		assert(&m_child_node.back() == node);
		Node *child_1 = node;

		// all nodes with child_1 or child_2 as their closest node
		// need to recalculate their closest node

		for (i = m_node.begin();i != m_node.end();++i)
		{
			if (i->m_closest == child_1 || i->m_closest == child_2)
			{
				i->m_closest = NULL;
				i->m_closest_dist = 9e99;

				typename NodeList::iterator j;

				for (j = m_node.begin();j != m_node.end();++j)
				{
					if (j != i)
					{
						double d = distance[i->m_val][j->m_val];
						
						if (d < i->m_closest_dist)
						{
							i->m_closest_dist = d;
							i->m_closest = &(*j);
						}
					}
				}

				if (m_node.size() != 1)
				{
					assert(i->m_closest != NULL);
				}
			}
		}

		//cout << "Node list size = " << m_node.size() << "\n";
		//cout << "Child list size = " << m_child_node.size() << "\n";

		// temporarily add the new (merged) node to the front of the
		// m_child_node list, and find the distance to all the nodes
		// in the m_node list.
		// 
		// (The reason for putting it in the m_child_node list is so
		// that it can be safely spliced into the m_node list later on).

		m_child_node.push_front(Node(-1, new_val));
		Node *new_node = &m_child_node.front();

		new_node->m_child_1 = child_1;
		new_node->m_child_2 = child_2;
		new_node->m_num = child_1->m_num + child_2->m_num;

		new_node->m_closest_dist = 9e99;

		for (i = m_node.begin();i != m_node.end();++i)
		{
			double d = distance[i->m_val][new_node->m_val];

			if (d < new_node->m_closest_dist)
			{
				new_node->m_closest_dist = d;
				new_node->m_closest = &(*i);
			}

			if (d < i->m_closest_dist)
			{
				i->m_closest_dist = d;
				i->m_closest = new_node;
			}
		}

		// set i to the position to insert the new node before

		if (new_node->m_closest_dist > threshold)
		{
			// no need to keep in sorted order; just append the node
			// to the end of the list
			i = m_node.end();
		}
		else
		{
			// find the correct (sorted) position for the new node

			for (i = m_node.begin();i != m_node.end();++i)
			{
				if (new_node->m_closest_dist < i->m_closest_dist)
				{
					break;
				}
			}
		}

#ifdef DEBUG
		cout << "Merging nodes "
			<< child_1 << " (#" << child_1->m_index << ") "
#ifdef CAN_PRINT_OBJECT
			<< child_1->m_val
#endif // CAN_PRINT_OBJECT
			<< " & "
			<< child_2 << " (#" << child_2->m_index << ") "
#ifdef CAN_PRINT_OBJECT
			<< child_2->m_val
#endif // CAN_PRINT_OBJECT
			<< " => "
#ifdef CAN_PRINT_OBJECT
			<< new_node->m_val
#endif // CAN_PRINT_OBJECT
			<< " dist "
			<< new_node->m_closest_dist
			<< " from "
			<< new_node->m_closest;

		if (new_node->m_closest != NULL)
		{
			cout << " (#"
				<< new_node->m_closest->m_index
				<< ")";
		}

		cout << "\n";

		cout << "\nAt this point, node list is:\n\n";
		dump_node_list(m_node);
		cout << "\nChild list is:\n\n";
		dump_node_list(m_child_node);
		cout << "\nNew node is " << new_node
			<< " (#" << new_node->m_index << ") dist "
			<< new_node->m_closest_dist
			<< "\n\n";
		cout << "INSERTING at i = " << &(*i) << " (#" << i->m_index << ")\n";
#endif // DEBUG

		// (splice() inserts before position i)
		m_node.splice(i, m_child_node, m_child_node.begin());

#ifdef DEBUG
		cout << "\n- After splicing:\n"
		 	 << "\n- Node list is:\n\n";
		dump_node_list(m_node);
		cout << "\n- Child list is:\n\n";
		dump_node_list(m_child_node);
#endif // DEBUG
	}

#ifdef DEBUG
	cout << "\n" << m_node.size() << " nodes left:\n\n";
	dump_node_list(m_node);

	cout << "\nChild list is:\n\n";
	dump_node_list(m_child_node);
#endif // DEBUG

	// create a vector of node pointers, then sort by number of children

	vector <NodePtr> node_vec;
	node_vec.resize(m_node.size());

	typename NodeList::iterator i;
	for (i = m_node.begin(), n = 0;i != m_node.end();++i, ++n)
	{
		node_vec[n].node_ptr = &(*i);
	}

	// sort by number of children (in reverse, so largest cluster comes first)
	sort(node_vec.begin(), node_vec.end());

	result.erase(result.begin(), result.end());
	result.resize(node_vec.size());
	result_centroids.erase(result_centroids.begin(), result_centroids.end());
	result_centroids.resize(node_vec.size());

	for (n = 0;n < (int) node_vec.size();n++)
	{
		result_centroids[n] = node_vec[n].node_ptr->m_val;
		find_child_indexes(node_vec[n].node_ptr, result[n]);
	}
}

template <class T>
void Cluster<T>::find_child_indexes(Node *node, vector<int> &child_vec)
{
	if (node->m_index != -1)
	{
		child_vec.push_back(node->m_index);
	}
	else
	{
		assert(node->m_child_1 != NULL);
		assert(node->m_child_2 != NULL);

		find_child_indexes(node->m_child_1, child_vec);
		find_child_indexes(node->m_child_2, child_vec);
	}
}

#endif // CLUSTER_H_INCLUDED
