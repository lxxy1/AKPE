#include <stdlib.h>
#include <cstdio>
#include <stdbool.h>
#include <string.h>
#include <string>
#include <vector>
#include <set>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <cassert>
#include <queue>
#include <cstdarg>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include "iterator"
#include <immintrin.h>

#define BIG 7572880
#define unfilled -1
using namespace std;
using namespace chrono;

typedef int int32;
typedef char char8;
typedef bool boolean;
typedef pair<int32, int32> upair;
typedef unsigned int uint32;

int K;
int T;
int EDGENUM;
int NODENUM;
int lnum = 1000000;
long long cost = 0;
int plexnum = 0;
vector<vector<int>> Lplex;
vector<vector<int>> Rplex;

class CuckooHash
{
private:
	/* data */
	int32 capacity;
	int32 mask;
	int32 size;
	int32 buff_size = sizeof(int32);
	int32 *hashtable;

	void rehash(int32 **_table)
	{
		int32 oldcapacity = capacity;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		capacity = (mask + 1) * buff_size;
		int32 *newhash = new int32[capacity];
		memset((newhash), unfilled, sizeof(int32) * capacity);
		for (int32 i = 0; i < oldcapacity; ++i)
		{
			if ((*_table)[i] != unfilled)
				insert((*_table)[i], &newhash);
		}
		swap((*_table), newhash);
		delete[] newhash;
	}
	void insert(const int32 &_u, int32 **_table)
	{

		int32 hs = hash1(_u);
		for (int32 i = 0; i < buff_size; ++i)
		{
			if ((*_table)[hs * buff_size + i] == unfilled)
			{
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}
		hs = hash2(_u);
		for (int32 i = 0; i < buff_size; ++i)
		{
			if ((*_table)[hs * buff_size + i] == unfilled)
			{
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}

		bool use_hash1 = true;
		int32 u = _u;
		for (int32 i = 0; i < mask; ++i)
		{
			int32 replaced;
			if (use_hash1)
				hs = hash1(u);
			else
				hs = hash2(u);
			int32 j = 0;
			for (; j < buff_size; ++j)
			{
				if ((*_table)[hs * buff_size + j] == unfilled)
					break;
			}
			if (buff_size == j)
			{
				replaced = move((*_table)[hs * buff_size]);
				j = 1;
				for (; j < buff_size; j++)
				{
					(*_table)[hs * buff_size + j - 1] =
						move((*_table)[hs * buff_size + j]);
				}
				(*_table)[hs * buff_size + j - 1] = u;
			}
			else
			{
				replaced = move((*_table)[hs * buff_size + j]);
				(*_table)[hs * buff_size + j] = u;
			}
			use_hash1 = hs == hash2(replaced);
			u = move(replaced);
			if (u == unfilled)
				return;
		}
		rehash(_table);
		insert(u, _table);
	}

	int32 hash1(const int32 &x) { return x & mask; }
	int32 hash2(const int32 &x) { return ~x & mask; }

public:
	CuckooHash(/* args */)
	{
		capacity = 0;
		hashtable = NULL;
		mask = 0;
		size = 0;
	}
	~CuckooHash()
	{
		if (hashtable)
			delete[] hashtable;
	}

	void reserve(int32 _size)
	{
		if (capacity >= _size)
			return;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		while (_size >= mask * buff_size)
			mask = (mask << 1) | 1;
		capacity = (mask + 1) * buff_size;
		if (hashtable)
			delete[] hashtable;
		hashtable = new int32[capacity];
		memset(hashtable, unfilled, sizeof(int32) * capacity);
	}

	void insert(const int32 &_u)
	{
		if (find(_u))
			return;
		insert(_u, &hashtable);
		size++;
	}

	bool find(const int32 &_u)
	{
		int32 hs1 = hash1(_u);
		int32 hs2 = hash2(_u);
		//cout << "buff_size=" << buff_size << endl;
		if(buff_size != 4 || sizeof(int32) != 4)
		{
			cout << "buff_size=" << buff_size << endl;
			cout << "sizeof(int32)=" << sizeof(int32) << endl;
		}

		assert(buff_size == 4 && sizeof(int32) == 4);
		__m128i cmp = _mm_set1_epi32(_u);
		__m128i b1 = _mm_load_si128((__m128i *)&hashtable[buff_size * hs1]);
		__m128i b2 = _mm_load_si128((__m128i *)&hashtable[buff_size * hs2]);
		__m128i flag = _mm_or_si128(_mm_cmpeq_epi32(cmp, b1), _mm_cmpeq_epi32(cmp, b2));

		return _mm_movemask_epi8(flag) != 0;
	}
	int32 getcapacity() { return capacity; }
	int32 getmask() { return mask; }
	int32 *gethashtable() { return hashtable; }
};

vector<CuckooHash> Pcuhash;
vector<CuckooHash> Ncuhash;
high_resolution_clock::time_point a_s = high_resolution_clock::now();
high_resolution_clock::time_point a_e = high_resolution_clock::now();
auto Dichromatictime = std::chrono::duration_cast<std::chrono::microseconds>(a_e - a_s).count();
auto Dichromatictwotime = std::chrono::duration_cast<std::chrono::microseconds>(a_e - a_s).count();
auto colortime = std::chrono::duration_cast<std::chrono::microseconds>(a_e - a_s).count();
auto colorchecktime = std::chrono::duration_cast<std::chrono::microseconds>(a_e - a_s).count();
auto colordefftime = std::chrono::duration_cast<std::chrono::microseconds>(a_e - a_s).count();
auto vneiPtime = std::chrono::duration_cast<std::chrono::microseconds>(a_e - a_s).count();
auto vneitime = std::chrono::duration_cast<std::chrono::microseconds>(a_e - a_s).count();

class edge
{
public:
	int n1;
	int n2;
	int sign;
};

class node
{
public:
	int rid;
	int Pid;
	int Nid;
	int Pdegree;
	int Ndegree;
	bool isExist;
};
class nodesetid
{
public:
	int id;
	int degree;
};
bool flag = false;
class Graph
{
public:
	edge *elist;  
	node *nlist;  
	int *Padjncy; 
	int *Nadjncy; 
	vector<nodesetid> nodeset;

	int maxpd = 0;
	int maxnd = 0;

	Graph(/* args */);
	virtual ~Graph();
	void mkGraph(string nverts);
	void showGraph();
	bool testadj(int a, int b, int c);
	bool isExistinG(int a);
	int vectorneighbor(int a, vector<int> &r);
	int vectorPneighbor(int a, vector<int> &r);
	int vectorNneighbor(int a, vector<int> &r);
	void Dichromatic(int nv, set<int> &Ln, set<int> &Rn);
	bool test_difference(vector<int> &a, int b, int c, bool f);
	void pivotselect(vector<int> &p, vector<int> &pn, vector<int> &c, int u, int f, vector<int> &r);

	void removefromneiadj(int v, int u, bool f);
	void VertexReduction(int k, int t);
	int ColorSort(vector<int> p, int &colorsetnum, vector<int> &colornode, bool f);
	int ColorCheck(int a, int colorsetnum, vector<int> &colornode, bool f);

	int gvectors_intersection(vector<int> &v1, int a, vector<int> &v);
	int gvectors_intersection(vector<int> &v1, int a, bool f, vector<int> &v);

    int gvectors_difference(vector<int> &v1, int a, vector<int> &v);
    int gvectors_difference(vector<int> &v1, int a, bool f, vector<int> &v);
};

Graph::Graph(/* args */)
{
}

Graph::~Graph()
{
	if (nlist != NULL)
	{
		delete nlist;
	}
	if (Padjncy != NULL)
	{
		delete Padjncy;
	}
	if (Nadjncy != NULL)
	{
		delete Nadjncy;
	}
}

bool comp(const edge &a, const edge &b)
{
	if (a.n1 < b.n1)
	{
		return true;
	}
	else if (a.n1 == b.n1 && a.n2 < b.n2)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool comp1(const nodesetid &a, const nodesetid &b)
{
	if (a.degree < b.degree)
	{
		return true;
	}
	else
	{
		return false;
	}
}

class Hasher
{
public:
	size_t operator()(std::vector<int> const &hnode) const
	{
		std::size_t seed = hnode.size();
		for (auto &i : hnode)
		{
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

class Equal
{
public:
	bool operator()(std::vector<int> const &hnode1, std::vector<int> const &hnode2) const
	{
		if (hnode1 == hnode2)
		{
			return true;
		}
		return false;
	}
};

int vectors_difference(vector<int> &v1, vector<int> &v2, vector<int> &v)
{
	v.clear();
	int nsize = 0;
	vector<int>::iterator i = v1.begin();
	vector<int>::iterator j = v2.begin();
	while (i != v1.end() && j != v2.end())
	{
		if ((*i) < (*j))
		{
			v.push_back((*i));
			i++;
			nsize++;
		}
		else if ((*i) > (*j))
		{
			j++;
		}
		else
		{
			i++;
			j++;
		}
	}
	for (; i != v1.end(); i++)
	{
		v.push_back((*i));
		nsize++;
	}
	return nsize;
}


int vectors_intersection(vector<int> &v1, vector<int> &v2, vector<int> &v)
{
	v.clear();
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v)); 
	return v.size();
}

int Graph::gvectors_difference(vector<int> &v1, int a, vector<int> &v)
{
    high_resolution_clock::time_point Colord_s = high_resolution_clock::now();
    v.clear();
    for (int i = 0; i < v1.size(); i++)
    {
        if (!testadj(v1[i], a, 0) && !testadj(v1[i], a, 1))
        {
            v.push_back(v1[i]);
        }
    }
    high_resolution_clock::time_point Colord_e = high_resolution_clock::now();
	colordefftime += std::chrono::duration_cast<std::chrono::microseconds>(Colord_e - Colord_s).count();
    return v.size();
}

int Graph::gvectors_difference(vector<int> &v1, int a, bool f, vector<int> &v)
{
    high_resolution_clock::time_point Colord_s = high_resolution_clock::now();
    v.clear();
    if(f)
    {
        for (int i = 0; i < v1.size(); i++)
        {
            if(!testadj(v1[i],a,0))
            {
                v.push_back(v1[i]);
            }
        }
    }
    else
    {
        for (int i = 0; i < v1.size(); i++)
        {
            if(!testadj(v1[i],a,1))
            {
                v.push_back(v1[i]);
            }
        }
    }
    high_resolution_clock::time_point Colord_e = high_resolution_clock::now();
	colordefftime += std::chrono::duration_cast<std::chrono::microseconds>(Colord_e - Colord_s).count();
    return v.size();
}

int Graph::gvectors_intersection(vector<int> &v1, int a, vector<int> &v)
{ 
	v.clear();
	set<int> vset;
	if (1)
	{
		int i = 0;
		vector<int>::iterator j = v1.begin();
		while (i < nlist[a].Pdegree && j != v1.end())
		{
			if (Padjncy[nlist[a].Pid + i] < (*j))
				i++;
			else if (Padjncy[nlist[a].Pid + i] > (*j))
				j++;
			else
			{
				vset.insert(Padjncy[nlist[a].Pid + i]);
				i++;
			}
		}
	}
	if (1)
	{
		int i = 0;
		vector<int>::iterator j = v1.begin();
		while (i < nlist[a].Ndegree && j != v1.end())
		{
			if (Nadjncy[nlist[a].Nid + i] < (*j))
				i++;
			else if (Nadjncy[nlist[a].Nid + i] > (*j))
				j++;
			else
			{
				vset.insert(Nadjncy[nlist[a].Nid + i]);
				i++;
			}
		}
	}
	v.assign(vset.begin(), vset.end());
	sort(v.begin(), v.end());

	return v.size();
}

int Graph::gvectors_intersection(vector<int> &v1, int a, bool f, vector<int> &v)
{ 
	v.clear();
	if (f)
	{
		int i = 0;
		vector<int>::iterator j = v1.begin();
		while (i < nlist[a].Pdegree && j != v1.end())
		{
			if (Padjncy[nlist[a].Pid + i] < (*j))
				i++;
			else if (Padjncy[nlist[a].Pid + i] > (*j))
				j++;
			else
			{
				v.push_back(Padjncy[nlist[a].Pid + i]);
				i++;
			}
		}
	}
	else
	{
		int i = 0;
		vector<int>::iterator j = v1.begin();
		while (i < nlist[a].Ndegree && j != v1.end())
		{
			if (Nadjncy[nlist[a].Nid + i] < (*j))
				i++;
			else if (Nadjncy[nlist[a].Nid + i] > (*j))
				j++;
			else
			{
				v.push_back(Nadjncy[nlist[a].Nid + i]);
				i++;
			}
		}
	}

	return v.size();
}

int vectors_union(vector<int> &v1, vector<int> &v2, vector<int> &v)
{
	v.clear();
	v.reserve(v1.size() + v2.size());
	set<int> nv;
	int nsize = 0;
	for (auto &i : v1)
	{
		nv.insert(i);
	}
	for (auto &j : v2)
	{
		nv.insert(j);
	}
	v.assign(nv.begin(), nv.end());

	sort(v.begin(), v.end());
	return v.size();
}

bool Graph::isExistinG(int a)
{
	return nlist[a].isExist;
}

void Graph::mkGraph(string nverts)
{
	ifstream filenverts;
	filenverts.open(nverts);
	int a1 = 0, a2 = 0, s = 0, en = 0;
	filenverts >> NODENUM;
	filenverts >> EDGENUM;

	elist = new edge[EDGENUM + 1];

	nlist = new node[NODENUM + 1];

    int positiveedge = 0;

	while (!filenverts.eof())
	{
		filenverts >> a1;
		filenverts >> a2;
		filenverts >> elist[en].sign;

        if(elist[en].sign==1)
        {
            positiveedge++;
        }
		if (a1 > a2)
		{
			elist[en].n1 = a2;
			elist[en].n2 = a1;
		}
		else if (a1 < a2)
		{
			elist[en].n1 = a1;
			elist[en].n2 = a2;
		}
		en++;
	}
    cout << "positiveedge=" << positiveedge << "negativeedge=" << EDGENUM - positiveedge << endl;
	sort(elist, elist + EDGENUM, comp);
	for (int i = 0; i < NODENUM + 1; i++)
	{
		nlist[i].rid = i;
		nlist[i].Pid = 0;
		nlist[i].Nid = 0;
		nlist[i].Pdegree = 0;
		nlist[i].Ndegree = 0;
		nlist[i].isExist = true;
	}
	int nown1 = 0, nown2 = 0;
	for (int i = 0; i < EDGENUM; i++)
	{
		if (nown1 != elist[i].n1 || nown2 != elist[i].n2)
		{
			if (elist[i].sign == 1)
			{
				nlist[elist[i].n1].Pdegree++;
				nlist[elist[i].n2].Pdegree++;
				maxpd = max(nlist[elist[i].n1].Pdegree, maxpd);
				maxpd = max(nlist[elist[i].n2].Pdegree, maxpd);
			}
			else
			{
				nlist[elist[i].n1].Ndegree++;
				nlist[elist[i].n2].Ndegree++;
				maxnd = max(nlist[elist[i].n1].Ndegree, maxnd);
				maxnd = max(nlist[elist[i].n2].Ndegree, maxnd);
			}
		}
		nown1 = elist[i].n1;
		nown2 = elist[i].n2;
	}
	Pcuhash.resize(NODENUM);
	Ncuhash.resize(NODENUM);
	for (int i = 1; i < NODENUM + 1; i++)
	{
		Pcuhash[i - 1].reserve(nlist[i - 1].Pdegree);
		Ncuhash[i - 1].reserve(nlist[i - 1].Ndegree);
		nlist[i].Pid = nlist[i - 1].Pid + nlist[i - 1].Pdegree;
		nlist[i].Nid = nlist[i - 1].Nid + nlist[i - 1].Ndegree;
	}
	Padjncy = new int[nlist[NODENUM].Pid + 1];
	Nadjncy = new int[nlist[NODENUM].Nid + 1];
	for (int i = 0; i < nlist[NODENUM].Pid + 1; i++)
	{
		Padjncy[i] = -1;
	}
	for (int i = 0; i < nlist[NODENUM].Nid + 1; i++)
	{
		Nadjncy[i] = -1;
	}
	nown1 = 0, nown2 = 0;
	for (int i = 0; i < EDGENUM; i++)
	{
		if (nown1 != elist[i].n1 || nown2 != elist[i].n2)
		{
			if (elist[i].sign == 1)
			{
				Pcuhash[elist[i].n1].insert(elist[i].n2);
				Pcuhash[elist[i].n2].insert(elist[i].n1);
				for (int j = nlist[elist[i].n1].Pid; j < nlist[elist[i].n1 + 1].Pid; j++)
				{
					if (Padjncy[j] == -1)
					{
						Padjncy[j] = elist[i].n2;
						break;
					}
				}
				for (int j = nlist[elist[i].n2].Pid; j < nlist[elist[i].n2 + 1].Pid; j++)
				{
					if (Padjncy[j] == -1)
					{
						Padjncy[j] = elist[i].n1;
						break;
					}
				}
			}
			else
			{
				Ncuhash[elist[i].n1].insert(elist[i].n2);
				Ncuhash[elist[i].n2].insert(elist[i].n1);
				for (int j = nlist[elist[i].n1].Nid; j < nlist[elist[i].n1 + 1].Nid; j++)
				{
					if (Nadjncy[j] == -1)
					{
						Nadjncy[j] = elist[i].n2;
						break;
					}
				}
				for (int j = nlist[elist[i].n2].Nid; j < nlist[elist[i].n2 + 1].Nid; j++)
				{
					if (Nadjncy[j] == -1)
					{
						Nadjncy[j] = elist[i].n1;
						break;
					}
				}
			}
		}
		nown1 = elist[i].n1;
		nown2 = elist[i].n2;
	}
}

void Graph::showGraph()
{
	for (int i = 0; i < NODENUM; i++)
	{
		cout << i << " Pneighbor" << endl;
		for (int j = 0; j < nlist[i].Pdegree; j++)
		{
			cout << Padjncy[nlist[i].Pid + j] << endl;
		}
	}
	for (int i = 0; i < NODENUM; i++)
	{
		cout << i << " Nneighbor" << endl;
		for (int j = 0; j < nlist[i].Ndegree; j++)
		{
			cout << Nadjncy[nlist[i].Nid + j] << endl;
		}
	}
}

bool Graph::testadj(int a, int b, int c) 
{
	int d;
	if (c == 0)
	{
		return Pcuhash[nodeset[a].id].find(nodeset[b].id);
	}
	else
	{
		return Ncuhash[nodeset[a].id].find(nodeset[b].id);
	}
}

int Graph::vectorneighbor(int a, vector<int> &r)
{
	high_resolution_clock::time_point vnei_s = high_resolution_clock::now();
	r.clear();
	for (int i = 0; i < nlist[a].Pdegree; i++)
	{
		r.push_back(Padjncy[nlist[a].Pid + i]);
	}
	for (int i = 0; i < nlist[a].Ndegree; i++)
	{
		r.push_back(Nadjncy[nlist[a].Nid + i]);
	}
	sort(r.begin(), r.end());
	high_resolution_clock::time_point vnei_e = high_resolution_clock::now();
	vneitime += std::chrono::duration_cast<std::chrono::microseconds>(vnei_e - vnei_s).count();
	return r.size();
}

int Graph::vectorPneighbor(int a, vector<int> &r)
{
	high_resolution_clock::time_point vneiP_s = high_resolution_clock::now();
	r.clear();
	for (int i = 0; i < nlist[a].Pdegree; i++)
	{
		r.push_back(Padjncy[nlist[a].Pid + i]);
	}
	high_resolution_clock::time_point vneiP_e = high_resolution_clock::now();
	vneiPtime += std::chrono::duration_cast<std::chrono::microseconds>(vneiP_e - vneiP_s).count();
	return r.size();
}

int Graph::vectorNneighbor(int a, vector<int> &r)
{
	r.clear();
	for (int i = 0; i < nlist[a].Ndegree; i++)
	{
		r.push_back(Nadjncy[nlist[a].Nid + i]);
	}
	return r.size();
}

void Graph::Dichromatic(int nv, set<int> &Ln, set<int> &Rn)
{
	vector<int> L;
	vector<vector<int>> LPN;
	LPN.resize(nlist[nv].Pdegree);
	vector<vector<int>> LNN;
	LNN.resize(nlist[nv].Pdegree);
	vector<int> R;
	vector<vector<int>> RPN;
	vector<vector<int>> RNN;
	RPN.resize(nlist[nv].Ndegree);
	RNN.resize(nlist[nv].Ndegree);
	for (int i = 0; i < nlist[nv].Pdegree; i++)
	{
		L.push_back(Padjncy[nlist[nv].Pid + i]);
	}
	for (int i = 0; i < nlist[nv].Ndegree; i++)
	{
		R.push_back(Nadjncy[nlist[nv].Nid + i]);
	}
	for (int i = 0; i < L.size(); i++)
	{
		int v = L[i];
		int k = 0, j = 0;
		while (j < nlist[v].Pdegree && k < L.size())
		{
			if (Padjncy[nlist[v].Pid + j] < L[k])
				j++;
			else if (Padjncy[nlist[v].Pid + j] > L[k])
			{
				k++;
			}
			else
			{
				LPN[i].push_back(k);
				j++;
			}
		}
		k = 0, j = 0;
		while (j < nlist[v].Ndegree && k < R.size())
		{
			if (Nadjncy[nlist[v].Nid + j] < R[k])
				j++;
			else if (Nadjncy[nlist[v].Nid + j] > R[k])
			{
				k++;
			}
			else
			{
				LNN[i].push_back(k);
				j++;
			}
		}
	}
	for (int i = 0; i < R.size(); i++)
	{
		int v = R[i];
		int k = 0, j = 0;
		while (j < nlist[v].Pdegree && k < R.size())
		{
			if (Padjncy[nlist[v].Pid + j] < R[k])
				j++;
			else if (Padjncy[nlist[v].Pid + j] > R[k])
			{
				k++;
			}
			else
			{
				RPN[i].push_back(k);
				j++;
			}
		}
		k = 0, j = 0;
		while (j < nlist[v].Ndegree && k < L.size())
		{
			if (Nadjncy[nlist[v].Nid + j] < L[k])
				j++;
			else if (Nadjncy[nlist[v].Nid + j] > L[k])
			{
				k++;
			}
			else
			{
				RNN[i].push_back(k);
				j++;
			}
		}
	}
	queue<int> qL;
	vector<bool> Lflag;
	Lflag.resize(L.size());
	vector<int> LPdegree;
	LPdegree.resize(L.size());
	vector<int> LNdegree;
	LNdegree.resize(L.size());
	for (int i = 0; i < L.size(); i++)
	{
		LPdegree[i] = LPN[i].size();
		LNdegree[i] = LNN[i].size();
		Lflag[i] = true;
		int a, b;
		a = LPN[i].size();
		b = LNN[i].size();
		if (a < T - 2 * K || a + b < 2 * T - 2 * K)
		{
			qL.push(i);
			Lflag[i] = false;
		}
	}

	queue<int> qR;
	vector<bool> Rflag;
	Rflag.resize(R.size());
	vector<int> RPdegree;
	RPdegree.resize(R.size());
	vector<int> RNdegree;
	RNdegree.resize(R.size());
	for (int i = 0; i < R.size(); i++)
	{
		RPdegree[i] = RPN[i].size();
		RNdegree[i] = RNN[i].size();
		Rflag[i] = true;
		int a, b;
		a = RPN[i].size();
		b = RNN[i].size();
		if (a < T - 2 * K || a + b < 2 * T - 2 * K)
		{
			qR.push(i);
			Rflag[i] = false;
		}
	}

	while (qL.size() != 0 || qR.size() != 0)
	{
		if (qL.size() != 0)
		{
			int ni = qL.front();
			for (int i = 0; i < LPN[ni].size(); i++)
			{
				if (Lflag[LPN[ni][i]] == false)
					continue;
				LPdegree[LPN[ni][i]]--;
				if (LPdegree[LPN[ni][i]] < T - 2 * K || LNdegree[LPN[ni][i]] + LPdegree[LPN[ni][i]] < 2 * T - 2 * K)
				{
					Lflag[LPN[ni][i]] = false;
					qL.push(LPN[ni][i]);
				}
			}
			for (int i = 0; i < LNN[ni].size(); i++)
			{
				if (Rflag[LNN[ni][i]] == false)
					continue;
				RNdegree[LNN[ni][i]]--;
				if (RNdegree[LNN[ni][i]] + RPdegree[LNN[ni][i]] < 2 * T - 2 * K)
				{
					Rflag[LNN[ni][i]] = false;
					qR.push(LNN[ni][i]);
				}
			}
			qL.pop();
		}
		else
		{
			int ni = qR.front();
			for (int i = 0; i < RPN[ni].size(); i++)
			{
				if (Rflag[RPN[ni][i]] == false)
					continue;
				RPdegree[RPN[ni][i]]--;
				if (RPdegree[RPN[ni][i]] < T - 2 * K || RNdegree[RPN[ni][i]] + RPdegree[RPN[ni][i]] < 2 * T - 2 * K)
				{
					Rflag[RPN[ni][i]] = false;
					qR.push(RPN[ni][i]);
				}
			}
			for (int i = 0; i < RNN[ni].size(); i++)
			{
				if (Lflag[RNN[ni][i]] == false)
					continue;
				LNdegree[RNN[ni][i]]--;
				if (LNdegree[RNN[ni][i]] + LPdegree[RNN[ni][i]] < 2 * T - 2 * K)
				{
					Lflag[RNN[ni][i]] = false;
					qL.push(RNN[ni][i]);
				}
			}
			qR.pop();
		}
	}
	high_resolution_clock::time_point Dichromatictwo_s = high_resolution_clock::now();

	vector<int> Lone;
	set<int> Ltwo;
	vector<int> Rone;
	set<int> Rtwo;
	for (int i = 0; i < L.size(); i++) // L
	{
		if (Lflag[i] == false)
			continue;
		Lone.push_back(L[i]);
		Ln.insert(L[i]);
		for (int j = 0; j < nlist[L[i]].Pdegree; j++) 
		{
			Ltwo.insert(Padjncy[nlist[L[i]].Pid + j]);
		}
		for (int j = 0; j < nlist[L[i]].Ndegree; j++) 
		{
			Rtwo.insert(Nadjncy[nlist[L[i]].Nid + j]);
		}
	}
	for (int i = 0; i < R.size(); i++) // R
	{
		if (Rflag[i] == false)
			continue;
		Rone.push_back(R[i]);
		Rn.insert(R[i]);
		for (int j = 0; j < nlist[R[i]].Pdegree; j++) 
		{
			Rtwo.insert(Padjncy[nlist[R[i]].Pid + j]);
		}
		for (int j = 0; j < nlist[R[i]].Ndegree; j++)
		{
			Ltwo.insert(Nadjncy[nlist[R[i]].Nid + j]);
		}
	}
	vector<int> Ltwov;
	Ltwov.assign(Ltwo.begin(), Ltwo.end());
	vector<int> Rtwov;
	Rtwov.assign(Rtwo.begin(), Rtwo.end());
	vector<int> nofnv;
	vectors_union(L, R, nofnv);
	vector<int> Ltwov1;
	vector<int> Rtwov1;
	vectors_difference(Ltwov, nofnv, Ltwov1);
	vectors_difference(Rtwov, nofnv, Rtwov1);
	vector<int> emptyset;
	for (int i = 0; i < Ltwov1.size(); i++)
	{
		vector<int> Pneiofi;
		vectorPneighbor(Ltwov1[i], Pneiofi);
		int a = vectors_intersection(Pneiofi, Lone, emptyset);
		if (a >= T - 2 * K)
		{
			vector<int> Nneiofi;
			vectorNneighbor(Ltwov1[i], Nneiofi);
			int b = vectors_intersection(Nneiofi, Rone, emptyset);
			if (a + b >= 2 * T - 2 * K)
			{
				Ln.insert(Ltwov1[i]);
			}
		}
	}
	for (int i = 0; i < Rtwov1.size(); i++)
	{
		vector<int> Pneiofi;
		vectorPneighbor(Rtwov1[i], Pneiofi);
		int a = vectors_intersection(Pneiofi, Rone, emptyset);
		if (a >= T - 2 * K)
		{
			vector<int> Nneiofi;
			vectorNneighbor(Rtwov1[i], Nneiofi);
			int b = vectors_intersection(Nneiofi, Lone, emptyset);
			if (a + b >= 2 * T - 2 * K)
			{
				Rn.insert(Rtwov1[i]);
			}
		}
	}
	high_resolution_clock::time_point Dichromatictwo_e = high_resolution_clock::now();
	Dichromatictwotime += std::chrono::duration_cast<std::chrono::microseconds>(Dichromatictwo_e - Dichromatictwo_s).count();
}

bool Graph::test_difference(vector<int> &a, int b, int c, bool f)
{
	vector<int> bn;
	if (f)
	{
		for (int i = 0; i < nlist[b].Pdegree; i++)
		{
			bn.push_back(Padjncy[nlist[b].Pid + i]);
		}
		for (int i = 0; i < nlist[c].Pdegree; i++)
		{
			bn.push_back(Padjncy[nlist[c].Pid + i]);
		}
	}
	else
	{
		for (int i = 0; i < nlist[b].Ndegree; i++)
		{
			bn.push_back(Nadjncy[nlist[b].Nid + i]);
		}
		for (int i = 0; i < nlist[c].Ndegree; i++)
		{
			bn.push_back(Nadjncy[nlist[c].Nid + i]);
		}
	}
	vector<int> mm;
	if (vectors_difference(a, bn, mm) == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Graph::pivotselect(vector<int> &p, vector<int> &pn, vector<int> &c, int u, int f, vector<int> &r)
{
	r.clear();
	vector<int> w;
	vectors_intersection(p, pn, w);
	for (vector<int>::iterator i = w.begin(); i != w.end(); i++)
	{
		if (u == (*i))
			continue;
		if (test_difference(c, u, (*i), f))
		{
			r.push_back((*i));
		}
	}
}

void showvec(string s, vector<int> v1)
{
    cout << s;
    for (auto &i : v1)
	{
		cout<<i<<" ";
	}
    cout << endl;
}

void KPlexEnumUtil(Graph *subg, vector<int> CL, vector<int> CR, vector<int> PL, vector<int> PR, vector<int> QL, vector<int> QR, vector<int> CLdegree, vector<int> CRdegree, int level)
{
	cost++;
	int f = 0;
	vector<int> PLnew;
	vector<int> PRnew;
	vector<int> QLnew;
	vector<int> QRnew;
	vector<int> emptyv;
	for (vector<int>::iterator i = PL.begin(); i != PL.end(); i++)
	{
		f = 0;
		for (int j = 0; j < CL.size(); j++)
		{
			if (subg->testadj((*i), CL[j], 1))
			{
				f = 1;
				break;
			}

			if (CLdegree[j] == CL.size() + CR.size() - K)
			{
				if (!subg->testadj((*i), CL[j], 0))
				{
					f = 1;
					break;
				}
			}
		}
		if (f == 1)
			continue;
		for (int k = 0; k < CR.size(); k++)
		{
			if (subg->testadj((*i), CR[k], 0))
			{
				f = 1;
				break;
			}

			if (CRdegree[k] == CL.size() + CR.size() - K)
			{
				if (!subg->testadj((*i), CR[k], 1))
				{
					f = 1;
					break;
				}
			}
		}
		if (f == 1)
			continue;
		if (subg->gvectors_intersection(CL, (*i), true, emptyv) + subg->gvectors_intersection(CR, (*i), false, emptyv) < CL.size() + CR.size() - K + 1)
			f = 1;
		if (f == 0)
		{
			PLnew.push_back((*i));
		}
	}
	for (vector<int>::iterator i = PR.begin(); i != PR.end(); i++)
	{
		f = 0;
		for (int j = 0; j < CL.size(); j++)
		{
			if (subg->testadj((*i), CL[j], 0))
			{
				f = 1;
				break;
			}

			if (CLdegree[j] == CL.size() + CR.size() - K)
			{
				if (!subg->testadj((*i), CL[j], 1))
				{
					f = 1;
					break;
				}
			}
		}
		if (f == 1)
			continue;
		for (int k = 0; k < CR.size(); k++)
		{
			if (subg->testadj((*i), CR[k], 1))
			{
				f = 1;
				break;
			}

			if (CRdegree[k] == CL.size() + CR.size() - K)
			{
				if (!subg->testadj((*i), CR[k], 0))
				{
					f = 1;
					break;
				}
			}
		}
		if (f == 1)
			continue;
		if (subg->gvectors_intersection(CL, (*i), false, emptyv) + subg->gvectors_intersection(CR, (*i), true, emptyv) < CL.size() + CR.size() - K + 1)
			f = 1;
		if (f == 0)
		{
			PRnew.push_back((*i));
		}
	}
	for (vector<int>::iterator i = QL.begin(); i != QL.end(); i++)
	{
		f = 0;
		for (int j = 0; j < CL.size(); j++)
		{
			if (subg->testadj((*i), CL[j], 1))
			{
				f = 1;
				break;
			}

			if (CLdegree[j] == CL.size() + CR.size() - K)
			{
				if (!subg->testadj((*i), CL[j], 0))
				{
					f = 1;
					break;
				}
			}
		}
		if (f == 1)
			continue;
		for (int k = 0; k < CR.size(); k++)
		{
			if (subg->testadj((*i), CR[k], 0))
			{
				f = 1;
				break;
			}

			if (CRdegree[k] == CL.size() + CR.size() - K)
			{
				if (!subg->testadj((*i), CR[k], 1))
				{
					f = 1;
					break;
				}
			}
		}
		if (f == 1)
			continue;
		if (subg->gvectors_intersection(CL, (*i), true, emptyv) + subg->gvectors_intersection(CR, (*i), false, emptyv) < CL.size() + CR.size() - K + 1)
			f = 1;
		if (f == 0)
		{
			QLnew.push_back((*i));
		}
	}
	for (vector<int>::iterator i = QR.begin(); i != QR.end(); i++)
	{
		f = 0;
		for (int j = 0; j < CL.size(); j++)
		{
			if (subg->testadj((*i), CL[j], 0))
			{
				f = 1;
				break;
			}

			if (CLdegree[j] == CL.size() + CR.size() - K)
			{
				if (!subg->testadj((*i), CL[j], 1))
				{
					f = 1;
					break;
				}
			}
		}
		if (f == 1)
			continue;
		for (int k = 0; k < CR.size(); k++)
		{
			if (subg->testadj((*i), CR[k], 1))
			{
				f = 1;
				break;
			}

			if (CRdegree[k] == CL.size() + CR.size() - K)
			{
				if (!subg->testadj((*i), CR[k], 0))
				{
					f = 1;
					break;
				}
			}
		}
		if (f == 1)
			continue;
		if (subg->gvectors_intersection(CL, (*i), false, emptyv) + subg->gvectors_intersection(CR, (*i), true, emptyv) < CL.size() + CR.size() - K + 1)
			f = 1;
		if (f == 0)
		{
			QRnew.push_back((*i));
		}
	}
	if (PLnew.size() == 0 && PRnew.size() == 0 && QLnew.size() == 0 && QRnew.size() == 0)
	{
		if (CL.size() >= T && CR.size() >= T)
		{
			cout << "\nplex" << endl;
			for (vector<int>::iterator i = CL.begin(); i != CL.end(); i++)
				cout << (*i) << " ";
			cout << endl;
			for (vector<int>::iterator i = CR.begin(); i != CR.end(); i++)
				cout << (*i) << " ";
        }
		return;
	}
	
	
	flag = !flag;
	vector<int> CLdegreenow;
	vector<int> CRdegreenow;
	if (flag)
	{
		if (PLnew.size() != 0)
		{
			for (vector<int>::iterator i = PLnew.begin(); PLnew.size()!=0; )
			{
				i= PLnew.begin();
				int id = 0;
				CLdegreenow = CLdegree;
				CRdegreenow = CRdegree;
				for (int i1 = 0; i1 < CL.size(); i1++)
				{
					if (subg->testadj(*i, CL[i1], 0))
					{
						CLdegreenow[i1]++;
						id++;
					}
				}
				for (int i2 = 0; i2 < CR.size(); i2++)
				{
					if (subg->testadj((*i), CR[i2], 1))
					{
						CRdegreenow[i2]++;
						id++;
					}
				}
				CL.push_back((*i));
				CLdegreenow.push_back(id);
				int ithis;
				ithis = (*i);
				PLnew.erase(i);
				KPlexEnumUtil(subg, CL, CR, PLnew, PRnew, QLnew, QRnew, CLdegreenow, CRdegreenow, level+1);
				CL.pop_back();
				CLdegreenow.clear();
				CRdegreenow.clear();
				QLnew.push_back(ithis);
				
			}
		}
		if (PRnew.size() != 0)
		{
			for (vector<int>::iterator i = PRnew.begin(); PRnew.size() != 0; )
			{
				i = PRnew.begin();
				int id = 0;
				CLdegreenow = CLdegree;
				CRdegreenow = CRdegree;
				for (int i1 = 0; i1 < CL.size(); i1++)
				{
					if (subg->testadj((*i), CL[i1], 1))
					{
						CLdegreenow[i1]++;
						id++;
					}
				}
				for (int i2 = 0; i2 < CR.size(); i2++)
				{
					if (subg->testadj((*i), CR[i2], 0))
					{
						CRdegreenow[i2]++;
						id++;
					}
				}
				CR.push_back((*i));
				CRdegreenow.push_back(id);
				int ithis;
				ithis = (*i);
				PRnew.erase(i);
				KPlexEnumUtil(subg, CL, CR, PLnew, PRnew, QLnew, QRnew, CLdegreenow, CRdegreenow, level+1);
				CR.pop_back();
				CLdegreenow.clear();
				CRdegreenow.clear();
				QRnew.push_back(ithis);
				
			}
		}
	}
	else
	{
		if (PRnew.size() != 0)
		{
			for (vector<int>::iterator i = PRnew.begin(); PRnew.size() != 0; )
			{
				i = PRnew.begin();
				int id = 0;
				CLdegreenow = CLdegree;
				CRdegreenow = CRdegree;
				for (int i1 = 0; i1 < CL.size(); i1++)
				{
					if (subg->testadj((*i), CL[i1], 1))
					{
						CLdegreenow[i1]++;
						id++;
					}
				}
				for (int i2 = 0; i2 < CR.size(); i2++)
				{
					if (subg->testadj((*i), CR[i2], 0))
					{
						CRdegreenow[i2]++;
						id++;
					}
				}
				CR.push_back((*i));
				CRdegreenow.push_back(id);
				int ithis;
				ithis = (*i);
				PRnew.erase(i);
				KPlexEnumUtil(subg, CL, CR, PLnew, PRnew, QLnew, QRnew, CLdegreenow, CRdegreenow, level+1);
				CR.pop_back();
				CLdegreenow.clear();
				CRdegreenow.clear();
				QRnew.push_back(ithis);
				
			}
		}
		if (PLnew.size() != 0)
		{
			for (vector<int>::iterator i = PLnew.begin(); PLnew.size()!=0 ; )
			{
				i = PLnew.begin();
				int id = 0;
				CLdegreenow = CLdegree;
				CRdegreenow = CRdegree;
				for (int i1 = 0; i1 < CL.size(); i1++)
				{
					if (subg->testadj((*i), CL[i1], 0))
					{
						CLdegreenow[i1]++;
						id++;
					}
				}
				for (int i2 = 0; i2 < CR.size(); i2++)
				{
					if (subg->testadj((*i), CR[i2], 1))
					{
						CRdegreenow[i2]++;
						id++;
					}
				}
				CL.push_back((*i));
				CLdegreenow.push_back(id);
				int ithis;
				ithis = (*i);
				PLnew.erase(i);
				KPlexEnumUtil(subg, CL, CR, PLnew, PRnew, QLnew, QRnew, CLdegreenow, CRdegreenow, level+1);
				CL.pop_back();
				CLdegreenow.clear();
				CRdegreenow.clear();
				QLnew.push_back(ithis);
				
			}
		}
	}
}

void KPlexEnum(Graph *subg, int k, int t)
{
	vector<int> CL;
	vector<int> CLdegree;
	vector<int> CR;
	vector<int> CRdegree;
	vector<int> PL;
	vector<int> PR;
	vector<int> QL;
	vector<int> QR;
	for (int i = 0; i < NODENUM; i++)
	{
		if (subg->isExistinG(i) == false)
			continue;
		CL.clear();
		CLdegree.clear();
		CR.clear();
		CRdegree.clear();
		PL.clear();
		PR.clear();
		QL.clear();
		QR.clear();
		CL.push_back(i);
		CLdegree.push_back(0);
		set<int> L;
		set<int> R;
		high_resolution_clock::time_point Dichromatic_s = high_resolution_clock::now();
		subg->Dichromatic(i, L, R);
		high_resolution_clock::time_point Dichromatic_e = high_resolution_clock::now();
		Dichromatictime += std::chrono::duration_cast<std::chrono::microseconds>(Dichromatic_e - Dichromatic_s).count();
		for (set<int>::iterator y = L.begin(); y != L.end(); y++)
		{
			if ((*y) > i)
				PL.push_back((*y));
			if ((*y) < i)
				QL.push_back((*y));
		}
        if (PL.size() + 1 < T)
			continue;
		for (set<int>::iterator y = R.begin(); y != R.end(); y++)
		{
			if ((*y) > i)
				PR.push_back((*y));
			if ((*y) < i)
				QR.push_back((*y));
		}
		if (PR.size() < T)
			continue;
		sort(PL.begin(), PL.end());
		sort(QL.begin(), QL.end());
		sort(PR.begin(), PR.end());
		sort(QR.begin(), QR.end());
		KPlexEnumUtil(subg, CL, CR, PL, PR, QL, QR, CLdegree, CRdegree,1);
	}
}

void Graph::removefromneiadj(int v, int u, bool f) 
{
	if (f)
	{
		for (int i = 0; i < nlist[v].Pdegree; i++)
		{
			if (u == Padjncy[nlist[v].Pid + i])
			{
				for (int j = i + 1; j < nlist[v].Pdegree; j++)
				{
					Padjncy[nlist[v].Pid + j - 1] = Padjncy[nlist[v].Pid + j];
				}
				nlist[v].Pdegree--;
				break;
			}
		}
	}
	else
	{
		for (int i = 0; i < nlist[v].Ndegree; i++)
		{
			if (u == Nadjncy[nlist[v].Nid + i])
			{
				for (int j = i + 1; j < nlist[v].Ndegree; j++)
				{
					Nadjncy[nlist[v].Nid + j - 1] = Nadjncy[nlist[v].Nid + j];
				}
				nlist[v].Ndegree--;
				break;
			}
		}
	}
}

void Graph::VertexReduction(int k, int t)
{
	queue<int> q;
	for (int i = 0; i < NODENUM; i++)
	{
		if (nlist[i].Pdegree < T - K || nlist[i].Ndegree < T - K + 1 || nlist[i].Pdegree + nlist[i].Ndegree < 2 * T - K)
		{
			q.push(i);
			nlist[i].isExist = false;
		}
	}
	int dnum = 0;
	while (q.size() != 0)
	{
		int ni = q.front();
		for (int i = 0; i < nlist[ni].Pdegree; i++)
		{
			int ti = Padjncy[nlist[ni].Pid + i];
			removefromneiadj(ti, ni, true);
			if (nlist[ti].isExist == true && (nlist[ti].Pdegree < T - K || nlist[ti].Pdegree + nlist[ti].Ndegree < 2 * T - K))
			{
				q.push(ti);
				nlist[ti].isExist = false;
			}
		}
		for (int i = 0; i < nlist[ni].Ndegree; i++)
		{
			int ti = Nadjncy[nlist[ni].Nid + i];
			removefromneiadj(ti, ni, false);
			if (nlist[ti].isExist == true && (nlist[ti].Ndegree < T - K + 1 || nlist[ti].Pdegree + nlist[ti].Ndegree < 2 * T - K))
			{
				q.push(ti);
				nlist[ti].isExist = false;
			}
		}
		q.pop();
		dnum++;
	}
	cout << "dnum=" << dnum << endl;
	for (int i = 0; i < NODENUM; i++)
	{
		if (nlist[i].isExist == true)
		{
			nodesetid nid;
			nid.id = i;
			nid.degree = min(nlist[i].Pdegree, nlist[i].Ndegree);
			nodeset.push_back(nid);
		}
	}
	sort(nodeset.begin(), nodeset.end(), comp1);
	vector<int> newid(NODENUM, -1); 

	for (int i = 0; i < nodeset.size(); i++)
	{
		newid[nodeset[i].id] = i; 
	}
	node *newnlist = new node[nodeset.size() + 1];
	for (int i = 0; i < nodeset.size(); i++)
	{
		newnlist[i].rid = nlist[nodeset[i].id].rid;
		newnlist[i].Pid = nlist[nodeset[i].id].Pid;
		newnlist[i].Nid = nlist[nodeset[i].id].Nid;
		newnlist[i].Pdegree = nlist[nodeset[i].id].Pdegree;
		newnlist[i].Ndegree = nlist[nodeset[i].id].Ndegree;
		newnlist[i].isExist = nlist[nodeset[i].id].isExist;
		vector<int> neiborsort;
		for (int j = 0; j < newnlist[i].Pdegree; j++)
		{
			neiborsort.push_back(newid[Padjncy[newnlist[i].Pid + j]]);
		}
		sort(neiborsort.begin(), neiborsort.end());
		for (int j = 0; j < newnlist[i].Pdegree; j++)
		{
			Padjncy[newnlist[i].Pid + j] = neiborsort[j];
		}
		neiborsort.clear();
		for (int j = 0; j < newnlist[i].Ndegree; j++)
		{
			neiborsort.push_back(newid[Nadjncy[newnlist[i].Nid + j]]);
		}
		sort(neiborsort.begin(), neiborsort.end());
		for (int j = 0; j < newnlist[i].Ndegree; j++)
		{
			Nadjncy[newnlist[i].Nid + j] = neiborsort[j];
		}
	}
	if (nlist != NULL)
	{
		delete nlist;
	}
	nlist = newnlist;
	NODENUM = nodeset.size();
}

void plexcheck(Graph *subg)
{
    cout <<"plexnum="<< Lplex.size() << endl;
    int f = 0;
    for (int i = 0; i < Lplex.size();i++)
    {
        int pd = 0, nd = 0;
        for (int j = 0; j < Lplex[i].size();j++)
        {
            for(int j1=0;j1 < Lplex[i].size();j1++)
            {
                if(j1==j)
                    continue;
                if(subg->testadj(Lplex[i][j],Lplex[i][j1],1))
                    f = 1;
                if(subg->testadj(Lplex[i][j],Lplex[i][j1],0))
                    pd++;
            }
            for (int a = 0; a < Rplex[i].size();a++)
            {
                if(subg->testadj(Lplex[i][j],Rplex[i][a],0))
                {
                    f = 1;
                }
                if(subg->testadj(Lplex[i][j],Rplex[i][a],1))
                {
                    nd++;
                }
            }
        }
        if(pd<Lplex[i].size()-K || pd+nd<Lplex[i].size()+Rplex[i].size()-K)
            f = 1;
        pd = 0, nd = 0;
        for (int j = 0; j < Rplex[i].size();j++)
        {
            for(int j1=0;j1 < Rplex[i].size();j1++)
            {
                if(j1==j)
                    continue;
                if(subg->testadj(Rplex[i][j],Rplex[i][j1],1))
                    f = 1;
                if(subg->testadj(Rplex[i][j],Rplex[i][j1],0))
                    pd++;
            }
            for (int a = 0; a < Lplex[i].size();a++)
            {
                if(subg->testadj(Rplex[i][j],Lplex[i][a],0))
                {
                    f = 1;
                }
                if(subg->testadj(Rplex[i][j],Lplex[i][a],1))
                {
                    nd++;
                }
            }
        }
        if(pd<Rplex[i].size()-K || pd+nd<Lplex[i].size()+Rplex[i].size()-K)
            f = 1;
        if(f==1)
        {
            cout << "wrong" << endl;
            showvec("CL", Lplex[i]);
            showvec("CR", Rplex[i]);
        }
    }
}
int main(int argc, char **argv)
{
	Graph *g = new Graph();
	cout << "Reading nverts from file " << argv[1] << endl;
	high_resolution_clock::time_point mkGraph_s = high_resolution_clock::now();
	g->mkGraph(argv[1]);
	T = atoi(argv[2]);
	K = atoi(argv[3]);
    cout << "T=" << T << "K=" << K << endl;
    high_resolution_clock::time_point mkGraph_e = high_resolution_clock::now();
	auto mkGraphtime = std::chrono::duration_cast<std::chrono::microseconds>(mkGraph_e - mkGraph_s).count();
	cout << "mkGraph time: " << mkGraphtime / 1e6 << endl;

	high_resolution_clock::time_point VertexReduction_s = high_resolution_clock::now();
	g->VertexReduction(K, T);
	high_resolution_clock::time_point VertexReduction_e = high_resolution_clock::now();
	auto VertexReductiontime = std::chrono::duration_cast<std::chrono::microseconds>(VertexReduction_e - VertexReduction_s).count();
	cout << "\nVertexReduction time: " << VertexReductiontime / 1e6 << endl;

	high_resolution_clock::time_point KPlexEnum_s = high_resolution_clock::now();
	KPlexEnum(g, K, T);
	high_resolution_clock::time_point KPlexEnum_e = high_resolution_clock::now();
	auto KPlexEnumtime = std::chrono::duration_cast<std::chrono::microseconds>(KPlexEnum_e - KPlexEnum_s).count();
	cout << "\nKPlexEnum time: " << KPlexEnumtime / 1e6 << endl;

	cout << "\nDichromatic time: " << Dichromatictime / 1e6 << endl;
	cout << "\nDichromatictwo time: " << Dichromatictwotime / 1e6 << endl;
	cout << "\ncolor time: " << colortime / 1e6 << endl;
	cout << "\ncolorcheck time: " << colorchecktime / 1e6 << endl;

	cout << "cost=" << cost << endl;
    //plexcheck(g);

    return 0;
}