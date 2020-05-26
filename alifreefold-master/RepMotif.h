/*
 * RepMotif.h
 *
 *  Created on: Jul 28, 2014
 *      Author: Séhi
 */

#ifndef REPMOTIF_H_
#define REPMOTIF_H_

#include <string>
#include <cctype>
#include <vector>
#include <stack>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <map>
#include <time.h>
#include <set>
#include <iterator>
#include <boost/graph/adjacency_list.hpp>
#include "math.h"
#include <unordered_set>
#include <unordered_map>
#include "Eigen/Dense"
#include "stdlib.h"

extern "C"{
#include  <ViennaRNA/data_structures.h>
#include  <ViennaRNA/subopt.h>
#include  <ViennaRNA/boltzmann_sampling.h>

};

extern int minNumberOfSubOpt;

using namespace std;

class RepNmotif {

private:
        int maxLevelNmotifs;
		string i_PathFileInput;
		bool isCircular;
		vector<string> headers,structures,sequences;
        vector< map<string,float> > nmotifsForEachStructure;
        map<string,float> nmotifsAllStructure;
        vector< map<string,vector<int>> > nmotifsForEachStructureWithPosNucOfnmotifs;

public:

	typedef boost::adjacency_list<
       boost::vecS, boost::vecS, boost::undirectedS> graph;

	typedef boost::graph_traits<graph>::vertex_descriptor vertex_t;
	typedef boost::graph_traits<graph>::vertex_iterator vertex_iterator_t;
	typedef boost::graph_traits<graph>::edge_iterator edge_iterator_t;
	typedef boost::graph_traits<graph>::adjacency_iterator adjacency_iterator_t;

	/*Motif Properties*/
	struct Motif
	{
 	 string motifName;
 	 string motifNameWithFeatLevel1;
 	 string motifNameWithFeatLevel2;
	 vector<int> posNuc;
	};

	RepNmotif();

	RepNmotif(vector<string>, int, bool);

	//PART2
	static bool isAlignmentFile(string);

	//RepNmotif(string,int,int);
	
	//RepNmotif(string,int);

	virtual ~RepNmotif();

	void readFileandParseSequenceStructure(string, vector<string>&, vector<string>&,
		vector<string>&, bool&, int&);

	void readFileandParseAlignmentStructure(string, vector<string>&, vector<string>&,
		vector<string>&, bool&);

    void parseSequenceStructure(string,string);

    void dbn2Grss(string,graph&);

    void specialCarac2SpecialMotifVect(string & ,bool&,vector<RepNmotif::Motif>&);

    void grss2motifVect(string,graph, vector<RepNmotif::Motif>&);

    void intersectMotifs(int, string, vector<RepNmotif::Motif>&, map<string,float>& , map<string, vector<int>>&, map<string,float>&);

    void extractFeatFromCurrentMotif(Motif currentMotif,map<string,float>& ,map<string, vector<int>>&,map<string, float> &);

	//getters
    const map<string, float>& getNmotifsAllStructure() const {
		return nmotifsAllStructure;
	}

    const vector<map<string, float> >& getNmotifsForEachStructure() const {
		return nmotifsForEachStructure;
	}

    const vector<map<string, vector<int> > >& getNmotifsForEachStructureWithPosNucOfnmotifs() const {
		return nmotifsForEachStructureWithPosNucOfnmotifs;
	}

	const vector<string>& getHeaders() const {
		return headers;
	}

	const vector<string>& getSequences() const {
		return sequences;
	}

	const vector<string>& getStructures() const {
		return structures;
	}

};

void initCurrentMotif(RepNmotif::Motif&);

void initNbUnpDegBpSideABL1L2L3symLoopPosNuc(unsigned int& , unsigned int& , unsigned int& , unsigned int& , unsigned int& , bool& , unsigned int& , unsigned int& , unsigned int& , unsigned int& , vector<int>& );

RepNmotif::vertex_iterator_t advanceVertex_iterator_t2adjacency_iterator_t(RepNmotif::vertex_iterator_t, RepNmotif::adjacency_iterator_t);

void parseTempLine(string,string,string&,string&);

string discretize(int);


#endif /* REPMOTIF_H_ */
