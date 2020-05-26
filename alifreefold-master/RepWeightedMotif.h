/*
 * RepWeightedMotif.h
 *
 *  Created on: Jul 28, 2014
 *      Author: Séhi
 */

#ifndef REPWEIGHTEDMOTIF_H_
#define REPWEIGHTEDMOTIF_H_

#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <math.h>       /* fabs */

#include "Eigen/Dense"
#include "Eigen/Core"

extern int minNumberOfSubOpt;

using namespace std;

class RepWeightedMotif {

    Eigen::ArrayXXf matALLMotifWeigthed;// matrice n-motif pondéré (25*n * m)
    Eigen::ArrayXXf matALLMotif; //matrice n-motif (25*n * m)
	Eigen::ArrayXXf distanceToEveryBody; //matrice (25*n * 25*n), distance de toutes les SS aux autres SS

    Eigen::RowVectorXi subOptIdx; //indice de début des SS dans la matrice des n-motif
    map<string,float> allStructureFeatureForWeigthOfMotifs; //association des poids et des motifs (utilisés pour la matrice n-motif)
    vector< map<string,float> > nmotifsForEachStructure; //occurence des motifs dans chaque SS (vecteur d'associations)
    vector<string> headers,structures,sequences; //structures SS des différentes séquences
	vector<string> headerOfRepresentatives1, sequenceOfRepresentatives1, structureOfRepresentatives1; //structures choisies parmi les 25 de chaque. Vecteur augmente en taille au fil des itérations
	vector<string> headerOfRepresentatives2, sequenceOfRepresentatives2, structureOfRepresentatives2; //structures choisies parmi les 25 de chaque. Vecteur augmente en taille au fil des itérations
    string headerOfRepresentative;
    string sequenceOfRepresentative;
    string structureOfRepresentative;
	Eigen::RowVectorXf centroid;
	Eigen::RowVectorXf groupCentroid;
    Eigen::RowVectorXf distToCentroid;
    Eigen::RowVectorXf allWeight; //version vecteur de "allStructureFeatureForWeigthOfMotifs"

	vector<string> nmotifsHeader;				//>> n-motifs header of the "w_nmRepOfSS.csv" file
	Eigen::ArrayXXf matAllNMotifsWeighted;		//>> matrix of all weighted n-motifs of the "w_nmRepOfSS.csv" file
	vector<float> rowAllWeights;				//>> row of all weights of the "weights.csv" file 
	vector<int> rowAllPositionsNMotifs;			//>> row of all positions of n-motifs in both "r_nmRepOfAS.csv" and "weights.csv"
	Eigen::ArrayXXf matDistances;				//>> matrix of all distances between alignment and subopt
	vector<string> closestSuboptHeaders;		//>> alignment headers closest to subopt headers
	vector<string> closestSuboptSequences;		//>> alignment sequences closest to subopt sequences 
	vector<string> closestSuboptStructures;		//>> alignment structures closest to subopt structures

public:
	//-- Constructor --
	RepWeightedMotif();
	RepWeightedMotif(vector<string>, map<string, float>, vector<map<string, float>>,
		vector<string>, vector<string>, vector<string>, int, int, bool, bool, bool);
	RepWeightedMotif(map<string, float>, vector<map<string, float>>, vector<string>,
		vector<string>, vector<string>, int, int, bool, bool);

	//-- Destructor --
	virtual ~RepWeightedMotif();
	
	//-- Methods --
	Eigen::ArrayXXf filterMotifs(Eigen::RowVectorXi, int);

	static bool isWNMotifFile(string& filePath);
	static bool isWeightsFile(string& filePath);
	static bool isSuboptFile(string& filePath);

	void readWeightedNMotifsSuboptStructures(string& filePath);
	void readWeightsNMotifs(string& filePath);
	void readSubopts(string& filePath);

	void findPositionsNMotifs(map<string,float>& nmotifsHeaderFirstRowOccurrencesAlignment, 
		vector<string>& nmotifsHeaderSubopt);
	void createMatWeightedNMotifsAlignment(Eigen::ArrayXXf& matOccurrencesNMotifAlignment, 
		vector<float>& rowAllWeightsSubopt, vector<int>& rowAllPositionsNMotifs);
	void createMatDistancesAlignment(Eigen::ArrayXXf& matAllNMotifsWeighted,
		Eigen::ArrayXXf& matWeightedNMotifsAlignment, vector<int>& rowAllPositionsNMotifs);
	void computeClosestAlignmentToSubopts(Eigen::ArrayXXf& matDistancesAlignmentSubopts, 
		vector<string>& suboptHeaders, vector<string>& suboptSequences, vector<string>& suboptStructures);

    void computeRepresentative(Eigen::ArrayXXf, vector<string>, vector<string>, vector<string>, 
		bool, bool);
	void computeClosestStructuresToCentroid(Eigen::RowVectorXf, Eigen::RowVectorXi, vector<string>, 
		vector<string>, vector<string>, bool, bool);
	void computeClosestStructuresToGroup(Eigen::ArrayXXf, Eigen::RowVectorXi, vector<string>, 
		vector<string>, vector<string>, bool, bool);

	//-- Getters --
    const map<string, float>& getAllStructureFeatureForWeigthOfMotifs() const {
		return allStructureFeatureForWeigthOfMotifs;
	}

	const Eigen::ArrayXXf& getMatALLMotif() const {
		return matALLMotif;
	}
    const Eigen::ArrayXXf& getMatAllMotifWeigthed() const {
		return matALLMotifWeigthed;
	}

	const Eigen::RowVectorXf& getAllWeight() const {
		return allWeight;
	}

    const Eigen::RowVectorXi& getSubOptIdx() const {
		return subOptIdx;
	}

    const string& getHeaderOfRepresentative() const {
		return headerOfRepresentative;
	}
    const string& getSequenceOfRepresentative() const {
		return sequenceOfRepresentative;
	}
    const string& getStructureOfRepresentative() const {
		return structureOfRepresentative;
	}

	const Eigen::RowVectorXf& getCentroid() const {
		return centroid;
	}
	const Eigen::RowVectorXf& getGroupCentroid() const {
		return groupCentroid;
	}
    const Eigen::RowVectorXf& getDistToCentroid() const {
		return distToCentroid;
	}
	const Eigen::ArrayXXf& getDistance() const {
		return distanceToEveryBody;
	}

	const Eigen::ArrayXXf& getMatDistancesAlignmentsSubopts() const {
		return matDistances;
	}

	const vector<string>& getHeaderOfRepresentatives1() const {
		return headerOfRepresentatives1;
	}
	const vector<string>& getSequenceOfRepresentatives1() const {
		return sequenceOfRepresentatives1;
	}
	const vector<string>& getStructureOfRepresentatives1() const {
		return structureOfRepresentatives1;
	}
	
	const vector<string>& getHeaderOfRepresentatives2() const {
		return headerOfRepresentatives2;
	}
	const vector<string>& getSequenceOfRepresentatives2() const {
		return sequenceOfRepresentatives2;
	}
	const vector<string>& getStructureOfRepresentatives2() const {
		return structureOfRepresentatives2;

	}

	const vector<string>& getClosestSuboptHeaders() const {
		return closestSuboptHeaders;
	}
	const vector<string>& getClosestSuboptSequences() const {
		return closestSuboptSequences;
	}
	const vector<string>& getClosestSuboptStructures() const {
		return closestSuboptStructures;
	}
};

#endif /* REPWEIGHTEDMOTIF_H_ */