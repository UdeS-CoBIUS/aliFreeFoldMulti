/*
 * RepWeightedMotif.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: Séhi
 */

#include "RepWeightedMotif.h"

/*********************************************************************************/

RepWeightedMotif::RepWeightedMotif() {
	// TODO Auto-generated constructor stub
}


/*********************************************************************************/
RepWeightedMotif::RepWeightedMotif(map<string,float> allStructureFeature, vector< map<string,float> > nmotifsForEachStructure,vector<string> headers,vector<string> sequences,vector<string> structures,int minNumberOfSubOpt,int computeConservationIndex, bool verbose, bool verbose2) {
	// TODO Auto-generated constructor stub
    this->allStructureFeatureForWeigthOfMotifs=allStructureFeature;
    this->nmotifsForEachStructure=nmotifsForEachStructure;
    //this->sequenceOfRepresentative="";
    //this->structureOfRepresentative="";
    this->headers=headers;
    this->sequences=sequences;
    this->structures=structures;

    unsigned int nrows=nmotifsForEachStructure.size();
    unsigned int k=0;
    Eigen::RowVectorXi suboptidx(nrows);
    int nbUniqueSeq=nrows/minNumberOfSubOpt;

	if(verbose2){cout << "\nsuboptidx"; cout << endl;}
	
    for (int i=0;i<nbUniqueSeq;i++)
    {
        for (int j=0;j<minNumberOfSubOpt;j++)
        {
            suboptidx(k)=i;
			if(verbose2){cout << suboptidx(k);cout << ", ";}
            k++;
        }
    }

	if(verbose || verbose2){cout << endl;if(verbose2){cout << endl;}}
	
    this->subOptIdx=suboptidx;

    this->matALLMotifWeigthed=filterMotifs(suboptidx,computeConservationIndex);

	//Version initiale
    RepWeightedMotif::computeRepresentative(matALLMotifWeigthed,headers,sequences,structures,verbose,verbose2);
	
	//Stratégie 1
	RepWeightedMotif::computeClosestStructuresToCentroid(distToCentroid, subOptIdx, headers, sequences, structures,verbose,verbose2);
	
	//Stratégie 2
	RepWeightedMotif::computeClosestStructuresToGroup(matALLMotifWeigthed, subOptIdx, headers, sequences, structures,verbose,verbose2);

   // RepWeightedMotif::computeRepresentative(this->matALLMotifWeigthed,this->headers,this->sequences,this->structures);

}

//-- TEST --
/*********************************************************************************/
RepWeightedMotif::RepWeightedMotif(vector<string> allInputFilePaths, 
	map<string,float> allStructureFeature, vector< map<string,float> > nmotifsForEachStructure, 
	vector<string> headers, vector<string> sequences, vector<string> structures, 
	int minNumberOfSubOpt, int computeConservationIndex, bool isPartialTask, bool verbose, 
	bool verbose2)
{
	this->allStructureFeatureForWeigthOfMotifs=allStructureFeature;
    this->nmotifsForEachStructure=nmotifsForEachStructure;

    unsigned int nrows=nmotifsForEachStructure.size();
    unsigned int k=0;
    Eigen::RowVectorXi suboptidx(nrows);
    int nbUniqueSeq=nrows/minNumberOfSubOpt;

	if(verbose2){cout << "\nsuboptidx"; cout << endl;}

    for (int i=0;i<nbUniqueSeq;i++)
    {
        for (int j=0;j<minNumberOfSubOpt;j++)
        {
            suboptidx(k)=i;
			if(verbose2){cout << suboptidx(k);cout << ", ";}
            k++;
        }
    }

	if(verbose || verbose2){cout << endl;if(verbose2){cout << endl;}}
	
    this->subOptIdx=suboptidx;
    this->matALLMotifWeigthed=filterMotifs(suboptidx,computeConservationIndex);

    vector<string>::iterator it;

    it = find_if(allInputFilePaths.begin(), allInputFilePaths.end(), isWNMotifFile);
    string wNMotifsFilePath = *it;
    it = find_if(allInputFilePaths.begin(), allInputFilePaths.end(), isWeightsFile);
    string weightsFilePath = *it;
	it = find_if(allInputFilePaths.begin(), allInputFilePaths.end(), isSuboptFile);
	string suboptFilePath = *it;

    readWeightedNMotifsSuboptStructures(wNMotifsFilePath);
	readWeightsNMotifs(weightsFilePath);
	readSubopts(suboptFilePath);

	// Find the n-motifs vector shared in both "r_nmRepOfAS.csv" and "w_nmRepOfSS.csv"
	findPositionsNMotifs(allStructureFeatureForWeigthOfMotifs, nmotifsHeader);

	// Create the weighted n-motifs matrix for the alignments
	createMatWeightedNMotifsAlignment(matALLMotif, rowAllWeights, rowAllPositionsNMotifs);

	// Calculate distances between "w_nmRepOfSS.csv" and "w_nmRepOfAS.csv"
	createMatDistancesAlignment(matAllNMotifsWeighted, matALLMotifWeigthed, rowAllPositionsNMotifs);

	// Compute closest alignment structures to subopt structures
	computeClosestAlignmentToSubopts(matDistances, this->headers, this->sequences, this->structures);

}

/*********************************************************************************/
/*
Inputs:
	- string filePath: File path in the -i parameter
Outputs:
	- None
Doing:
	+ Return true if the file path contains the string "w_nmRepOfSS.csv"
	+ Return false else
*/
bool RepWeightedMotif::isWNMotifFile(string& filePath) {
	if (filePath.find("w_nmRepOfSS.csv") != string::npos)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*********************************************************************************/
/*
Inputs:
	- string filePath: File path in the -i parameter
Outputs:
	- None
Doing:
	+ Return true if the file path contains the string "weights.csv"
	+ Return false else 
*/
bool RepWeightedMotif::isWeightsFile(string& filePath) {
	if (filePath.find("weights.csv") != string::npos)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*********************************************************************************/
/*
Inputs:
	- string filePath: File path in the -i parameter
Outputs:
	- None
Doing:
	+ Return true if the file path contains the string "subopt.db"
	+ Return false else
*/
bool RepWeightedMotif::isSuboptFile(string& filePath) {
	if (filePath.find("subopt.db") != string::npos)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*********************************************************************************/
/*
Inputs:
	- string filePath: "w_nmRepOfSS.csv" file path
Outputs:
	- None
Doing:
	+ Read the "w_nmRepOfSS.csv" file 
	+ Push in vector the name of n-motifs from the "w_nmRepOfSS.csv" file 
	+ Push in vector the name of sequences (headers) from the "w_nmRepOfSS.csv" file 
	+ Push in matrix the weighted n-motif values for each sequences from the
		"w_nmRepOfSS.csv" file
*/
void RepWeightedMotif::readWeightedNMotifsSuboptStructures(string& filePath) {
	vector<vector<float>> tempMatAllNMotifsWeighted;
	string line;
	int lineCounter = 0;

	ifstream wNMotifsFile(filePath.c_str());

	while (getline(wNMotifsFile, line)) {
		stringstream lineStream(line);

		if (lineCounter == 0)
		{
			string nmotif;

			while(getline(lineStream, nmotif, ','))
			{
				nmotifsHeader.push_back(nmotif);
			}
		}
		else 
		{
			string value;
			string::size_type sizeString;
			vector<float> rowAllNMotifsWeighted;
			int valueCounter = 0;

			while (getline(lineStream, value, ','))
			{
				if (valueCounter == 0)
				{
					headers.push_back(value);
				}
				else 
				{
					rowAllNMotifsWeighted.push_back(stof(value, &sizeString));
				}
				valueCounter += 1;
			}
			tempMatAllNMotifsWeighted.push_back(rowAllNMotifsWeighted);
			
		}
		lineCounter += 1;
	}
	int nRows = tempMatAllNMotifsWeighted.size();
	int nCols = nmotifsHeader.size();

	Eigen::ArrayXXf matAllNMotifsWeighted = Eigen::ArrayXXf::Zero(nRows, nCols);

	for (unsigned int row = 0; row < nRows; row++)
	{
		for (unsigned int col = 0; col < nCols; col++)
		{
			matAllNMotifsWeighted(row, col) = tempMatAllNMotifsWeighted[row][col];
		}
	}
	this->matAllNMotifsWeighted = matAllNMotifsWeighted;
}

/*********************************************************************************/
/*
Inputs:
	- string filePath: "weights.csv" file path 
Outputs:
	- None
Doing:
	+ Read the "weights.csv" file 
	+ Verify the name of n-motifs from the "weights.csv" file are the same n-motifs 
		from the "w_nmRepOfSS.csv" file
	+ Push in vetor the weights of n-motifs from the "weights.csv" file
*/
void RepWeightedMotif::readWeightsNMotifs(string& filePath) {	
	string line;
	int lineCounter = 0;

	ifstream weightsNMotifsFile(filePath.c_str());

	while (getline(weightsNMotifsFile, line)) {
		stringstream lineStream(line);

		if (lineCounter == 0) 
		{
			string nmotif;
			int nmotifCounter = 0;

			while(getline(lineStream, nmotif, ','))
			{
				if (nmotifsHeader[nmotifCounter + 1] != nmotif)
				{
					cout << "Error" << endl;
					cout << "The n-motifs header are different between the files" << endl;
					exit(EXIT_FAILURE);
				}
				nmotifCounter += 1;
			}
		}
		else 
		{
			string weight;
			string::size_type sizeString;

			while (getline(lineStream, weight, ','))
			{
				rowAllWeights.push_back(stof(weight, &sizeString));
			}
		}
		lineCounter += 1;
	}
}

/*********************************************************************************/
/*
Inputs:
	- string filePath: "subopt.db" file path 
Outputs:
	- None
Doing:
	+ Read the "subopt.db" file 
	+ Push in vector the name of sequences (headers) of subopt headers from the 
		"subopt.db" file
	+ Push in vector the sequences of subopt sequences from the "subopt.db" file
	+ Push in vector the structures of subopt structures from the "subopt.db" file
*/
void RepWeightedMotif::readSubopts(string& filePath) {
	vector<string> headers;
	vector<string> sequences;
	vector<string> structures;
	string line = "";
	string tempHeader = "";
	string tempSequence = "";
	string tempStructure = "";

	ifstream suboptsFile(filePath.c_str());

	while (getline(suboptsFile, line))
	{
		if (line[0] == '>')
		{
			if (tempSequence != "" and tempStructure != "")
			{
				headers.push_back(tempHeader);
				sequences.push_back(tempSequence);
				structures.push_back(tempStructure);
			}
			tempHeader = line;
			tempHeader.erase(0, 1);
			tempHeader.erase(remove(tempHeader.begin(), tempHeader.end(), ' '),
				tempHeader.end());

			tempHeader.erase(remove_if(tempHeader.begin(), tempHeader.end(),
				(int(*)(int))isspace), tempHeader.end());
		}
		else
		{
			if (line.length() > 1)
			{
				if (line[0] == '.' or line[0] == '(' or line[0] == ')')
				{
					tempStructure = line;
				}
				else
				{
					tempSequence = line;
				}
			}
		}
	}
	headers.push_back(tempHeader);
	sequences.push_back(tempSequence);
	structures.push_back(tempStructure);

	this->headers = headers;
	this->sequences = sequences;
	this->structures = structures;
}

/*********************************************************************************/
/*
Inputs:
	- map<string, float> nmotifsHeaderFirstRowOccurrencesAlignment: n-motifs header 
		and the first row of occurences of n-motifs from the "r_nmRepOfAS.csv" file
	- vector<string> nmotifsHeaderSubopt: n-motifs header from the "w_nmRepOfSS.csv" 
		file
Outputs:
	- None
Doing:
	+ Find the n-motifs from the "r_nmRepOfAS.csv" in the n-motifs from the "w_nmRepOfSS.csv"
		file
	+ Push in vector the positions of the shared n-motifs from the "w_nmRepOfSS.csv" file
	+ Push in vector the positions of the non-shared n-motifs from the "w_nmRepOfSS.csv" 
		file as -1
*/
void RepWeightedMotif::findPositionsNMotifs(map<string, float>& nmotifsHeaderFirstRowOccurrencesAlignment,
	vector<string>& nmotifsHeaderSubopt) {
	map<string, float>::iterator mapIt;
	vector<string>::iterator vectorIt;
	int positionNMotif;

	for (mapIt = nmotifsHeaderFirstRowOccurrencesAlignment.begin();
		mapIt != nmotifsHeaderFirstRowOccurrencesAlignment.end(); mapIt++)
	{
		vectorIt = find(nmotifsHeaderSubopt.begin(), nmotifsHeaderSubopt.end(), mapIt->first);

		if (vectorIt - nmotifsHeaderSubopt.begin() != nmotifsHeaderSubopt.size())
		{
			positionNMotif = vectorIt - nmotifsHeaderSubopt.begin();
		}
		else
		{
			positionNMotif = -1;
		}
		rowAllPositionsNMotifs.push_back(positionNMotif);
	}
}

/*********************************************************************************/
/*
Inputs:
	- ArrayXXf matOccurrencesNMotifAlignment: matrix of the occurrences of the n-motifs
		for each alignment structures from the "r_nmRepOfAS.csv" file
	- vector<float> rowAllWeightsSubopt:  vector of all the weights of the subopts from 
		the "weights.csv" file
	- vector<int> rowAllPositionsNMotifs: vector of the positions of the shared n-motifs 
		from the "w_nmRepOfSS.csv" file
Outputs:
	- None
Doing:
	+ Create a matrix of zeros m-alignment sequences by n-n-motifs (m-rows x n-cols)
	+ Calculate the weighted n-motif for each element of the matrix 
		(weighted = weight x occurrence)
	+ Fill the matrix with weighted n-motif values for all the m-alignment sequences
*/
void RepWeightedMotif::createMatWeightedNMotifsAlignment(Eigen::ArrayXXf& matOccurrencesNMotifAlignment,
	vector<float>& rowAllWeightsSubopt, vector<int>& rowAllPositionsNMotifs) {
	Eigen::ArrayXXf matWeightedNMotifsAlignment = Eigen::ArrayXXf::Zero(matOccurrencesNMotifAlignment.rows(),
		matOccurrencesNMotifAlignment.cols());
	int positionValue;
	float weightAtPosition;
	float occurrenceAtPosition;

	for (unsigned int rowNumber = 0; rowNumber < matWeightedNMotifsAlignment.rows(); rowNumber++)
	{
		for (unsigned int colNumber = 0; colNumber < matWeightedNMotifsAlignment.cols(); colNumber++)
		{
			positionValue = rowAllPositionsNMotifs[colNumber] - 1;

			if (positionValue < 0)
			{
				weightAtPosition = 0.f;
			}
			else
			{
				weightAtPosition = rowAllWeightsSubopt[positionValue];
			}
			occurrenceAtPosition = matOccurrencesNMotifAlignment(rowNumber, colNumber);

			matWeightedNMotifsAlignment(rowNumber, colNumber) = weightAtPosition * occurrenceAtPosition;
		}
	}
	this->matALLMotifWeigthed = matWeightedNMotifsAlignment;
}

/*********************************************************************************/
/*
Inputs:
	- ArrayXXf matWeightedNMotifsSubopt: matrix of the weighted n-motifs from the 
		"w_nmRepOfSS.csv" file
	- ArrayXXf matWeightedNMotifsAlignment: matrix of the weighted n-motifs from the 
		"w_nmRepOfAS.csv" file
	- vector<int> rowAllPositionsNMotifs: vector of the positions of the shared n-motifs 
		from the "w_nmRepOfSS.csv" file
Outputs:
	- None
Doing:
	+ Create a matrix of zeros m-alignment sequences by n-subopt sequences (rows x cols)
	+ Create a temp row vector containing a row of the m-alignment sequences weighted 
		n-motif values 
	+ Create a temp row vector containing a row of the n-subopt sequences weighted n-motif 
		values using the vector of the positions of the shared n-moifs 
	+ Calculate the euclidean distance between the 2 temp row vectors 
	+ Fill the matrix with the euclidean distances for all the m-alignment sequences
*/
void RepWeightedMotif::createMatDistancesAlignment(Eigen::ArrayXXf& matWeightedNMotifsSubopt,
	Eigen::ArrayXXf& matWeightedNMotifsAlignment, vector<int>& rowAllPositionsNMotifs) {
	int nCols = matWeightedNMotifsSubopt.rows();
	int nRows = matWeightedNMotifsAlignment.rows();
	int row;
	int positionValue;
	float weightedNMotifAtPosition;

	assert(nCols / 25 == nRows);

	matDistances = Eigen::ArrayXXf::Zero(nRows, nCols);

	for (int col = 0; col < nCols; col++)
	{
		row = col / 25;

		Eigen::RowVectorXf tempRowWeightedNMotifsAlignment = matWeightedNMotifsAlignment.row(row).array();
		Eigen::RowVectorXf tempRowWeightedNMotifsSubopt = Eigen::RowVectorXf::Zero(1, rowAllPositionsNMotifs.size());

		for (unsigned int i = 0; i < rowAllPositionsNMotifs.size(); i++)
		{
			positionValue = rowAllPositionsNMotifs[i] - 1;

			if (positionValue < 0)
			{
				weightedNMotifAtPosition = 0.f;
			}
			else
			{
				weightedNMotifAtPosition = matWeightedNMotifsSubopt(col, positionValue);
			}		
			tempRowWeightedNMotifsSubopt(0, i) = weightedNMotifAtPosition;
		}
		matDistances(row, col) = (tempRowWeightedNMotifsAlignment - tempRowWeightedNMotifsSubopt).matrix().lpNorm<2>();
	}
}

/*********************************************************************************/
/*
Inputs:
	- ArrayXXf matDistancesAlignmentSubopts: matrix of the euclidean distances between
		the alignment structures and the related subopt structures 
	- vector<string> suboptHeaders: vector of the name of the structures (headers) from 
		the "subopt.db" file
	- vector<string> suboptSequences: vector of the sequences from the "subopt.db" file
	- vector<string> suboptStructures: vector of the structures from the "subopt.db" 
		file
Outputs:
	- None
Doing:
	+ Select 25 non-zero columns for each row
	+ Find the index of the minimum of the row
	+ Push in vector the name of the sequence (header) of the closest subopt related
	+ Push in vector the sequence of the closest subopt related
	+ Push in vector the structure of the closest subopt related
*/
void RepWeightedMotif::computeClosestAlignmentToSubopts(Eigen::ArrayXXf& matDistancesAlignmentSubopts,
	vector<string>& suboptHeaders, vector<string>& suboptSequences, vector<string>& suboptStructures) {
	int nRows = matDistancesAlignmentSubopts.rows();
	int firstCol = 0;
	string tempSuboptHeader;
	string tempSuboptSequence;
	string tempSuboptStructure;

	for (int row = 0; row < nRows; row++)
	{
		Eigen::RowVectorXf tempRowDistances = matDistancesAlignmentSubopts.block(row, firstCol, 1, 25);

		Eigen::ArrayXXf::Index minDistanceIndex;
		tempRowDistances.array().minCoeff(&minDistanceIndex);

		minDistanceIndex = (row * 25) + minDistanceIndex;

		closestSuboptHeaders.push_back(suboptHeaders[minDistanceIndex]);
		closestSuboptSequences.push_back(suboptSequences[minDistanceIndex]);
		closestSuboptStructures.push_back(suboptStructures[minDistanceIndex]);

		firstCol += 25;
	}
}

/*********************************************************************************/
RepWeightedMotif::~RepWeightedMotif() {
	// TODO Auto-generated destructor stub
}


/*******************************************************************************************************/
Eigen::ArrayXXf RepWeightedMotif::filterMotifs(Eigen::RowVectorXi suboptidx,int computeConservationIndex){


    Eigen::ArrayXXf matS_nmotifs = Eigen::ArrayXXf::Zero(nmotifsForEachStructure.size(), allStructureFeatureForWeigthOfMotifs.size());
    std::map<string, float>::iterator it;
    int ind = 0;

    for (unsigned int i = 0; i < nmotifsForEachStructure.size(); i++)
    {
        for (std::map<string, float>::iterator it2 =nmotifsForEachStructure[i].begin();it2 != nmotifsForEachStructure[i].end(); ++it2)
        {
            it = allStructureFeatureForWeigthOfMotifs.find(it2->first);
            ind = distance(allStructureFeatureForWeigthOfMotifs.begin(), it);
            matS_nmotifs(i, ind) = it2->second;
        }
    }

    this->matALLMotif=matS_nmotifs;

    Eigen::ArrayXXf matS_weightdNmotifs = Eigen::ArrayXXf::Zero(nmotifsForEachStructure.size(), allStructureFeatureForWeigthOfMotifs.size());
    Eigen::RowVectorXf allWeight(matS_nmotifs.cols());

    if (computeConservationIndex==1)
    {

        vector<int> countNmotifVect;
        vector<int> uniquecountNmotifVect;

        Eigen::RowVectorXf vectFreq;
        Eigen::RowVectorXf vectFreqLog;
        Eigen::RowVectorXf currentCountCol;

        float currentSum;
        float entCount;
        float entLabel;
        float entY;
        float entXY;

        float currentWeight;
        Eigen::RowVectorXf allEntCount(matS_nmotifs.cols());
        Eigen::RowVectorXf allEntLabel(matS_nmotifs.cols());

        vector<int> currentLabelY;
        vector<int> currentUniqueLabelY;

        vector<int> currentLabelXY;
        vector<int> currentUniqueLabelXY;

        vector<int> currentLabelCol;
        vector<int> currentUniqueLabelCol;

        vector<int> currentCol;
        vector<int> currentUniqueCol;

        int tempTest;
        for (unsigned int j = 0; j < matS_nmotifs.cols(); j++)
        {
                for (unsigned int i = 0; i < matS_nmotifs.matrix().col(j).size(); i++)
                {
                    if (matS_nmotifs(i,j)>0)
                    {currentLabelCol.push_back(suboptidx(i)); }
                    currentCol.push_back(int(matS_nmotifs(i,j)));
                    currentLabelY.push_back(suboptidx(i));

                    tempTest= stoi(std::to_string(int(matS_nmotifs(i,j)))+std::to_string(suboptidx(i)));
                    currentLabelXY.push_back(tempTest);


                }

                //Entropy Count
                 std::sort (currentCol.begin(), currentCol.end());
                currentUniqueCol=currentCol;
                std::vector<int>::iterator it;
                it = std::unique (currentUniqueCol.begin(), currentUniqueCol.end());
                currentUniqueCol.resize( std::distance(currentUniqueCol.begin(),it));

                entCount=0;
                for (unsigned int i = 0; i < currentUniqueCol.size(); i++)
                {
                    currentSum=0;
                    for (unsigned int k = 0; k < currentCol.size(); k++)
                        {
                            if(currentCol[k]==currentUniqueCol[i])
                            {currentSum++;}
                        }
                    entCount=entCount+(currentSum/currentCol.size())*log(currentSum/currentCol.size());
                }
                entCount=entCount*-1;
                currentCol.clear();

                currentWeight=1/exp(entCount);//inverse diversity-->nmfold1

                allWeight(j)=currentWeight;
                matS_weightdNmotifs.col(j)=(matS_nmotifs.col(j))*currentWeight;

        }

    }

    else
    {

        Eigen::RowVectorXf allWeight=Eigen::RowVectorXf::Ones(matS_nmotifs.cols());
        matS_weightdNmotifs=matS_nmotifs;
    }

    this->allWeight=allWeight;

   return matS_weightdNmotifs;
}


/*******************************************************************************************************/
void RepWeightedMotif::computeRepresentative(Eigen::ArrayXXf matALLMotifWeigthed,vector<string> headers,vector<string> sequences,vector<string> structures, bool verbose, bool verbose2){

	if(verbose || verbose2){cout << "Compute Initial : " << endl << endl;}
    //Eigen::VectorXf allNorm;
    int n=matALLMotifWeigthed.rows();
	int i,j;

    //Compute Centroid
    Eigen::RowVectorXf centroid;
    centroid=matALLMotifWeigthed.colwise().mean();

	//compute distance between centroids and structures
    Eigen::RowVectorXf distToCentroid(n);
	Eigen::ArrayXXf distance = Eigen::ArrayXXf::Zero(n,n); //taille (25*n * 25*n)

	for (i=0; i<n;i++){
	    distToCentroid(i) = (centroid.array() - matALLMotifWeigthed.row(i).array()).matrix().lpNorm<2>();
    }
	
	//compute distance to everybody to everybody
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++) {
			distance(i,j) = (matALLMotifWeigthed.row(i).array() - matALLMotifWeigthed.row(j).array()).matrix().lpNorm<2>();
		}
	}

    //Get structure close to the centroid
    Eigen::ArrayXXf::Index minRow, minCol;
    distToCentroid.array().minCoeff(&minRow, &minCol);

    this->headerOfRepresentative=headers[minCol];
    this->sequenceOfRepresentative=sequences[minCol];
    this->structureOfRepresentative=structures[minCol];

	this->centroid=centroid.array();
	this->distanceToEveryBody=distance;
    this->distToCentroid=distToCentroid;

}


/*******************************************************************************************************/
void RepWeightedMotif::computeClosestStructuresToCentroid(Eigen::RowVectorXf distToCentroid, Eigen::RowVectorXi subOptIdx, vector<string> headers, vector<string> sequences, vector<string> structures, bool verbose, bool verbose2){
	//Stratégie 1 : plus proche de chaque du Centroid.
	//on calcule pour chaque groupe de SS la SS la plus proche du Centroid.
	if(verbose || verbose2){cout << "Compute Strategy 1 : " <<endl;}
	int k;
	int n = subOptIdx.size()/minNumberOfSubOpt;
	int m = minNumberOfSubOpt;
	Eigen::ArrayXXf StructuresFromSequence;
	Eigen::ArrayXXf::Index minRow, minCol;
	
	vector<string> h;
	vector<string> st;
	vector<string> seq;
	
	for(k=0; k<n; k++){
		//Pour tous les groupes de SS
		if(verbose2){cout << "\tn : " << n << ", "; cout << "m : " << m << ", "; cout << "k : " << k << ", "; cout << endl;}
		StructuresFromSequence = distToCentroid.block(0,m*k,1,m); //extraire les 25 valeurs depuis l'indice subOptIdx(k)
		StructuresFromSequence.array().minCoeff(&minRow, &minCol); //extraire l'indice du plus proche du centroid
		//on le stocke dans les représentants.
		h.push_back(headers[m*k + minCol]); 
		seq.push_back(sequences[m*k + minCol]);
		st.push_back(structures[m*k + minCol]);
	}
	if(verbose || verbose2){cout << endl;}
	
	this->headerOfRepresentatives1=h;
	this->structureOfRepresentatives1=st;
	this->sequenceOfRepresentatives1=seq;
}


/*******************************************************************************************************/
void RepWeightedMotif::computeClosestStructuresToGroup(Eigen::ArrayXXf matALLMotifWeigthed, Eigen::RowVectorXi subOptIdx, vector<string> headers, vector<string> sequences, vector<string> structures, bool verbose, bool verbose2){
	/**Stratégie 2 : plus proche du groupe des précédents résultats (stratégie itérative)
	//début : le groupe est formé uniquement du Centroid
	//Tant qu'on n'a pas n structures trouvées :
	//	1 - on trouve le plus proche du centroid du groupe
	//	2 - on l'ajoute au groupe
	//	3 - on élimine les autres résultats de la même séquence
	//	4 - on calcule le nouveau centroid du groupe
	//sortie : les n structures trouvées**/
	
	if(verbose || verbose2){cout << "Compute Strategy 2 : " << endl;}
	int i,j,k,glength;
	int n = subOptIdx.size()/minNumberOfSubOpt; //nombre d'itérations (fixé par le nombre de séquence en entrée) 
	int N = matALLMotifWeigthed.rows();
	int m = matALLMotifWeigthed.cols();
	int indexMin;
	Eigen::NoChange_t NoChange;
	
	vector<string> h;
	vector<string> st;
	vector<string> seq;
	
	Eigen::ArrayXXf groupe(1,m); 							//taille n+1 car on y inclus le centroid global initialement.
	Eigen::ArrayXXf::Index minRow, minCol;					//index des minimums
    Eigen::RowVectorXf tempCentroid;						//centroid initial
    tempCentroid=this->centroid;
	
	//association des indices des SS conservées pour le calcul dans "matALLMotifWeigthed"
	//Initialement, on conserve tous les index
	//Ainsi, "matALLMotifWeigthed.row(groupeMotifWeigthed(i))" donne la bonne structure de matALLMotifWeigthed associé à l'index "i" dans groupeMotifWeigthed
	Eigen::RowVectorXi groupeMotifWeigthed = Eigen::RowVectorXi::LinSpaced(N,0,N-1);
	Eigen::RowVectorXi groupeMotifWeigthedTemp = Eigen::RowVectorXi::LinSpaced(N,0,N-1); //copie temporaire
	
	if(verbose2){cout << "groupeMotifWeigthed: " << groupeMotifWeigthed.rows() << " rows, " << groupeMotifWeigthed.cols() << " cols" << endl;
				 cout << "groupe: " << groupe.rows() << " rows, " << groupe.cols() << " cols" << endl;}
	
	//Initialisation : le groupe de comparaison n'est composée que du Centroid global
	groupe.row(0) = tempCentroid.array();
	Eigen::RowVectorXf groupCentroid = tempCentroid;
	
	if(verbose || verbose2){cout << endl;}
	//Itérations
	for(k=1; k<n+1; k++){
		if(verbose || verbose2){cout << "iteration " << k << " : " <<endl;}
		
		//step 1
		if(verbose || verbose2){cout << "Step 1... ";}
		glength = groupeMotifWeigthed.cols(); 		//nombre de structures non traitées restantes
		Eigen::RowVectorXf distToCentroid(glength);	//distance au groupe déjà construit (initialement le centroid global seul)
		if(verbose || verbose2){cout << "glength=" << glength << ". compute distToCentroid... ";}
		for (i=0; i<glength;i++){
			distToCentroid(i) = (groupCentroid.array() - matALLMotifWeigthed.row(groupeMotifWeigthed(i)).array()).matrix().lpNorm<2>();
		}
		if(verbose || verbose2){cout << "done... extract representative : ";}
		distToCentroid.array().minCoeff(&minRow, &minCol);
		indexMin = groupeMotifWeigthed(minCol); //index réel de la structure dans "matALLMotifWeigthed"
		if(verbose || verbose2){cout << "Structure " << minCol << " ==> Distance = " << distToCentroid(minCol) << endl;}
		
		//step 2
		if(verbose || verbose2){cout << "Step 2... pushing... ";}
		//On push la structure trouvée dans les optimales
		h.push_back(headers[indexMin]); 
		seq.push_back(sequences[indexMin]);
		st.push_back(structures[indexMin]);
		if(verbose || verbose2){cout << "done... ";
								cout << "resize groupe... ";}
		//On l'ajoute également dans le nouveau groupe auquel se rapprocher						
		groupe.resize(k+1,NoChange);
		groupe.row(k) = matALLMotifWeigthed.row(indexMin).array();
		if(verbose || verbose2){cout << "OK : "  << endl;}
		if(verbose2){cout << "\t- groupe: " << groupe.rows() << " rows, " << groupe.cols() << " cols" << endl;}
		
		//step 3
		//On retire les 25 structures de la séquence représentée de la liste de celles restantes
		if(verbose || verbose2){cout << "Step 3... "; cout << "resize groupeMotifWeigthed... ";}
		groupeMotifWeigthed.resize(glength-25); //redimensionnement du vecteur des structures de travail courantes
		if(verbose || verbose2){cout << "OK : ";}
		if(verbose2){cout << endl << "\t- groupeMotifWeigthed: " << groupeMotifWeigthed.rows() << " rows, " << groupeMotifWeigthed.cols() << " cols" << endl;
					 cout << "eliminate result structure... " << endl;
					 cout << "\t--> Index of finding Structure : " << indexMin << ". For Sequence " << subOptIdx(indexMin) << endl;}
		m=groupeMotifWeigthedTemp.cols();		//taille du vecteur précédent sur lequel itérer
		j=0;									//index courant dans le nouveau vecteur (pour éviter le dépassement d'indice)
		for(i=0;i<m;i++){
			if(verbose2){cout << "\t- Case of Structure " << groupeMotifWeigthedTemp(i) << " of the Sequence " << subOptIdx(groupeMotifWeigthedTemp(i));}
			if(j < glength-25 && subOptIdx(groupeMotifWeigthedTemp(i)) != subOptIdx(indexMin)){
				if(verbose2){cout << "\t -- False" << endl;}
				groupeMotifWeigthed(j) = groupeMotifWeigthedTemp(i); //Structure conservée pour l'itération suivante
				j++;												 //Case suivante de "groupeMotifWeigthed"
			}
			if(verbose2 && subOptIdx(groupeMotifWeigthedTemp(i)) == subOptIdx(indexMin)){cout << "\t -- True" << endl;}
		}
		//Construction de la copie du vecteur d'indice
		if(verbose || verbose2){cout << "done... resize groupeMotifWeigthedTemp... ";}
		groupeMotifWeigthedTemp.resize(glength-25); //redimensionnement du vecteur
		if(verbose || verbose2){cout << "OK. " << endl;}
		if(verbose2){cout << "\t- groupeMotifWeigthedTemp: " << groupeMotifWeigthedTemp.rows() << " rows, " << groupeMotifWeigthedTemp.cols() << " cols" << endl;}
		m=groupeMotifWeigthedTemp.cols();
		for(i=0;i<m;i++){groupeMotifWeigthedTemp(i) = groupeMotifWeigthed(i);} //copie du vecteur
		
		//step 4
		if(verbose || verbose2){cout << "Step 4... compute new centroid... ";}
		groupCentroid = groupe.colwise().mean(); //nouveau centroid du groupe
		if(verbose || verbose2){cout << "done." << endl << endl;}
	}
	
	//écriture des résultats
	this->groupCentroid=groupCentroid;
	this->headerOfRepresentatives2=h;
	this->structureOfRepresentatives2=st;
	this->sequenceOfRepresentatives2=seq;
}