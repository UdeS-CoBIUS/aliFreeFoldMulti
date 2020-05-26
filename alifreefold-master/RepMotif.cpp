/*
 * RepMotif.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: Séhi
 */

#include "RepMotif.h"

RepNmotif::RepNmotif() {
	// TODO Auto-generated constructor stub
}

RepNmotif::RepNmotif(vector<string> allInputFilePaths, int maxLevelNmotifs, bool isFullTask) {
	// TODO Auto-generated constructor stub
    this->maxLevelNmotifs = maxLevelNmotifs;
	this->isCircular = false;

	if (isFullTask)
	{
		readFileandParseSequenceStructure(allInputFilePaths[0], headers, structures, 
			sequences, isCircular, minNumberOfSubOpt);
    	if (headers.size()<2)
    	{
    		throw invalid_argument("Not enougth structures to compare. Enter more than one RNA secondary structure in your input file.");
    	}
	}
	else
	{
		// Find in the -i parameter the file path of the "alignement.db" file 
		vector<string>::iterator it;
		it = find_if(allInputFilePaths.begin(), allInputFilePaths.end(), isAlignmentFile);
		string alignmentFilePath = *it;

		readFileandParseAlignmentStructure(alignmentFilePath, headers, structures,
			sequences, isCircular);
	}
}

RepNmotif::~RepNmotif() {
	// TODO Auto-generated destructor stub
}

/******************************************************************************/
void RepNmotif::readFileandParseSequenceStructure(string filePath,
						  vector<string>& headers, vector<string>& structures,
						  vector<string>& sequences, bool& isCircular,int & minNumberOfSubOpt) {
  string line, tempLine, tempLineHeader, tempSequence, tempStructure;
  ifstream myfile(filePath.c_str());
  int nb_object = 0;
  
  bool first = true;
  tempLine = "";
  
  
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.uniq_ML = 1;
  
  while (getline(myfile, line)) {
    if (line[0] == '>') {
      
      if (tempLine != "") {
	
	vrna_fold_compound_t *vc = vrna_fold_compound(tempLine.c_str(), NULL, 
		VRNA_OPTION_DEFAULT);
	/* call zuker subopt function */
	vrna_subopt_solution_t * vcs = vrna_subopt_zuker(vc);
	/* call subopt function */
	//vrna_subopt_solution_t * vcs = vrna_subopt(vc,1000,1,NULL);
	
	for (int j=0;j<minNumberOfSubOpt;j++)
	  {
	    tempStructure=vcs[j].structure;
	    
	    //cout<<tempStructure<<" : "<<vcs[j].energy<<endl;
	    tempSequence=tempLine;
	    structures.push_back(tempStructure);
	    
	    sequences.push_back(tempSequence);
	    parseSequenceStructure(tempStructure,tempSequence);
	  }
	vrna_fold_compound_free(vc);
	
      }
      tempLineHeader = line;
      
      tempLineHeader.erase(0, 1);
      tempLineHeader.erase(remove(tempLineHeader.begin(), tempLineHeader.end(), ' '), 
      	tempLineHeader.end());
      if ((tempLineHeader.at(0) == 'c')
	  && (tempLineHeader.at(1) == '_')) {
	isCircular = true;
	//tempLineHeader.erase(0, 2);
      }else{isCircular = false; }
      tempLineHeader.erase(
			   std::remove_if(tempLineHeader.begin(), tempLineHeader.end(),
					  (int (*)(int))isspace), tempLineHeader.end());
      
      for (int j=0;j<minNumberOfSubOpt;j++)
	{
	  headers.push_back(tempLineHeader);
	}
      
      // cout<<tempLineHeader<<endl;
      if (first == false) {
	nb_object = nb_object + 1;
      }
      
      first = false;
      tempLine = "";
    } else {
      if (line.length() > 1) {
	tempLine += line;
      }
    }
  }
  
  vrna_fold_compound_t *vc = vrna_fold_compound(tempLine.c_str(), NULL, VRNA_OPTION_DEFAULT);
  /* call subopt function */
  vrna_subopt_solution_t * vcs = vrna_subopt_zuker(vc);
  
  for (int j=0;j<minNumberOfSubOpt;j++)
    {
      tempStructure=vcs[j].structure;
      tempSequence=tempLine;
      structures.push_back(tempStructure);
      sequences.push_back(tempSequence);
      parseSequenceStructure(tempStructure,tempSequence);
    }
  vrna_fold_compound_free(vc);
  
  
  free(vcs);
  myfile.close();
  
}

/******************************************************************************/
/*
Inputs:
	- string filePath: File path in the -i parameter
Outputs:
	- None
Doing:
	+ Return true if the file path contains the string "alignement.db"
	+ Return false else
*/
bool RepNmotif::isAlignmentFile(string filePath) {
	if (filePath.find("alignment.db") != string::npos)
	{
		return true;
	}
	else
	{
		return false;
	}
}


/******************************************************************************/
void RepNmotif::readFileandParseAlignmentStructure(string filePath, vector<string>& headers, 
	vector<string>& structures, vector<string>& sequences, bool& isCircular) {
	ifstream myfile(filePath.c_str());

    string tempHead = "";
    string tempSeq = "";
    string tempStruct = "";
    string line = "";
    string tempLineHeader = "";
    
    while (getline(myfile, line)) {
    	if (line[0] == '>') {
    		if (tempSeq != "" and tempStruct != "") 
    		{
	  			structures.push_back(tempStruct);  
	  			sequences.push_back(tempSeq);
	  			parseSequenceStructure(tempStruct,tempSeq);
	  			headers.push_back(tempHead);
	  			tempSeq = "";
	  			tempStruct = ""; 
	  		}

	  		string tempLineHeader = line;
			tempLineHeader.erase(0, 1);
			tempLineHeader.erase(remove(tempLineHeader.begin(), tempLineHeader.end(), ' '), 
				tempLineHeader.end());

			if ((tempLineHeader.at(0) == 'c') && (tempLineHeader.at(1) == '_')) 
			{
	  			isCircular = true;
	  			//tempLineHeader.erase(0, 2)
	  		}
	  		else
	  		{
	  			isCircular = false; 
	  		}

	  		tempLineHeader.erase(
			     std::remove_if(tempLineHeader.begin(), tempLineHeader.end(),
					    (int (*)(int))isspace), tempLineHeader.end()); 
			tempHead =tempLineHeader;
		}
		else
		{
			if (line.length() > 1) 
			{
				if(line[0] == '.' or line[0] == '(' or line[0] == ')')
				{
					tempStruct += line;
				}
				else
				{
					tempSeq += line;
				}
			}
		}
	}

    structures.push_back(tempStruct);  
    sequences.push_back(tempSeq);
    headers.push_back(tempHead);
    parseSequenceStructure(tempStruct,tempSeq);
}

/******************************************************************************/
void RepNmotif::parseSequenceStructure(string currentStructure,string currentSequence) {
	graph currentGrss;
	vector<RepNmotif::Motif> currentMotifVect;
    map<string, float> currentStructureFeature;
	map<string, vector<int>> currentStructureFeatureWithPosNuc;

	bool onlyG4motifinStruct = false;

        //Parse Structure
        /*special dot and bracket --> vect of motif*/
        specialCarac2SpecialMotifVect(currentStructure ,
                onlyG4motifinStruct, currentMotifVect);

        if (onlyG4motifinStruct == false)
        /*convert structure to shapiro and compute statistics on rna secondary structure*/
        {
            /*dot and bracket --> rss graph*/
            dbn2Grss(currentStructure, currentGrss);

            /*rss graph --> vect of motif*/
            grss2motifVect(currentStructure, currentGrss, currentMotifVect);
        }
        for (int i = 0; i <= maxLevelNmotifs; i++)
        {intersectMotifs(i, currentStructure, currentMotifVect,
                    currentStructureFeature, currentStructureFeatureWithPosNuc,nmotifsAllStructure);
        }


	nmotifsForEachStructure.push_back(currentStructureFeature);

	nmotifsForEachStructureWithPosNucOfnmotifs.push_back(
			currentStructureFeatureWithPosNuc);

}

/*********************************************************************************/
void RepNmotif::specialCarac2SpecialMotifVect(string& currentStructure, bool& onlyG4motifinStruct,
		vector<RepNmotif::Motif>& motifVect) {

	std::ostringstream motifFeatTemp;
	int countTemp;
	multimap<char, int> mymm;
	multimap<char, int>::iterator it;

	/*current motif and parameters*/
	Motif motifCurrent;
	initCurrentMotif(motifCurrent);

	string specialCaracPseudoknotsLeft = "[{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	string specialCaracPseudoknotsRight = "]}>abcdefghijklmnopqrstuvwxyz";
	//add here special caracters
	string allspecialCarac =
			"[]{}<>ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+";

	/*Pseudoknots*/
	if (currentStructure.find_first_of(specialCaracPseudoknotsLeft)
			!= string::npos) {
		//voir algorithm de knuth et pratt
		for (unsigned int i = 0; i < currentStructure.size(); i++) {
			mymm.insert(std::make_pair(currentStructure.at(i), i));
		}

		//std::multiset<int>::key_compare mycomp = mymm.key_comp();
		for (unsigned int i = 0; i < specialCaracPseudoknotsLeft.size(); i++) {
			countTemp = mymm.count(specialCaracPseudoknotsLeft.at(i));
			if (countTemp != 0) {
				for (it =mymm.equal_range(specialCaracPseudoknotsLeft.at(i)).first;it!= mymm.equal_range(specialCaracPseudoknotsLeft.at(i)).second;++it) {
					motifCurrent.posNuc.push_back((*it).second);
				}
				for (it =mymm.equal_range(specialCaracPseudoknotsRight.at(i)).first;it!= mymm.equal_range(specialCaracPseudoknotsRight.at(i)).second;++it)
					motifCurrent.posNuc.push_back((*it).second);

				motifFeatTemp << "P";
				motifCurrent.motifName = motifFeatTemp.str();
				motifFeatTemp << "_" << countTemp;//*modif
				motifCurrent.motifNameWithFeatLevel1 = motifFeatTemp.str();
				motifFeatTemp.str("");

				motifVect.push_back(motifCurrent);
				initCurrentMotif(motifCurrent);
			}
		}
	}

	/*Gquadruplexes*/
	vector<int> posG4;
	vector<int> posNuc;
	int countG ;
	int nbloop = 1;
	countG =0;

	unsigned int countCar = 0;

	std::size_t firstPlus = currentStructure.find_first_of("+");
	//parse g4
	if (firstPlus != string::npos) {
		for (unsigned int i = 0; i < currentStructure.size(); i++) {
			if (currentStructure.at(i) == '+') {
				posG4.push_back(i);
			}
			if ((currentStructure.at(i) == '+')	|| (currentStructure.at(i) == '.')) {
				countCar++;
			}
		}
		if (countCar == currentStructure.size()) {
			onlyG4motifinStruct = true;
		}

		std::size_t lastPlus = currentStructure.find_last_of("+");

		if (onlyG4motifinStruct == true) {
			if(firstPlus>0)
			{
				motifFeatTemp << "E5";
				motifCurrent.motifName = motifFeatTemp.str();
				motifFeatTemp << "_" << firstPlus;
				motifCurrent.motifNameWithFeatLevel1 = motifFeatTemp.str();
				motifFeatTemp.str("");
				for (unsigned int i = 0; i <= firstPlus; i++) {
					posNuc.push_back(i);
                }
				motifCurrent.posNuc = posNuc;
				motifVect.push_back(motifCurrent);
				initCurrentMotif(motifCurrent);
				posNuc.clear();
            }else{
                motifFeatTemp << "E5";
                motifCurrent.motifName = motifFeatTemp.str();
                motifFeatTemp << "_" << 0;
                motifCurrent.motifNameWithFeatLevel1 = motifFeatTemp.str();
                motifFeatTemp.str("");
                posNuc.push_back(0);
                motifCurrent.posNuc = posNuc;
                motifVect.push_back(motifCurrent);
                initCurrentMotif(motifCurrent);
                posNuc.clear();
                 }

			int temp=currentStructure.length() - lastPlus - 1;
			if(temp>0)
			{
				motifFeatTemp << "E3";
				motifCurrent.motifName = motifFeatTemp.str();
				motifFeatTemp << "_" << currentStructure.length() - lastPlus - 1;
				motifCurrent.motifNameWithFeatLevel1 = motifFeatTemp.str();
				motifFeatTemp.str("");
				for (unsigned int i = lastPlus; i < currentStructure.length();
						i++) {
					posNuc.push_back(i);
				}
                motifCurrent.posNuc = posNuc;
				motifVect.push_back(motifCurrent);
				initCurrentMotif(motifCurrent);
				posNuc.clear();
            }else{
                motifFeatTemp << "E3";
                motifCurrent.motifName = motifFeatTemp.str();
                motifFeatTemp << "_" << 0;
                motifCurrent.motifNameWithFeatLevel1 = motifFeatTemp.str();
                motifFeatTemp.str("");
                posNuc.push_back(currentStructure.length()- 1);
                motifCurrent.posNuc = posNuc;
                motifVect.push_back(motifCurrent);
                initCurrentMotif(motifCurrent);
                posNuc.clear();
                }
		}

		//unsigned int j;
		for (unsigned int i = 0; i < posG4.size(); i++) {
			if (i + 1 < posG4.size()) {
				if (posG4[i + 1] - 1 == posG4[i]) {
					posNuc.push_back(posG4[i]);
					countG++;
				} else if (nbloop <= 3) {
					posNuc.push_back(posG4[i]);
					nbloop++;
					countG++;
				} else if (nbloop == 4) {
					countG++;
					posNuc.push_back(posG4[i]);
                    motifFeatTemp << "G4";
					motifCurrent.motifName = motifFeatTemp.str();
					motifFeatTemp << "_" << countG / 4;
					motifCurrent.motifNameWithFeatLevel1 = motifFeatTemp.str();
                    motifFeatTemp.str("");
					motifCurrent.posNuc = posNuc;
					motifVect.push_back(motifCurrent);
					initCurrentMotif(motifCurrent);
					posNuc.clear();
					nbloop = 1;
					countG = 0;
				}
			} else {
				countG++;
				posNuc.push_back(posG4[i]);
                motifFeatTemp << "G4";
				motifCurrent.motifName = motifFeatTemp.str();
				motifFeatTemp << "_" << countG / 4;
				motifCurrent.motifNameWithFeatLevel1 = motifFeatTemp.str();
                motifFeatTemp.str("");
				motifCurrent.posNuc = posNuc;
				motifVect.push_back(motifCurrent);
				initCurrentMotif(motifCurrent);
				posNuc.clear();
				nbloop = 1;
				countG = 0;
			}
		}
	}
	/*replace all special caracters by . except '(', ')', '.'*/
	for (unsigned int k = 0; k < allspecialCarac.size(); k++) {
		replace(currentStructure.begin(), currentStructure.end(),
				allspecialCarac[k], '.');
	}

}

/*********************************************************************************/
void RepNmotif::dbn2Grss(string currentStructure,
		RepNmotif::graph& currentGrss) {

	unsigned int i, j;
	vector<bool> rightPaired(currentStructure.size());

	// inspired from dot2ct source code of RNAstructure software
	// Parse the dots and brackets.
	for (i = 0; i < currentStructure.size(); i++) {
		if (currentStructure.at(i) == '(') {
			rightPaired[i] = true;
		} else {
			rightPaired[i] = false;
		}
		boost::add_vertex(currentGrss);
	}

	for (i = 0; i < currentStructure.length(); i++) {
		if (currentStructure.at(i) == ')') {
			// Find a 3' paired nuc, then find the 5' partner.
			j = i - 1;
			while ((rightPaired[j] == false) && j >= 0) {
				j--;
			}
			rightPaired[j] = false;
			boost::add_edge(i, j, currentGrss);
		}
	}
}

/******************************************************************************/
void RepNmotif::grss2motifVect(string currentStructure,
        graph currentGrss, vector<RepNmotif::Motif>& motifVect)
{

	pair<vertex_iterator_t, vertex_iterator_t> itVertex = boost::vertices(
			currentGrss);
	vertex_iterator_t itVertex_previouscurrent, itVertex_aftercurrent;
	adjacency_iterator_t neighbourIt, neighbourEnd, neighbourItLoop,
			neighbourEndLoop;
	adjacency_iterator_t neighbourIt_itVertex_previouscurrent,
			neighbourEnd_itVertex_previouscurrent;
	adjacency_iterator_t neighbourIt_itVertex_aftercurrent,
			neighbourEnd_itVertex_aftercurrent;
	ostringstream motifFeatTemp;


    map<int, vector<int>> mapPosVSIndexMotifVect;//


	//circular RNA
	bool ispresentE5, ispresentE3, ispresentMwithE5E3, caseStemStemLinkedAtE5E3;
	ispresentE5 = ispresentE3 = ispresentMwithE5E3 = caseStemStemLinkedAtE5E3 =
			false;

	unsigned int nbUnpairedNucCircE5, nbUnpairedNucCircE3, nbUnpairedNucCircM,
			nbUnpairedNucCircE5E3, degreeTemp, lastStemInd;
	vector<int> posNucCirc;
	nbUnpairedNucCircE5 = nbUnpairedNucCircE3 = nbUnpairedNucCircM =
			nbUnpairedNucCircE5E3 = degreeTemp = lastStemInd = 0;
	// circular RNA

	int neighbourItValue_current, neighbourItValue_previouscurrent;
	unsigned int neighbourItValue_aftercurrent;
	bool inloop, newStem = true;
	inloop = true;
	newStem = false;

	vector<bool> visitedPosVect(currentStructure.size());
    vector<bool> visitedPosVectSpecialCaseMLoopWithoutNuc(currentStructure.size());

	/*current motif and parameters*/
	Motif motifCurrent;
	unsigned int nbUnpairedNuc, degree, nbBp, nbNucSidea, nbNucSideb, nbQuartet,
			nbNucL1, nbNucL2, nbNucL3;
	bool symLoop, parseE5, parseE3;
	parseE5 = parseE3 = false;
    vector<int> posNuc,posNucLoop;
    list<int> posNuc1,posNucLoop1;
	initCurrentMotif(motifCurrent);
	initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc, degree, nbBp,
			nbNucSidea, nbNucSideb, symLoop, nbQuartet, nbNucL1, nbNucL2,
			nbNucL3, posNuc);

	for (unsigned i = 0; i < currentStructure.size(); i++) {
		visitedPosVect[i] = false;
        visitedPosVectSpecialCaseMLoopWithoutNuc[i]=false;
	}


	/*Motif Extraction*/
    for (vertex_iterator_t itVertex_current = itVertex.first; itVertex_current != itVertex.second; ++itVertex_current)
    {
        if (visitedPosVect[*itVertex_current] == false)
        {
			visitedPosVect[*itVertex_current] = true;
            boost::tie(neighbourIt, neighbourEnd) = boost::adjacent_vertices(*itVertex_current, currentGrss);

            //Parse 5' and 3' ends
            if ((itVertex_current == itVertex.first))
            {
				//Loop: Free Nuc E5'
				vertex_iterator_t itE5 = itVertex_current;
				while ((*itE5 < *itVertex.second) && (degree < 1)) {
                    boost::tie(neighbourItLoop, neighbourEndLoop) =	boost::adjacent_vertices(*itE5, currentGrss);
					if (neighbourItLoop != neighbourEndLoop) {
						//posNuc.push_back(*itE5);
						degree++;
					} else {
						posNuc.push_back(*itE5);
						visitedPosVect[*itE5] = true;
						nbUnpairedNuc++;
						itE5++;
					}
				}
                posNuc.push_back(*itE5);

                //if ((nbUnpairedNuc >= 0))	//Eloops:freeEnd 5'
                //{
                    if (isCircular == false)
                    {
						motifFeatTemp << "E5";
						motifCurrent.motifName = motifFeatTemp.str();
						motifFeatTemp << "_" << nbUnpairedNuc;
                        motifCurrent.motifNameWithFeatLevel1 =	motifFeatTemp.str();
						motifFeatTemp.str("");
                        if(nbUnpairedNuc>0){motifCurrent.posNuc = posNuc;}
                        else{posNuc.push_back(0);
                            motifCurrent.posNuc=posNuc; }
                        motifVect.push_back(motifCurrent);
                    }
                    else if (nbUnpairedNuc>0)
                    {
						ispresentE5 = true;
						nbUnpairedNucCircE5 = nbUnpairedNuc;
						posNucCirc.insert(posNucCirc.end(), posNuc.begin(),
								posNuc.end());
					}

					initCurrentMotif(motifCurrent);
                    parseE5  = true;
                 //}
                 initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc, degree,
                        nbBp, nbNucSidea, nbNucSideb, symLoop, nbQuartet,
                        nbNucL1, nbNucL2, nbNucL3, posNuc);

				//Loop: Free Nuc E3'
				vertex_iterator_t itE3 = itVertex.second;
                --itE3;
                while ((*itE3 > *itVertex.first) && (degree < 1))
                {
                    boost::tie(neighbourItLoop, neighbourEndLoop) = boost::adjacent_vertices(*itE3, currentGrss);
                    if (neighbourItLoop != neighbourEndLoop)
                    {
						//posNuc.push_back(*itE3);
						degree++;
                    }
                    else
                    {
                        posNuc1.push_front(*itE3);
						visitedPosVect[*itE3] = true;
						nbUnpairedNuc++;
						itE3--;
					}
                }
                posNuc1.push_front(*itE3);
                posNuc.insert(posNuc.end(), posNuc1.begin(),posNuc1.end());
                posNuc1.clear();

                //if (nbUnpairedNuc >= 0)	//Eloops:freeEnd 3'
                //{
                    if (isCircular == false)
                    {
						motifFeatTemp << "E3";
						motifCurrent.motifName = motifFeatTemp.str();
						motifFeatTemp << "_" << nbUnpairedNuc;
						motifCurrent.motifNameWithFeatLevel1 =
								motifFeatTemp.str();
						motifFeatTemp.str("");
                        if(nbUnpairedNuc>0){motifCurrent.posNuc = posNuc;}
                        else {int temp=currentStructure.length()- 1;
                            posNuc.push_back(temp);
                            motifCurrent.posNuc=posNuc; }

                        motifVect.push_back(motifCurrent);
                    }
                    else if(nbUnpairedNuc>0){
						ispresentE3 = true;
						nbUnpairedNucCircE3 = nbUnpairedNuc;
						posNucCirc.insert(posNucCirc.end(), posNuc.begin(),
								posNuc.end());
					}

					initCurrentMotif(motifCurrent);
					parseE3 = true;
                //}
				initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc, degree,
						nbBp, nbNucSidea, nbNucSideb, symLoop, nbQuartet,
						nbNucL1, nbNucL2, nbNucL3, posNuc);

				//Component
				//itE5--;itE3++;
				if (parseE5) {
					itE5--;
				}
				if (parseE3) {
					itE3++;
				}

				vertex_iterator_t itComponent = itE5;
                itComponent++;
                while ((*itComponent < *itE3)) {
                    boost::tie(neighbourItLoop, neighbourEndLoop) = boost::adjacent_vertices(*itComponent, currentGrss);
                    if (neighbourItLoop != neighbourEndLoop)
                    {
						posNuc.push_back(*itComponent);
                        posNuc.push_back(*neighbourItLoop);
                        itComponent = advanceVertex_iterator_t2adjacency_iterator_t(itComponent, neighbourItLoop);
						degree++;
                    }
                    else
                    {
						posNuc.push_back(*itComponent);
						visitedPosVect[*itComponent] = true;
						nbUnpairedNuc++;
					}
					itComponent++;
                }

				//Component
				if (degree > 1) {

                    if (isCircular == false)
                    {
						motifFeatTemp << "M";
						motifCurrent.motifName = motifFeatTemp.str();
						motifFeatTemp << "_" << nbUnpairedNuc;
                        motifCurrent.motifNameWithFeatLevel1 =	motifFeatTemp.str();
                        motifFeatTemp.str("");
						motifCurrent.posNuc = posNuc;
						motifVect.push_back(motifCurrent);
                    }
                    else
                    {
						ispresentMwithE5E3 = true;
						nbUnpairedNucCircM = nbUnpairedNuc;
						degreeTemp = degree;
						posNucCirc.insert(posNucCirc.end(), posNuc.begin(),
								posNuc.end());
					}

					initCurrentMotif(motifCurrent);
					initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
							degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
							nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
				}
                else {
                    initCurrentMotif(motifCurrent);
                    initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
                            degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
                            nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
                    posNuc.push_back(*itE5);
                    itE3--;
                    posNuc1.push_front(*itE3);
                }

				//Circular
                if (isCircular == true)
                {
					nbUnpairedNucCircE5E3 = nbUnpairedNucCircE5
							+ nbUnpairedNucCircE3;
                    sort (posNucCirc.begin(), posNucCirc.end());
					if (ispresentMwithE5E3 == false) {
						if ((ispresentE5 == true) || (ispresentE3 == true)) {
							motifFeatTemp << "H";
							motifCurrent.motifName = motifFeatTemp.str();
							motifFeatTemp << "_" << nbUnpairedNucCircE5E3;
							motifCurrent.motifNameWithFeatLevel1 =
									motifFeatTemp.str();
                            motifFeatTemp.str("");
							motifCurrent.posNuc = posNucCirc;
                            motifVect.push_back(motifCurrent);
							initCurrentMotif(motifCurrent);
                            initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
                                    degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
                                    nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
						}
					} else {
                        //if ((degreeTemp == 2) && (nbUnpairedNucCircM == 0)) {
                        if ((degreeTemp == 2) && (nbUnpairedNucCircM == 0)) {
                            caseStemStemLinkedAtE5E3 = true;
						} else if (degreeTemp == 2) {
							if ((nbUnpairedNucCircM == 0)
									|| (nbUnpairedNucCircE5E3 == 0)) {
								motifFeatTemp << "B";
								motifCurrent.motifName = motifFeatTemp.str();
								motifFeatTemp << "_"
										<< nbUnpairedNucCircE5E3+nbUnpairedNucCircM;
								motifCurrent.motifNameWithFeatLevel1 =
										motifFeatTemp.str();
                                motifFeatTemp.str("");
								motifCurrent.posNuc = posNucCirc;
								motifVect.push_back(motifCurrent);
								initCurrentMotif(motifCurrent);
                                initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
                                        degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
                                        nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
							} else {

								motifFeatTemp << "I";
								motifCurrent.motifName = motifFeatTemp.str();
								if (nbUnpairedNucCircE5E3
										== nbUnpairedNucCircM) {
									motifFeatTemp << "_" << true;
								} else {
									motifFeatTemp << "_" << false;
								}
								motifCurrent.motifNameWithFeatLevel1 =
										motifFeatTemp.str();
                                motifFeatTemp.str("");
                                motifCurrent.posNuc = posNucCirc;
								motifVect.push_back(motifCurrent);
								initCurrentMotif(motifCurrent);
                                initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
                                        degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
                                        nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
							}

						} else if (degreeTemp > 2) {
							motifFeatTemp << "M";
							motifCurrent.motifName = motifFeatTemp.str();
                            motifFeatTemp << "_" <<  nbUnpairedNucCircM;
                            motifCurrent.motifNameWithFeatLevel1 = motifFeatTemp.str();
                            motifFeatTemp.str("");
                            motifCurrent.posNuc = posNucCirc;
                            motifVect.push_back(motifCurrent);

                            //cout<<motifCurrent.motifNameWithFeatLevel1<<endl;
                            //for (int l=0; l<posNucCirc.size();l++)
                            //{cout<<";"<<posNucCirc[l]<<";";} cout<<endl;

                            initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
                                    degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
                                    nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
						}

					}
				}

				if (neighbourIt != neighbourEnd) {
					inloop = false;
				} else {
					continue;
				}
			}

            if (itVertex_current != itVertex.first)
            {
				itVertex_previouscurrent = itVertex_current;
				itVertex_previouscurrent--;
				itVertex_aftercurrent = itVertex_current;
				itVertex_aftercurrent++;
                boost::tie(neighbourIt_itVertex_previouscurrent,neighbourEnd_itVertex_previouscurrent) =boost::adjacent_vertices(*itVertex_previouscurrent,	currentGrss);
                boost::tie(neighbourIt_itVertex_aftercurrent,neighbourEnd_itVertex_aftercurrent) =boost::adjacent_vertices(*itVertex_aftercurrent,currentGrss);
                if (neighbourIt != neighbourEnd)
                {
					neighbourItValue_current = *neighbourIt;
				}

                if ((neighbourIt_itVertex_previouscurrent!= neighbourEnd_itVertex_previouscurrent)&& (neighbourIt != neighbourEnd))	//neighbor
                {
                    neighbourItValue_previouscurrent =*neighbourIt_itVertex_previouscurrent;
                    if (neighbourItValue_current + 1 == neighbourItValue_previouscurrent)
                    {
                        inloop = false;

                    } else
                    {
						newStem = true;
                        inloop = false;
                    }
                }
                else if (neighbourIt != neighbourEnd)
                {
					inloop = false;
                }
                else
                {
					inloop = true;
				}

                if ((neighbourIt_itVertex_aftercurrent!= neighbourEnd_itVertex_aftercurrent) && (neighbourIt != neighbourEnd))
                    {
                        neighbourItValue_aftercurrent = *neighbourIt_itVertex_aftercurrent;
                        if (*itVertex_current == neighbourItValue_aftercurrent)
                        {
                            newStem = true;
                            inloop = false;
                            nbBp++;
                        }
                    }
            }
            else
            {
                if (neighbourIt != neighbourEnd)
                {
					inloop = false;
                }
                else
                {
					inloop = true;
				}
			}

            if ((newStem == true))
			{
                if (nbBp > 0) { //Stem
                    posNuc.push_back(*itVertex_current);
                    posNuc1.push_front(*neighbourIt);
                    //posNuc.pop_back();//cause some bugs
                    //posNuc1.pop_front();//cause some bugs
                    visitedPosVect[*neighbourIt] = true;

                    motifFeatTemp << "S";
					motifCurrent.motifName = motifFeatTemp.str();
					motifFeatTemp << "_" << nbBp;
					motifCurrent.motifNameWithFeatLevel1 = motifFeatTemp.str();
                    motifFeatTemp.str("");
                    posNuc.insert( posNuc.end(), posNuc1.begin(), posNuc1.end() );
                    posNuc1.clear();
					motifCurrent.posNuc = posNuc;

					motifVect.push_back(motifCurrent);
					initCurrentMotif(motifCurrent);
					initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
							degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
							nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);

                    //Special multiple loop without nucleotide
                    vertex_iterator_t it = itVertex_current;
                    while ( (*it < *itVertex.second)&& (visitedPosVectSpecialCaseMLoopWithoutNuc[*it]==false))
                        {
                            visitedPosVectSpecialCaseMLoopWithoutNuc[*it] = true;
                            boost::tie(neighbourItLoop, neighbourEndLoop) =  boost::adjacent_vertices(*it, currentGrss);
                            if (neighbourItLoop != neighbourEndLoop)
                            {
                                posNuc1.push_back(*it);
                                posNuc1.push_back(*neighbourItLoop);
                                it = advanceVertex_iterator_t2adjacency_iterator_t(it, neighbourItLoop);
                                degree++;
                            }
                            else
                            {nbUnpairedNuc++;break;}
                            it++;
                        }

                        if (nbUnpairedNuc==0 && degree>2)
                        {
                            posNuc1.push_front(posNuc1.back());
                            posNuc1.pop_back();
                            posNuc.insert( posNuc.end(), posNuc1.begin(), posNuc1.end() );
                            posNuc1.clear();
                            motifFeatTemp << "M";
                            motifCurrent.motifName = motifFeatTemp.str();
                            motifFeatTemp << "_" << nbUnpairedNuc;
                            motifCurrent.motifNameWithFeatLevel1 = motifFeatTemp.str();
                            motifFeatTemp.str("");
                            motifCurrent.posNuc = posNuc;
                            motifVect.push_back(motifCurrent);
                            initCurrentMotif(motifCurrent);
                            initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
                                    degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
                                    nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
                        }

				}
                newStem=false;
            }

            if ((inloop == false))
            {

                nbBp++;
                posNuc.push_back(*itVertex_current);
                posNuc1.push_front(*neighbourIt);
                visitedPosVect[*neighbourIt] = true;

            }

            else //no neighbor
            {
                if (nbBp != 0) //Stem
                    {
                    motifFeatTemp << "S";
					motifCurrent.motifName = motifFeatTemp.str();
					motifFeatTemp << "_" << nbBp;
					motifCurrent.motifNameWithFeatLevel1 = motifFeatTemp.str();
                    motifFeatTemp.str("");
                    posNuc.insert( posNuc.end(), posNuc1.begin(), posNuc1.end() );
                    posNuc1.clear();
					motifCurrent.posNuc = posNuc;


					motifVect.push_back(motifCurrent);

                    if (find(posNuc.begin(), posNuc.end(),currentStructure.length() - 1) != posNuc.end())
                    {lastStemInd = motifVect.size() - 1;}
					initCurrentMotif(motifCurrent);
					initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
							degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
							nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
                    }
                if (nbBp == 0) //Loops
                    {
					vertex_iterator_t it = itVertex_current;
                    vertex_iterator_t itTemp = itVertex_current;
                    itTemp--;
                    posNuc.insert( posNuc.end(), posNuc1.begin(), posNuc1.end());
                    posNuc.push_back(*itTemp);
                    posNuc.push_back(*itVertex_current);
                    it++;
					nbUnpairedNuc++;
					nbNucSidea = nbUnpairedNuc;

                    while ((*it != *itVertex_current) && (*it < *itVertex.second))
                        {
                            boost::tie(neighbourItLoop, neighbourEndLoop) =
                                    boost::adjacent_vertices(*it, currentGrss);
                            if (neighbourItLoop != neighbourEndLoop)
                            {
                                posNuc.push_back(*it);
                                posNuc.push_back(*neighbourItLoop);
                                it = advanceVertex_iterator_t2adjacency_iterator_t(
                                        it, neighbourItLoop);
                                degree++;
                            } else

                            {
                                posNuc.push_back(*it);
                                visitedPosVect[*it] = true;
                                nbUnpairedNuc++;
                                if (degree == 0) {
                                    nbNucSidea = nbUnpairedNuc;
                                } else if (degree == 1) {
                                    nbNucSideb = nbUnpairedNuc - nbNucSidea;
                                }

                            }
                            it++;
                        }

                        posNuc.pop_back();

                    if (degree == 1) //Hloops
                        {
						motifFeatTemp << "H";
						motifCurrent.motifName = motifFeatTemp.str();
						motifFeatTemp << "_" << nbUnpairedNuc;
						motifCurrent.motifNameWithFeatLevel1 =
								motifFeatTemp.str();
                        motifFeatTemp.str("");
						motifCurrent.posNuc = posNuc;
						motifVect.push_back(motifCurrent);
						initCurrentMotif(motifCurrent);
						initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
								degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
								nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
                        }
                    else if ((degree == 2) && (nbNucSideb != 0) && (nbNucSidea != 0)) //Iloops
                        {

                            if (nbNucSidea == nbNucSideb) {
                                symLoop = true;
                            }

                            motifFeatTemp << "I";
                            motifCurrent.motifName = motifFeatTemp.str();
                            motifFeatTemp << "_" << symLoop;
                            motifCurrent.motifNameWithFeatLevel1 =
                                    motifFeatTemp.str();
                            if (nbNucSidea <= nbNucSideb) {
                                motifFeatTemp << "_" << nbNucSidea << "x"
                                        << nbNucSideb;
                            } else {
                                motifFeatTemp << "_" << nbNucSideb << "x"
                                        << nbNucSidea;
                            }
                            motifCurrent.motifNameWithFeatLevel2 =
                                    motifFeatTemp.str();
                            motifFeatTemp.str("");

                            motifCurrent.posNuc = posNuc;
                            motifVect.push_back(motifCurrent);
                            initCurrentMotif(motifCurrent);

                            initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
                                    degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
                                    nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
                        }
                    else if ((degree == 2)) //Bulge
                        {
                            motifFeatTemp << "B";
                            motifCurrent.motifName = motifFeatTemp.str();
                            motifFeatTemp << "_" << nbUnpairedNuc;
                            motifCurrent.motifNameWithFeatLevel1 =
                                    motifFeatTemp.str();
                            motifFeatTemp.str("");

                            if (posNuc.front()>posNuc.back())//case where posNuc are not sorted
                            {
                                posNuc1.push_front(posNuc.back());
                                posNuc.pop_back();
                                posNuc1.push_front(posNuc.back());
                                posNuc.pop_back();
                                posNuc1.insert( posNuc1.end(), posNuc.begin(), posNuc.end() );
                                posNuc.clear();
                                posNuc.insert( posNuc.begin(), posNuc1.begin(), posNuc1.end() );
                                posNuc1.clear();
                            }

                            motifCurrent.posNuc = posNuc;
                            motifVect.push_back(motifCurrent);
                            initCurrentMotif(motifCurrent);
                            initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
                                    degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
                                    nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
                        }

                    else if ((degree > 2)) // Mloops
                        {
                            motifFeatTemp << "M";
                            motifCurrent.motifName = motifFeatTemp.str();
                            motifFeatTemp << "_" << nbUnpairedNuc;
                            motifCurrent.motifNameWithFeatLevel1 =
                                    motifFeatTemp.str();
                            motifFeatTemp.str("");

                            if (posNuc.front()>posNuc.back())//case where posNuc are not sorted
                            {
                                posNuc1.push_front(posNuc.back());
                                posNuc.pop_back();
                                posNuc1.push_front(posNuc.back());
                                posNuc.pop_back();
                                posNuc1.insert( posNuc1.end(), posNuc.begin(), posNuc.end() );
                                posNuc.clear();
                                posNuc.insert( posNuc.begin(), posNuc1.begin(), posNuc1.end() );
                                posNuc1.clear();
                            }

                            motifCurrent.posNuc = posNuc;
                            motifVect.push_back(motifCurrent);
                            initCurrentMotif(motifCurrent);
                            initNbUnpDegBpSideABL1L2L3symLoopPosNuc(nbUnpairedNuc,
                                    degree, nbBp, nbNucSidea, nbNucSideb, symLoop,
                                    nbQuartet, nbNucL1, nbNucL2, nbNucL3, posNuc);
                        }
                    }
            }
		}
	}

    //Circular RNa case where no unpaired nucleotide and 2 stem at the starting and end of dot n bracket
	Motif motifFirstStem, motifLastStem, motifFirstLastStem;
	if (caseStemStemLinkedAtE5E3 == true) {

		motifFirstStem = motifVect.at(0);
		motifLastStem = motifVect.at(lastStemInd);

		motifVect.erase(motifVect.begin() + lastStemInd); //remove last stem in vector of motifs

		motifVect.erase(motifVect.begin()); //remove first stem in vector of motifs

        motifFirstLastStem = motifFirstStem;
		motifFirstLastStem.posNuc.insert(motifFirstLastStem.posNuc.end(),
				motifLastStem.posNuc.begin(), motifLastStem.posNuc.end());
        motifFeatTemp << "S" << motifFirstLastStem.posNuc.size() / 2;
		motifFirstLastStem.motifNameWithFeatLevel1 = motifFeatTemp.str();
		motifFeatTemp.str("");

		motifVect.push_back(motifFirstLastStem);
        caseStemStemLinkedAtE5E3=false;
	}

}


/******************************************************************************/
void RepNmotif::intersectMotifs(int currentNumberOfNeighbor,
        string currentStructure, vector<RepNmotif::Motif>& motifVect,
        map<string, float>& currentStructureFeature,
        map<string, vector<int>>& currentStructureFeatureWithPosNuc, map<string, float> & nmotifsAllStructure)
{

    unsigned int j;
    std::vector<int>::iterator it,itcurrentMotif;
	Motif currentMotif;
	Motif currentMotifwithFeatlevel2;
	vector<RepNmotif::Motif> motifVect1;
	vector<string> subStringVectForSortedFeat, subStringVectForSortedFeat1;

    int sizeStructure = currentStructure.size() + 1000;

    //cout<<currentStructure<<endl;

	if (currentNumberOfNeighbor == 0) {
		for (unsigned int i = 0; i < motifVect.size(); i++) {

			extractFeatFromCurrentMotif(motifVect[i], currentStructureFeature,
                    currentStructureFeatureWithPosNuc,nmotifsAllStructure);
            /*if(currentNumberOfNeighbor==0){
            cout<<motifVect[i].motifName<<'[';
            for (unsigned int l = 0; l < motifVect[i].posNuc.size(); l++) {
              cout<<motifVect[i].posNuc[l]<<',';
            }
            cout<<']'<<endl;}*/
		}

    } else if (currentNumberOfNeighbor >= 1) {

        for (unsigned int i = 0; i < motifVect.size(); i++) {

            /*if(currentNumberOfNeighbor==1){
            cout<<motifVect[i].motifName<<'[';
            for (unsigned int l = 0; l < motifVect[i].posNuc.size(); l++) {
              cout<<motifVect[i].posNuc[l]<<',';
            }
            cout<<']'<<endl;}*/

			//skip the current motif which cover all the nucleotides of the structure
            if(motifVect[i].posNuc.size()>=currentStructure.size())
            {continue;}

			bool temp = false;
			if ((motifVect[i].motifName.compare("H") == 0)
					|| (motifVect[i].motifName.compare("I") == 0)
					|| (motifVect[i].motifName.compare("M") == 0)
					|| (motifVect[i].motifName.compare("B") == 0)) {
				temp = true;
			}

			j = 0;

            currentMotif.motifName.append(motifVect[i].motifName+"[");
            currentMotif.motifNameWithFeatLevel1.append(motifVect[i].motifNameWithFeatLevel1+"[");
			currentMotif.posNuc = motifVect[i].posNuc;

            bool haveNeighbor=false;


			while (j < motifVect.size()) {


                    if ((i != j))
                    {
                        if ((temp == true)	&& (((motifVect[j].motifName.compare("H") == 0)	|| (motifVect[j].motifName.compare("I") == 0)
                                        || (motifVect[j].motifName.compare("M") == 0)|| (motifVect[j].motifName.compare("B") == 0)))) {
                            j++;
                            continue;
                        }

                        if ((motifVect[i].posNuc.back()<motifVect[j].posNuc.front())||(motifVect[j].posNuc.back()<motifVect[i].posNuc.front()))
                        {
                            j++; continue;
                        }

                        vector<int> vInterSect(	motifVect[i].posNuc.size()	+ motifVect[j].posNuc.size());
                        it = set_intersection(motifVect[i].posNuc.begin(),motifVect[i].posNuc.end(),motifVect[j].posNuc.begin(),motifVect[j].posNuc.end(), vInterSect.begin());
                        vInterSect.resize(it - vInterSect.begin());
                        if (!(vInterSect.empty())) {
                            vector<int> vUnion(sizeStructure);
                            subStringVectForSortedFeat.push_back(motifVect[j].motifName);
                            subStringVectForSortedFeat1.push_back(motifVect[j].motifNameWithFeatLevel1);
                            it = std::set_union(currentMotif.posNuc.begin(),currentMotif.posNuc.end(),motifVect[j].posNuc.begin(),motifVect[j].posNuc.end(), vUnion.begin());
                            vUnion.resize(it - vUnion.begin());
                            currentMotif.posNuc = vUnion;
                            vUnion.clear();
                            haveNeighbor=true;
                        }

                        vInterSect.clear();
                    }

				j++;
			}

            if (haveNeighbor==false)
            {
                initCurrentMotif(currentMotif);
                continue;}

			std::sort(subStringVectForSortedFeat.begin(),subStringVectForSortedFeat.end());
			std::sort(subStringVectForSortedFeat1.begin(),subStringVectForSortedFeat1.end());

            //use append
			for (unsigned int j = 0; j < subStringVectForSortedFeat.size();j++) {
				currentMotif.motifName += subStringVectForSortedFeat[j];
				currentMotif.motifNameWithFeatLevel1 +=	subStringVectForSortedFeat1[j];
			}
            currentMotif.motifName.append(+"]");
			currentMotif.motifNameWithFeatLevel1.append("]");

			motifVect1.push_back(currentMotif);
            extractFeatFromCurrentMotif(currentMotif, currentStructureFeature,currentStructureFeatureWithPosNuc,nmotifsAllStructure);

			initCurrentMotif(currentMotif);
			subStringVectForSortedFeat.clear();
			subStringVectForSortedFeat1.clear();
		}
		motifVect = motifVect1;

    }
}

/******************************************************************************/
void RepNmotif::extractFeatFromCurrentMotif(Motif currentMotif,
        map<string, float>& currentStructureFeature,
        map<string, vector<int>>& currentStructureFeatureWithPosNuc, map<string, float> & nmotifsAllStructure) {

    if (!currentMotif.motifName.empty())
    {
		if (!((currentMotif.motifName.compare("H[S]") == 0)
				|| (currentMotif.motifName.compare("I[SS]") == 0)
				|| (currentMotif.motifName.compare("B[SS]") == 0)
				|| (currentMotif.motifName.compare("E5[S]") == 0)
                || (currentMotif.motifName.compare("E3[S]") == 0)))
        {
			currentStructureFeature[currentMotif.motifName]++;
            currentStructureFeatureWithPosNuc[currentMotif.motifName].insert(currentStructureFeatureWithPosNuc[currentMotif.motifName].end(),currentMotif.posNuc.begin(), currentMotif.posNuc.end());

            nmotifsAllStructure[currentMotif.motifName]++;
		}
	}

    if (!currentMotif.motifNameWithFeatLevel1.empty())
    {
		currentStructureFeature[currentMotif.motifNameWithFeatLevel1]++;
        currentStructureFeatureWithPosNuc[currentMotif.motifNameWithFeatLevel1].insert(currentStructureFeatureWithPosNuc[currentMotif.motifNameWithFeatLevel1].end(),currentMotif.posNuc.begin(), currentMotif.posNuc.end());
        nmotifsAllStructure[currentMotif.motifNameWithFeatLevel1]++;
	}

}

/******************************************************************************/
void initNbUnpDegBpSideABL1L2L3symLoopPosNuc(unsigned int& nbUnpairedNuc,
		unsigned int& degree, unsigned int& nbBp, unsigned int& nbNucSidea,
		unsigned int& nbNucSideb, bool& symLoop, unsigned int& nbQuartet,
		unsigned int& nbNucL1, unsigned int& nbNucL2, unsigned int& nbNucL3,
		vector<int>& posNuc) {

	nbUnpairedNuc = degree = nbBp = nbNucSidea = nbQuartet = nbNucSideb =
			nbNucL1 = nbNucL2 = nbNucL3 = 0;
	symLoop = false;
	posNuc.clear();
}

/*********************************************************************************/
void initCurrentMotif(RepNmotif::Motif& motifCurrent) {

	motifCurrent.posNuc.clear();
	motifCurrent.motifName.clear();
	motifCurrent.motifNameWithFeatLevel1.clear();
	motifCurrent.motifNameWithFeatLevel2.clear();
}

/******************************************************************************/
RepNmotif::vertex_iterator_t advanceVertex_iterator_t2adjacency_iterator_t(
		RepNmotif::vertex_iterator_t itbegin,
		RepNmotif::adjacency_iterator_t itend) {

	RepNmotif::vertex_iterator_t it;
     unsigned int dist = 0;
	if (*itend > *itbegin) {

        dist = *itend - *itbegin;
        for (unsigned int i = 0; i < dist; i++) {
            itbegin++;
        }
	} else {
		dist = *itbegin - *itend;

        for (unsigned int i = 0; i < dist; i++) {
            --itbegin;
        }
	}

	return it = itbegin;
}
