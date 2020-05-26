#!/usr/bin/python3
#-*- coding: utf-8 -*-

"""

``align_seq_struct.py`` **module description**

This module receives as input a ref RNA structure (input fasta file-name1 in sys.argv[1]) and a set of RNA sequences (multifasta file-name2 in sys.argv[2]) and computes a set of structures for all the RNA in file-name2 (the result is written in the standard output). sys.argv[3] is mandatory but unused. sys.argv[4] is also mandatory and should be "yes" (accounts for Wobble base pairs) or "no" (does not account for Wobble base pairs).

Note:
A stem in a RNA structure is represented as an 3-tuple [s,e,l] 
such that s, e, l are respectively the start position, end position 
and number of base pairs of the stem in the RNA sequence.

An array of stem is an array of 3-tuples representing stems.

.. moduleauthor:: Aida Ouangraoua

2017
update june 2018

"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys
from math import*

#local import
#from path import *
import numpy as np


MINSEP = 3 #minimum length of hairpin loop

#Variables de tests
#test1 = outputs+"setOf_5_g2intron_1/"
#test2 = outputs+"setOf_5_RF00002+5_8S_rRNA_16/"
#test3 = outputs+"setOf_5_RF00037+IRE_18/"
#test4 = outputs+"setOf_5_RF00163+Hammerhead_1_21/"
#test5 = outputs+"setOf_5_RF00168+Lysine_11/"
#test6 = outputs+"setOf_all_g2intron_all/"
statistiques = dict()
option = 1#valeur possible : 1 (first stem), 2(best stem), 3(3-best stems)

def externStems(refStems,refSegment):
    start,end = refSegment
    startStem,endStem = start,start
    refReturned = []
    for stem in refStems:
        if((start <= stem[0] and stem[1] <= end) and (stem[0] >= endStem)):
            refReturned.append(stem)
            starStem,endStem = stem[0],stem[1]
    return refReturned

def charToInt(c):
    """Function that maps nucleotides to integers

    Args:
        c (char) : a nucleotide char

    Returns:
        code (int): unique integer corresponding to nucleotide c

    """
    code = 0
    if c == 'A':
        code = 0
    elif c == 'U':
        code = 1
    elif c == 'G':
        code = 2
    elif c == 'C':
        code = 3
    else:
        code = 0
    return code

def computeMaxStems(seq,initStems,k):
    """Function that computes maximum size stems of a RNA sequence 
       given initial stems

    Args:
        seq (nucleotide string): A RNA sequence
        initStems (array of stems) : An array of all stems of length k in seq
        k (int) : length of all stems contained in initStems

    Returns:
        stems (array of stems) : An array of all maximum size stems of seq

    """
    initStems.sort()
    stems = initStems
    i = 0
    # For each stem stems[i] of length k
    while(i < len(stems)):
        s,e,k = stems[i]
        j=1

        # While extendable
        while([s+j,e-j,k] in stems[i+1:]):
            # Extend it
            stems[i] = [s,e,k+j]
            # Remove included stem
            stems.remove([s+j,e-j,k])
            j+=1
        i+=1
    return stems
    
def computeRefStems(refSeq,refStr):
    """Function that computes the stems of a RNA secondary structure

    Args:
        refSeq (string of nucleotides) : A RNA sequence
        refStr (string of {.,(,)}): A secondary structure for refSeq

    Returns:
        stems (array of stems) : An array of all stems of refStr

    """
    bracketPile = []
    initStems = []
    # Compute all base pairs (stems of length 1)
    for i in range(len(refStr)):
        if(refStr[i]=='('):
            bracketPile.append(i)
        if((refStr[i]==')')):
            initStems.append([bracketPile[-1],i,1])
            del bracketPile[-1]
    # Compute maximum size stems
    stems = computeMaxStems(refSeq,initStems,1)
    return stems

def computeKmerCode(seq,k):
    """Function that computes all positions of k-mers in a RNA sequence

    Args:
        seq (string of nucleotides) : A RNA sequence
        k (int): size of k-mer to consider

    Returns:
        kmerCode (dict) : A map that associates to each k-mer, its list of 
        strat position in seq

    Note:
        Example: there are 12 2-mers : AA, AC, AG, AU, CA, CC, CG, CU, ...  
        k-mer are identified in kmerCode by a unique code
    """
    
    kmerCode = {}
    for i in range(len(seq)-(k-1)):
        kmer = seq[i:i+k]
        code = sum([charToInt(c)*4**j for (j,c) in enumerate(kmer)])
        if( not (code in kmerCode)):
            kmerCode[code] = []
        kmerCode[code].append(i)
    return kmerCode

    
def computeStemsk(rnaSeq,k):
    """Function that computes all stems of length k  in a RNA sequence

    Args:
        rnaSeq (string of nucleotides): A RNA sequence
        k (int): size of stems to consider

    Returns:
        initStems (array of stems): An array containing all stems of 
        length k in rnaSeq

    """
    
    lengthRnaSeq = len(rnaSeq)
    forwardCode = computeKmerCode(rnaSeq,k)
    reverseRnaSeq = str(Seq(rnaSeq).reverse_complement())
    reverseCode = computeKmerCode(reverseRnaSeq,k)
    initStems = []    
    for code in forwardCode.keys():
        if(code in reverseCode.keys()):
            for i in forwardCode[code]:
                for j in reverseCode[code]:
                    if(i < lengthRnaSeq-1-j-2*(k-1)):
                        initStems.append([i,lengthRnaSeq-1-j,k])
                        
    return initStems

def computeRnaMaxStems(rnaSeq,withWobble):
    """Function that computes all maximum size stems in a RNA sequence

    Args:
        rnaSeq (string of nucleotides): A RNA sequence
        withWobble (string): if "yes" accounts for wobble base pairs (G-U)

    Returns:
        rnaMaxStems (array of stems): An array containing all maximum size
        stems of rnaSeq
        wobbleRnaMaxStems (array of stems): An array empty or containing all 
        maximum size stems of rnaSeq while accounting for wobble base pairs

    """
    initStems = computeStemsk(rnaSeq,3)
    rnaMaxStems = computeMaxStems(rnaSeq,initStems,3)
    i = len(rnaMaxStems)-1
    while(i >= 0):
        s,e,l = rnaMaxStems[i]
        if(s+(l-1) >= e-(l-1)-MINSEP):
            rnaMaxStems[i][2] = l-1
        i-=1
    wobbleRnaMaxStems = []
    if(withWobble=="yes"):
        for i in range(len(rnaMaxStems)):
            s,e,l = rnaMaxStems[i]
            initL = l
            j = 1
            while(s-j >= 0 and e+j < len(rnaSeq) and abs(charToInt(rnaSeq[s-j]) - charToInt(rnaSeq[e+j])) == 1):
                j+=1
                
            s,e,l = s-(j-1),e+(j-1),l+(j-1)
            j = 1
            while(s+(l-1)+j < e-(l-1)-j-MINSEP and abs(charToInt(rnaSeq[s+(l-1)+j]) - charToInt(rnaSeq[e-(l-1)-j])) == 1):
                j+=1
            if(l+(j-1) > initL):
                wobbleRnaMaxStems.append([s,e,l+(j-1)])

    return rnaMaxStems, wobbleRnaMaxStems

################################################################
#A modifier : formule de calcul de la distance entre deux stems#
################################################################

def computeStemDistance(stemInRef,stemInRna,refSeq,rnaSeq,refSegment, rnaSegment):
    """Function that computes a distance between two stems from two RNA

    Args:
        stemInRef (stem (3-tuple)): A ref stem
        stemInRna (stem (3-tuple)): A rna stem
        refSegment (2-tuple): start and end position of a ref segment containing stemInRef
        rnaSegment (2-tuple): start and end position of a rna segment containing stemInRna

    Returns:
        dist (float): distance between stems

    """

    startRef,endRef = refSegment
    startRna, endRna = rnaSegment
    diff_dist_start = abs((stemInRef[0]-startRef)-(stemInRna[0]-startRna))# difference in distance to start
    diff_dist_end = abs((endRef - stemInRef[1])-(endRna-stemInRna[1]))# difference in distance to end
    diff_loop_length = abs((stemInRef[1]-stemInRef[0]+1-2*(stemInRef[2]))-(stemInRna[1]-stemInRna[0]+1-2*(stemInRna[2])))# difference in loop length
    diff_length = abs(stemInRef[2]-stemInRna[2])# difference in number of base pairs
    diff_seq = 0#compute_distance_seq(refSeq,rnaSeq)
    all_dists = [diff_dist_start,diff_dist_end,diff_loop_length,diff_length,diff_seq]
    all_dists.sort()
    #dist = sqrt(all_dists[0]**2 + all_dists[1]**2 + all_dists[2]**2 + all_dists[3]**2)##euclidian distance
    
    return  all_dists 

def computeClosestRnaStems(refStems,rnaMaxStems, refSeq, rnaSeq, refSegment, rnaSegment, method="first"):
    """Function that computes a distance between two stems from two RNA

    Args:
        refStems (array of stems): An array of ref stems contained in refSegment
        rnaMaxStems (stem (3-tuple)): An array of rna stems contained in rnaSegment
        refSegment (2-tuple): start and end position of a ref segment containing the ref stems
        rnaSegment (2-tuple): start and end position of a rna segment containing the rna stems
        method (string in {"First","Best","Last","Worst","3"})

    Returns:
        closestRnaStems (array of stems): An array of compatible RNA stems corresponding to ref stems
        distance (int) : the total distance between each closertRnaStems' distance and their refStem

    """

    startRef,endRef = refSegment
    startRna, endRna = rnaSegment
    closestRnaStemsRetuned = [];

    # Compute set of ref stems in segment
    refStemsInSegment = externStems(refStems,refSegment)
    # refStemsInSegment = []                #ancien code
    # for stem in refStems:
    #     if (startRef <= stem[0] and stem[1]<=endRef):
    #         refStemsInSegment.append(stem)
    
    # Compute set of rna stems in segment
    rnaStemsInSegment = []
    for stem in rnaMaxStems:
        if (startRna <= stem[0] and stem[1]<=endRna):
            rnaStemsInSegment.append(stem)

    distance = 0
    if(len(refStemsInSegment) > 0 and len(rnaStemsInSegment) > 0):
        #listes des stems
        listdist = []
        listcorrStemInRna = []
        corrStemInRnaFinal = []
        
        #calcule de l'exactitude des stems
        for stemInRef in refStemsInSegment:
            # Start by left-side stem
            corrStemInRna = []
    
            # Distance between stemInRef and any stemInRna 
            # dist = [[computeStemDistance(stemInRef,stemInRna,refSegment, rnaSegment),stemInRna] for stemInRna in rnaStemsInSegment]   #old version          
            distBrut = np.array([computeStemDistance(stemInRef,stemInRna,refSeq,rnaSeq,refSegment, rnaSegment) for stemInRna in rnaStemsInSegment])
            #print(distBrut[0:2]) #code de test
            #dist = [[sqrt(distBrut[i,0]**2 + distBrut[i,1]**2 + distBrut[i,2]**2 + distBrut[i,3]**2 + distBrut[i,4]**2),rnaStemsInSegment[i]] for i in range(len(rnaStemsInSegment))]   #distance-5
            dist = [[sqrt(distBrut[i,0]**2 + distBrut[i,1]**2 + distBrut[i,2]**2 + distBrut[i,3]**2),rnaStemsInSegment[i]] for i in range(len(rnaStemsInSegment))]                     #distance-4
            #print(dist[0:2]) #code de test
            
            dist.sort()
            
            #ajout dans une liste
            if(method == "worst"):
                listdist.append(dist[-1][0])
                listcorrStemInRna.append(dist[-1][1])
            else:
                listdist.append(dist[0][0])
                listcorrStemInRna.append(dist[0][1])
    
        # Choose the closest stem in RNA
        if(method == "first" or method == "last"):
            arraycorrStemInRna = listcorrStemInRna #trie par ordre de stem
            arraycorrStemInRef = refStemsInSegment #trie par ordre de stem
            arraydist = listdist                   #trie par ordre de stem
        else:
            arraydist = np.sort(np.array(listdist))                                          #trie par meilleur stem
            arraycorrStemInRna = np.array(listcorrStemInRna)[np.argsort(np.array(listdist))] #trie par meilleur stem
            arraycorrStemInRef = np.array(refStemsInSegment)[np.argsort(np.array(listdist))] #trie par meilleur stem
        
        # Choose the n closest stems and loop recursively
        if(method == "3"):
            n = min(3,len(arraydist))
        else:
            n = 1
        distances = [arraydist[i] for i in range(n)]
        for i in range(n):
            if(method == "last" or method == "worst"):
                corrStemInRna = arraycorrStemInRna[-1]
                stemInRef = arraycorrStemInRef[-1]
            else:
                corrStemInRna = arraycorrStemInRna[i]
                stemInRef = arraycorrStemInRef[i]

            # Call computeClosestRnaStems recursively (inside loop and on the right side)
            stemsBefore,distB = computeClosestRnaStems(refStems,rnaMaxStems,refSeq,rnaSeq,[startRef, stemInRef[0]], [startRna,corrStemInRna[0]],method)
            stemsMiddle,distM = computeClosestRnaStems(refStems,rnaMaxStems,refSeq,rnaSeq,[stemInRef[0]+stemInRef[2],stemInRef[1]-stemInRef[2]],[corrStemInRna[0]+corrStemInRna[2],corrStemInRna[1]-corrStemInRna[2]],method)
            stemsAfter,distA = computeClosestRnaStems(refStems,rnaMaxStems,refSeq,rnaSeq,[stemInRef[1]+1, endRef], [corrStemInRna[1]+1,endRna],method)
            closestRnaStems = stemsBefore + [corrStemInRna] + stemsMiddle + stemsAfter
            #closestRnaStems = [corrStemInRna] + stemsMiddle + stemsAfter
            
            # Append the final structure into a list (size n) of structures
            distances[i] += distB + distM + distA
            corrStemInRnaFinal.append(closestRnaStems)
        
        # return the closest structures
        distStructures = np.array(distances)
        closestRnaStemsRetuned = corrStemInRnaFinal[np.argmin(distStructures)]
        
    return closestRnaStemsRetuned,distance

def extendStems(rnaSeq,closestRnaStr,closestRnaStems):
    """Function that extends stems for a RNA

    Args:
        rnaSeq (string of nucleotides): A RNA sequence
        closestRnaStr (string of {.,(,)}): the secondary structure corresponding to initial Rna stems
        closestRnaStems (array of stems): An array of initial compatible RNA stems 

    Returns:
        closestRnaStems (array of stems): An array of final compatible RNA stems 
        closestRnaStr (string of {.,(,)}): the secondary structure corresponding to final Rna stems

    """
    closestRnaStr_ = list(closestRnaStr)
    # For each stem
    for k in range(len(closestRnaStems)):
        # Extend by exterior as far as there are Watson-Crick pairs
        s,e,l = closestRnaStems[k]
        i = 1
        while(s-i >= 0 and e+i < len(rnaSeq) and abs(charToInt(rnaSeq[s-i])- charToInt(rnaSeq[e+i]))==1 and closestRnaStr_[s-i]=='.' and closestRnaStr_[e+i]=='.'):
            closestRnaStr_[s-i] = '('
            closestRnaStr_[e+i] = ')'
            i += 1
        closestRnaStems[k] = [s-(i-1),e+(i-1),l+(i-1)]
        
        # Extend by interior as far as there are Watson-Crick pairs
        s,e,l = closestRnaStems[k]
        i = 1
        while(s+(l-1)+i <  e-(l-1)-i-MINSEP  and abs(charToInt(rnaSeq[s+(l-1)+i])- charToInt(rnaSeq[e-(l-1)-i]))==1 and closestRnaStr_[s+(l-1)+i]=='.' and closestRnaStr_[e-(l-1)-i]=='.'):
            closestRnaStr_[s+(l-1)+i] = '('
            closestRnaStr_[e-(l-1)-i] = ')'
            i += 1
        closestRnaStems[k] = [s,e,l+(i-1)]
    closestRnaStr = "".join(closestRnaStr_)
    return closestRnaStems, closestRnaStr


def main(refFile,rnaFile,withWobble,method="first",oldRun = False,extendRef = True,extendRNA = True, verbose = True):
    """
    Inputs:
        - refFile (string) : fichier de la structure de référence
        - rnaFile (string) : fichier des structures sous-optimales de la séquence à aligner
        - withWobble (bool) :  prise en compte de paires wobbles (/!\ inactif dans cette version /!\)
        - method (string) : méthode d'alignement choisie. Peut valoir
            - First : on aligne en commençant par la première stem
            - Last  : on aligne en commençant par la dernière stem
            - Best  : on aligne en commençant par la meilleure stem
            - Worst : on aligne en commençant par la pire stem
            - 3     : on aligne en parcourant les 3-meilleures stems
        - oldRun (bool) : joue l'ancienne version du script
        - extendRef (bool) : étend les stems de la référence
        - extendRNA (bool) : étend les stems des structures sous-optimales
        - verbose (bool) : imprime des informations sur la console
    Doing :
        aligne selon la technique "stem-extérieures en priorité"
    """
    
    #refFile = sys.argv[1] #resultat d'alifreefold
    #rnaFile = sys.argv[2] #fichier fasta des sequences 
    #withWobble = sys.argv[4] #avec ou sans wobble (G-U)
    
    #global option #code de test
    output = []
    
    if(verbose):
        print("Option : \nwithWobble : {}\nextendRef : {}\nextendRNA : {}\nAncien Programme : {}\n".format(withWobble,extendRef,extendRNA,oldRun))
    
    refLines = open(refFile,"r")

    ref = []
    for record in SeqIO.parse(refFile, "fasta"):
        ref.append([record.id,str(record.seq)])

    refLen = min([ref[0][1].find(i) for i in {".","("}]) #(len(ref[0][1]))/2
    refId = ref[0][0]
    refSeq = str(ref[0][1][:refLen])
    refStr = str(ref[0][1][refLen:])
    #refSeq,refStr = InvertStructSeq(refSeq,refStr) #inversion de la structure

    refStems = computeRefStems(refSeq,refStr)

    if(extendRef):
        refStems, refStr = extendStems(refSeq,refStr,refStems)
    else:
        refStr = "".join(refStr)
        
    if(verbose):
        print(">"+refId)
        print(refSeq)
        print(refStr)
        print(refStems)
        print("####\n")

    if(oldRun):
        rna = []
        for record in SeqIO.parse(rnaFile, "fasta"):
            rna.append([record.id,str(record.seq)])
    else:
        rna = {}
        # print("extract rna subopt : ",end="") #code de test
        for record in SeqIO.parse(rnaFile, "fasta"):
            old = rna.get(record.id, list())
            old.append(str(record.seq))
            rna[record.id] = old
        #     print(".",end="") #code de test
        # print("\nCompute Stems :\n") #code de test
        

    if(not oldRun):
        rna = rna.items()
        
    for rnaId,rnaSeqs in rna:
        if(oldRun):
            rnaSeq = rnaSeqs
            
        rnaMaxStems = []
        
        if(oldRun):
            rnaLen = min([rnaSeq.find(i) for i in {".","("}]) #(len(ref[0][1]))/2
            rnaStr = str(rnaSeq[rnaLen:])
            rnaSeq = str(rnaSeq[:rnaLen])
            #rnaSeq,rnaStr = InvertStructSeq(rnaSeq,rnaStr)  #inversion de la structure
            rnaMaxStems, wobbleRnaMaxStems = computeRnaMaxStems(rnaSeq,withWobble)
            rnaMaxStems += wobbleRnaMaxStems
        else:
            p1 = []
            for rnaSeq in rnaSeqs:
                rnaLen = min([rnaSeq.find(i) for i in {".","("}]) #(len(ref[0][1]))/2
                rnaStr = str(rnaSeq[rnaLen:])
                rnaSeq = str(rnaSeq[:rnaLen])
                #rnaSeq,rnaStr = InvertStructSeq(rnaSeq,rnaStr)  #inversion de la structure
                a = computeRefStems(rnaSeq,rnaStr)
                p1.extend(a) 
            for el in p1:
                if el in rnaMaxStems:
                    pass
                else:
                    rnaMaxStems.append(el)
        
        lenRnaSeq = len(rnaSeq)            
        lenRefSeq = len(refSeq)
        
        # Compute set of closest RNA stems
        # option = 1 #code de test
        closestRnaStems = computeClosestRnaStems(list(refStems),rnaMaxStems,refSeq,rnaSeq,[0,lenRefSeq],[0,lenRnaSeq],method)[0]
        # print closestRnaStems #code de test

        # Compute structure corresponding to closest RNA stems
        closestRnaStr = list('.'*len(rnaSeq))
        for stem in closestRnaStems:
            s,e,l = stem
            for i in range(l):
                closestRnaStr[s+i] = '('
                closestRnaStr[e-i] = ')'

        # Extend stems with all possible Watson-Crick pairs
        if(extendRNA):
            closestRnaStems, closestRnaStr = extendStems(rnaSeq,closestRnaStr,closestRnaStems)
        else:
            closestRnaStr = "".join(closestRnaStr)
        if(verbose):
            print(">"+rnaId)
            cent = ""
            dix = ""
            zero = ""
            for i in range(len(rnaSeq)):
                cent+= str(i//100)
                dix += str((i%100)//10)
                zero += str((i%100)%10)
            print(cent)
            print(dix)
            print(zero)
            print(rnaSeq)
            print(closestRnaStr)
            print(closestRnaStems)
            print(refStems)
            if(withWobble):
                print("")
            print("#####\n")
        #rnaSeq,closestRnaStr = InvertStructSeq(rnaSeq,closestRnaStr)  #inversion de la structure
        output.append((rnaId,rnaSeq,closestRnaStr))
        
    return output
    #écrire dans un fichier sous le même format que dans "alifreefold_result"
    
#main()
