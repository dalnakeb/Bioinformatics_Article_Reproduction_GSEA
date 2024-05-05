import argparse
import math
import random
import time
from copy import copy

import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt

import json
from icecream import ic





   # def _rankGeneExpr(self, geneExprValues, geneNames, classOfDistinction):
"""
        Rank the gene expression set according to the metric chosen to rank every line, we choose the class of
        distinction to be evaluated and sort accordingly in descending order.

        :param rankedGeneExpr: dictionary (gene name, list of expression values)
        :param geneNames: geneNames
        :return: ranked numpy list according to the metric described above [gene name, expression values, metric],
        """
      #  rankedGeneExprValues = np.copy(geneExprValues)
       # print(rankedGeneExprValues[0])
"""
        # Correlation of genes with phenotype 1
        GroupASize = self._phenotypesCount[1][0]
        GroupBSize = self._phenotypesCount[1][1]

        # Extract expression values for both groups
        GroupA = rankedGeneExprValues[:, :GroupASize]
        GroupB = rankedGeneExprValues[:, GroupBSize:]

        # Perform t-test on each gene to calculate p-value
        t_statistics, p_values = ttest_ind(GroupA, GroupB, axis=1)

        # Calculate fold change
        foldChange = (np.sum(GroupB, axis=1) / GroupBSize) - (np.sum(GroupA, axis=1) / GroupASize)
        foldChange[foldChange > 0] = np.log2(foldChange[foldChange > 0])
        foldChange[foldChange < 0] = -np.log2(np.abs(foldChange[foldChange < 0]))
        foldChange[foldChange == 0] = 0

        # Calculate differential gene expression
        rankedDGE = foldChange * (-np.log10(p_values))

        # Sort rankedGeneExpr by differential gene expression in descendant order

        permutations = rankedDGE.argsort()[::-1]
        rankedDGE = rankedDGE[permutations]
        rankedGeneExprValues = rankedGeneExprValues[permutations]
        rankedGeneNames = np.array(geneNames)[permutations]"""
      #  return 0,0,0
#        return rankedGeneExprValues, rankedGeneNames, rankedDGE

"""def _enrichmentScore(self, rankedDGE, rankedGeneNames, pathwayGeneSets, p):
        pathwaysNum = min(1000, len(pathwayGeneSets))
        ERScorePerPathwayValues = np.zeros(len(list(pathwayGeneSets.items())[0:pathwaysNum]))
        ERScorePerPathwayNames = np.array([None for _ in range((len(list(pathwayGeneSets.items())[0:pathwaysNum])))])
        ERScorePerPathwayIndex = 0

        for pathwayGeneSetName, geneSet in list(pathwayGeneSets.items())[0:pathwaysNum]:
            runningSum = 0
            minER = math.inf
            maxER = -math.inf

            ERScorePerPathwayNames[ERScorePerPathwayIndex] = pathwayGeneSetName
            # if ERScorePerPathwayIndex % 75 == 0:
            #    ic(ERScorePerPathwayIndex)
            # Extract gene indices for the current pathway
            pathwayGeneMatch = np.isin(rankedGeneNames, geneSet)
            pathwayGeneIndices = np.where(pathwayGeneMatch)[0]
            nonPathwayGeneIndices = np.where(~pathwayGeneMatch)[0]

            # Calculate running sum for pathway genes
            pathwayRunningSum = np.sum(np.abs(rankedDGE[pathwayGeneIndices]) ** p)

            # Calculate running sum for non-pathway genes
            nonPathwayRunningSum = np.sum(-rankedDGE[nonPathwayGeneIndices] ** p)

            # Update minER and maxER based on running sums
            minER = min(minER, pathwayRunningSum, nonPathwayRunningSum)
            maxER = max(maxER, pathwayRunningSum, nonPathwayRunningSum)

            # Append minER or maxER to ERScorePerPathway
            if abs(minER) > abs(maxER):
                ERScorePerPathwayValues[ERScorePerPathwayIndex] = minER
            else:
                ERScorePerPathwayValues[ERScorePerPathwayIndex] = maxER

            ERScorePerPathwayIndex += 1

        return ERScorePerPathwayValues, ERScorePerPathwayNames

    def _shuffleGeneExprValues(self, randomGeneExprValues):
        np.apply_along_axis(np.random.shuffle, axis=1, arr=randomGeneExprValues)

    def _estimateSignificance(self, ERScorePerPathwayValues, geneExprValues, geneNames, pathwayGeneSets, p,
                              permutationsNum):
        ERScorePerPermutationPerPathway = np.zeros((permutationsNum, len(ERScorePerPathwayValues)))
        randomGeneExprValues = np.copy(geneExprValues)
        for perm in range(permutationsNum):
            self._shuffleGeneExprValues(randomGeneExprValues)
            ic(perm)
            rankedGeneExprValues, rankedGeneNames, rankedDGE = self._rankGeneExpr(randomGeneExprValues, geneNames)
            ERScorePerPermutationPerPathway[perm], _ = self._enrichmentScore(rankedDGE, rankedGeneNames,
                                                                             pathwayGeneSets, p)

        PValuePerPathway = np.zeros(len(ERScorePerPathwayValues))

        for pathwayERScoreIndex in range(len(ERScorePerPathwayValues)):
            ERNullHypothesisDistribution = ERScorePerPermutationPerPathway[:, pathwayERScoreIndex]
            if ERScorePerPathwayValues[pathwayERScoreIndex] >= 0:
                PValuePerPathway[pathwayERScoreIndex] = len(ERNullHypothesisDistribution[
                                                                ERNullHypothesisDistribution >= ERScorePerPathwayValues[
                                                                    pathwayERScoreIndex]]) / len(
                    ERNullHypothesisDistribution)
            else:
                PValuePerPathway[pathwayERScoreIndex] = len(ERNullHypothesisDistribution[
                                                                ERNullHypothesisDistribution <= ERScorePerPathwayValues[
                                                                    pathwayERScoreIndex]]) / len(
                    ERNullHypothesisDistribution)

        return PValuePerPathway, ERScorePerPermutationPerPathway"""

"""def _getCorrectedPValuePerPathway(self,pathwayGeneSets , ERScorePerPathwayNames, ERScorePerPathwayValues, ERScorePerPermutationPerPathway):
        correctedPValuePerPathway = np.zeros(len(ERScorePerPathwayValues))
        geneSetLengthPerPathway = np.array([len(pathwayGeneSets[pathwayName]) for pathwayName in ERScorePerPathwayNames])
        NERScorePerPathwayValues = ERScorePerPathwayValues/geneSetLengthPerPathway
        NERScorePerPermutationPerPathway = ERScorePerPermutationPerPathway/geneSetLengthPerPathway

        for pathwayERScoreIndex in range(len(NERScorePerPathwayValues)):
            NERNullHypothesisDistribution = ERScorePerPermutationPerPathway[:, pathwayERScoreIndex]
            if NERScorePerPathwayValues[pathwayERScoreIndex] >= 0:
                correctedPValuePerPathway[pathwayERScoreIndex] = len(NERNullHypothesisDistribution[
                                                                NERNullHypothesisDistribution >= NERScorePerPathwayValues[
                                                                    pathwayERScoreIndex]]) / len(
                    NERNullHypothesisDistribution)
            else:
                correctedPValuePerPathway[pathwayERScoreIndex] = len(NERNullHypothesisDistribution[
                                                                NERNullHypothesisDistribution <= NERScorePerPathwayValues[
                                                                    pathwayERScoreIndex]]) / len(
                    NERNullHypothesisDistribution)

        return correctedPValuePerPathway"""

"""    def getERPerPathwayWithPValue(self, classOfDistinction = 1):
        # Rank in correlation to phenotypes
        rankedGeneExprValues, rankedGeneNames, rankedDGE = self._rankGeneExpr(self.geneExprValues,
                                                                              self.geneNames, classOfDistinction)

        # Enrichment score
       # ERScorePerPathwayValues, ERScorePerPathwayNames = self._enrichmentScore(rankedDGE,
       #                                                                         rankedGeneNames,
        #                                                                        self.pathwayGeneSets,
        #                                                                        p=1)
        # Estimating significance
       # PValuePerPathway, ERScorePerPermutationPerPathway = self._estimateSignificance(ERScorePerPathwayValues,
        #                                                                               self.geneExprValues,
        #                                                                               self.geneNames,
        #                                                                               self.pathwayGeneSets,
        #                                                                               p=1,
        #                                                                               permutationsNum=10)
       # correctedPValuePerPathway = self._getCorrectedPValuePerPathway(self.pathwayGeneSets, ERScorePerPathwayNames,
       #                                                                ERScorePerPathwayValues, ERScorePerPermutationPerPathway)
        return 0,0,0,0
#        return PValuePerPathway, ERScorePerPathwayValues, correctedPValuePerPathway, list(self.pathwayGeneSets.keys())"""


def fetchSDict(SDictFilePath):
    """
    Extracts gene sets from a json file given the filepath as an input. Gene sets with 15 genes and more are chosen
    :param SDictFilePath: filepath to a json file containing the gene sets
    :return: a dictionary of {geneSet name: set of genes}
    """
    with open(SDictFilePath, 'r') as file:
        geneSetFileDataC1 = file.read()
    geneSetsJSON = json.loads(geneSetFileDataC1)

    SDict = {}
    for SName, items in geneSetsJSON.items():
        if len(items["geneSymbols"]) >= 15:
            S = set(items["geneSymbols"])
            SDict[SName] = S
    print(f"{len(SDict)} pathways genes sets have been fetched")
    return SDict


def fetchD(DFilePath):
    """
    Extract the list of gene expressions from a gct file given as input the file-path.
    :param DFilePath: file-path to the gct file.
    :return: a list of lists of [gene name, list of expression values from different samples]
    """
    with open(DFilePath, 'r') as file:
        # Read the lines of the file
        lines = file.readlines()

        # Extract gene expression data
        D = []
        DNames = []
        for line in lines[3:]:
            columns = line.split('\t')
            DName = columns[0].strip("\t")
            if "_at" not in DName:
                expressionValues = list(map(float, columns[2:]))
                if "///" in DName:
                    DName = DName.split(" ")[0]

                DNames.append(DName)
                D.append(expressionValues)

    print(f"{len(D)} genes have been fetched")
    return np.array(D), DNames


def fetchP(PFilePath):
    with open(PFilePath, 'r') as file:
        lines = file.readlines()
        phenotypes = lines[2].strip().split()
        print(f"{len(phenotypes)} phenotypes have been fetched")
        return phenotypes


def sortD(D, DNames, P, C):
    D = copy(D)
    DNames = copy(DNames)
    #  Extract indices of the chosen phenotype
    uniqueP = list(set(P))
    if C == 0:
        indexP = P.index(uniqueP[1])
        PIndexes = np.arange(0, indexP)
    else:
        indexP = P.index(uniqueP[1])
        PIndexes = np.arange(indexP, len(P))
    #  Sort according to one prob for each gene of the chosen phenotype
    sortedD = []
    for i in range(len(D)):
        sortedD.append([D[i], D[i][random.choice(PIndexes)], DNames[i]])
    sortedD = sorted(sortedD, key=lambda x: x[1], reverse=True)
    #  Extract names and values for gene expressions
    LNames = []
    L = []
    LProbs = []
    for i in range(len(sortedD)):
        LNames.append(sortedD[i][2])
        L.append(sorted(sortedD[i][0], reverse=True))
        LProbs.append(sortedD[i][1])
    return np.array(L), np.array(LProbs), np.array(LNames)

def ER(S, LProbs, LNames, p):
    NH = len(S)
    N = len(LProbs)
    NR = np.sum(abs(LProbs[np.where(np.isin(LNames, S))])**p)
    ERS = 0
    absMaxDev = -math.inf

    oldPHits = 0
    oldPMisses = 0
    LInS = np.isin(LNames, S)
    for i in range(len(L)):
        PHits = oldPHits
        PMisses = oldPMisses

        if LInS[i]:
            PHits += (abs(LProbs[i])**p) / NR
        else:
            PMisses += 1/(N-NH)

        oldPHits = PHits
        oldPMisses = PMisses
        if absMaxDev < abs(PHits-PMisses):
            absMaxDev = abs(PHits-PMisses)
            ERS = PHits-PMisses

    return ERS


def parseArgs():
    parser = argparse.ArgumentParser(
        description="Process pathway gene sets, gene expressions dataset, and phenotypes file paths.")
    parser.add_argument("SDictFilePath", type=str, help="Path to the pathway gene sets JSON file")
    parser.add_argument("LFilePath", type=str, help="Path to the gene expressions dataset file")
    parser.add_argument("PFilePath", type=str, help="Path to the phenotypes file")
    args = parser.parse_args()

    return parser, args


if __name__ == '__main__':
    # Input files path
    parser, args = parseArgs()

    # Fetch data from files
    SDict = fetchSDict(args.SDictFilePath)  # SDict: Dict {pathway name: gene sets values}
    D, DNames = fetchD(args.DFilePath)  # D: list of lists of gene expression values of k samples. DName: list of names of every gene in the respective order
    P = fetchP(args.PFilePath)  # P: list of phenotypes (p1,..,p1n,p21,...,p2m)

    # choose a class of distinction and a pathway to evaluate
    pathway = "pathway"
    S = np.array(list(SDict["chr14q32"]))
    C = 0  # class of distinction
    L, LProbs, LNames = sortD(D, DNames, P, C)
    p = 1  # weight of the step
    ERS = ER(S, LProbs, LNames, p)
    print(ERS)

