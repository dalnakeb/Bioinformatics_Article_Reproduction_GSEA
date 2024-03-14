import argparse
import copy
import math
import random

import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt

import json
from icecream import ic


class GSEA:
    def __init__(self, pathwaysGeneSetsFilePath, geneExpressionsDatasetFilePath, phenotypesFilePath):
        # Fetch data from files
        self.pathwayGeneSets = self._fetchPathwayGeneSetsFromJSON(pathwaysGeneSetsFilePath)
        self.geneExprValues, self.geneNames = self._fetchGeneExprFromGCTFile(geneExpressionsDatasetFilePath)
        self.phenotypes = self._fetchPhenotypes(phenotypesFilePath)
        self._phenotypesCount = np.unique(self.phenotypes, return_counts=True)

        self.ERPerPathway = None

    def _fetchPathwayGeneSetsFromJSON(self, geneSetsFilePath):
        """
        Extracts gene sets from a json file given the filepath as an input
        :param geneSetsFilePath: filepath to a json file containing the gene sets
        :return: a dictionary of {geneSet name: set of genes}
        """
        with open(geneSetsFilePath, 'r') as file:
            geneSetFileDataC1 = file.read()
        geneSetsJSON = json.loads(geneSetFileDataC1)

        geneSets = {}
        for geneSet, items in geneSetsJSON.items():
            getSetName = geneSet
            genesList = set(items["geneSymbols"])
            geneSets[getSetName] = genesList

        return geneSets

    def _fetchGeneExprFromGCTFile(self, genesExprFilePath):
        """
        Extract the list of gene expressions from a gct file given as input the file-path.
        :param genesExprFilePath: file-path to the gct file.
        :return: a list of lists of [gene name, list of expression values from different samples]
        """
        with open(genesExprFilePath, 'r') as file:
            # Read the lines of the file
            lines = file.readlines()

            # Extract gene expression data
            geneExprValues = []
            geneNames = []
            for line in lines[3:]:
                columns = line.split('\t')
                geneName = columns[0].strip("\t")
                if "_at" not in geneName:
                    expressionValues = list(map(float, columns[2:]))
                    if "///" in geneName:
                        geneName = geneName.split(" ")[0]

                    geneNames.append(geneName)
                    geneExprValues.append(expressionValues)

        return np.array(geneExprValues), geneNames

    def _fetchPhenotypes(self, phenotypesFilePath):
        with open(phenotypesFilePath, 'r') as file:
            lines = file.readlines()
            phenotypes = lines[2].strip().split()
            return phenotypes

    def _rankGeneExpr(self, geneExprValues, geneNames):
        """
        Rank the gene expression set according to the metric chosen to rank every line, we average the gene expression
        values for every phenotype calculate the FC in from phenotype 1 to 2 and weighting it by the p-value and rank in
        descendant order.

        :param rankedGeneExpr: dictionary (gene name, list of expression values)
        :param geneNames: geneNames
        :return: ranked numpy list according to the metric described above [gene name, expression values, metric],
        """
        rankedGeneExprValues = np.copy(geneExprValues)

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
        rankedGeneNames = np.array(geneNames)[permutations]

        return rankedGeneExprValues, rankedGeneNames, rankedDGE

    def _enrichmentScore(self, rankedDGE, rankedGeneNames, pathwayGeneSets, p):
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

        return PValuePerPathway, ERScorePerPermutationPerPathway

    def _getCorrectedPValuePerPathway(self,pathwayGeneSets , ERScorePerPathwayNames, ERScorePerPathwayValues, ERScorePerPermutationPerPathway):
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

        return correctedPValuePerPathway

    def getERPerPathwayWithPValue(self):
        # Rank in correlation to phenotypes
        rankedGeneExprValues, rankedGeneNames, rankedDGE = self._rankGeneExpr(self.geneExprValues,
                                                                              self.geneNames)

        # Enrichment score
        ERScorePerPathwayValues, ERScorePerPathwayNames = self._enrichmentScore(rankedDGE,
                                                                                rankedGeneNames,
                                                                                self.pathwayGeneSets,
                                                                                p=1)
        # Estimating significance
        PValuePerPathway, ERScorePerPermutationPerPathway = self._estimateSignificance(ERScorePerPathwayValues,
                                                                                       self.geneExprValues,
                                                                                       self.geneNames,
                                                                                       self.pathwayGeneSets,
                                                                                       p=1,
                                                                                       permutationsNum=10)
        correctedPValuePerPathway = self._getCorrectedPValuePerPathway(self.pathwayGeneSets, ERScorePerPathwayNames,
                                                                       ERScorePerPathwayValues, ERScorePerPermutationPerPathway)

        return PValuePerPathway, ERScorePerPathwayValues, correctedPValuePerPathway, list(self.pathwayGeneSets.keys())




def parseArgs():
    parser = argparse.ArgumentParser(
        description="Process pathway gene sets, gene expressions dataset, and phenotypes file paths.")
    parser.add_argument("pathwaysGeneSetsFilePath", type=str, help="Path to the pathway gene sets JSON file")
    parser.add_argument("geneExpressionsDatasetFilePath", type=str, help="Path to the gene expressions dataset file")
    parser.add_argument("phenotypesFilePath", type=str, help="Path to the phenotypes file")
    args = parser.parse_args()

    return parser, args


if __name__ == '__main__':
    # Input files path
    parser, args = parseArgs()

    # Fetch data from files
    gsea = GSEA(args.pathwaysGeneSetsFilePath, args.geneExpressionsDatasetFilePath, args.phenotypesFilePath)
    PValuesPerPathway, ERScorePerPathwayValues, correctedPValuePerPathway, pathwaysNames = gsea.getERPerPathwayWithPValue()
    ic(PValuesPerPathway[0])
    ic(ERScorePerPathwayValues[0])
    ic(correctedPValuePerPathway[0])
    ic(pathwaysNames[0])

