import copy
import math
import random

import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt

import json
from icecream import ic


def enrichmentScore(rankedGeneExprValues, rankedGeneNames, pathwayGeneSets, p):
    ERScorePerPathway = []
    ERScorePerPathwayIndex = 0

    for pathwayGeneSetName, geneSet in list(pathwayGeneSets.items())[0:200]:
        runningSum = 0
        minER = math.inf
        maxER = -math.inf
        ERScorePerPathway.append([pathwayGeneSetName])
        #if ERScorePerPathwayIndex % 75 == 0:
        #    ic(ERScorePerPathwayIndex)
        # Extract gene indices for the current pathway
        pathwayGeneMatch = np.isin(rankedGeneNames, geneSet)
        pathwayGeneIndices = np.where(pathwayGeneMatch)[0]
        nonPathwayGeneIndices = np.where(~pathwayGeneMatch)[0]

        # Calculate running sum for pathway genes
        pathwayRunningSum = np.sum(np.abs(rankedGeneExprValues[pathwayGeneIndices, -1]) ** p)

        # Calculate running sum for non-pathway genes
        nonPathwayRunningSum = np.sum(-rankedGeneExprValues[nonPathwayGeneIndices, -1] ** p)

        # Update minER and maxER based on running sums
        minER = min(minER, pathwayRunningSum, nonPathwayRunningSum)
        maxER = max(maxER, pathwayRunningSum, nonPathwayRunningSum)

        # Append minER or maxER to ERScorePerPathway
        if abs(minER) > abs(maxER):
            ERScorePerPathway[ERScorePerPathwayIndex].append(minER)
        else:
            ERScorePerPathway[ERScorePerPathwayIndex].append(maxER)

        ERScorePerPathwayIndex += 1

        """runningSum = 0
        minER = math.inf
        maxER = -math.inf
        ERScorePerPathway.append([pathwayGeneSetName])

        for geneIndex in range(len(rankedGeneExprValues)):
            if rankedGeneNames[geneIndex] in pathwayGeneSets[pathwayGeneSetName]:
                runningSum += abs(rankedGeneExprValues[geneIndex][-1]**p)
            else:
                runningSum -= rankedGeneExprValues[geneIndex][-1]**p

            minER = min(minER, runningSum)
            maxER = max(maxER, runningSum)

        if abs(minER) > abs(maxER):
            ERScorePerPathway[ERScorePerPathwayIndex].append(minER)
        else:
            ERScorePerPathway[ERScorePerPathwayIndex].append(maxER)
        ERScorePerPathwayIndex += 1"""

    return ERScorePerPathway


def rankGeneExpr(geneExprValues, geneNames, phenotypesCount):
    """
    Rank the gene expression set according to the metric chosen to rank every line, we average the gene expression
    values for every phenotype calculate the FC in from phenotype 1 to 2 and weighting it by the p-value and rank in
    descendant order.

    :param rankedGeneExpr: dictionary (gene name, list of expression values)
    :param phenotypes: a list of phenotypes indicating the phenotype of an expression in a certain position
    :return: ranked list according to the metric described above [gene name, expression values, metric],
    """
    rankedGeneExprValues = np.copy(geneExprValues)
    rankedGeneNames = np.copy(geneNames)
    # Correlation of genes with phenotype 1
    GroupASize = phenotypesCount[1][0]
    GroupBSize = phenotypesCount[1][1]

    # Extract expression values for both groups
    GroupA = rankedGeneExprValues[:, :GroupASize]
    GroupB = rankedGeneExprValues[:, GroupBSize:]

    # Perform t-test on each gene
    t_statistics, p_values = ttest_ind(GroupA, GroupB, axis=1)

    # Calculate fold change
    foldChange = (np.sum(GroupB, axis=1) / GroupBSize) - (np.sum(GroupA, axis=1) / GroupASize)
    foldChange[foldChange > 0] = np.log2(foldChange[foldChange > 0])
    foldChange[foldChange < 0] = -np.log2(np.abs(foldChange[foldChange < 0]))
    foldChange[foldChange == 0] = 0

    # Calculate differential gene expression
    differentialGeneExpr = foldChange * (-np.log10(p_values))

    # Append differential gene expression to rankedGeneExpr
    rankedGeneExprValues = np.column_stack((rankedGeneExprValues, differentialGeneExpr))

    # Sort rankedGeneExpr by differential gene expression
    permutations = rankedGeneExprValues[:, -1].argsort()[::-1]
    rankedGeneExprValues = rankedGeneExprValues[permutations]
    rankedGeneNames = np.array(geneNames)[permutations]

    """# Order every gene samples by ascendant value of expression (and the phenotype list)
    phenotypesPerGene = []
    for i in range(len(rankedGeneExpr)):
        permutationOrder = np.argsort(rankedGeneExpr[i][1])
        rankedGeneExpr[i][1] = rankedGeneExpr[i][1][permutationOrder]
        phenotypesPerGene.append([])
        phenotypesPerGene[i].append(phenotypes[permutationOrder])"""

    return rankedGeneExprValues, rankedGeneNames


def fetchGeneExprFromGCTFile(genesExprFilePath):
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


def fetchPathwayGeneSetsFromJSON(geneSetsFilePath):
    """
    Extracts gene sets from a json file given the filepath as an input
    :param geneSetsFilePath: filepath to a json file containing the gene sets
    :return: a dictionary of [geneSet name, set of genes]
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


def shuffleGeneExpr(randomGeneExpr):
    random.shuffle(randomGeneExpr)


def estimateSignificance(enrichmentScorePerPathway, geneExprValues, geneNames, phenotypesCount, pathwayGeneSets, p):
    ERNullDistributionPerPathway = []
    randomGeneExprValues = copy.deepcopy(geneExprValues)
    for i in range(10):
        shuffleGeneExpr(randomGeneExprValues)
        ic(i)
        rankedGeneExprValues, rankedGeneNames = rankGeneExpr(randomGeneExprValues, geneNames, phenotypesCount)
        ERNullDistributionPerPathway.append(enrichmentScore(rankedGeneExprValues, rankedGeneNames, pathwayGeneSets, p))

    for pathwayGeneSetIndex in range(len(ERNullDistributionPerPathway[0])):
        pathwayERScores = []
        for erScoreIndex in range(len(ERNullDistributionPerPathway)):
            pathwayERScores.append(ERNullDistributionPerPathway[erScoreIndex][pathwayGeneSetIndex])
        pathwayERScores = np.array(pathwayERScores)
        if enrichmentScorePerPathway[pathwayGeneSetIndex][1] >= 0:
            p_value = 0
            for pathwayERScore in pathwayERScores:
                if float(pathwayERScore[1]) >= enrichmentScorePerPathway[pathwayGeneSetIndex][1]:
                    p_value += 1
            p_value /= len(pathwayERScores)
            enrichmentScorePerPathway[pathwayGeneSetIndex].append(p_value)
        else:
            p_value = 0
            for pathwayERScore in pathwayERScores:
                if float(pathwayERScore[1]) <= enrichmentScorePerPathway[pathwayGeneSetIndex][1]:
                    p_value += 1
            p_value /= len(pathwayERScores)
            enrichmentScorePerPathway[pathwayGeneSetIndex].append(p_value)

    return enrichmentScorePerPathway


def fetchPhenotypes(phenotypesFilePath):
    with open(phenotypesFilePath, 'r') as file:
        lines = file.readlines()
        phenotypes = lines[2].strip().split()
        return phenotypes

if __name__ == '__main__':
    # Input files path
    pathwaysGeneSetsFilePath = "Gene Expressions Datasets/Pathways Gene Sets/c2.all.v2023.2.Hs.json"
    geneExpressionsDatasetFilePath = "Gene Expressions Datasets/Lymphoblastoid Gender/Gender_collapsed_symbols.gct"
    phenotypesFilePath = "Gene Expressions Datasets/Lymphoblastoid Gender/Gender.cls"

    # Fetch data from files
    pathwayGeneSets = fetchPathwayGeneSetsFromJSON(pathwaysGeneSetsFilePath)
    geneExprValues, geneNames = fetchGeneExprFromGCTFile(geneExpressionsDatasetFilePath)
    phenotypes = fetchPhenotypes(phenotypesFilePath)

    # Rank in correlation to phenotypes
    phenotypesCount = np.unique(phenotypes, return_counts=True)
    rankedGeneExprValues, rankedGeneNames = rankGeneExpr(geneExprValues, geneNames, phenotypesCount)

    # Enrichment score
    p = 1
    enrichmentScorePerPathway = enrichmentScore(rankedGeneExprValues, rankedGeneNames, pathwayGeneSets, p)
    # Estimating significance
    enrichmentScorePerPathwayWithPValue = estimateSignificance(enrichmentScorePerPathway, geneExprValues, geneNames, phenotypesCount, pathwayGeneSets, p)

    # MHT
