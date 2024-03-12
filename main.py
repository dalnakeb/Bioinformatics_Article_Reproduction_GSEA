import math

import numpy as np
import json
from icecream import ic
def enrichmentScore(rankedGeneExpr, pathwayGeneSets, p):
    ERScorePerPathway = []

    for pathwayGeneSetIndex in range(len(pathwayGeneSets)//10):
        runningSum = 0
        minER = math.inf
        maxER = -math.inf
        ERScorePerPathway.append([pathwayGeneSets[pathwayGeneSetIndex][0]])
        for geneIndex in range(len(rankedGeneExpr)):
            if rankedGeneExpr[geneIndex][0] in pathwayGeneSets[pathwayGeneSetIndex][1]:
                runningSum += abs(rankedGeneExpr[geneIndex][2]**p)
            else:
                runningSum -= rankedGeneExpr[geneIndex][2]**p

            minER = min(minER, runningSum)
            maxER = max(maxER, runningSum)

        if abs(minER) > abs(maxER):
            ERScorePerPathway[pathwayGeneSetIndex].append(minER)
        else:
            ERScorePerPathway[pathwayGeneSetIndex].append(maxER)
    ic(sorted(ERScorePerPathway)[-1])

def rankGeneExpr(geneExpr, phenotypes):
    """
    Rank the gene expression set first according to the metric chosen to rank every line, it's simply by log of
    subtracting the averages of phenotype 2 from 1 (change fold) and ordering the list accordingly in descendant order.
    Then order every line by value and then by order of correlation to the phenotype.

    :param geneExpr: dictionary (gene name, list of expression values)
    :param phenotypes: a list of phenotypes indicating the phenotype of an expression in a certain position
    :return: ranked list according to the metric described above [gene name, expression values, metric],
    permuted phenotype lists that correlates with the expression values.
    """

    # Correlation of genes with phenotype 1
    phenotypes = np.array(phenotypes)
    phenotypesCount = np.unique(phenotypes, return_counts=True)
    for i in range(len(geneExpr)):
        differentialGeneExpr = np.sum(geneExpr[i][1][0:phenotypesCount[1][0]])/phenotypesCount[1][0] - \
                               np.sum(geneExpr[i][1][phenotypesCount[1][1]:])/phenotypesCount[1][1]
        if differentialGeneExpr > 0:
            geneExpr[i].append(math.log(differentialGeneExpr, 2))  # log(phenotype 1-2) fold change
        elif differentialGeneExpr < 0:
            geneExpr[i].append(-math.log(abs(differentialGeneExpr), 2))
        else:
            geneExpr[i].append(0)

    geneExpr.sort(key=lambda x: x[2], reverse=True)

    # Order every gene samples by ascendant value of expression (and the phenotype list)
    phenotypesPerGene = []
    for i in range(len(geneExpr)):
        permutationOrder = np.argsort(geneExpr[i][1])
        geneExpr[i][1] = geneExpr[i][1][permutationOrder]
        phenotypesPerGene.append([])
        phenotypesPerGene[i].append(phenotypes[permutationOrder])
    return geneExpr, phenotypesPerGene


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
        geneExpr = []
        index = 0
        for line in lines[3:]:
            columns = line.split('\t')
            geneName = columns[0].strip("\t")
            if "_at" not in geneName:
                expressionValues = list(map(float, columns[2:]))
                if "///" in geneName:
                    geneName = geneName.split(" ")[0]

                geneExpr.append([])
                geneExpr[index].append(geneName)
                geneExpr[index].append(np.array(expressionValues))
                index += 1
    return geneExpr


def fetchPathwayGeneSetsFromJSON(geneSetsFilePath):
    """
    Extracts gene sets from a json file given the filepath as an input
    :param geneSetsFilePath: filepath to a json file containing the gene sets
    :return: a list of lists of [geneSet name, list of genes]
    """
    with open(geneSetsFilePath, 'r') as file:
        geneSetFileDataC1 = file.read()
    geneSetsJSON = json.loads(geneSetFileDataC1)

    geneSets = []
    for geneSet, items in geneSetsJSON.items():
        getSetName = items["systematicName"]
        genesList = items["geneSymbols"]
        geneSets.append([getSetName, genesList])

    return geneSets


def estimateSignificance(enrichmentScore):
    pass


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
    geneExpr = fetchGeneExprFromGCTFile(geneExpressionsDatasetFilePath)
    phenotypes = fetchPhenotypes(phenotypesFilePath)

    #Rank in correlation to phenotypes
    rankedGeneExpr, phenotypesPerGene = rankGeneExpr(geneExpr, phenotypes)
    p = 1
    enrichmentScore = enrichmentScore(rankedGeneExpr, pathwayGeneSets, p)
    # Estimating significance
    significanceEstimation = estimateSignificance(enrichmentScore)

    # MHT
