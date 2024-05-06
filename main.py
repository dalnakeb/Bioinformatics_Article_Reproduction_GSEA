import argparse
import math
import random
import time
from copy import copy

import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import json
from icecream import ic
import plotly.graph_objs as go

def parseArgs():
    parser = argparse.ArgumentParser(
        description="Process pathway gene sets, gene expressions dataset, and phenotypes file paths.")
    parser.add_argument("SDictFilePath", type=str, help="Path to the pathway gene sets JSON file")
    parser.add_argument("DFilePath", type=str, help="Path to the gene expressions dataset file")
    parser.add_argument("PFilePath", type=str, help="Path to the phenotypes file")
    args = parser.parse_args()

    return parser, args


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
    return np.array(D), np.array(DNames)


def fetchP(PFilePath):
    with open(PFilePath, 'r') as file:
        lines = file.readlines()
        phenotypes = lines[2].strip().split()
        print(f"{len(phenotypes)} phenotypes have been fetched")
        return np.array(sorted(phenotypes))


def sortD(D, DNames, P, C):
    #  Extract indices of the chosen phenotype
    uniqueP = list(set(P))
    PIndices = np.where(P == sorted(uniqueP)[C])[0]

   #  Sort according to one prob for each gene of the chosen phenotype
    DProbes = D[:, np.random.choice(PIndices)]
    indices = np.argsort(DProbes) 
    L = D[indices][::-1]
    LProbes = DProbes[indices][::-1]
    LNames = DNames[indices][::-1]

    return L, LProbes, LNames


def computeES(S, LProbes, LNames, p, withRunningSum=False):
    NH = len(S)
    N = len(LProbes)
    LInS = np.isin(LNames, S)
    NR = np.sum(abs(LProbes[np.where(LInS)])**p)
    ES = 0

    absMaxDev = -math.inf
    PHits = 0
    PMisses = 0
    if withRunningSum:
        runningSum = np.zeros(N)
    for i in range(N):
        if LInS[i]:
            PHits += ((LProbes[i])**p) / NR
        else:
            PMisses += 1/(N-NH)

        if withRunningSum:
            runningSum[i] = PHits-PMisses

        if absMaxDev < abs(PHits-PMisses):
            absMaxDev = abs(PHits-PMisses)
            ES = PHits-PMisses

    if withRunningSum:
        return ES, runningSum

    return ES


def shuffleD(D, P, PIndices):
    PIndices = np.random.permutation(PIndices)
    P = P[PIndices]
    D = D[:, PIndices]
    return D, P, PIndices


def computePValue(ES, S, D, DNames, P, C, p):
    D = np.array(copy(D))
    P = np.array(copy(P))
    PIndices = list(range(len(P)))

    permutations = 1000
    ESNullDistribution = np.zeros(permutations)
    for i in range(permutations):
        D, P, PIndices = shuffleD(D, P, PIndices)
        L, LProbes, LNames = sortD(D, DNames, P, C)
        ESNullDistribution[i] = computeES(S, LProbes, LNames, p)

        if i % 100 == 0:
            print(i)

    if ES >= 0:
        PValue = len([0 <= ESNullValue <= ES for ESNullValue in ESNullDistribution]) / len(ESNullDistribution[ESNullDistribution>=0])
    else:
        PValue = len([ESNullValue >= ES for ESNullValue in ESNullDistribution]) / len(ESNullDistribution[ESNullDistribution<0])

    return PValue, np.sort(ESNullDistribution)


def plotRunningSum(runningSum, xLabel, yLabel, title):
    plt.plot(runningSum)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.title(title)
    plt.grid(True)
    plt.show()


def normalizeES(ES, ESNullDistribution):
    meanES = np.mean(ESNullDistribution)
    NES = ES/meanES
    NESNullDistribution = ESNullDistribution/meanES

    return NES, np.sort(NESNullDistribution)


def computeQValue(NES, NESNullDistribution):
    if NES >= 0:
        QValue = len([0 <= NESNullValue <= NES for NESNullValue in NESNullDistribution]) / len(NESNullDistribution[NESNullDistribution>=0])
    else:
        QValue = len([NESNullValue >= NES for NESNullValue in NESNullDistribution]) / len(NESNullDistribution[NESNullDistribution<0])

    return QValue


if __name__ == '__main__':
    # Input
    #parser, args = parseArgs()
    SDictFilePath = "Gene Expressions Datasets\Pathways Gene Sets\c1.all.v2023.2.Hs.json"
    DFilePath = "Gene Expressions Datasets\Lymphoblastoid Gender\Gender_collapsed_symbols.gct"
    PFilePath = "Gene Expressions Datasets\Lymphoblastoid Gender\Gender.cls"
    C = 1  # Class of distinction
    pathway = "chrYp11"  # pathway gene set name
    p = 1  # Weight of step

    # Fetch data from files
    SDict = fetchSDict(SDictFilePath)  # Dict {pathway name: gene set values}
    D, DNames = fetchD(DFilePath)  # Gene expressions set of multiple samples each (values and names)
    P = fetchP(PFilePath)  # Phenotypes of samples in their respected order
    S = np.array(list(SDict[pathway]))  # Pathway gene set values

    # Rank the gene set expression values.
    L, LProbes, LNames = sortD(D, DNames, P, C)  # Ranked gene set expression values, chosen probe and gene names

    # Compute the ES
    ES, runningSum = computeES(S, LProbes, LNames, p, withRunningSum=True)
    ic(ES)
    plotRunningSum(runningSum, "Gene", "Value", "Running sum")

    #  Estimate significance by computing the p-value
    PValue, ESNullDistribution = computePValue(ES, S, D, DNames, P, C, p)

    #  Multi-hypothesis testing
    #  Normalize ES and ERNull
    NES, NESNullDistribution = normalizeES(ES, ESNullDistribution)
    ic(NES)
    #  Compute q-value
    QValue = computeQValue(NES, NESNullDistribution)
