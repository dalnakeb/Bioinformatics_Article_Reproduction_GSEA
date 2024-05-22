import math
import random
import numpy as np
import matplotlib.pyplot as plt
import json
import pandas as pd

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
    Extract the list of genes expressions from a gct file given as input the file-path.
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
    """
    :param PFilePath: filepath to the phenotype list (indicating the position of a phenotype sample in the genes
      expressions
    :return: list of phenotypes
    """
    with open(PFilePath, 'r') as file:
        lines = file.readlines()
        phenotypes = lines[2].strip().split()
        print(f"{len(phenotypes)} phenotypes have been fetched")
        return np.array(sorted(phenotypes))


def sortD(D, DNames, P, C):
    """
    This method ranks the gene set expression and the list of gene names according to a certain metric which is the
      Fold change between two randomly picked samples from both phenotypes (in decreasing order). The class of
      distinction determines which phenotype is of interest to compute the fold change from the other phenotype to the
      one of interest.
    :param D: unranked gene set expression values
    :param DNames: unranked genes names in the respective order of the gene set expressions
    :param P: list of phenotypes indicating the phenotype of a sample at a given position in the gene set expressions
    :param C: class of distinction index in the alphabetic order
    :return: (sorted gene set expressions values, sorted fold changes for every gene, sorted gene names)
    """
    #  Extract indices of the chosen phenotype
    uniqueP = list(set(P))  # a set of unique phenotypes
    PIndices = np.where(P == np.sort(uniqueP)[C])  # indices of the phenotype of interest
    notPIndices = np.where(P == np.sort(uniqueP)[(C + 1) % 2])  # indices of the other phenotype

    DFC =  np.mean(D[:, notPIndices[0]], axis=1) / np.mean(D[:, PIndices[0]], axis=1)  # Unranked list of fold changes
    indices = np.argsort(DFC)
    L = D[indices][::-1]  # sorted gene set expressions
    LFC = DFC[indices][::-1]  # sorted list fold changes of
    LNames = DNames[indices][::-1]  # sorted list of gene names

    return L, LFC, LNames


def computeES(S, LFC, LNames, p, plotRS=False, pathway=None, CName=None, LE=False):
    """
    Computes the enrichment score of a pathway given a gene set expressions and a class of distinction. a high positive
      value of the enrichment score indicates an enrichment of the pathway. Whereas a low negative value indicates an
      under-representation of the pathway. Values close to 0 indicate a randomly distributed genes in the gene set
      expressions
    :param S: pathway gene set
    :param LFC: ranked list of fold changes of genes
    :param LNames: ranked list of gene names
    :param p:  weight of step in the running sum
    :param plotRS: boolean to plot the running sum or not
    :param pathway: pathway name
    :param CName: name of phenotype of interest
    :param LE: boolean to print the leading edge of the pathway gene set
    :return: ES: enrichment score
    """
    NH = len(S)  # number of genes in the pathway
    N = len(LFC)  # number of genes in the genes expressions set
    NMinusNH = N-NH  # number of genes that are not in the pathway

    LInS = np.isin(LNames, S)  # list of booleans indicating if a gene is present in the pathway for a given index
    NR = np.sum(abs(LFC[np.where(LInS)])**p)  # sum of fold changes of genes in the pathway
    ES = 0

    absMaxDev = -math.inf
    PHits = 0
    PMisses = 0
    if plotRS:
        runningSum = np.zeros(N)

    for i in range(N):
        if LInS[i]:  # if gene i is present in pathway
            PHits += (abs(LFC[i])**p) / NR
        else:
            PMisses += 1/NMinusNH

        if plotRS:
            runningSum[i] = PHits-PMisses

        if absMaxDev < abs(PHits-PMisses):  # if a greater maximum deviation from 0 is found
            absMaxDev = abs(PHits-PMisses)
            ES = PHits-PMisses

    if plotRS:
        plotRunningSum(runningSum, "Position in L", "Running Enrichment Score", f"{pathway}", ES, CName)
    if LE:
        print(f"Leading Edge: {np.intersect1d(LNames, S)[0:5]}")

    return ES

def computePValue(ES, S, D, DNames, P, C, p):
    """
    Compute the p-value of the enrichment score using the permutation method
    :param ES: Enrichment score
    :param S: Pathway gene set
    :param D: Unranked gene set expressions
    :param DNames: Unranked list of genes names
    :param P: List of phenotypes
    :param C: Class of distinction
    :param p: weight of the running sum step
    :return: p-value
    """
    print(f"Computing P-value...")
    permutations = 1000
    ESNullDistribution = np.zeros(permutations)
    for i in range(permutations):
        PShuffled = np.random.permutation(P)
        L, LFC, LNames = sortD(D, DNames, PShuffled, C)
        ESNullDistribution[i] = computeES(S, LFC, LNames, p)

    if ES >= 0:
        PValue = np.sum((ES <= ESNullDistribution)) / np.sum((0 <= ESNullDistribution))
    else:
        PValue = np.sum((ESNullDistribution <= ES)) / np.sum((ESNullDistribution <= 0))

    return PValue


def plotRunningSum(runningSum, xLabel, yLabel, title, ES, CName):
    plt.plot(runningSum)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.title(title + f"| {CName}")
    plt.grid(True)
    plt.axvline(x=np.where(runningSum == ES)[0], color = "r")
    plt.axhline(y=ES, color="black", linestyle='--')

    plt.annotate(f'ES: {round(ES,4)}',
                 xy=(np.where(runningSum == ES)[0], ES),
                 xytext=(np.where(runningSum == ES)[0]*1.2, ES/1.1),  # Text position
                 arrowprops=dict(facecolor='black', shrink=0.01))

    plt.show()


def plotPathwaysRanking(pathways, C, p, LFC, LNames, P, SDict):
    """
    Plot the ranking of enrichment score for a given list of pathways. A P-value threshold determines the color of the
    bar (red if P-value >0.2, green otherwise)
    :param pathways: list of pathways
    :param C: Class of distinction
    :param p: Weight of step for running sum
    :param LFC: Ranked list of fold changes of genes
    :param LNames: Ranked list of gene names according to fold changes
    :param P: List of phenotypes
    :param SDict: Dictionary containing all the pathways and their associated genes
    :return: None
    """
    PValues = []
    ESs = []
    for pathway in pathways:
        print(pathway)
        S = np.array(list(SDict[pathway]))  # Pathway gene set
        random.seed(2)
        CName = np.sort(list(set(P)))[C]  # name of phenotype of interest (class of distinction)

        # Compute the ES [optional: plot the running sum graph]
        ES = computeES(S, LFC, LNames, p, pathway=pathway, CName=CName, LE=True)
        ESs.append(ES)

        #  Estimate significance by computing the p-value
        PValue = computePValue(ES, S, D, DNames, P, C, p)
        PValues.append(PValue)

    data = {
        "Pathways": pathways,
        "ESs": ESs,
        "P adj <0.2": [PValue < 0.2 for PValue in PValues]
    }

    df = pd.DataFrame(data)

    # Set the color based on the significance
    colors = df["P adj <0.2"].map({True: "green", False: "red"})

    plt.figure(figsize=(10, 8))
    plt.barh(df["Pathways"], df["ESs"], color=colors)

    red = plt.Line2D([0], [0], color='red', lw=4, label='False')
    green = plt.Line2D([0], [0], color='green', lw=4, label='True')
    plt.legend(handles=[red, green], title='P adj <0.2')

    plt.title(f"Gene Set Enrichment: {CName}")
    plt.xlabel("ES")
    plt.ylabel("Pathway")

    plt.show()


if __name__ == '__main__':
    # Input
    SDictFilePath = r"gene_Expressions_Datasets\pathways_Gene_Sets\c2.all.v2023.2.Hs.json"  # Genes expression file path
#    DFilePath = r"gene_Expressions_Datasets\lymphoblastoid_Gender\Gender_collapsed_symbols.gct"  # Pathways get sets filepath
 #   PFilePath = r"gene_Expressions_Datasets\lymphoblastoid_Gender\Gender.cls"  # Phenotypes filepath

    DFilePath = r"gene_Expressions_Datasets\diabetes\Diabetes_collapsed_symbols.gct"  # Pathways get sets filepath
    PFilePath = r"gene_Expressions_Datasets\diabetes\Diabetes.cls"  # Phenotypes filepath

    C = 0  # Class of distinction (in the alphabetic order)
    p = 1  # Weight of step of the running sum (0 for Kolmogorov–Smirnov statistic)
    plotAllPathways = True

    # Fetch data from files
    SDict = fetchSDict(SDictFilePath)  # Dict {pathway name: gene set values}
    D, DNames = fetchD(DFilePath)  # Gene expressions set of multiple samples for each phenotype (values and names)
    P = fetchP(PFilePath)  # Phenotypes of samples in their respected order

    random.seed(2)
    # Rank the genes set expression according to their fold change.
    L, LFC, LNames = sortD(D, DNames, P, C)  # Ranked gene set expression values, Fold changes and gene names

    if not plotAllPathways:
        pathway = "chrXp11"  # Pathway gene set name
        S = np.array(list(SDict[pathway]))  # Pathway gene set
        CName = np.sort(list(set(P)))[C]  # name of phenotype of interest (class of distinction)
        print(f"Class of distinction: {CName}")
        # Compute the ES [optional: plot the running sum graph]
        ES = computeES(S, LFC, LNames, p, plotRS=True, pathway=pathway, CName=CName, LE=True)
        print(f"ES: {ES}")

        #  Estimate significance by computing the p-value
        PValue = computePValue(ES, S, D, DNames, P, C, p)
        print(f"PValue: {PValue}")
    else:
        # Plot the ranking of the list of given pathways
        #pathways = ["chrYp11", "chrXp21", "chrYq11", "chrXq13", "chrXp11", "chrXp22", "chrXq11", "chrXq12"]
        pathways = ["REACTOME_INSULIN_PROCESSING", "REACTOME_INSULIN_RECEPTOR_RECYCLING", "REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE",
                    "REACTOME_REGULATION_OF_INSULIN_SECRETION"]

        plotPathwaysRanking(pathways, C, p, LFC, LNames, P, SDict)

