##import neccesary libraries
import pandas as pd
import pyranges as pr
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import operator
import scipy
from scipy.cluster.hierarchy import dendrogram, leaves_list
import scipy.spatial.distance as ssd
import os, glob
import sys
import seaborn as sns
import statistics
import math
import argparse

##Primary function for associating ATAC and GWAS
def gaiaAssociation(atacLocation, gwasLocation, chromosomeSize, outputLocation, uniqueCount=0, lociCutoff=0):
    
    ##Bring in our atac files for each cell type
    dataFrameList = []
    cellNames = []
    for filename in glob.glob(atacLocation + "/" + '*.txt'):
        with open(os.path.join(os.getcwd(), filename), 'r') as f: # open in readonly mode
            loopFrame = pd.read_csv(filename, sep="\t")
            if "chr1" in loopFrame.columns[0] or "Chr1" in loopFrame.columns[0]:
                columnFrame = loopFrame.columns.to_frame().T
                loopFrame = pd.concat([columnFrame,loopFrame], ignore_index=True)
            loopFrame = loopFrame.rename({loopFrame.columns[0] : 'Chromosome', loopFrame.columns[1] : 'Start', loopFrame.columns[2] : "End"}, axis=1)
            
            loopFrame["Start"] = pd.to_numeric(loopFrame["Start"])
            loopFrame["End"] = pd.to_numeric(loopFrame["End"])
            
            loopFrame['Chromosome'] = loopFrame['Chromosome'].str.strip()
            loopFrame['Chromosome'] = loopFrame['Chromosome'].astype(str)
            
            dataFrameList.append(loopFrame)
        
        ##Save cell type names by formatting the file names
        s = filename[:-4]
        cellNames.append(s[s.rindex('/')+1:])
    
    ##do some light formating to ensure all atac-seq have the same labeling
    for i in range(len(dataFrameList)):
        if "seqnames" in dataFrameList[i].columns:
            dataFrameList[i] = dataFrameList[i].rename(columns={"seqnames": "chr"})
        if "chr" in dataFrameList[i].columns:
            dataFrameList[i] = dataFrameList[i].rename(columns={"chr": "Chromosome"})
        if "start" in dataFrameList[i].columns:
            dataFrameList[i] = dataFrameList[i].rename(columns={"start": "Start"})
        if "end" in dataFrameList[i].columns:
            dataFrameList[i] = dataFrameList[i].rename(columns={"end": "End"})
        
        ##ensure there is a size column
        dataFrameList[i]["size"] = dataFrameList[i]["End"] - dataFrameList[i]["Start"]
        dataFrameList[i]["Size"] = dataFrameList[i]["End"] - dataFrameList[i]["Start"]
            
    ##create an empty matrix to store our associations in
    matrix1 = np.zeros((len(dataFrameList),len(dataFrameList)))

    print("Formatting ATAC-seq: ")
    
    ##create list to store our pyRanges objects in, merge them to ensure regions arent double represented
    prRanges = []
    for i in range(len(dataFrameList)):
        
        ##this particular bit of circular coding creates a pyrange, merges itself to remove redundancy, turns it back into a dataframe, creates a size column, and then turns it back into a pyrange. This is because pyrange doesnt have a function to add a size column and merging removes the size column, this tedium is their fault, but it isn't super time intensive
        tempRange = pr.PyRanges(dataFrameList[i])
        loopRange = tempRange.merge(strand=False)
        df = pd.DataFrame(list(zip(loopRange.Chromosome, loopRange.Start, loopRange.End)),columns =['Chromosome', 'Start', 'End'])
        df["size"] = df["End"] - df["Start"]
        testRange = pr.PyRanges(df)
        prRanges.append(testRange)
        
    print("Comparing ATAC-seq between cell types: ")
    
    if uniqueCount != 0:
        
        ## Create unique sets of ATAC based on a cutoff value of number of overlaps
        overlapCutoff = int(uniqueCount)

        prRangesUnique = [] 
        prFrames = []

        print("Finding Unique ATAC peaks between cell types:")
        for j in range(len(prRanges)):
            prFrames.append(prRanges[j].df)

        for j in range(len(prFrames)):

            if( j != len(prFrames) - 1):
                print("Cell Type " + str(j+1), end = " - ")
                testFrame = pd.concat(prFrames[:j] + prFrames[j+1 :])
            else:
                print("Cell Type " + str(j+1), end = "\n") 
                testFrame = pd.concat(prFrames[:j])


            allButMain = pr.PyRanges(testFrame)
            loopRange = prRanges[j].coverage(allButMain)
            mainFrame = loopRange.df
            mainFrame = mainFrame[mainFrame["NumberOverlaps"] <= overlapCutoff]
            mainFrame = mainFrame.drop(columns=['NumberOverlaps', 'FractionOverlaps'])
            mainRange = pr.PyRanges(mainFrame)
            prRangesUnique.append(mainRange)

        prRanges = prRangesUnique
    
    ##this function loops through each cell type and compares it to every cell type with a lower index number than itself. the comparison is a simple overlap calculation, it then adds these values to its respective matrix
    for i in range(len(dataFrameList)):
        if( i != len(dataFrameList) - 1):
            print("Cell Type " + str(i+1), end = " - ") 
        else:
            print("Cell Type " + str(i+1), end = "\n") 
        
        for j in range(len(dataFrameList)):
            
            if i > j:
                
                loopRange = prRanges[i].intersect(prRanges[j])
                if loopRange.df.empty:
                    matrix1[i,j] = 1
                    matrix1[j,i] = 1
                else:
                    overlap = sum(loopRange.End) - sum(loopRange.Start)
                    sizeA = sum(prRanges[i].size)
                    sizeB = sum(prRanges[j].size)

                    weightValue = 1- ((overlap)/ (sizeA) * (overlap)/(sizeB))

                    matrix1[i,j] = weightValue

                    matrix1[j,i] = weightValue
    
    ##Create the first figure, the dendrogram of ATAC-seq datasets
    distArray = ssd.squareform(matrix1) 
    plt.figure(figsize=(10, 16))
    Z = scipy.cluster.hierarchy.linkage(distArray, method='single', metric='euclidean')
    dn = dendrogram(Z, labels=cellNames, orientation = "left")
    plt.savefig(outputLocation + '/dendrogram.pdf', bbox_inches = "tight")
    
    ##Bring in our GWAS sets
    gwasFrame = 0
    for filename in glob.glob(gwasLocation + "/" + '*.tsv'):
        with open(os.path.join(os.getcwd(), filename), 'r') as f: # open in readonly mode
            loopFrame = pd.read_csv(filename, sep="\t", low_memory=False)
            if not isinstance(gwasFrame, pd.DataFrame):
                gwasFrame = loopFrame
            else:
                gwasFrame = gwasFrame.append(loopFrame)
                
    if "Start" in gwasFrame.columns:
        gwasFrame["CHR_POS"] = gwasFrame["Start"].astype(int)
    if "Reference_name" in gwasFrame.columns:
        gwasFrame["Chromosome"] = gwasFrame["Reference_name"]
    if "Chromosome" in gwasFrame.columns:
        gwasFrame["CHR_ID"] = gwasFrame["Chromosome"]
        
    
    ##format our GWAS frame by getting rid of NAN values and then getting rid of studies with non-numerical chromosome positions, if you want to fix these chromosome positions manually, you can
    gwasNoNAN = gwasFrame[(gwasFrame["CHR_POS"].notna()) & (gwasFrame["DISEASE/TRAIT"].notna())]
    
    gwasNoNAN = gwasNoNAN.reset_index(drop=True)
    
    
    floatIndexes = []
    nonNumericIndexes = []
    for index, row in gwasNoNAN.iterrows():
        if isinstance(gwasNoNAN.at[index, "CHR_POS"], float):
            floatIndexes.append(index)
            gwasNoNAN.at[index, "CHR_POS"] = int(gwasNoNAN.at[index, "CHR_POS"])
        elif isinstance(gwasNoNAN.at[index, "CHR_POS"], str):
            if gwasNoNAN.at[index, "CHR_POS"].isdigit():
                gwasNoNAN.at[index, "CHR_POS"] = int(gwasNoNAN.at[index, "CHR_POS"])
            else:
                nonNumericIndexes.append(index)
        elif not np.issubdtype(gwasNoNAN.at[index, "CHR_POS"], int):
            nonNumericIndexes.append(index)
        
    gwasFormatted = gwasNoNAN.drop(index=nonNumericIndexes)
    
    ###
    
    ###Temporary Formating to check against only loci with a count greater than our loci cutoff
    
    ###
    
    gwasFormatted = gwasFormatted.rename(columns={"DISEASE/TRAIT": "DT"})
    
    vc = gwasFormatted.DT.value_counts()
    highCount = list(vc[vc > int(lociCutoff)].index)
    gwasFormatted = gwasFormatted[(gwasFormatted["DT"].isin(highCount))]
    
    
    ###
    
    ###
    
    ###
    
    ##Create columns for turning into pyranges
    gwasFormatted["Start"] = gwasFormatted["CHR_POS"]
    if "End" not in gwasFrame.columns:
        gwasFormatted["End"] = gwasFormatted["CHR_POS"]
    if "Chromosome" not in gwasFrame.columns:
        gwasFormatted["Chromosome"] = "chr" + gwasFormatted["CHR_ID"].astype(str)
    
    ##save pyranges for each disease/trait
    gwasPyranges = []
    for item in set(gwasFormatted["DT"]):
        
        loopFrame = gwasFormatted[(gwasFormatted["DT"] == item)]
        
        tempRange = pr.PyRanges(loopFrame)
        loopRange = tempRange
        df = pd.DataFrame(list(zip(loopRange.Chromosome, loopRange.Start, loopRange.End)),columns =['Chromosome', 'Start', 'End'])
        df["size"] = df["End"] - df["Start"]
        testRange = pr.PyRanges(df)
        gwasPyranges.append(testRange)
        
        
    ##Create a reordered list of our ATAC-sets based on how they showed up in the dendrogram
    reorderPyranges = []
    reorderDataframeList = []
    for item in reversed(leaves_list(Z)):
        reorderPyranges.append(prRanges[item])
        reorderDataframeList.append(dataFrameList[item])
    
    
    heatmapMatrix = np.zeros((len(reorderPyranges),len(gwasPyranges)))
    
    gwasNames = list(set(gwasFormatted["DT"]))
    
    
    
    ## create the chromosome ratios for every atac-seq set
    
    chromSize = pd.read_csv(chromosomeSize, thousands=',')  
    chromSize = chromSize.sort_values(by=['Chromosome'])
    ATACratios = []
    accpetableChrom = list(chromSize['Chromosome'])
    for item in reorderPyranges:

        itemLoop = item.df
        itemTemp = itemLoop[itemLoop['Chromosome'].isin(accpetableChrom)]
        itemTemp = itemTemp.rename(columns={"size": "Size"})
        itemTemp = itemTemp.rename(columns={"Size": "sizeATAC"})
        
        dff = itemTemp.groupby(["Chromosome"]).sizeATAC.sum().reset_index()
        dff = dff.merge(chromSize, left_on='Chromosome', right_on='Chromosome')
        dff = dff.rename(columns={"size in bp": "sizeChrom"})
        dff["ratio"] = dff["sizeATAC"]/dff["sizeChrom"]
        
        dff["ratio"] = dff["sizeATAC"]/dff["sizeChrom"]
        
        ATACratios.append(dff)
        
        
    ##Create weight matrix
    print("Creating GWAS Weight Comparisons: ")
    for i in range(len(gwasPyranges)):

        tempCountFrame = pd.DataFrame(gwasPyranges[i].df.Chromosome.value_counts())
        tempCountFrame = tempCountFrame.reset_index(level=0)
        tempCountFrame = tempCountFrame.sort_values(by=['index'])
        tempCountFrame = tempCountFrame.rename(columns={"Chromosome": "count", "index": "Chromosome"})

        if( i != len(gwasPyranges) - 1):
            print("GWAS: " + gwasNames[i], end = " - ") 
        else:
            print("GWAS: " + gwasNames[i], end = "\n")

        for j in range(len(reorderPyranges)):

            finalCountFrame = tempCountFrame.merge(ATACratios[j], left_on='Chromosome', right_on='Chromosome')
            expectedCount = sum(finalCountFrame["count"] * finalCountFrame["ratio"])



            loopRange = reorderPyranges[j].coverage(gwasPyranges[i])

            overlaps = sum(loopRange.NumberOverlaps)

            if expectedCount != 0:
                zscoreVal = statistics.NormalDist(mu=expectedCount, sigma=math.sqrt(expectedCount)).zscore(overlaps)
            else:
                zscoreVal = statistics.NormalDist(mu=expectedCount, sigma=1).zscore(overlaps)

            heatmapMatrix[j,i] = zscoreVal

         
        
    print("Highest Z-score: " + str(heatmapMatrix.max()))
        
    if(len(gwasPyranges) > 1 ):

        gwasWeightMatrix = np.zeros((len(gwasPyranges),len(gwasPyranges)))

        ## Creating gwas dendrogram
        print("Creating GWAS dendrogram: ")
        for i in range(len(gwasPyranges)):
            if( i != len(gwasPyranges) - 1):
                print("GWAS Source " + str(i+1), end = " - ") 
            else:
                print("GWAS Source " + str(i+1), end = "\n") 

            for j in range(len(gwasPyranges)):
                if i > j:

                    weightValue = np.linalg.norm(heatmapMatrix[:,i] - heatmapMatrix[:,j])

                    gwasWeightMatrix[i,j] = weightValue

                    gwasWeightMatrix[j,i] = weightValue
                
        ##Create the second figure, the dendrogram of GWAS datasets

        distArray = ssd.squareform(gwasWeightMatrix) 
        plt.figure(figsize=(16, 10))
        Z2 = scipy.cluster.hierarchy.linkage(distArray, method='single', metric='euclidean')
        dn2 = dendrogram(Z2, labels=gwasNames, orientation = "top", leaf_rotation = 90)
        plt.savefig(outputLocation + '/gwasdendrogram.pdf', bbox_inches = "tight")
    
        
    
        newHeatMatrix = heatmapMatrix[:,leaves_list(Z2)]
        ax = sns.heatmap(newHeatMatrix, cmap="YlGnBu")
        plt.savefig(outputLocation + '/heatmap.pdf', bbox_inches = "tight")
    
    else:
        newHeatMatrix = heatmapMatrix
        ax = sns.heatmap(newHeatMatrix, cmap="YlGnBu",annot=True)
        plt.savefig(outputLocation + '/heatmap.pdf', bbox_inches = "tight")
    
    
    figure, axis = plt.subplots(2, 2, figsize=(40, 30))
    sns.heatmap(newHeatMatrix, cmap="YlGnBu", cbar=False, ax=axis[1, 1])
    sns.heatmap(newHeatMatrix, cmap="YlGnBu", annot=True, ax=axis[0, 0])
    if(len(gwasPyranges) > 1 ):
        dendrogram(Z2, labels=gwasNames, orientation = "top", leaf_rotation = 90, ax=axis[0, 1])
    else:
        axis[0, 1].axis('off')
    dendrogram(Z,  labels=cellNames, orientation = "left", ax=axis[1, 0])
    axis[1, 1].margins(x=0)
    plt.savefig(outputLocation + '/complete.pdf', bbox_inches = "tight")
    

def main():
    
    ## parse the given arguments to see if the neccesary libraries have been given and 
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--ATACuniqueness',  default=0, required = False, type=int,
                        help='a cutoff value for ATAC uniqueness (e.g. if given 12, then any atac peak found in more than 12 atac sets will be removed from all of them) - by default uniqueness is not considered')
    parser.add_argument('-l', '--lociCutoff',default=0, required = True, type=int,
                        help='a loci cutoff value, will only consider loci groups (phenotypes or cohorts) with more loci than this cutoff value - by default this cutoff is 0')
    parser.add_argument('-a', '--atac', required = True, 
                        help='the folder location of the atac bed files stored in .txt format')
    parser.add_argument('-g', '--loci', required = True, 
                        help='the folder location of the loci files stored in .tsv format')
    parser.add_argument('-c', '--chrom', required = True, 
                        help='the folder location of the chromsome size files stored in .csv format')
    parser.add_argument('-o', '--output', required = True, 
                        help='the folder location you want the results to be output inxs')

    #(atacLocation, gwasLocation, chromosomeSize, outputLocation):

    args = parser.parse_args()
    
    if not os.path.exists(args.atac):
        print("ATAC folder cannot be found")
    if not os.path.exists(args.loci):
        print("Loci folder cannot be found")
    if not os.path.exists(args.chrom):
        print("Chrom size file cannot be found")
    if not os.path.exists(args.output):
        print("Output folder cannot be found")

    gaiaAssociation(args.atac, args.loci, args.chrom, args.output, args.ATACuniqueness, args.lociCutoff)
