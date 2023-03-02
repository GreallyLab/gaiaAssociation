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
def gaiaAssociation(atacLocation, gwasLocation, chromosomeSize, outputLocation, uniqueCount=0, lociCutoff=0, lociSelection = 0, subsettingRegion=0, saveOverlap = False, zscoreRange =0):
        
    ## ensure all the given locations and files are accesible and real
    if not os.path.exists(atacLocation):
        sys.exit("ATAC folder cannot be found")
    if not os.path.exists(gwasLocation):
        sys.exit("Loci folder cannot be found")
    if not os.path.exists(chromosomeSize):
        sys.exit("Chrom size file cannot be found")
    if not os.path.exists(outputLocation):
        print("Output folder cannot be found, attempting to build it:")
        os.mkdir(outputLocation)
        if not os.path.exists(outputLocation):
            sys.exit("Output folder cannot be built")
        else:
            print("Output folder succesfully built:")
    
    ##Bring in our atac files for each cell type
    dataFrameList = []
    cellNames = []
    for filename in glob.glob(atacLocation + "/" + '*.txt'):
        
        loopFrame = pd.read_csv(filename, sep="\t")
        
        ##this checks if there was no header, if the first column name is chr1, it shifts columns down into a row
        if "chr1" in loopFrame.columns[0] or "Chr1" in loopFrame.columns[0]:
            columnFrame = loopFrame.columns.to_frame().T
            loopFrame = pd.concat([columnFrame,loopFrame], ignore_index=True)
        
        ##Gives first three columns the names 'Chromosome', 'Start', "End"
        loopFrame = loopFrame.rename({loopFrame.columns[0] : 'Chromosome', loopFrame.columns[1] : 'Start', loopFrame.columns[2] : "End"}, axis=1)
        
        ##ensure start and stop are numeric values
        loopFrame["Start"] = pd.to_numeric(loopFrame["Start"])
        loopFrame["End"] = pd.to_numeric(loopFrame["End"])
        
        ##ensure the chromosome id doesnt have white space and is the string datatype
        loopFrame['Chromosome'] = loopFrame['Chromosome'].str.strip()
        loopFrame['Chromosome'] = loopFrame['Chromosome'].astype(str)
        
        ##add this dataframe to our list of atac-sets
        dataFrameList.append(loopFrame)
    
        ##Save cell type names by formatting the file names
        s = filename[:-4]
        cellNames.append(s[s.rindex('/')+1:])
    
    if(len(dataFrameList)==0):
        sys.exit("No ATAC files found")
    
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
        
        ##ensure there is a size column and a Size column, this redundancy exists because of different package dependencies request for formatting
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
        df["Size"] = df["End"] - df["Start"]
        testRange = pr.PyRanges(df)
        prRanges.append(testRange)
    
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
        
        
    ## If a merging region file is included then subset our ATAC sets using this region set
    if subsettingRegion != 0:
    
        print("Subsetting ATAC peaks based on given region file:")
    
        geneLocation = subsettingRegion
        geneFrame = pd.read_csv(geneLocation, sep="\t")
    
        geneFrame = geneFrame.rename(columns={"Chr": "Chromosome", "start": "Start", "end" : "End"})

        geneRange = pr.PyRanges(geneFrame)
        
        if "seqnames" in geneFrame.columns:
            geneFrame = geneFrame.rename(columns={"seqnames": "chr"})
        if "chr" in geneFrame.columns:
            geneFrame = geneFrame.rename(columns={"chr": "Chromosome"})
        if "start" in geneFrame.columns:
            geneFrame = geneFrame.rename(columns={"start": "Start"})
        if "end" in geneFrame.columns:
            geneFrame = geneFrame.rename(columns={"end": "End"})
    
        overlappedForms = []
        overlappedRanges = []
        for item in prRanges:
            loopRange = geneRange.intersect(item)
            loopRange = loopRange.merge(strand=False)
            loopFrame = loopRange.df
    
            loopFrame["Size"] = loopFrame["End"] - loopFrame["Start"]
            
            loopRange = pr.PyRanges(loopFrame)
            
            overlappedForms.append(loopFrame)
            overlappedRanges.append(loopRange)
            
        prRanges = overlappedRanges
        dataFrameList = overlappedForms
        
        
    print("Comparing ATAC-seq between cell types: ")
    
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
                    sizeA = sum(prRanges[i].Size)
                    sizeB = sum(prRanges[j].Size)

                    weightValue = 1- ((overlap)/ (sizeA) * (overlap)/(sizeB))

                    matrix1[i,j] = weightValue

                    matrix1[j,i] = weightValue
    
    
    ## If the user only included one atac-seq we need to skip steps which would fail with only atac, this includes making a dendrogram using scipy
    if len(dataFrameList) > 1:
        ##Create the first figure, the dendrogram of ATAC-seq datasets
        distArray = ssd.squareform(matrix1)
        plt.figure(figsize=(10, 16))
        Z = scipy.cluster.hierarchy.linkage(distArray, method='single', metric='euclidean')
        dn = dendrogram(Z, labels=cellNames, orientation = "left")
        plt.savefig(outputLocation + '/atac_dendrogram.pdf', bbox_inches = "tight")
    else:
        print("You have only included one ATAC-seq set, this will omit the ATAC dendrogram")
    
    print("Formatting loci:")
    
    ##Bring in our Loci sets
    gwasFrame = 0
    
    ##add GWAS in tsv format
    for filename in glob.glob(gwasLocation + "/" + '*.tsv'):
        loopFrame = pd.read_csv(filename, sep="\t", low_memory=False)
        
        ## give a backup id by file count, in case they dont have an identity column
        s = filename[:-4]
        loopFrame["fileId"] =  (s[s.rindex('/')+1:])
        
        ##These if statements are to account for various formatting mistakes
        if "DISEASE/TRAIT" not in loopFrame.columns:
            loopFrame["DISEASE/TRAIT"] = loopFrame["fileId"]
        if "Start" in loopFrame.columns:
            loopFrame["CHR_POS"] = loopFrame["Start"].astype(int)
        if "start" in loopFrame.columns:
            loopFrame["CHR_POS"] = loopFrame["start"].astype(int)
        if "Reference_name" in loopFrame.columns:
            loopFrame["Chromosome"] = loopFrame["Reference_name"]
        if "Chromosome" in loopFrame.columns:
            loopFrame["CHR_ID"] = loopFrame["Chromosome"]
        if "chr" in loopFrame.columns:
            loopFrame["CHR_ID"] = loopFrame["chr"]
            loopFrame["Chromosome"] = loopFrame["chr"]
        
        if not isinstance(gwasFrame, pd.DataFrame):
            gwasFrame = loopFrame
        else:
            gwasFrame = pd.concat([gwasFrame, loopFrame], axis=0, ignore_index=True)
            
    ## add GWAS in csv format
    for filename in glob.glob(gwasLocation + "/" + '*.csv'):
    
            loopFrame = pd.read_csv(filename, low_memory=False)
            
            ## give a backup id by file count, in case they dont have an identity column
            s = filename[:-4]
            loopFrame["fileId"] =  (s[s.rindex('/')+1:])
            
            ##These if statements are to account for various formatting mistakes
            if "DISEASE/TRAIT" not in loopFrame.columns:
                loopFrame["DISEASE/TRAIT"] = loopFrame["fileId"]
            if "Start" in loopFrame.columns:
                loopFrame["CHR_POS"] = loopFrame["Start"].astype(int)
            if "start" in loopFrame.columns:
                loopFrame["CHR_POS"] = loopFrame["start"].astype(int)
            if "Reference_name" in loopFrame.columns:
                loopFrame["Chromosome"] = loopFrame["Reference_name"]
            if "Chromosome" in loopFrame.columns:
                loopFrame["CHR_ID"] = loopFrame["Chromosome"]
            if "chr" in loopFrame.columns:
                loopFrame["CHR_ID"] = loopFrame["chr"]
                loopFrame["Chromosome"] = loopFrame["chr"]
            
            if not isinstance(gwasFrame, pd.DataFrame):
                gwasFrame = loopFrame
            else:
                gwasFrame = pd.concat([gwasFrame, loopFrame],axis=0, ignore_index=True)
                
    if not isinstance(gwasFrame, pd.DataFrame):
        sys.exit("No Loci files found")
    
    ##format our Loci frame by getting rid of NAN values and then getting rid of studies with non-numerical chromosome positions, if you want to fix these chromosome positions manually, you can
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
    
    ### rename the Disease/TRait column to more easily callable DT column
    gwasFormatted = gwasFormatted.rename(columns={"DISEASE/TRAIT": "DT"})
    
    ## subset loci based on cutoff value
    vc = gwasFormatted.DT.value_counts()
    if lociCutoff != 0:
        print("Subsetting loci based on user-defined count:")
    highCount = list(vc[vc > int(lociCutoff)].index)
    gwasFormatted = gwasFormatted[(gwasFormatted["DT"].isin(highCount))]
    
    
    ## Take only the user selected loci if the user has given a list
    if lociSelection != 0:
        print("Subsetting loci based on user-defined list:")
        lociTypeList = []
        with open(lociSelection) as f:
            for line in f:
                lociTypeList.append(line)
                lociTypeList = [loci.rstrip('\n') for loci in lociTypeList]
    
        gwasFormatted = gwasFormatted[(gwasFormatted["DT"].isin(lociTypeList))]
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
        df["Size"] = df["End"] - df["Start"]
        testRange = pr.PyRanges(df)
        gwasPyranges.append(testRange)
        
        
    ##Create a reordered list of our ATAC-sets based on how they showed up in the dendrogram. This step will be skipped if there is only one ATAC-seq set
    reorderPyranges = []
    reorderDataframeList = []
    reorderCellNames = []
    if len(dataFrameList) > 1:
        ## using the leaves from our scipy Z set, resort your pyranges and dataframes to be in the correct order
        for item in reversed(leaves_list(Z)):
            reorderPyranges.append(prRanges[item])
            reorderDataframeList.append(dataFrameList[item])
            reorderCellNames.append(cellNames[item])
    else:
        reorderPyranges = prRanges
        reorderDataframeList = dataFrameList
    
    
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
        itemTemp = itemTemp.rename(columns={"Size": "sizeATAC"})
        
        dff = itemTemp.groupby(["Chromosome"]).sizeATAC.sum().reset_index()
        dff = dff.merge(chromSize, left_on='Chromosome', right_on='Chromosome')
        dff = dff.rename(columns={"size in bp": "sizeChrom"})
        dff["ratio"] = dff["sizeATAC"]/dff["sizeChrom"]
        
        dff["ratio"] = dff["sizeATAC"]/dff["sizeChrom"]
        
        ATACratios.append(dff)
        
    ##if saveOverlap flag is sent, then make sure we have a subfolder to save them to
    if saveOverlap != False:
        if not os.path.exists(outputLocation + '/Overlaps'):
            os.mkdir(outputLocation + '/Overlaps')
        
    ##Create weight matrix
    print("Calculating Loci Z-scores: ")
    for i in range(len(gwasPyranges)):

        ## create a temporary dataframe where we store the number of counts in this loci set per chromosome
        tempCountFrame = pd.DataFrame(gwasPyranges[i].df.Chromosome.value_counts())
        tempCountFrame = tempCountFrame.reset_index(level=0)
        tempCountFrame = tempCountFrame.sort_values(by=['index'])
        tempCountFrame = tempCountFrame.rename(columns={"Chromosome": "count", "index": "Chromosome"})

        ##print a label for each itteration
        if( i != len(gwasPyranges) - 1):
            print("Loci: " + gwasNames[i], end = " - ")
        else:
            print("Loci: " + gwasNames[i], end = "\n")

        ## Compare this Loci set against every atac set
        for j in range(len(reorderPyranges)):

            ## Combine out our counts with the atac ratios we made earlier
            finalCountFrame = tempCountFrame.merge(ATACratios[j], left_on='Chromosome', right_on='Chromosome')
            
            ## find how many overlaps we would expect based purely on chance, that is: percentage of chromosome covered in open regions * number of counts in that chromosome. Do this for every chromosome and then sum it
            expectedCount = sum(finalCountFrame["count"] * finalCountFrame["ratio"])

            ## find how many overlaps there actually were by comparing our loci set to this atac set
            loopRange = reorderPyranges[j].coverage(gwasPyranges[i])

            ## get the sum of those overlaps
            overlaps = sum(loopRange.NumberOverlaps)
            
            ##if saveOverlap flag is sent, then save every single individual overlap file
            if saveOverlap != False:
                saveRange = gwasPyranges[i].coverage(reorderPyranges[j])
                saveFrame = saveRange.df
                saveFrame = saveFrame[saveFrame["NumberOverlaps"]!=0]
                saveFrame.to_csv(outputLocation + '/Overlaps/' + reorderCellNames[j] + '_' + gwasNames[i] + '.txt', sep='\t', index=False)
                

            ##we generate a z-score based on a normal distribution centered on expected count with a variance equal to sqrt(mean)
            if expectedCount != 0:
                zscoreVal = statistics.NormalDist(mu=expectedCount, sigma=math.sqrt(expectedCount)).zscore(overlaps)
                
            ##if there are 0 expected overlaps (theoretically impossible edge case, but, rounding errors in computer may make it happen) make a normal distribution centered on 0 with a standard deviation of 1
            else:
                zscoreVal = statistics.NormalDist(mu=expectedCount, sigma=1).zscore(overlaps)

            ## add our z score to the matrix
            heatmapMatrix[j,i] = zscoreVal

         
    try:
        print("Highest Z-score: " + str(heatmapMatrix.max()))
    except ValueError:
        pass
        
    if zscoreRange!=0:
        print("Subsetting Loci by the Z-score Range Cutoff:")
        columnsToRemove = []
        for i in range(len(gwasPyranges)):
            if np.ptp(heatmapMatrix[:,i]) < zscoreRange:
                columnsToRemove.append(i)
        
        heatmapMatrix = np.delete(heatmapMatrix, columnsToRemove, 1)
        print("Original Count of Loci Groups: " + str(len(gwasPyranges)))
        print("Amount of Loci Groups removed using Z-Score Filter: " + str(len(columnsToRemove)))
        for index in sorted(columnsToRemove, reverse=True):
            del gwasPyranges[index]
            del gwasNames[index]
    ## If more than one loci set is included, make a dendrogram relating them. This is done by using their z-scores with each ATAC-seq set as a vector, and then comparing the distances between these vectors. This is a simple way of comparing their similiarity.
    if(len(gwasPyranges) > 1 ):

        gwasWeightMatrix = np.zeros((len(gwasPyranges),len(gwasPyranges)))

        ## Creating loci dendrogram
        print("Creating Loci dendrogram: ")
        for i in range(len(gwasPyranges)):
            if( i != len(gwasPyranges) - 1):
                print("Loci Source " + str(i+1), end = " - ")
            else:
                print("Loci Source " + str(i+1), end = "\n")

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
        plt.savefig(outputLocation + '/loci_dendrogram.pdf', bbox_inches = "tight")
    
        
    
        newHeatMatrix = heatmapMatrix[:,leaves_list(Z2)]
        ax = sns.heatmap(newHeatMatrix, cmap="YlGnBu")
        plt.savefig(outputLocation + '/zscore_heatmap.pdf', bbox_inches = "tight")
    
    else:
        
        print("You have only included one loci set, this will omit the loci dendrogram")
    
        newHeatMatrix = heatmapMatrix
        ax = sns.heatmap(newHeatMatrix, cmap="YlGnBu",annot=True)
        plt.savefig(outputLocation + '/zscore_heatmap.pdf', bbox_inches = "tight")
    
    ## Create the final figure
    figure, axis = plt.subplots(2, 2, figsize=(40, 30))
    sns.heatmap(newHeatMatrix, cmap="YlGnBu", cbar=False, ax=axis[1, 1])
    sns.heatmap(newHeatMatrix, cmap="YlGnBu", annot=True, ax=axis[0, 0])
    if(len(gwasPyranges) > 1 ):
        dendrogram(Z2, labels=gwasNames, orientation = "top", leaf_rotation = 90, ax=axis[0, 1])
    else:
        axis[0, 1].text(0.1, 0.5, gwasNames[0], horizontalalignment='center', verticalalignment='bottom', transform=axis[0, 1].transAxes)
        axis[0, 1].grid(False)
        axis[0, 1].axis('off')
    if len(dataFrameList) > 1:
        dendrogram(Z,  labels=cellNames, orientation = "left", ax=axis[1, 0])
    else:
        axis[1, 0].text(0.1, 0.5, cellNames[0], horizontalalignment='right', verticalalignment='center', transform=axis[1, 0].transAxes)
        axis[1, 0].grid(False)
        axis[1, 0].axis('off')
    axis[1, 1].margins(x=0)
    plt.savefig(outputLocation + '/complete_figure.pdf', bbox_inches = "tight")
    
    ##save the matrix of z-scores in case they want to do something else with the formatting
    ##function for sorting by indexes from w3resource
    def sort_by_indexes(lst, indexes, reverse=False):
      return [val for (_, val) in sorted(zip(indexes, lst), key=lambda x: \
              x[0], reverse=reverse)]
    
    newHeatFrame = pd.DataFrame(newHeatMatrix)
    newHeatFrame.columns = sort_by_indexes(gwasNames, leaves_list(Z2))
    newHeatFrame.index = reorderCellNames

    newHeatFrame.to_csv(outputLocation + '/zscore_matrix.txt',sep='\t')

    
def main():
    
    ## parse the given arguments to see if the neccesary libraries have been given and 
    parser = argparse.ArgumentParser(prog ='gaia',
                                     description ='Compare ATAC-seq data to loci.')
    
    parser.add_argument('-a', '--atac', required = True, 
                        help='the folder location of the atac bed files stored in .txt format')
    parser.add_argument('-g', '--loci', required = True, 
                        help='the folder location of the loci files stored in .tsv format')
    parser.add_argument('-c', '--chrom', required = True, 
                        help='the location of the chromsome size file stored in a .csv format')
    parser.add_argument('-o', '--output', required = True, 
                        help='the folder location you want the results to be output into')
    parser.add_argument('-u', '--ATACuniqueness',  default=0, required = False, type=int,
                        help='a cutoff value for ATAC uniqueness (e.g. if given 12, then any atac peak found in more than 12 atac sets will be removed from all of them) - by default uniqueness is not considered')
    parser.add_argument('-l', '--lociCutoff',default=0, required = False, type=int,
                        help='a loci cutoff value, will only consider loci groups (phenotypes or cohorts) with more loci than this cutoff value - by default this cutoff is 0')
    parser.add_argument('-s', '--specificLoci',default=0, required = False,
                        help='a txt file with the specific loci phenotype you would like to use. This can be very helpful if using a large loci set with with many phenotypes, and you want to sort by more than just loci count.' )
    parser.add_argument('-m', '--maskRegion',default=0, required = False,
                        help='a txt file containing a set of regions that you want to subset each ATAC region by, for example a set of regions around the TSSs of genes of interest. This will reduce the ATAC regions to just those that overlap with this set of regions' )
    parser.add_argument('-p', action='store_true',
                        help='a flag which will save the overlapping LOCI for every single cell type for every single loci group as an independant file in a subfolder in the output, this is by default false' )
    parser.add_argument('-z', '--zcoreRangeCutoff',  default=0, required = False, type=int,
                        help='a cutoff value for selecting which loci groups to keep, after calculating the zscore enrichment for each loci for each given atac-seq set, only loci groups with a zscore range (max-min) equal or greater to this cutoff will be kept. This can be useful for subsetting large groups of loci data into only those which have interesting insights, or unique cell enrichments.')

    args = parser.parse_args()

    gaiaAssociation(args.atac, args.loci, args.chrom, args.output, args.ATACuniqueness, args.lociCutoff, args.specificLoci, args.maskRegion, args.p, args.zcoreRangeCutoff)
