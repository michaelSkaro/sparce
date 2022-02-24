import os
import re
import sys
import math
import seaborn as sns
import glob
import io
import math
import os
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import where
import sklearn
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.feature_selection import RFE, SelectFromModel, SelectKBest, chi2
from sklearn.linear_model import LogisticRegression


class bin_manipulations:
    def __inti__(self, genomeAnnotationFile, nFeatures, binSize, nSamples):
        self.genomeAnnotationFile = genomeAnnotationFile
        self.nFeatures = nFeatures
        self.binSize = binSize
        self.nSamples = nSamples
        pass

    def miniBin(genomeAnnotationFile, binSize):
        '''
        :param genomeAnnotationFile: gene annotation file
        :type genomeAnnotationFile: string
        :param binSize: bin size
        :type binSize: int
        :return: binned gene annotation file
        :rtype: pandas dataframe


        Objective: Using a slideing winodw iterate over the gene annotation file and bin the genes based on the bin size. 
        If the number of featues in the bin is greater than the number of samples, then the bin is broken into mini bins of size binSize/10. 
        

        '''



        gff_df = pd.read_table(
            genomeAnnotationFile,
            delimiter="\t",
        ).fillna(0)
        gff_df["bin"] = gff_df["start"].apply(lambda x: math.floor(x/binSize))
        gff_df["bin"] = gff_df["bin"].astype(int)
        gff_df = gff_df.sort_values(["chr", "bin"])
        gff_df.index = range(1, len(gff_df) + 1)

        # iterate over the bins and break them into mini bins if the number of features is greater than the number of samples

        for chr in gff_df["chr"].unique():
            df_chr = gff_df[gff_df["chr"] == chr]
            df_chr = df_chr.sort_values(["start"])
            df_chr.index = range(1, len(df_chr) + 1)
            df_chr_binned = pd.DataFrame()
            for bin in df_chr["bin"].unique():
                df_chr_bin = df_chr[df_chr["bin"] == bin]
                df_chr_bin = df_chr_bin.sort_values(["start"])
                df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                if len(df_chr_bin) > nFeatures:
                    miniBins = prepare_transcriptomics_file_for_SPARCE.chunkIt(df_chr_bin, 10)
                    for miniBin in miniBins:
                        miniBin = miniBin.sort_values(["start"])
                        miniBin.index = range(1, len(miniBin) + 1)
                        df_chr_binned = df_chr_binned.append(miniBin)
                else:
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
            gff_df = gff_df.append(df_chr_binned)
        gff_df.index = range(1, len(gff_df) + 1)
        return gff_df


    def expand_bin(gff_df, binSize):
        '''
        :param gff_df: gene annotation file
        :type gff_df: pandas dataframe
        :param binSize: bin size
        :type binSize: int
        :return: expanded gene annotation file
        :rtype: pandas dataframe


        Objective: if a bin has significantly less features than the average number of features in a bin, expand the binSize for that bin
        until it captures approxiamtely equal features to the average number of values in a bin. Adjust the start and end of the next bin accordingly.
        '''

        gff_df["bin"] = gff_df["start"].apply(lambda x: math.floor(x/binSize))
        gff_df["bin"] = gff_df["bin"].astype(int)
        gff_df = gff_df.sort_values(["chr", "bin"])
        gff_df.index = range(1, len(gff_df) + 1)

        for chr in gff_df["chr"].unique():
            df_chr = gff_df[gff_df["chr"] == chr]
            df_chr = df_chr.sort_values(["start"])
            df_chr.index = range(1, len(df_chr) + 1)
            df_chr_binned = pd.DataFrame()
            for bin in df_chr["bin"].unique():
                df_chr_bin = df_chr[df_chr["bin"] == bin]
                df_chr_bin = df_chr_bin.sort_values(["start"])
                df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                if len(df_chr_bin) < (nFeatures/10):
                    df_chr_bin = df_chr_bin.sort_values(["start"])
                    df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                    df_chr_bin["bin"] = df_chr_bin["start"].apply(lambda x: math.floor(x/binSize))
                    df_chr_bin["bin"] = df_chr_bin["bin"].astype(int)
                    df_chr_bin = df_chr_bin.sort_values(["chr", "bin"])
                    df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
                else:
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
            gff_df = gff_df.append(df_chr_binned)
        gff_df.index = range(1, len(gff_df) + 1)
        
        return gff_df

    def merge_last_bin_on_chromosome_to_avoid_thin_bins(gff_df, binSize):
        '''
        :param gff_df: gene annotation file
        :type gff_df: pandas dataframe
        :param binSize: bin size
        :type binSize: int
        :return: merged gene annotation file
        :rtype: pandas dataframe


        Objective: while iterating over the bins look ahead to the next bin, if the next bin is the last bin on the chromosome merge the last bin with the current bin
        '''
            
        gff_df["bin"] = gff_df["start"].apply(lambda x: math.floor(x/binSize))
        gff_df["bin"] = gff_df["bin"].astype(int)
        gff_df = gff_df.sort_values(["chr", "bin"])
        gff_df.index = range(1, len(gff_df) + 1)

        for chr in gff_df["chr"].unique():
            df_chr = gff_df[gff_df["chr"] == chr]
            df_chr = df_chr.sort_values(["start"])
            df_chr.index = range(1, len(df_chr) + 1)
            df_chr_binned = pd.DataFrame()
            for bin in df_chr["bin"].unique():
                df_chr_bin = df_chr[df_chr["bin"] == bin]
                df_chr_bin = df_chr_bin.sort_values(["start"])
                df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                if bin == df_chr["bin"].unique().max():
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
                else:
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
            gff_df = gff_df.append(df_chr_binned)
        gff_df.index = range(1, len(gff_df) + 1)
        
        return gff_df
    
    def dynamic_bin_adjustment(gff_df, binSize, nFeatures, nSamples):
        '''
        :param gff_df: gene annotation file
        :type gff_df: pandas dataframe
        :param binSize: bin size
        :type binSize: int
        :param nFeatures: number of features
        :type nFeatures: int
        :param nSamples: number of samples
        :type nSamples: int
        :return: adjusted gene annotation file
        :rtype: pandas dataframe


        Objective: while iterating over the bins calculate the average number of features left in the remaining bins, 
        adjust the binSize to capture approximately equal numbers of features in each bin
        '''
            
        gff_df["bin"] = gff_df["start"].apply(lambda x: math.floor(x/binSize))
        gff_df["bin"] = gff_df["bin"].astype(int)
        gff_df = gff_df.sort_values(["chr", "bin"])
        gff_df.index = range(1, len(gff_df) + 1)

        for chr in gff_df["chr"].unique():
            df_chr = gff_df[gff_df["chr"] == chr]
            df_chr = df_chr.sort_values(["start"])
            df_chr.index = range(1, len(df_chr) + 1)
            df_chr_binned = pd.DataFrame()
            for bin in df_chr["bin"].unique():
                df_chr_bin = df_chr[df_chr["bin"] == bin]
                df_chr_bin = df_chr_bin.sort_values(["start"])
                df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                if bin == df_chr["bin"].unique().max():
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
                else:
                    nFeatures_remaining = nFeatures - len(df_chr_bin)
                    nSamples_remaining = nSamples - len(df_chr_bin)
                    nFeatures_remaining_per_bin = nFeatures_remaining/df_chr["bin"].unique().max()
                    nSamples_remaining_per_bin = nSamples_remaining/df_chr["bin"].unique().max()
                    if nFeatures_remaining_per_bin > nSamples_remaining_per_bin:
                        binSize = binSize + 1
                    elif nFeatures_remaining_per_bin < nSamples_remaining_per_bin:
                        binSize = binSize - 1
                    else:
                        pass
                    df_chr_bin["bin"] = df_chr_bin["start"].apply(lambda x: math.floor(x/binSize))
                    df_chr_bin["bin"] = df_chr_bin["bin"].astype(int)
                    df_chr_bin = df_chr_bin.sort_values(["chr", "bin"])
                    df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
            gff_df = gff_df.append(df_chr_binned)
        gff_df.index = range(1, len(gff_df) + 1)

        return gff_df, binSize
