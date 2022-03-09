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

    def filter_for_genomic_variants(gff_df, binSize, vcf_file):
        '''
        :param gff_df: gene annotation file
        :type gff_df: pandas dataframe
        :param binSize: bin size
        :type binSize: int
        :param vcf_file: vcf file
        :type vcf_file: str
        :return: filtered vcf file
        :rtype: pandas dataframe
        :return: filtered gene annotation file
        :rtype: pandas dataframe
        Objective: iterate over the gff_df and the vcf_file to assign a bin to each feature in the gff_df and vcf_file. Remove features in the vcf_file that do not land in the genes in the gff_df.
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

        vcf_df = pd.read_csv(vcf_file, sep="\t", header=None)
        vcf_df.columns = ["chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "format"]
        vcf_df["bin"] = vcf_df["pos"].apply(lambda x: math.floor(x/binSize))
        vcf_df["bin"] = vcf_df["bin"].astype(int)
        vcf_df = vcf_df.sort_values(["chr", "bin"])
        vcf_df.index = range(1, len(vcf_df) + 1)

        for chr in vcf_df["chr"].unique():
            df_chr = vcf_df[vcf_df["chr"] == chr]
            df_chr = df_chr.sort_values(["pos"])
            df_chr.index = range(1, len(df_chr) + 1)
            df_chr_binned = pd.DataFrame()
            for bin in df_chr["bin"].unique():
                df_chr_bin = df_chr[df_chr["bin"] == bin]
                df_chr_bin = df_chr_bin.sort_values(["pos"])
                df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                if bin == df_chr["bin"].unique().max():
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
                else:
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
            vcf_df = vcf_df.append(df_chr_binned)
        vcf_df.index = range(1, len(vcf_df) + 1)

        for index, row in gff_df.iterrows():
            chr = row["chr"]
            start = row["start"]
            bin = row["bin"]
            df_chr = vcf_df[vcf_df["chr"] == chr]
            df_chr = df_chr.sort_values(["bin"])
            df_chr.index = range(1, len(df_chr) + 1)
            df_chr_binned = pd.DataFrame()
            for bin in df_chr["bin"].unique():
                df_chr_bin = df_chr[df_chr["bin"] == bin]
                df_chr_bin = df_chr_bin.sort_values(["pos"])
                df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                if bin == df_chr["bin"].unique().max():
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
                else:
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
            vcf_df = vcf_df.append(df_chr_binned)
        vcf_df.index = range(1, len(vcf_df) + 1)

        vcf_df = vcf_df[vcf_df["chr"] == gff_df["chr"]]
        vcf_df = vcf_df.sort_values(["pos"])
        vcf_df.index = range(1, len(vcf_df) + 1)
        vcf_df = vcf_df.drop_duplicates(subset=["chr", "pos"])
        vcf_df.index = range(1, len(vcf_df) + 1)

        gff_df = gff_df[gff_df["chr"] == vcf_df["chr"]]
        gff_df = gff_df.sort_values(["start"])
        gff_df.index = range(1, len(gff_df) + 1)
        gff_df = gff_df.drop_duplicates(subset=["chr", "start"])

        
        return gff_df, vcf_df

    def weight_variants_in_vcf_file(vcf_file, nFeatures, binSize):
        '''
        :param vcf_file: vcf file
        :type vcf_file: str
        :return: weighted vcf file
        :param nFeatures: number of features
        :type nFeatures: int
        :param binSize: bin size
        :type binSize: int
        
        Objective: Use the filtered vcf file to weight the variants in the vcf file according to the size of the difference between the start and end position of the variant. 
        Re-assign variants into bins with approximately equal nFeatures in each bin and approxmiately equal variant weight. Use the nFeatures to determine the bin size.
        '''

        # read the vcf file
        vcf_df = pd.read_csv(vcf_file, sep="\t", header=None,usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8], skiprows=9)
        vcf_df.columns = ["chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "format"]
        vcf_df["bin"] = vcf_df["pos"].apply(lambda x: math.floor(x/binSize))
        vcf_df["bin"] = vcf_df["bin"].astype(int)
        vcf_df = vcf_df.sort_values(["chr", "bin"])
        vcf_df.index = range(1, len(vcf_df) + 1)

        # assign variants into bins with approximately equal nFeatures in each bin and approxmiately equal variant weight
        for chr in vcf_df["chr"].unique():
            df_chr = vcf_df[vcf_df["chr"] == chr]
            df_chr = df_chr.sort_values(["pos"])
            df_chr.index = range(1, len(df_chr) + 1)
            df_chr_binned = pd.DataFrame()
            for bin in df_chr["bin"].unique():
                df_chr_bin = df_chr[df_chr["bin"] == bin]
                df_chr_bin = df_chr_bin.sort_values(["pos"])
                df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                if bin == df_chr["bin"].unique().max():
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
                else:
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
            vcf_df = vcf_df.append(df_chr_binned)
        vcf_df.index = range(1, len(vcf_df) + 1)

        # assign a weight to each variant based on the size of the difference between the start and end position of the variant
        vcf_df["weight"] = vcf_df["pos"].apply(lambda x: math.floor(x/binSize))
        vcf_df["weight"] = vcf_df["weight"].astype(int)
        vcf_df = vcf_df.sort_values(["chr", "weight"])
        vcf_df.index = range(1, len(vcf_df) + 1)

        # re-assign variants into bins with approximately equal nFeatures in each bin and approxmiately equal variant weight
        for chr in vcf_df["chr"].unique():
            df_chr = vcf_df[vcf_df["chr"] == chr]
            df_chr = df_chr.sort_values(["weight"])
            df_chr.index = range(1, len(df_chr) + 1)
            df_chr_binned = pd.DataFrame()
            for bin in df_chr["weight"].unique():
                df_chr_bin = df_chr[df_chr["weight"] == bin]
                df_chr_bin = df_chr_bin.sort_values(["pos"])
                df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                if bin == df_chr["weight"].unique().max():
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
                else:
                    df_chr_binned = df_chr_binned.append(df_chr_bin)
            vcf_df = vcf_df.append(df_chr_binned)
        vcf_df.index = range(1, len(vcf_df) + 1)

        return vcf_df

    def match_variants_in_vcf_to_bins_in_gff(vcf_df, gff_df):
        '''
        :param vcf_df: weighted vcf file
        :type vcf_df: pd.DataFrame
        :param gff_df: gff file
        :type gff_df: pd.DataFrame
        :return: matched features
        :rtype: pd.DataFrame
        Objective: Use the weighted vcf file to match the variants in the vcf file to the features in the gff_df. Return a merged pandas dataframe with 
        the matched features and the corresponding variants.
        '''

        # match variants in vcf file to features in gff_df
        vcf_df = vcf_df[vcf_df["chr"] == gff_df["chr"]]
        vcf_df = vcf_df.sort_values(["pos"])
        vcf_df.index = range(1, len(vcf_df) + 1)
        vcf_df = vcf_df.drop_duplicates(subset=["chr", "pos"])
        vcf_df.index = range(1, len(vcf_df) + 1)

        gff_df = gff_df[gff_df["chr"] == vcf_df["chr"]]
        gff_df = gff_df.sort_values(["start"])
        gff_df.index = range(1, len(gff_df) + 1)
        gff_df = gff_df.drop_duplicates(subset=["chr", "start"])

        gff_df["bin"] = gff_df["start"].apply(lambda x: math.floor(x/binSize))
        gff_df["bin"] = gff_df["bin"].astype(int)
        gff_df = gff_df.sort_values(["chr", "bin"])
        gff_df.index = range(1, len(gff_df) + 1)

        vcf_df["bin"] = vcf_df["pos"].apply(lambda x: math.floor(x/binSize))
        vcf_df["bin"] = vcf_df["bin"].astype(int)
        vcf_df = vcf_df.sort_values(["chr", "bin"])
        vcf_df.index = range(1, len(vcf_df) + 1)

        # merge the two dataframes on the bin column
        merged_df = vcf_df.merge(gff_df, how="inner", on=["chr", "bin"])
        merged_df.index = range(1, len(vcf_df) + 1)

        # which column number of the merged df "bin" column
        bin_col = merged_df.columns.get_loc("bin")




        return merged_df, bin_col

    def transpose_vcf_into_tensor(vcf_file, genotype_classes):
        '''
        :param vcf_file: vcf file
        :type vcf_file: pd.DataFrame
        :param genotype_classes: list of genotype classes
        :type genotype_classes: list
        :return: tensor
        :rtype: torch.Tensor
        :return: transposed vcf file
        :rtype: pd.DataFrame
        Objective: Transpose the vcf file into a tensor. With the columns of the vcf file being the variants and rows being the genontypes
        and the values being the calls of the variant in the genontype.
        '''

        # read in the first 8 columns of the vcf file, skip the header
        vcf_df = pd.read_csv(vcf_file, sep="\t", header=None, usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8], skiprows=9)
        vcf_df.columns = ["chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "format"]
        vcf_df.index = range(1, len(vcf_df) + 1)

        # index the rows of the vcf files
        vcf_df["index"] = vcf_df.index

        # read in the columns of the vcf files after the first 8 columns
        genotypes_mat = pd.read_csv(vcf_file, sep="\t", header=None, usecols=range(8, len(vcf_df.columns)), skiprows=9)

        # index the rows of the genotypes matrix
        genotypes_mat["index"] = genotypes_mat.index

        # transpose the genotypes matrix into a dataframe with the columns being the genotypes and the rows being the variants
        genotypes_df = genotypes_mat.transpose()
        genotypes_df.columns = genotype_classes
        genotypes_df.index = range(1, len(genotypes_df) + 1)

        # transpose the vcf dataframe into a dataframe with the columns being the variants and the rows being the index
        vcf_df = vcf_df.transpose()
        vcf_df.columns = genotype_classes
        vcf_df.index = range(1, len(vcf_df) + 1)


        # read in the genotype classes file
        genotype_classes_df = pd.read_csv(genotype_classes, sep="\t", header=None)
        genotype_classes_df.columns = ["genotype"]
        genotype_classes_df.index = range(1, len(genotype_classes_df) + 1)

        # create a tensor in which the x axis is the genotype class, the y axis is the variant and the z axis is the genotype call
        tensor = torch.zeros(len(genotype_classes_df), len(vcf_df), len(genotypes_df))

        # fill the tensor with the genotype calls for each variant
        for i in range(1, len(genotypes_df) + 1):
            for j in range(1, len(vcf_df) + 1):
                for k in range(1, len(genotype_classes_df) + 1):
                    tensor[k - 1, j - 1, i - 1] = genotypes_df.iloc[i - 1, j - 1]
        
        return tensor, vcf_df

    def extract_dataframe_from_tensor_based_on_bin(vcf_df, tensor, bin_col):
        '''
        :param vcf_df: vcf file
        :type vcf_df: pd.DataFrame
        :param tensor: tensor
        :type tensor: torch.Tensor
        :param bin_col: column number of the bin column
        :type bin_col: int
        :return: dataframe
        :rtype: pd.DataFrame
        Objective: Extract the dataframe from the tensor based on the bin column.
        '''

        # extract the dataframe from the tensor based on the bin column
        dataframe = tensor[:, vcf_df.iloc[:, bin_col].values - 1, :]
        dataframe = dataframe.transpose(0, 2, 1)

        return dataframe

    def extract_dataframe_from_tensor_based_on_variant(vcf_df, tensor, variant_col):
        '''
        :param vcf_df: vcf file
        :type vcf_df: pd.DataFrame
        :param tensor: tensor
        :type tensor: torch.Tensor
        :param variant_col: column number of the variant column
        :type variant_col: int
        :return: dataframe
        :rtype: pd.DataFrame
        Objective: Extract the dataframe from the tensor based on the variant column.
        '''

        # extract the dataframe from the tensor based on the variant column
        dataframe = tensor[:, :, vcf_df.iloc[:, variant_col].values - 1]

        return dataframe

    def extract_dataframe_from_tensor_based_on_genotype(vcf_df, tensor, genotype_col):
        '''
        :param vcf_df: vcf file
        :type vcf_df: pd.DataFrame
        :param tensor: tensor
        :type tensor: torch.Tensor
        :param genotype_col: column number of the genotype column
        :type genotype_col: int
        :return: dataframe
        :rtype: pd.DataFrame
        Objective: Extract the dataframe from the tensor based on the genotype column.
        '''

        # extract the dataframe from the tensor based on the genotype column
        dataframe = tensor[:, :, vcf_df.iloc[:, genotype_col].values - 1]

        return dataframe