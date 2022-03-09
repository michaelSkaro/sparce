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


class parseVariantAnnotationFile:
    def __inti__(self, variatnAnnotationFile):
        self.variatnAnnotationFile = variatnAnnotationFile
        pass

    def parse_bedfile(variantAnnotationFile):

        '''
        :param variantAnnotationFile: variant annotation file
        :type variantAnnotationFile: string
        :return: variant annotation file
        :rtype: pandas dataframe
        format: chr start end name score strand
        Objective: parse a bed file that contains variants from an organism of interest into apandas dataframe. Bin variants into a sliding window of size 100kb.
        '''
        # read in the bed file
        df = pd.read_csv(variantAnnotationFile, sep="\t", header=None, comment="#")
        df.columns = ["chr", "start", "end", "name", "score", "strand"]
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        df["score"] = df["score"].astype(int)
        df["strand"] = df["strand"].astype(str)
        df = df.sort_values(["chr", "start", "end"])
        df.index = range(1, len(df) + 1)

        # bin the variants into a sliding window of size 100kb
        df["bin"] = df["start"].apply(lambda x: math.ceil(x / 100000))


        return df

    def count_variants_in_bins(df):
        '''
        :param df: variant annotation file
        :type df: pandas dataframe
        :return: variant annotation file
        :rtype: pandas dataframe
        format: chr start end name score strand
        Objective: count the number of variants in each bin.
        '''
        # count the number of variants in each bin
        df["bin_count"] = df.groupby(["bin"])["bin"].transform("count")

        return df
        
        
    def count_variants_in_bins_per_chromosome(df):
        '''
        :param df: variant annotation file
        :type df: pandas dataframe
        :return: variant annotation file
        :rtype: pandas dataframe
        format: chr start end name score strand
        Objective: count the number of variants in each bin.
        '''
        # count the number of variants in each bin
        
        df["bin_count"] = df.groupby(["bin"])["bin"].transform("count")

        # count the number of variants in each chromosome
        df["chr_count"] = df.groupby(["chr"])["chr"].transform("count")



        return df
    

    def count_variants_in_bins_per_chromosome_per_strand(df):
        '''
        :param df: variant annotation file
        :type df: pandas dataframe
        :return: variant annotation file
        :rtype: pandas dataframe
        format: chr start end name score strand
        Objective: count the number of variants in each bin.
        '''
        # count the number of variants in each bin
        df["bin_count"] = df.groupby(["bin"])["bin"].transform("count")

        # count the number of variants in each chromosome
        df["chr_count"] = df.groupby(["chr"])["chr"].transform("count")

        df["strand_count"] = df.groupby(["strand"])["strand"].transform("count")

        return df
    
    def count_variants_in_bins_per_chromosome_per_strand_per_score(df):
        '''
        :param df: variant annotation file
        :type df: pandas dataframe
        :return: variant annotation file
        :rtype: pandas dataframe
        
        format: chr start end name score strand
        Objective: count the number of variants in each bin.
        '''
        # count the number of variants in each bin
        df["bin_count"] = df.groupby(["bin"])["bin"].transform("count")

        # count the number of variants in each chromosome
        df["chr_count"] = df.groupby(["chr"])["chr"].transform("count")

        df["strand_count"] = df.groupby(["strand"])["strand"].transform("count")

        df["score_count"] = df.groupby(["score"])["score"].transform("count")

        return df
    
    def count_variants_in_bins_per_chromosome_per_strand_per_score_per_name(df):
        '''
        :param df: variant annotation file
        :type df: pandas dataframe
        :return: variant annotation file
        :rtype: pandas dataframe
        format: chr start end name score strand
        Objective: count the number of variants in each bin.
        '''
        # count the number of variants in each bin
        df["bin_count"] = df.groupby(["bin"])["bin"].transform("count")

        # count the number of variants in each chromosome
        df["chr_count"] = df.groupby(["chr"])["chr"].transform("count")

        # count the number of variants in each chromosome
        df["strand_count"] = df.groupby(["strand"])["strand"].transform("count")

        df["score_count"] = df.groupby(["score"])["score"].transform("count")

        df["name_count"] = df.groupby(["name"])["name"].transform("count")

        return df
    
    def count_variants_in_bins_per_chromosome_per_strand_per_score_per_name_per_start(df):
        '''
        :param df: variant annotation file
        :type df: pandas dataframe
        :return: variant annotation file
        :rtype: pandas dataframe
        
        format: chr start end name score strand
        Objective: count the number of variants in each bin.
        '''
        # count the number of variants in each bin
        df["bin_count"] = df.groupby(["bin"])["bin"].transform("count")

        # count the number of variants in each chromosome
        df["chr_count"] = df.groupby(["chr"])["chr"].transform("count")

        df["strand_count"] = df.groupby(["strand"])["strand"].transform("count")

        df["score_count"] = df.groupby(["score"])["score"].transform("count")

        df["name_count"] = df.groupby(["name"])["name"].transform("count")
        df["start_count"] = df.groupby(["start"])["start"].transform("count")

        return df
    
    def count_variants_in_bins_per_chromosome_per_strand_per_score_per_name_per_start_per_end(df):
        '''
        :param df: variant annotation file
        :type df: pandas dataframe
        :return: variant annotation file
        :rtype: pandas dataframe
        format: chr start end name score strand
        Objective: count the number of variants in each bin.
        '''
        # count the number of variants in each bin
        df["bin_count"] = df.groupby(["bin"])["bin"].transform("count")
        # count the number of variants in each chromosome
        df["chr_count"] = df.groupby(["chr"])["chr"].transform("count")

        df["strand_count"] = df.groupby(["strand"])["strand"].transform("count")

        df["score_count"] = df.groupby(["score"])["score"].transform("count")

        df["name_count"] = df.groupby(["name"])["name"].transform("count")

        df["start_count"] = df.groupby(["start"])["start"].transform("count")

        df["end_count"] = df.groupby(["end"])["end"].transform("count")

        return df
    
    def plot_variant_density_in_bins_per_chromosome(df):
        '''
        :param df: variant annotation file
        :type df: pandas dataframe
        :return: variant annotation file
        :rtype: pandas dataframe
        :plot: variant density in bins per chromosome
        Objective: Sum the variants in the each sliding window and vsiaulize the density of the variants in each of the bins.
        '''
        # count the number of variants in each bin
        df["bin_count"] = df.groupby(["bin"])["bin"].transform("count")

        # count the number of variants in each chromosome
        df["chr_count"] = df.groupby(["chr"])["chr"].transform("count")

        # plot the density of variants in each bin
        plt.figure(figsize=(20,10))
        sns.distplot(df["bin_count"], bins=100, kde=False)
        plt.xlabel("Number of variants in each bin")
        plt.ylabel("Number of bins")
        plt.title("Density of variants in each bin")
        plt.savefig("density_of_variants_in_each_bin.png")
        plt.show()

        return df


    def associate_variants_with_genes(genomeAnnotationFile, annotatedBedFile):
        '''
        :param genomeAnnotationFile: genome annotation file
        :type genomeAnnotationFile: string
        :param annotatedBedFile: annotated bed file
        :type annotatedBedFile: string
        :return: variant annotation file
        :rtype: pandas dataframe
    
        input: genomeAnnotationFile and annotatedBedFile
        Output: dataframe
        Objective: for each variant in the annotatedBedFile find any gene within 100kb of the variant using the genomeAnnotaitonFile
        '''
        # read the genome annotation file
        genomeAnnotation = pd.read_csv(genomeAnnotationFile, sep="\t",  header=None, comment="#")

        # read the annotated bed file
        annotatedBed = pd.read_csv(annotatedBedFile, sep="\t",header=None, comment="#")

        # extract the chromosome, start and end of the variant
        annotatedBed["chr"] = annotatedBed["chr"].astype(str)
        annotatedBed["start"] = annotatedBed["start"].astype(int)
        annotatedBed["end"] = annotatedBed["end"].astype(int)

        # extract the chromosome, start and end of the gene
        genomeAnnotation["chr"] = genomeAnnotation["chr"].astype(str)
        genomeAnnotation["start"] = genomeAnnotation["start"].astype(int)
        genomeAnnotation["end"] = genomeAnnotation["end"].astype(int)

        # iterate over the variants in the annotated bed file and find genes within 100kb of the variant
        for index, row in annotatedBed.iterrows():
            # extract the chromosome, start and end of the variant
            chr = row["chr"]
            start = row["start"]
            end = row["end"]

            # extract the chromosome, start and end of the gene
            genes = genomeAnnotation[(genomeAnnotation["chr"] == chr) & (genomeAnnotation["start"] <= start + 100000) & (genomeAnnotation["end"] >= start - 100000)]

            # if there are genes within 100kb of the variant, add the gene name to the annotated bed file
            if genes.shape[0] > 0:
                annotatedBed.loc[index, "gene"] = genes["name"].values[0]
            else:
                annotatedBed.loc[index, "gene"] = "NA"
        
        return annotatedBed
    
    
    def parse_vcf(variantAnnotationFile):
        '''
        :param variantAnnotationFile: variant annotation file
        :type variantAnnotationFile: string
        :return: variant annotation file
        :rtype: pandas dataframe
        format: chr start end name score strand
        '''
        df = pd.read_csv(variantAnnotationFile, sep="\t", header=None, comment="#")
        df.columns = ["chr", "start", "end", "name", "score", "strand"]
        df["chr"] = df["chr"].astype(str)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        df["score"] = df["score"].astype(int)
        df["strand"] = df["strand"].astype(str)

        df = df.sort_values(["chr", "start", "end"])
        df.index = range(1, len(df) + 1)

        return df