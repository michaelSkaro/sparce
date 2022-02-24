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


class parse_gene_annotation:
    def __init__(self, geneAnnotationFile):
        self.geneAnnotationFile = geneAnnotationFile
        pass

    def parse_gfffile(geneAnnotationFile):
        '''
        :param geneAnnotationFile: gene annotation file
        :type geneAnnotationFile: string
        :return: gene annotation file
        :rtype: pandas dataframe


        format: chr start end name score strand
        Objective: parse organism gff file into pandas dataframe
        '''
        df = pd.read_csv(geneAnnotationFile, sep="\t", header=None, comment =  "#")
        df.columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        df["chr"] = df["chr"].astype(str)
        df["source"] = df["source"].astype(str)
        df["feature"] = df["feature"].astype(str)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        df["score"] = df["score"].astype(int)
        df["strand"] = df["strand"].astype(str)
        df["frame"] = df["frame"].astype(str)
        df["attribute"] = df["attribute"].astype(str)
        

        df = df.sort_values(["chr", "start", "end"])
        df.index = range(1, len(df) + 1)

        return df

    def bin_genome(df, binSize):
        '''
        :param df: gene annotation file
        :type df: pandas dataframe
        :param binSize: bin size
        :type binSize: int
        :return: gene annotation file
        :rtype: pandas dataframe
        

        format: chr start end name score strand

        Objective: bin features of the genome into bins accoring to sliding window equal to the binSize

        '''
        df["bin"] = df["start"].apply(lambda x: math.floor(x/binSize))
        df["bin"] = df["bin"].astype(int)
        df = df.sort_values(["chr", "bin"])
        df.index = range(1, len(df) + 1)

        

        return df
