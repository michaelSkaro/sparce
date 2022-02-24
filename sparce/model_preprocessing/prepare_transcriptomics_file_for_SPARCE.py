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

class prepare_transcriptomics_file_for_SPARCE:
    def __inti__(self, transcriptomics_file):
            self.transcriptomics_file = transcriptomics_file
            pass
    
    def _parse_transcriptomics_file(transcriptomics_file):
        '''
        :param transcriptomics_file: transcriptomics file
        :type transcriptomics_file: string
        :return: transcriptomics file
        :rtype: pandas dataframe
        

        Objective: parse the transcriptomics file into pandas dataframe
        '''
        
        df = pd.read_table(
            transcriptomics_file,
            delimiter=",",
        ).fillna(0)

        if "barcode" in df.columns:
            df = df.drop("barcode", axis=1)
        if "__no_feature" in df.columns:
            df = df.drop("__no_feature", axis=1)
        if "Unnamed: 0" in df.columns:
            df = df.drop("Unnamed: 0", axis=1)
        if "__ambiguous" in df.columns:
            df = df.drop("__ambiguous", axis=1)
        if "__too_Low_aQual" in df.columns:
            df = df.drop("__too_low_aQual", axis=1)
        if "__not_aligned" in df.columns:
            df = df.drop("__not_aligned", axis=1)
        if "__alignment_not_unique" in df.columns:
            df = df.drop("__alignment_not_unique", axis=1)
        
        return df
    
    def _split(df):
        '''
        :param df: transcriptomics file
        :type df: pandas dataframe
        :return: transcriptomics file
        :rtype: pandas dataframe


        Objective: split the transcriptomics file into two dataframes, one for the transcript expression and one for the target classes

        
        '''
        X = df.drop("target", axis=1)
        y = df["target"]
        
        return X, y
    

    
    def _balance_data_in_classes(X, y):
        '''
        :param X: transcriptomics file
        :type X: pandas dataframe
        :param y: target classes
        :type y: pandas dataframe
        :return: balanced transcriptomics file
        :rtype: pandas dataframe

        Objective: balance the data in the classes using the sklearn oversample and SMOTE methods
        
        '''
        from imblearn.over_sampling import SMOTE
        from imblearn.pipeline import Pipeline
        from imblearn.under_sampling import RandomUnderSampler, TomekLinks

        counter = Counter(y_train)
        oversample = SMOTE()
        Xsm, ysm = oversample.fit_resample(X, y)
        Xsm = Xsm.astype(int)

        return Xsm, ysm

    def _chunkIt(seq, num):
        """
        :param seq: total length of sequence
        :type seq: int
        :param num: number of chunks
        :type num: int
        :return: tuples of sequence bins
        :rtype: list of tuples

        Objective: chunk the data into bins of size num

        """
        avg = len(seq) / float(num)
        out = []
        last = 0.0

        while last < len(seq):
            out.append(seq[int(last) : int(last + avg)])
            last += avg

        return out
    
    def _bin_transcript_features_based_on_gene_annotation(df, gff_df, binSize):
        '''
        :param df: transcriptomics file
        :type df: pandas dataframe
        :param gff_df: gene annotation file
        :type gff_df: pandas dataframe
        :param binSize: bin size
        :type binSize: int
        :return: binned transcriptomics file
        :rtype: pandas dataframe



        input: df, gff_df, binSize
        output: df_binned

        Objective: bin transcript features based on gene annotation

        '''
        df["bin"] = df["start"].apply(lambda x: math.floor(x/binSize))
        df["bin"] = df["bin"].astype(int)
        df = df.sort_values(["chr", "bin"])
        df.index = range(1, len(df) + 1)
        df_binned = pd.DataFrame()
        for chr in gff_df["chr"].unique():
            df_chr = df[df["chr"] == chr]
            df_chr = df_chr.sort_values(["start"])
            df_chr.index = range(1, len(df_chr) + 1)
            df_chr_binned = pd.DataFrame()
            for bin in df_chr["bin"].unique():
                df_chr_bin = df_chr[df_chr["bin"] == bin]
                df_chr_bin = df_chr_bin.sort_values(["start"])
                df_chr_bin.index = range(1, len(df_chr_bin) + 1)
                df_chr_binned = df_chr_binned.append(df_chr_bin)
            df_binned = df_binned.append(df_chr_binned)
        df_binned.index = range(1, len(df_binned) + 1)
        return df_binned
