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

#from sparce import fs.args_parse

class feature_selection:
   def __init__(self, X, y, nFeatures,nJobs):
       self.X = X
       self.y = y
       self.nFeatures = nFeatures
       self.nJobs = nJobs
       pass
   
   
   def chiSQ_selector(X, y, nFeatures):
       '''
       input: X, y, nFeatures
       output: chiSQ_support, feature_list

       Objective: Use the chi-squared test to select features.
       Reference: https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.chi2.html

       '''
       
       feature_list = X.columns
       chiSQ = SelectKBest(chi2, k= nFeatures)
       chiSQ.fit(X, y.ravel())
       chiSQ_support = chiSQ.get_support()

       return chiSQ_support, feature_list
   
   
   def rfR(X, y, nFeatures,nJobs):
       '''
       input: X, y, nFeatures, nJobs
       output: rfR_support

       Objective: Use the RandomForestRegressor to select features.
       Reference: https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html


       '''


       rfR_features = SelectFromModel(
           RandomForestRegressor(n_estimators=100, n_jobs = nJobs), max_features=nFeatures
       )
       rfR_features.fit(X, y.ravel())
       rfR_support = rfR_features.get_support()

       return rfR_support
   
   
   def recursiveFeatureSelection(X, y, nFeatures,nJobs):
       '''
       input: X, y, nFeatures, nJobs
       output: rfe_support

       Objective: Use the Recursive Feature Elimination to select features.
       Reference: https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.RFE.html

       '''



       rfe_selector = RFE(
           estimator=LogisticRegression(n_jobs = nJobs),
           n_features_to_select=nFeatures,
           step= math.celi(len(X.columns)/10),
           verbose=0,
       )
       rfe_selector.fit(X, y.ravel())
       rfe_support = rfe_selector.get_support()

       return rfe_support
   
   
   
   def lassoR(X, y, nFeatures,nJobs):
       '''
       input: X, y, nFeatures, nJobs
       output: lassoR_support

       Objective: Use the Lasso Regression to select features.
       Reference: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Lasso.html
       
       '''



       lR_selector = SelectFromModel(
           LogisticRegression(penalty="l2", n_jobs= nJobs), max_features=nFeatures
       )
       lR_selector.fit(X, y.ravel())
       lR_support = lR_selector.get_support()

       return lR_support
   
   
   
   def rfC(X, y, nFeatures,nJobs):
       '''
       input: X, y, nFeatures, nJobs
       output: rfC_support

       Objective: Use the RandomForestClassifier to select features.
       Reference: https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html
       
       '''



       rfC_features = SelectFromModel(
           RandomForestClassifier(n_estimators=100, n_jobs=nJobs), max_features=nFeatures
       )
       rfC_features.fit(X, y.ravel())
       rfC_support = rfC_features.get_support()

       return rfC_support
   
   # cross validation
   
   def cross_validate_feature_selection(feature_list, chiSQ_support, rfe_support, lR_support, rfC_support,rfR_support):
        '''
        input: feature_list, chiSQ_support, rfe_support, lR_support, rfC_support,rfR_support
        output: df

        Objective: Use the cross validation to select features.
        Reference: https://towardsdatascience.com/feature-selection-techniques-in-machine-learning-with-python-f24e7da3f36e
        
        '''
        df = pd.DataFrame(
            {
                "Feature": feature_list,
                "chi2": chiSQ_support,
                "RFE": rfe_support,
                "Logistics": lR_support,
                "RandomForestClassifier": rfC_support,
                "RandomForstRegression": rfR_support,
            }
        )
        df["Total"] = np.sum(df, axis=1)
        df = df.sort_values(["Total", "Feature"], ascending=False)
        df.index = range(1, len(df) + 1)
        
        return df

   def grade_features(X, y, nFeatures, n_jobs):
       '''
       input: X, y, nFeatures, n_jobs
       output: CV

       Objective: Call each function in the feature selection class, grade each set of features with each of the algorithms.
       Sum the total boolen supports for each feature and then sort the features by the sum.
       
       Reference: https://towardsdatascience.com/feature-selection-techniques-in-machine-learning-with-python-f24e7da3f36e

       '''
       chiSQ_support, feature_list = feature_selection.chiSQ_selector(X, y, nFeatures=nFeatures)
       rfe_support = feature_selection.recursiveFeatureSelection(X, y, nFeatures=nFeatures, nJobs = nJobs)
       lR_support = feature_selection.lassoR(X, y, nFeatures=nFeatures, nJobs = nJobs)
       rfC_support = feature_selection.rfC(X, y, nFeatures=nFeatures, nJobs = nJobs)
       rfR_support = feature_selection.rfR(X, y, nFeatures=nFeatures, nJobs = nJobs)

       CV = feature_selection.cross_validate_feature_selection(
           feature_list,
           chiSQ-support,
           rfe_support,
           lR_support,
           rfC_support,
           rfR_support,
       )
       CV = CV[1:nFeatures]

       return CV

class parseVariantAnnotationFile:
    def __inti__(self, variatnAnnotationFile):
        self.variatnAnnotationFile = variatnAnnotationFile
        pass

    def parse_bedfile(variantAnnotationFile):

        '''
        input: bedfile
        output: dataframe

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
        input: annotated bedfile
        output: dataframe

        format: chr start end name score strand

        Objective: count the number of variants in each bin.

        '''
        # count the number of variants in each bin
        df["bin_count"] = df.groupby(["bin"])["bin"].transform("count")

        return df
        
        
    def count_variants_in_bins_per_chromosome(df):
        '''
        input: annotated bedfile
        output: dataframe

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
        input: annotated bedfile
        output: dataframe

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
        input: annotated bedfile
        output: dataframe

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
        input: annotated bedfile
        output: dataframe

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
        input: annotated bedfile
        output: dataframe

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
        input: annotated bedfile
        output: dataframe

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
        input: annotated bed file
        ouput: visualized plot

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
        input: vcf file
        output: dataframe

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

#--------------------------------#

class parse_gene_annotation:
    def __init__(self, geneAnnotationFile):
        self.geneAnnotationFile = geneAnnotationFile
        pass

    def parse_gfffile(geneAnnotationFile):

        '''
        input: gfffile
        output: dataframe

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
        input: parsed gff file
        output: dataframe

        format: chr start end name score strand

        Objective: bin features of the genome into bins accoring to sliding window equal to the binSize

        '''
        df["bin"] = df["start"].apply(lambda x: math.floor(x/binSize))
        df["bin"] = df["bin"].astype(int)
        df = df.sort_values(["chr", "bin"])
        df.index = range(1, len(df) + 1)

        

        return df


#--------------------------------#
class prepare_transcriptomics_file_for_SPARCE:
    def __inti__(self, transcriptomics_file):
            self.transcriptomics_file = transcriptomics_file
            pass
    
    def _parse_transcriptomics_file(transcriptomics_file):
        '''
        input: transcriptomics_file
        output: parsed_transcriptomics_file

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
        input: parsed transcriptomics file
        output: X, y

        Objective: split the transcriptomics file into two dataframes, one for the transcript expression and one for the target classes

        
        '''
        X = df.drop("target", axis=1)
        y = df["target"]
        
        return X, y
    

    
    def _balance_data_in_classes(X, y):
        '''
        input: X, y
        output: X_balanced, y_balanced

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
        input: seq, num
        output: list of tuples for bin start and end points

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

#--------------------------------#
'''
def main(cla_args=sys.argv[1:]):
   #--------------------------------#
   args = fs.args_parse(parse_cla(cla_args))
   X = args['x']
   y = args['y']
   nFeatures = args['f']
   nJobs = args['j']
   

if __name__ == '__main__':
    main()
'''
#--------------------------------#


