#### Statistical Processing of attributes via Recursive Cross Elimination




*SPARCE*

The sparce software is a statistical machine learning software that automates
feature seleciton in genomics data files. The software was originally outiftted
for general use in genetics, transcirptomics, methylomics and ATAC-seq data.

Installation

```{python}
conda create -n sparce pip
conda activate sparce
```

```{python}
pip install sparce
```


***HOW TO RUN***

```{python}
'''
Run inside script
'''


import sparce
from sparce import feature_selection
import pandas as pd
from sklearn.preprocessing import OrdinalEncoder

def preprocess(file): 
  X = pd.read_csv('file')
  enc = OrdinalEncoder()
  enc.fit(X['a column in X'])
  X['a column in X'] = enc.transform(X['a column in X'])
  y = X['a column in X']
  X = X.drop('a column in X', axis = 1)
  
  return X,y

X, y = preprocess(file)

nFeatures = 5
nJobs = 10

CV = feature_selection.grade_features(X = X, y = y, nFeatures = nFeatures , nJobs = nJobs)

```


# CLI

Clone the repository and re-invoke the main function.
import args_parse into the sparce.py
Ready to run in the cli

```console

python sparce.py -x <file> -y <target> -nFeatures <int> -nJobs <int>

conda deactivate sparce

```

sparce assumptions

The data is in tidy format where (Features x samples) with a column labeled "target"
The features are continuous attributes in a classificaiton problem
The classes are mutually exclusive
nFeatures > nSamples, you are attempting to reduce the dimensionality of the problem to produce nSamples > nFeatures





