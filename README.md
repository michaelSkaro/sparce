<h1 align="center">Statistical Preprocessing of Attributes via Cross-validated Recursive Elimination</h1>

*sparce*

The sparce software is a machine learning based software for automated feature seleciton in genomics data files. The software was originially outfitted to select features in transcriptomic profiling data but has now been outfitted for general use in genetics, transcriptomics, methylomics and ATAC-seq data.

### Installation

We utilized the python programming language to build ths package. We have uploaded this package onto the public python development forum pypi. To install this package:


```
conda create -n sparce pip

conda activate sparce
```

```
pip install sparce
```

### Usage

```python
'''
Run inside script
'''


import sparce
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

CV = sparce.sparce.feature_selection.grade_features(X = X, y = y, nFeatures = nFeatures , n_jobs = nJobs)


```

### CLI
1.  Clone the repository and re-invoke the main function. 
2.  import args_parse into the sparce.py 
3.  Ready to run in the cli

```
python sparce.py -x <file> -y <target> -nFeatures <int> -nJobs <int>
```

```
conda deactivate sparce
```

### sparce assumptions
- The data is in tidy format where (Features x samples) with a column labeled "target"
- The features are continuous attributes in a classificaiton problem
- The classes are mutually exclusive
- nFeatures > nSamples,  you are attempting to reduce the dimensionality of the problem to produce nSamples > nFeatures


### sparce documentaitons:

RTD build *coming soon*


### Descritption of variables
  - X = The descriptive variables measuring features of the target variable y. X should be in long format as standard. Where variables are columns and instances are rows. 
  - y = The target variable. This variable is expected to be a numpy array of length = len(X['A_column']). The y variable should be encoded using the sklearn OHE encoding or ordinally encoded. 
  - nFeatures = Top n features toreturn for the selection algorithms. The objetive of this parameter is to limit the recursive feature selection from running indefinitely. 
  - nJobs = compute resources to dedicate to the job at hand. RFE is slow, be prepared to get a cup of coffee. This can be shut off if the user wishes. 

### Help python3 sparce.py --help 
