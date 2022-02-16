<h1 align="center">Statistical Preprocessing of Attributes via Cross-validated Recursive Elimination</h1>

*sparce*

The sparce software is a machine learning based software for automated feature seleciton in genomics data files. The software was originially outfitted to select features in transcriptomic profiling data but has now been outfitted for general use in genetics, transcriptomics, methylomics and ATAC-seq data.

### Installation

We utilized the python programming language to build ths package. We have uploaded this package onto the public python development forum pypi. To install this package:

```
pip install -i sparce
```

### Usage

```python
import sparce as fs
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
nFeatures = 10

CV = fs.grade_features(X = X, y = y, nFeatures = nFeatures , nJobs = nJobs)


```

### sparce assumptions
- The data is in tidy format where (Features x samples) with a column labeled "target"
- The features are continuous attributes in a classificaiton problem
- The classes are mutuallyexclusive
- There are more samples than features


### sparce documentaitons:

RTD build *coming soon*


### Descritption of variables
  - X = The descriptive variables measuring features of the target variable y. X should be in long format as standard. Where variables are columns and instances are rows. 
  - y = The target variable. This variable is expected to be a numpy array of length = len(X['A_column']). The y variable should be encoded using the sklearn OHE encoding or ordinally encoded. 
  - nFeatures = Top n features toreturn for the selection algorithms. The objetive of this parameter is to limit the recursive feature selection from running indefinitely. 
  - nJobs = compute resources to dedicate to the job at hand. RFE is slow, be prepared to get a cup of coffee. This can be shut off if the user wishes. 

### Help python3 sparce.py --help 
