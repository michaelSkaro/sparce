# sparce
A software resource to split a genomics, methylomics or trasncriptomics file into bins and select features within the bins using recursive feature elimination and machine learning


## Feature selection for ML workflows
This python package will automate feature selection in ML projects in the Arnold lab. The method was used as a feature selection method in my first chapter but seems to service the lab. I will keep this package updated for the projects in the coming months and as a vehicle for developing my introduction to my thesis, a review of feature selection methods in ML and which ones I used for the feature selection methods in CH1. 

In the coming weeks I will add a plethora of information describing the use of the package and the algorithms that are working under the hood.

For now we will settle for a simple Read me to get us started on 0.0.1. 

As a general description this package is intended for the use of novice users looking for general feature selection in an automated fashion. This will not replace your own featue analysis but can help users make a good first step eliminating redundant features and begin looking for strong signals in their feature columns. 

### Installation

We utilized the python programming language to build ths package. We have uploaded this package onto the public python development forum pypi. To install this package:

```
pip install -i https://test.pypi.org/simple/ sparce
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

### Descritption of variables
  - X = The descriptive variables measuring features of the target variable y. X should be in long format as standard. Where variables are columns and instances are rows. 
  - y = The target variable. This variable is expected to be a numpy array of length = len(X['A_column']). The y variable should be encoded using the sklearn OHE encoding or ordinally encoded. 
  - nFeatures = Top n features toreturn for the selection algorithms. The objetive of this parameter is to limit the recursive feature selection from running indefinitely. 
  - nJobs = compute resources to dedicate to the job at hand. RFE is slow, be prepared to get a cup of coffee. This can be shut off if the user wishes. 

### Help python3 sparce.py --help 
