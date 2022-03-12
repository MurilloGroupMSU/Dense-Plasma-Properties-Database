
____
## DeSilva Data
____

This database contains electrical conductivity data from exploding wire experiments conducted by DeSilva. The original publication is given as a PDF in this repo. Specifially, the database contains elements:
* Al,
* Fe,
* Ni,
* Cu,
* W.

The data is in the form of {Z_{nuc}, $\rho$, T, $\sigma$} with units {none, g/cc, K, S/m}. The numerical values in the CSV files are transformed such that the base-10 logarithm of the temperature and conducvity are given. 


### Modified Lee-More Model (MLM)

In the original paper [Murillo, Frontiers in Physics, 2022] that led to the creation of this database, a modification of the Lee-More conductivity model was used. The MLM has been written as an easy to use Python library and is included here with the data. To use this library:
1. install the library in the apprpriate path on your computer
2. `import mlm`
3. `conds = mlm(Z, $\rho$, T)`


### Files Included 

I have organized the files in different ways and have provided a Jupyter notebook that allows one to explore the data and perhaps allow the user to create customized files. These files are:
* `MLM.py` -- modified Lee-More Model
* `DeSilvaKatsouros.pdf` -- original DeSilva paper
* `View_DeSilva_data.ipynb` -- Jupyter notebook for exploring the data
* `DeSilva_combined.csv` -- data for all elements
* `Cu.csv`, `W.csv`, `Fe.csv`, `Al.csv` -- individual files for each element


### Exploratory Data Analysis

Here is what the data looks like:

![](image.png)


____
Thanks to Prof. Stephens at Texas Tech and Dr. Hansen at Sandia for the data. If you have data beyond what is here, please send it to me! 


