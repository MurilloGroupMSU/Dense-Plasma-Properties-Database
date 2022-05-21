
____
## DeSilva Data
____

This database contains electrical conductivity data from exploding wire experiments conducted by DeSilva and co-workers. The original publications are given as a PDF in this repo. Specifially, the database contains elements:
* Al,
* Fe,
* Ni,
* Cu,
* W.

The data is in the form of {Z_{nuc}, $\rho$, T, $\sigma$} with units {none, g/cc, K, S/m}. The numerical values in the CSV files are transformed such that the base-10 logarithm of the temperature and conductivity are given. (DeSilva and Vunni include other elements, which are not given here because they are in terms of internal energy rather than temperature.)


### Files Included 

I have organized the files in different ways and have provided a Jupyter notebook that allows one to explore the data and perhaps allow the user to create customized files. These files are:
* `DeSilvaKatsouros.pdf` -- original DeSilva paper
* `View_DeSilva_data.ipynb` -- Jupyter notebook for exploring the data
* `DeSilva_combined.csv` -- data for all elements
* `Cu_cond.csv`, `W_cond.csv`, `Fe_cond.csv`, `Al_cond.csv` -- individual files for each element


### Exploratory Data Analysis

Here is what the data looks like:

![](combined_plot.png)



### Modified Lee-More Model (MLM)

This data was originally collected here as part of the paper,
* _Data-driven Electrical Conductivities of Dense Plasmas_, Michael S. Murillo, Frontiers in Physics, 2022,
which also included the MLM, a modfied version of the Lee-More conductivity model, more of which can be found [here](https://murillogroupmsu.github.io/Modified-Lee-More-Transport/).

____
Thanks to Prof. Stephens at Texas Tech and Dr. Hansen at Sandia for the data. If you have electrical conductivity data beyond what is here, please send it to me!

Special thanks to Prof. DeSilva for many useful discussions regarding his data. 



