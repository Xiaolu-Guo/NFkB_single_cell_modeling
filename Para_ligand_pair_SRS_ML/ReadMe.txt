Para_pair_SRS_ML.ipynb

[to run the code]
put the codes and corresponding data files under the same filefolder, use Jupyter noteboot to run Para_pair_SRS_ML.ipynb


[version requirement]
Python version:
3.11.5 (main, Sep 11 2023, 08:19:27) [Clang 14.0.6 ]
csv version info: 1.0
Pandas version: 2.1.4
numpy version info: 1.26.0rc1
sklearn version info: 1.3.2
xgboost version info: 2.0.3

[codes description]
1. using parameter as input, the SRS score for pairwise ligand as output, to train the XGBoost regressor 
2. printed results are the R square of the regression model, and after thresholding (SRS>1 as distinguishing, SRS<=1 as confusion) the SRS score, the accuracy. 

