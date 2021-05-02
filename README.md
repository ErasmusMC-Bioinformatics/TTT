# TTT project
## Introduction
Here we give the R script of the implementation of the prediction models described in the study “The nasal microbiota in the development of a multi-parametric prediction model to differentiate bacterial versus viral infections in lower respiratory tract infections”. In this study, we investigated whether adding nasal cavity microbiota profiles on top of clinical variables increased the predictive value of a machine learning classifier developed to distinguish between bacterial and viral infection in patients with lower respiratory tract infections. 

A cohort of 242 patients were included in the internal evaluation phase according to the date of recruitment to the TAILORED-Treatment study. This cohort was used to compare the prediction performances of the classifiers using clinical variables alone, as well as classifiers using both clinical and microbiota variables (Internal evaluation phase). In the expanded cohort (51 extra patients), 5-fold cross-validation (CV) analysis was conducted to evaluate the contribution of clinical and microbiota variables to prediction performance (Cross-validation phase).

## Data
Clinical data and microbiota data of 293 patients recruited in the TAILORED-Treatment study were used. The clinical data were generated using standard laboratory methods according to the published TAILORED-Treatment clinical trial protocol [[1]](#1). All clinical information as listed in the eCRFs was used by the expert panel to determine the etiology of infection. Microbiota data was generated from nasal cavity swab samples by 16S rRNA gene sequencing.

## Model
Random Forests classification method was implemented using R package *randomForestSRC*.

In the internal evaluation phase, 242 patients were included to compare the prediction performances of classifiers using clinical eCRF variables alone (CC and CE) and classifiers using both eCRF and microbiota variables (CEM). More specifically, we first ranked the variable importance of the eCRF variables and microbiota variables separately by function vimp in the initial cohort. Then the ranked eCRF variables were included incrementally into the Random Forests model followed by the ranked microbiota variables. Each time 1000 forests were grown and the best forest was kept. This was iterated 10 times to evaluate the robustness of the model. 

In the expanded cohort with 51 extra patients, 5-fold cross-validation was conducted to evaluate the contribution of the clinical eCRF and microbiota variables to the prediction performance. The total cohort of 293 patients was randomly split into 5 subsets, or so called ‘folds’. Each time one fold containing 58 or 59 patients was taken as the test set and the remaining were part of the training set. Since the two classes were unbalanced, the majority class label ‘viral infection’ was subsampled to reach equal prior probability as class ‘bacterial infection’ in each training set. Unlike the internal evaluation phase where the eCRF and microbiota variables were ranked separately and added sequentially (i.e. first the eCRF variables and then the microbiota ones), here all variables were ranked simultaneously based on the training set at hand using function *vimp* and added into the model accordingly. In both phases, the performance of the models was assessed by calculating the area under the receiver operating characteristic curve (AUC) for overall performance and the percentage of correctly predicted cases. In the 5-fold CV, this was performed for both training set (using ‘out-of-bag’ predictions) and test set at hand. 

## Reference
<a id="1">[1]</a> 
van Houten CB, Oved K, Eden E, Cohen A, Engelhard D, Boers S, et al. 
Observational multi-centre, prospective study to characterize novel pathogen-and host-related factors in hospitalized patients with lower respiratory tract infections and/or sepsis - the "TAILORED-Treatment" study. 
BMC Infect Dis. 2018;18(1):377.
