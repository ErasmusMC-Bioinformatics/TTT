## TTT project

This is an R Markdown document of TTT project. See the details in the manuscript *The nasal microbiota in the development of a multi-parametric prediction model to differentiate bacterial versus viral infections in lower respiratory tract infections*

    # set working directory
    # setwd("D:/TTT") # set working directory

    # install.packages("randomForestSRC")
    library(randomForestSRC)

    ## 
    ##  randomForestSRC 2.11.0 
    ##  
    ##  Type rfsrc.news() to see new features, changes, and bug fixes. 
    ## 

## Internal Evaluation

In the internal evaluation phase, we performed feature reduction based
on 242 patients in the initial TAILORED-Treatment cohort (see main
text). In total, 29 clinical variables and 7 microbiota variables
(Table\_S3.xlsx) entered the classifier modelling phase to build CC, CE
and CEM –&gt; CC: Classifier using CRP only in the initial cohort. CE:
Classifiers using two or more eCRF variables (incl. CRP) in the initial
cohort. CEM: Classifiers using all input eCRF variables (incl. CRP) and
at least one microbiota in the initial cohort.

Rank the 29 eCRF variables according to their variable importance
    
    Train_eCRF <- read.delim("Input/TTT_eCRF_242Train21CatFeature_8NumFeature.txt",sep="\t",header=T,row.names=1)
    character_vars <- lapply(Train_eCRF, class) == "character"
    Train_eCRF[, character_vars] <- lapply(Train_eCRF[, character_vars], as.factor)

    err_lowest = 1
    for (i in 1:1000) {
      #print(i)
      Train.temp <- rfsrc(Class ~ ., Train_eCRF, ntree = 500,
                          mtry = NULL, importance = "none", block.size = 1, ensemble = "oob",
                          bootstrap = "by.root", samptype = "swr",
                          na.action = "na.impute", nimpute = 1,
                          proximity = "inbag", forest.wt = "inbag", forest = TRUE)
      
      err <- Train.temp$err.rate
      if (err[nrow(err),1] < err_lowest) {
        err_lowest <- err[nrow(err),1]
        Train.classifier <- Train.temp
      }
    }

    VI_blocksize1 <- vimp(Train.classifier, importance = "permute", block.size = 1)$importance  # Obtain the eCRF variable importance
    write.csv(as.data.frame(VI_blocksize1),file="Result/TTT_eCRF_242Train21CatFeature_8NumFeature_VIMP1.csv")

Similarly, rank the 7 nasal cavity microbiota variables according to their variable importance

    Train_micro <- read.delim("Input/TTT_microbiome_242Train_7genera.txt",sep="\t",header=T,row.names=1)
    Train_micro$Class <- as.factor(Train_micro$Class)

    err_lowest = 1
    for (i in 1:1000) {
      #print(i)
      Train.temp <- rfsrc(Class ~ ., Train_micro, ntree = 500,
                          mtry = NULL, importance = "none", block.size = 1, ensemble = "oob",
                          bootstrap = "by.root", samptype = "swr",
                          na.action = "na.impute", nimpute = 1,
                          proximity = "inbag", forest.wt = "inbag", forest = TRUE)
      
      err <- Train.temp$err.rate
      if (err[nrow(err),1] < err_lowest) {
        err_lowest <- err[nrow(err),1]
        Train.classifier <- Train.temp
      }
    }

    VI_blocksize1 <- vimp(Train.classifier, importance = "permute", block.size = 1)$importance  # Obtain the microbiome variable importance
    write.csv(as.data.frame(VI_blocksize1),file="Result/TTT_micro_242Train_7genera_VIMP1.csv")

Training set of eCRF + microbiome features, sorted by the VIMP in the eCRF and microbiome obtained above

    Train_data <- read.delim("Input/TTT_242Train_eCRF29Feature_SortByVIMP1_7genera_SortByVIMP1.txt",sep="\t",header=T,row.names=1)
    character_vars <- lapply(Train_data, class) == "character"
    Train_data[, character_vars] <- lapply(Train_data[, character_vars], as.factor)

    R = 10 #Number of repetitions
    results = c("Number of top features","Iteration","AUC_train","Accuracy_bacterial_train","Accuracy_viral_train")

Build CC, CE and CEM in the initial cohort, using top 1 (CRP), 29 (eCRF)
and all 36 (eCRF+microbiom) sorted features

    for (N in 2:37) {    #Number of features (=N-1) to be included in the classifier
      print(N-1)
      for (r in 1:R) { #Iteration
        # Build the classifier and select the best one
        err_lowest = 1
        for (i in 1:1000) {
          Train.temp <- rfsrc(Class ~ ., Train_data[,1:N], ntree = 500,   
                              mtry = NULL, importance = "none", block.size = 1, ensemble = "oob",
                              bootstrap = "by.root", samptype = "swr",
                              na.action = "na.impute", nimpute = 1,
                              proximity = "inbag", forest.wt = "inbag", forest = TRUE)
          
          err <- Train.temp$err.rate
          if (mean(err[,1]) < err_lowest) {
            err_lowest <- mean(err[,1])
            Train.classifier <- Train.temp
          }
        }
        #Retrieve the results
        Confusion_train <- get.confusion(Train_data$Class,Train.classifier$predicted.oob)
        AUC_train <- get.auc(Train_data$Class,Train.classifier$predicted.oob)
        results = rbind(results,c(N-1,r,AUC_train,Confusion_train[1,1]/sum(Train_data$Class=="bacterial"),Confusion_train[2,2]/sum(Train_data$Class=="viral"))) 
      }
    }

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8
    ## [1] 9
    ## [1] 10
    ## [1] 11
    ## [1] 12
    ## [1] 13
    ## [1] 14
    ## [1] 15
    ## [1] 16
    ## [1] 17
    ## [1] 18
    ## [1] 19
    ## [1] 20
    ## [1] 21
    ## [1] 22
    ## [1] 23
    ## [1] 24
    ## [1] 25
    ## [1] 26
    ## [1] 27
    ## [1] 28
    ## [1] 29
    ## [1] 30
    ## [1] 31
    ## [1] 32
    ## [1] 33
    ## [1] 34
    ## [1] 35
    ## [1] 36

    write.csv(results,file="Result/TTT_InitialCohort_CC_CE_CEM.csv")

## Cross Validation

In the expanded cohort with 51 extra patients, 5-fold cross-validation
was conducted to evaluate the contribution of the eCRF and microbiota
variables to the prediction performance (see main text) of CC*, CEM* and
CCM\* –&gt; CC*: Classifier using CRP only in the 5-fold CV of the
expanded cohort. CEM*: Classifiers using two or more variables
(regardless eCRF or microbiota) in the 5-fold CV of the expanded cohort.
CCM\*: Classifiers using CRP and all input microbiota variables in the
5-fold CV of the expanded cohort.

1.The total cohort of 293 patients was randomly split into 5 subsets, or
so called ‘folds’. 

2.Since the two classes were unbalanced, the majority
class label ‘viral infection’ was subsampled to reach equal prior
probability as class ‘bacterial infection’ in each training set. 

3.For each training fold, the 36 input features were ranked simultaneously
according to their variable importance, as described above

    for (F in 1:5) {    #Fold
      print(paste("Fold =",F))
      Train_fold <- read.delim(paste("Input/TTT_EqualPrior_5foldCVtrain",F,".txt",sep = ""),sep="\t",header=T,row.names=1)
      character_vars <- lapply(Train_fold, class) == "character"
      Train_fold[, character_vars] <- lapply(Train_fold[, character_vars], as.factor)

      # Use all features simultaneously to get the VIMP
      N = 37
      err_lowest = 1
      for (i in 1:1000) {
        Train.temp <- rfsrc(Class ~ ., Train_fold[,1:N], ntree = 500,
                          mtry = NULL, importance = "none", block.size = 1, ensemble = "oob",
                          bootstrap = "by.root", samptype = "swr",
                          na.action = "na.impute", nimpute = 1,
                          proximity = "inbag", forest.wt = "inbag", forest = TRUE)
      
        err <- Train.temp$err.rate
        if (mean(err[,1]) < err_lowest) {
          err_lowest <- mean(err[,1])
          Train.classifier <- Train.temp
        }
      }

      VI_blocksize1 <- vimp(Train.classifier, importance = "permute", block.size = 1)$importance
      write.csv(as.data.frame(VI_blocksize1),file=paste("Result/TTT_EqualPrior_5foldCVtrain",F,"_VIMP1.csv",sep = ""))
    }

    ## [1] "Fold = 1"
    ## [1] "Fold = 2"
    ## [1] "Fold = 3"
    ## [1] "Fold = 4"
    ## [1] "Fold = 5"

4.Sort the features in the training and test folds according to the VIMP
in the corresponding training fold 

5.Now load the sorted training folds
and test folds to build CC\* and CEM\*

    results = c("Fold","Number of top features","AUC_train","Accuracy_bacterial_train","Accuracy_viral_train","AUC_test","Accuracy_bacterial_test","Accuracy_viral_test")

    for (F in 1:5) {    #Fold
      for (N in 2:37) {    #Number of features (=N-1) to be included in the classifier
        print(paste("Fold =",F, "; #Featue =", N-1))
        Train_fold <- read.delim(paste("Input/TTT_EqualPrior_5foldCVtrain",F,"_SortByEqualPriorVIMP1.txt",sep = ""),sep="\t",header=T,row.names=1)
        character_vars <- lapply(Train_fold, class) == "character"
        Train_fold[, character_vars] <- lapply(Train_fold[, character_vars], as.factor)
        
        Test_fold <- read.delim(paste("Input/TTT_5foldCVtest",F,"_SortByEqualPriorVIMP1.txt",sep = ""),sep="\t",header=T,row.names=1)
        character_vars <- lapply(Test_fold, class) == "character"
        Test_fold[, character_vars] <- lapply(Test_fold[, character_vars], as.factor)
        
        # Build the classifier many times (=i) and select the best one
        err_lowest = 1
        for (i in 1:1000) {
          Train.temp <- rfsrc(Class ~ ., Train_fold[,1:N], ntree = 500,
                              mtry = NULL, importance = "none", block.size = 1, ensemble = "oob",
                              bootstrap = "by.root", samptype = "swr",
                              na.action = "na.impute", nimpute = 1,
                              proximity = "inbag", forest.wt = "inbag", forest = TRUE)
          
          err <- Train.temp$err.rate
          if (mean(err[,1]) < err_lowest) {
            err_lowest <- mean(err[,1])
            Train.classifier <- Train.temp
          }
        }
        
        # Deploy the best classifier on the test fold
        Test.pred <-predict.rfsrc(Train.classifier,
                                  Test_fold[,1:N],
                                  na.action = "na.impute",
                                  proximity = TRUE,forest.wt = TRUE)
        
        #Retrieve the results
        Confusion_train <- get.confusion(Train_fold$Class,Train.classifier$predicted.oob)
        AUC_train <- get.auc(Train_fold$Class,Train.classifier$predicted.oob)
        Confusion_test <- get.confusion(Test_fold$Class,Test.pred$predicted)
        AUC_test <- get.auc(Test_fold$Class,Test.pred$predicted)
        
        results = rbind(results,c(F,N-1,AUC_train,Confusion_train[1,1]/sum(Train_fold$Class=="bacterial"),Confusion_train[2,2]/sum(Train_fold$Class=="viral"),AUC_test,Confusion_test[1,1]/sum(Test_fold$Class=="bacterial"),Confusion_test[2,2]/sum(Test_fold$Class=="viral"))) 
      }
    }

    ## [1] "Fold = 1 ; #Featue = 1"
    ## [1] "Fold = 1 ; #Featue = 2"
    ## [1] "Fold = 1 ; #Featue = 3"
    ## [1] "Fold = 1 ; #Featue = 4"
    ## [1] "Fold = 1 ; #Featue = 5"
    ## [1] "Fold = 1 ; #Featue = 6"
    ## [1] "Fold = 1 ; #Featue = 7"
    ## [1] "Fold = 1 ; #Featue = 8"
    ## [1] "Fold = 1 ; #Featue = 9"
    ## [1] "Fold = 1 ; #Featue = 10"
    ## [1] "Fold = 1 ; #Featue = 11"
    ## [1] "Fold = 1 ; #Featue = 12"
    ## [1] "Fold = 1 ; #Featue = 13"
    ## [1] "Fold = 1 ; #Featue = 14"
    ## [1] "Fold = 1 ; #Featue = 15"
    ## [1] "Fold = 1 ; #Featue = 16"
    ## [1] "Fold = 1 ; #Featue = 17"
    ## [1] "Fold = 1 ; #Featue = 18"
    ## [1] "Fold = 1 ; #Featue = 19"
    ## [1] "Fold = 1 ; #Featue = 20"
    ## [1] "Fold = 1 ; #Featue = 21"
    ## [1] "Fold = 1 ; #Featue = 22"
    ## [1] "Fold = 1 ; #Featue = 23"
    ## [1] "Fold = 1 ; #Featue = 24"
    ## [1] "Fold = 1 ; #Featue = 25"
    ## [1] "Fold = 1 ; #Featue = 26"
    ## [1] "Fold = 1 ; #Featue = 27"
    ## [1] "Fold = 1 ; #Featue = 28"
    ## [1] "Fold = 1 ; #Featue = 29"
    ## [1] "Fold = 1 ; #Featue = 30"
    ## [1] "Fold = 1 ; #Featue = 31"
    ## [1] "Fold = 1 ; #Featue = 32"
    ## [1] "Fold = 1 ; #Featue = 33"
    ## [1] "Fold = 1 ; #Featue = 34"
    ## [1] "Fold = 1 ; #Featue = 35"
    ## [1] "Fold = 1 ; #Featue = 36"
    ## [1] "Fold = 2 ; #Featue = 1"
    ## [1] "Fold = 2 ; #Featue = 2"
    ## [1] "Fold = 2 ; #Featue = 3"
    ## [1] "Fold = 2 ; #Featue = 4"
    ## [1] "Fold = 2 ; #Featue = 5"
    ## [1] "Fold = 2 ; #Featue = 6"
    ## [1] "Fold = 2 ; #Featue = 7"
    ## [1] "Fold = 2 ; #Featue = 8"
    ## [1] "Fold = 2 ; #Featue = 9"
    ## [1] "Fold = 2 ; #Featue = 10"
    ## [1] "Fold = 2 ; #Featue = 11"
    ## [1] "Fold = 2 ; #Featue = 12"
    ## [1] "Fold = 2 ; #Featue = 13"
    ## [1] "Fold = 2 ; #Featue = 14"
    ## [1] "Fold = 2 ; #Featue = 15"
    ## [1] "Fold = 2 ; #Featue = 16"
    ## [1] "Fold = 2 ; #Featue = 17"
    ## [1] "Fold = 2 ; #Featue = 18"
    ## [1] "Fold = 2 ; #Featue = 19"
    ## [1] "Fold = 2 ; #Featue = 20"
    ## [1] "Fold = 2 ; #Featue = 21"
    ## [1] "Fold = 2 ; #Featue = 22"
    ## [1] "Fold = 2 ; #Featue = 23"
    ## [1] "Fold = 2 ; #Featue = 24"
    ## [1] "Fold = 2 ; #Featue = 25"
    ## [1] "Fold = 2 ; #Featue = 26"
    ## [1] "Fold = 2 ; #Featue = 27"
    ## [1] "Fold = 2 ; #Featue = 28"
    ## [1] "Fold = 2 ; #Featue = 29"
    ## [1] "Fold = 2 ; #Featue = 30"
    ## [1] "Fold = 2 ; #Featue = 31"
    ## [1] "Fold = 2 ; #Featue = 32"
    ## [1] "Fold = 2 ; #Featue = 33"
    ## [1] "Fold = 2 ; #Featue = 34"
    ## [1] "Fold = 2 ; #Featue = 35"
    ## [1] "Fold = 2 ; #Featue = 36"
    ## [1] "Fold = 3 ; #Featue = 1"
    ## [1] "Fold = 3 ; #Featue = 2"
    ## [1] "Fold = 3 ; #Featue = 3"
    ## [1] "Fold = 3 ; #Featue = 4"
    ## [1] "Fold = 3 ; #Featue = 5"
    ## [1] "Fold = 3 ; #Featue = 6"
    ## [1] "Fold = 3 ; #Featue = 7"
    ## [1] "Fold = 3 ; #Featue = 8"
    ## [1] "Fold = 3 ; #Featue = 9"
    ## [1] "Fold = 3 ; #Featue = 10"
    ## [1] "Fold = 3 ; #Featue = 11"
    ## [1] "Fold = 3 ; #Featue = 12"
    ## [1] "Fold = 3 ; #Featue = 13"
    ## [1] "Fold = 3 ; #Featue = 14"
    ## [1] "Fold = 3 ; #Featue = 15"
    ## [1] "Fold = 3 ; #Featue = 16"
    ## [1] "Fold = 3 ; #Featue = 17"
    ## [1] "Fold = 3 ; #Featue = 18"
    ## [1] "Fold = 3 ; #Featue = 19"
    ## [1] "Fold = 3 ; #Featue = 20"
    ## [1] "Fold = 3 ; #Featue = 21"
    ## [1] "Fold = 3 ; #Featue = 22"
    ## [1] "Fold = 3 ; #Featue = 23"
    ## [1] "Fold = 3 ; #Featue = 24"
    ## [1] "Fold = 3 ; #Featue = 25"
    ## [1] "Fold = 3 ; #Featue = 26"
    ## [1] "Fold = 3 ; #Featue = 27"
    ## [1] "Fold = 3 ; #Featue = 28"
    ## [1] "Fold = 3 ; #Featue = 29"
    ## [1] "Fold = 3 ; #Featue = 30"
    ## [1] "Fold = 3 ; #Featue = 31"
    ## [1] "Fold = 3 ; #Featue = 32"
    ## [1] "Fold = 3 ; #Featue = 33"
    ## [1] "Fold = 3 ; #Featue = 34"
    ## [1] "Fold = 3 ; #Featue = 35"
    ## [1] "Fold = 3 ; #Featue = 36"
    ## [1] "Fold = 4 ; #Featue = 1"
    ## [1] "Fold = 4 ; #Featue = 2"
    ## [1] "Fold = 4 ; #Featue = 3"
    ## [1] "Fold = 4 ; #Featue = 4"
    ## [1] "Fold = 4 ; #Featue = 5"
    ## [1] "Fold = 4 ; #Featue = 6"
    ## [1] "Fold = 4 ; #Featue = 7"
    ## [1] "Fold = 4 ; #Featue = 8"
    ## [1] "Fold = 4 ; #Featue = 9"
    ## [1] "Fold = 4 ; #Featue = 10"
    ## [1] "Fold = 4 ; #Featue = 11"
    ## [1] "Fold = 4 ; #Featue = 12"
    ## [1] "Fold = 4 ; #Featue = 13"
    ## [1] "Fold = 4 ; #Featue = 14"
    ## [1] "Fold = 4 ; #Featue = 15"
    ## [1] "Fold = 4 ; #Featue = 16"
    ## [1] "Fold = 4 ; #Featue = 17"
    ## [1] "Fold = 4 ; #Featue = 18"
    ## [1] "Fold = 4 ; #Featue = 19"
    ## [1] "Fold = 4 ; #Featue = 20"
    ## [1] "Fold = 4 ; #Featue = 21"
    ## [1] "Fold = 4 ; #Featue = 22"
    ## [1] "Fold = 4 ; #Featue = 23"
    ## [1] "Fold = 4 ; #Featue = 24"
    ## [1] "Fold = 4 ; #Featue = 25"
    ## [1] "Fold = 4 ; #Featue = 26"
    ## [1] "Fold = 4 ; #Featue = 27"
    ## [1] "Fold = 4 ; #Featue = 28"
    ## [1] "Fold = 4 ; #Featue = 29"
    ## [1] "Fold = 4 ; #Featue = 30"
    ## [1] "Fold = 4 ; #Featue = 31"
    ## [1] "Fold = 4 ; #Featue = 32"
    ## [1] "Fold = 4 ; #Featue = 33"
    ## [1] "Fold = 4 ; #Featue = 34"
    ## [1] "Fold = 4 ; #Featue = 35"
    ## [1] "Fold = 4 ; #Featue = 36"
    ## [1] "Fold = 5 ; #Featue = 1"
    ## [1] "Fold = 5 ; #Featue = 2"
    ## [1] "Fold = 5 ; #Featue = 3"
    ## [1] "Fold = 5 ; #Featue = 4"
    ## [1] "Fold = 5 ; #Featue = 5"
    ## [1] "Fold = 5 ; #Featue = 6"
    ## [1] "Fold = 5 ; #Featue = 7"
    ## [1] "Fold = 5 ; #Featue = 8"
    ## [1] "Fold = 5 ; #Featue = 9"
    ## [1] "Fold = 5 ; #Featue = 10"
    ## [1] "Fold = 5 ; #Featue = 11"
    ## [1] "Fold = 5 ; #Featue = 12"
    ## [1] "Fold = 5 ; #Featue = 13"
    ## [1] "Fold = 5 ; #Featue = 14"
    ## [1] "Fold = 5 ; #Featue = 15"
    ## [1] "Fold = 5 ; #Featue = 16"
    ## [1] "Fold = 5 ; #Featue = 17"
    ## [1] "Fold = 5 ; #Featue = 18"
    ## [1] "Fold = 5 ; #Featue = 19"
    ## [1] "Fold = 5 ; #Featue = 20"
    ## [1] "Fold = 5 ; #Featue = 21"
    ## [1] "Fold = 5 ; #Featue = 22"
    ## [1] "Fold = 5 ; #Featue = 23"
    ## [1] "Fold = 5 ; #Featue = 24"
    ## [1] "Fold = 5 ; #Featue = 25"
    ## [1] "Fold = 5 ; #Featue = 26"
    ## [1] "Fold = 5 ; #Featue = 27"
    ## [1] "Fold = 5 ; #Featue = 28"
    ## [1] "Fold = 5 ; #Featue = 29"
    ## [1] "Fold = 5 ; #Featue = 30"
    ## [1] "Fold = 5 ; #Featue = 31"
    ## [1] "Fold = 5 ; #Featue = 32"
    ## [1] "Fold = 5 ; #Featue = 33"
    ## [1] "Fold = 5 ; #Featue = 34"
    ## [1] "Fold = 5 ; #Featue = 35"
    ## [1] "Fold = 5 ; #Featue = 36"

    write.csv(results,file="Result/TTT_CrossValidation_CC_CEM.csv")

6.Now load the sorted training folds and test folds containing only CRP
and 7 genera to build CCM\*

    results = c("Fold","Number of features","AUC_train","Accuracy_bacterial_train","Accuracy_viral_train","AUC_test","Accuracy_bacterial_test","Accuracy_viral_test")

    for (F in 1:5) {    #Fold
        print(paste("Fold =",F))
        Train_fold <- read.delim(paste("Input/TTT_EqualPrior_5foldCVtrain",F,"_SortByEqualPriorVIMP1_CRPand7genera.txt",sep = ""),sep="\t",header=T,row.names=1)
        Train_fold$Class <- as.factor(Train_fold$Class)
        
        Test_fold <- read.delim(paste("Input/TTT_5foldCVtest",F,"_SortByEqualPriorVIMP1_CRPand7genera.txt",sep = ""),sep="\t",header=T,row.names=1)
        Test_fold$Class <- as.factor(Test_fold$Class)
        
        # Build the classifier many times (=i) and select the best one
        err_lowest = 1
        for (i in 1:1000) {
          Train.temp <- rfsrc(Class ~ ., Train_fold, ntree = 500,
                              mtry = NULL, importance = "none", block.size = 1, ensemble = "oob",
                              bootstrap = "by.root", samptype = "swr",
                              na.action = "na.impute", nimpute = 1,
                              proximity = "inbag", forest.wt = "inbag", forest = TRUE)
          
          err <- Train.temp$err.rate
          if (mean(err[,1]) < err_lowest) {
            err_lowest <- mean(err[,1])
            Train.classifier <- Train.temp
          }
        }
        
        # Deploy the best classifier on the test fold
        Test.pred <-predict.rfsrc(Train.classifier,
                                  Test_fold,
                                  na.action = "na.impute",
                                  proximity = TRUE,forest.wt = TRUE)
        
        #Retrieve the results
        Confusion_train <- get.confusion(Train_fold$Class,Train.classifier$predicted.oob)
        AUC_train <- get.auc(Train_fold$Class,Train.classifier$predicted.oob)
        Confusion_test <- get.confusion(Test_fold$Class,Test.pred$predicted)
        AUC_test <- get.auc(Test_fold$Class,Test.pred$predicted)
        
        results = rbind(results,c(F,8,AUC_train,Confusion_train[1,1]/sum(Train_fold$Class=="bacterial"),Confusion_train[2,2]/sum(Train_fold$Class=="viral"),AUC_test,Confusion_test[1,1]/sum(Test_fold$Class=="bacterial"),Confusion_test[2,2]/sum(Test_fold$Class=="viral"))) 
      
    }

    ## [1] "Fold = 1"
    ## [1] "Fold = 2"
    ## [1] "Fold = 3"
    ## [1] "Fold = 4"
    ## [1] "Fold = 5"

    write.csv(results,file="Result/TTT_CrossValidation_CCM.csv")
