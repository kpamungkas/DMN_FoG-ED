# DMN_FoG-ED
Default Mode Network of Freezing of Gait and Executive Dysfunction

This project is my master's thesis of Translational Neuroscience at the University of Würzburg, Germany.

In this project, I am investigating the default mode network role in freezing of gait episode and its association with executive dysfunction in Parkinson's diseased patients who received deep brain stimulation therapy.

This project is divided into two parts: 1) DBS-induced freezing of gait changes classification using default mode network connectivity profile, 2) DBS-induced executive dysfunction classification using previously FoG-trained models.

----
The first part of this project is designed to test default mode network models in solving deep brain stimulation-induced freezing of gait classification. Here I use volume of tissue activated form 51 PD patients with subthalamic nucleus DBS as input data. Prior to data analysis, pre-processing on individual data is done. Pre-processing pipeline is described as follows:
1) VTAs-based resting-state functional connectivity analysis on LeadDBS software using normative scans (74 PD patients; PPMI)
2) Network masking (default mode network, basal ganglia, and motor) of whole-brain resting state functional connectivity map (RSFC) based on AAL 3 brain atlas (MATLAB script: Preproc_Masking.m)
3) Extraction of features i.e., mean network connectivity and mean connectivity of network main hubs (MATLAB script: Preproc_FeaturesExtraction_DMN.m)

Afterwards, masked RSFC maps and their features are used for 3D convolutional neural network model and ordinal logistic regression model training, respectively. 
1) Ordinal logistic regression model (5-fold cross-validation)
    R script: OLR_Final_Revised.R
3) 3D convolutional neural network model (training/validation split validation)
    Google Colab script: 3DCNN_Final_DMN.ipynb (change path of files for Basal Ganglia and Motor models)

Models are assessed by their in-sample accuracy and balanced accuracy.


----
The second part of this project is designed to find DMN role in the association between freezing of gait episode and executive dysfunction. A third cohort (10 PD patients with STN-DBS) experiencing DBS-induced changes in executive dysfunction is utilized to test the models from the first part of the project. Models (OLR vs 3D CNN) are assessed by their accuracy and balanced accuracy.

R script: OLR_Final_Revised.R (Part 2)
Google Colab script: ED_3DCNN_Final_DMN.ipynb (change path of files for Basal Ganglia and Motor models)

----
Kristia Pamungkas
Schweinfurt, 21.06.21
