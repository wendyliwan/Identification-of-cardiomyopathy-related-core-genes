# Identification-of-cardiomyopathy-related-core-genes
---
## Identification of cardiomyopathy-related core genes through human metabolic networks and expression data
### Background:
Cardiomyopathy is a complex type of cardiomyopathy. It is very important to accurately distinguish the specific types of diseases in the treatment of cardiomyopathy.    
Dilated cardiomyopathy (DCM) and ischemic cardiomyopathy (ICM) are two common and indistinguishable types of cardiomyopathy. This is a comprehensive method based on   
a human metabolic network using cardiomyopathy expression data was proposed to identify cardiomyopathy-related core genes that could distinguish normal, DCM and ICM   
samples using cardiomyopathy expression profile data based on a human metabolic network.    
### Method：    
First, according to the differentially expressed genes between different states (DCM/ICM and normal, or DCM and ICM) of samples, three sets of initial modules were   
obtained from the human metabolic network. Two permutation tests were used to evaluate the significance of the Pearson correlation coefficient difference score of the  
initial modules, and three candidate modules were screened out. Then, a cardiomyopathy risk module that was significantly related to DCM and ICM was determined  
according to the significance of the module score based on Markov random field (MRFms). Finally, based on the shortest path between cardiomyopathy known genes, 13 core   
genes related to cardiomyopathy were identified.
The expression profile data can be obtained from the GEO Database, and the Metabolic network from the Virtual Metabolic Human Database.  
### Usage:  
This project uses R.  
Detection of candidate modules: The Pearson permutation test(degree conserved).R and Pearson permutation test(size conserved) files are code that runs two permutation  
tests to assess the significance of the Pearson correlation coefficient difference score of the initial modules.    

Identification of cardiomyopathy risk modules: MRFms.R  is a code that calculates MRFms of candidate modules  based on Markov model. The MRFms permutation test(size conserved).R   
and MRFms permutation test(degree conserved).R are programs for evaluating the significance of MRFms of candidate modules through two permutation tests. 
