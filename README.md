# HCMHS
HCMHS: High-Order SNP Interactions Detection Based on Hierarchical Clustering and Multi-Task Harmony Search Algorithm

# main codes:
The experimental results were attained by running matlab codes.

# Abstract
The interaction between SNPs plays a key role in revealing the genetic mechanisms of complex diseases. However, as the order of SNP interactions increases, the number of SNP combinations increases exponentially, presenting a serious combinatorial explosion problem. Although swarm intelligence algorithms can alleviate the combinatorial explosion problem by optimizing search paths and have been widely used in SNP interaction detection, the performance of swarm intelligence algorithms is greatly affected by initialization and search direction. To mitigate these issues, we propose a high-order SNP interaction detection algorithm based on hierarchical clustering and multi-task harmony search (abbreviated as HCMHS). In harmony memory initialization, considering that similar SNPs are more likely to form pathogenic combinations, hierarchical clustering is first used to cluster SNPs into different clusters, and then based on these clusters, harmony memory is initialized. In  harmony search, considering the correlation between different orders of SNP interactions, the search directions for SNPs of various orders are optimized through a multi-task framework. To verify the performance of HCMHS, we carried out experiments on 58 simulated datasets and one real dataset, and HCMHS achieved the best results compared to nine advanced algorithms.

# Keywords:
High-Order SNP Interactions, Harmony Search Algorithm, Hierarchical Clustering


code list:
(A modified version of HCMHS)
  ============
     |-- cengci.m  :   the code for hierarchical clustering
    
     |-- dataexample.txt : a sample file of the dataset
   
     |-- Gtest_score.m : the code for g-test test
   
     |-- HS_2021_multiTask_UnifiedCoding2021.m : the score for multitask search
   
     |-- multi_criteriaEvaluationFuns2021.m: the score for multiple fitness function 
   
     |-- test.m  : the score for program startup code
