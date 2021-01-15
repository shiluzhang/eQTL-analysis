# Expression QTL analysis in expression genetic data



## 1. Introduction
Identifying the genetic loci that contribute to variation in quantitative traits in nature is an important problem to biologists. Studying the effects of such quantitative trait loci (QTLs) can help understand the biological basis and evolution of these traits. Expression quantitative trait loci (eQTL) refer to the genomic locations that influence variation in gene expression levels (mRNA abundances) (Tian et al., 2015). Expression QTL near the genomic location of the influenced gene are called local eQTL, while eQTL that are far away from the influenced gene are called trans-eQTL. Measuring gene expression in disease-relevant tissues in QTL experiments is of great interest because gene expressions mapping to the same location may indicate the existence of regulators that cause the disease (Tian et al., 2015). In expression genetics studies, genome-wide gene expression is assayed along with genotypes at genetic markers to identify eQTL with profound effect on gene expression (Broman et al., 2015). Some of the expression quantitative trait loci showed effects in multiple tissues, whereas some were specific to a single tissue (Tian et al., 2015). 

To identify the features of local expression quantitative trait loci, where a gene’s mRNA abundance is strongly associated with genotype near its genomic location, we focus on an F2 intercross experiment between diabetes-resistant C57BL/6J (abbreviated B6 or B) and diabetes-susceptible BTBR T+ tf/J (abbreviated BTBR or R) mouse lines, with gene expression microarray data on six tissues (adipose, gastrocnemius muscle, hypothalamus, pancreatic islets, kidney, and liver). F2 mice were genotyped with the 5K GeneChip (Affymetrix). All mice were genetically obese through introgression of the leptin mutation Lep-ob/ob (Tian et al., 2015). Leptin regulates appetite by signaling to the brain that the animal has had enough to eat. Since the Lep-ob/ob mouse cannot produce leptin, its food intake is uncontrollable by this mechanism. 

Combining the gene expression microarray data on six tissues with the genotype data, I am interested in the following questions: 
1. Are there additive allele effects or dominance effects in the quantitative trait loci analysis? 
2. Are the effects of expression QTL different in the two sexes? 
3. Do some pairs of tissues have similar effects of expression QTL? If it is expression QTL in one tissue, is it more likely to be expression QTL in other tissues?

