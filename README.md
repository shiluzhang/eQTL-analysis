# Expression QTL analysis in expression genetic data



## 1. Introduction
Identifying the genetic loci that contribute to variation in quantitative traits in nature is an important problem to biologists. Studying the effects of such quantitative trait loci (QTLs) can help understand the biological basis and evolution of these traits. Expression quantitative trait loci (eQTL) refer to the genomic locations that influence variation in gene expression levels (mRNA abundances) (Tian et al., 2015). Expression QTL near the genomic location of the influenced gene are called local eQTL, while eQTL that are far away from the influenced gene are called trans-eQTL. Measuring gene expression in disease-relevant tissues in QTL experiments is of great interest because gene expressions mapping to the same location may indicate the existence of regulators that cause the disease (Tian et al., 2015). In expression genetics studies, genome-wide gene expression is assayed along with genotypes at genetic markers to identify eQTL with profound effect on gene expression (Broman et al., 2015). Some of the expression quantitative trait loci showed effects in multiple tissues, whereas some were specific to a single tissue (Tian et al., 2015). 

To identify the features of local expression quantitative trait loci, where a gene’s mRNA abundance is strongly associated with genotype near its genomic location, we focus on an F2 intercross experiment between diabetes-resistant C57BL/6J (abbreviated B6 or B) and diabetes-susceptible BTBR T+ tf/J (abbreviated BTBR or R) mouse lines, with gene expression microarray data on six tissues (adipose, gastrocnemius muscle, hypothalamus, pancreatic islets, kidney, and liver). F2 mice were genotyped with the 5K GeneChip (Affymetrix). All mice were genetically obese through introgression of the leptin mutation Lep-ob/ob (Tian et al., 2015). Leptin regulates appetite by signaling to the brain that the animal has had enough to eat. Since the Lep-ob/ob mouse cannot produce leptin, its food intake is uncontrollable by this mechanism. 

Combining the gene expression microarray data on six tissues with the genotype data, I am interested in the following questions: 
1. Are there additive allele effects or dominance effects in the quantitative trait loci analysis? 
2. Are the effects of expression QTL different in the two sexes? 
3. Do some pairs of tissues have similar effects of expression QTL? If it is expression QTL in one tissue, is it more likely to be expression QTL in other tissues?


## 2. Methods:
### 2.1 Genotyping data:
To identify genes and pathways that cause obesity-induced type II diabetes, two different mouse lines: C57BL/6J (abbreviated B6 or B) and BTBR T+ tf/J (abbreviated BTBR or R) mice were used to conduct an intercross. B6 mice are resistant to diabetes, while BTBR mice are susceptible to diabetes. BTBR females were crossed to B6 males to generate F1 heterozygotes, and F1 parents were crossed to generate F2 mice. All mice were genetically obese through introgression of the leptin mutation Lepob/ob. F2 mice were genotyped with the 5K GeneChip (Affymetrix). Totally, there were 519 F2 mice genotyped at 2057 informative markers. There were three phenotypes in F2 mice, BB (homozygous B6), BR (heterozygous) and RR (homozygous BTBR). 

Below is an example of part of the genotyping data. Mouse3051 is one of the 519 F2 mice. rs13475697, rs3681603, rs13475703, and rs13475710 are four of the 2,057 informative markers. The first row denotes which of the chromosomes each marker is on. The second row is the position of the marker on that chromosome in centimorgan (cM). (A centiMorgan (cM) is a unit of recombinant frequency which is used to measure genetic distance.)

| MouseNum  | Sex | rs13475697 | rs3681603 | rs13475703	| rs13475710 |
| ------------- | ------------- |------------- | ------------- |------------- | ------------- |
|   |   | 1 | 1 | 1 | 1 | 
|   |   | 1.6449 | 1.648 | 1.86 | 2.078| 
| Mouse3051  | Male  | RR | RR | RR | RR | 
| Mouse3551  | Male  | BR | BR |BR | BR |

Table 1: An example of part of the genotyping data.

### 2.2 Gene expression data:
Gene expression was assayed with custom two-color, ink-jet microarrays manufactured by Agilent Technologies (Palo Alto, CA). RNA preparations were performed at Rosetta Inpharmatics (Merck & Co.). Six tissues from F2 mice were considered for expression profiling: adipose, gastrocnemius muscle, hypothalamus, pancreatic islets, liver, and kidney. Tissue-specific messenger RNA (mRNA) pools for each tissue were used for the reference channel, and gene expression was quantified as the ratio of the mean log10 intensity (mlratio). In the datasets, there were 519 mice with gene expression data on at least one of the six tissues (487 mice for adipose, 490 mice for gastrocnemius muscle, 369 mice for hypothalamus, 491 mice for pancreatic islets, 474 mice for kidney, and 483 mice for liver). The microarray included 37,827 probes, and we focused on the 36,364 probes with known location on the autosomes.

Below is an example of part of the gene expression data. Mouse3051 is one of the 519 F2 mice. In the header row, 497628, 497629, 497630, 497632 and 497637 are some of the 36,364 probes. For each row, the values are the ratio of the mean log10 intensity (mlratio).

|MouseNum|	497628|	497629|	497630|	497632|	497637|
| ------------- | ------------- |------------- | ------------- |------------- | ------------- |
|Mouse3051|	0.5912|	0.02161|	-0.08041|	-0.00478|	0.2248|
|Mouse3551|	0.626|	0.000572|	-0.06688|	-0.2142| -0.07748|	
					
Table 2: An example of part of the expression data.

### 2.3 QTL analysis
We calculated multipoint genotype probabilities at all genetic markers and at a set of pseudomarkers inserted into marker intervals (Table 3). The positions of pseudomarkers were placed at evenly spaced locations between markers, with a maximum spacing of 0.5 cM between adjacent markers or pseudomarkers (Broman et al., 2015). The conditional genotype probabilities were calculated using a hidden Markov model assuming a genotyping error rate of 0.2%, and with genetic distances converted to recombination fractions with the Carter-Falconer map function (CARTER & FALCONER, 1951).
|	|PBB	|PBR	|PRR|
| ------------- | ------------- |------------- | ------------- |
|Mouse3051	|1.0000	|0.0000	|0.0000|
|Mouse3551	|0.0000	|0.0000	|1.0000|
Table 3: An example for genotype probabilities at genetic markers and pseudomarkers. 

We considered each of the six tissues individually and focused on the 36,364 probes with known genomic location on the autosomes, and then identified the nearest marker or pseudomarker to the location of the probe, and using Haley-Knott regression (Haley & Knott, 1992) with sex included as an interactive covariate to estimate the association between genotype at that location and the gene expression of that probe. 

### 2.4 Regression model
Let B and R denote the two alleles in the cross, and let μ ̂_BB, μ ̂_BR and μ ̂_RR denote the average gene expression levels for genotypes BB, BR, and RR, respectively. Additive effect (a) provides a measure of the degree of change in the expression level that occurs with the substitution of R allele for B allele. The estimated additive effect is half the difference in the expression level between the two homozygotes BB and RR: 
```math
a^2+b^2=c^2
```

$a ̂=(μ ̂_RR-μ ̂_BB)/2$
The dominance effect (d) is the deviation of the heterozygote BR from the midpoint of the two homozygotes (BB, RR), estimated as:
d ̂=μ ̂_BR-(μ ̂_RR+μ ̂_BB)/2
For a given QTL model based on Haley-Knott regression, we have 
Y=Xβ+ϵ
where Y is an n×1 vector of expression levels, with n as the number of F2 individuals. 
X is an n×3 matrix of covariates: 
X=■(( 1 ⃑&((P_RR-P_BB))/2&P_BR-((P_RR+P_BB))/2) )
P_RR, P_BBand P_BRare conditional genotype probabilities for genotypes RR, BB and BR given observed multipoint marker genotype data (Table 3).
We did a linear regression that regress expression levels Y on probability matrix X. Thus, the estimated coefficient β ̂_1 is the estimated additive effect a ̂ and β ̂_2 is the estimated dominance effect d ̂. 
The ratio of d/a provides a measure of the degree of dominance. 
d/a={█(0 (pure additivity or no dominance)@1 (complate dominace)@-1   (complete recessive))┤

We plot estimated dominance effect d ̂ vs estimated additive effect a ̂ for all expression traits in six tissues. If all the data points are distributed around the d=0 line, the effect is additive. If all the data points are distributed around the diagonal line d=a, R is dominant. If all the data points are distributed around the diagonal line d=-a, B is dominant.

