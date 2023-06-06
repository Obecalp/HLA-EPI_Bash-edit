# HLA-EPI_Bash-edit
Re-edition of HLA-EPI in Bash

#Introduction \n

Epitopes are parts of the proteins corresponding to the allele that are particulary immunogen, meaning many an allele can contain many epitopes. \n
The polymorph part of an epitope is called a Eplet, which are 1-5 amino-acid.

It has been observed that the more there were epitopes' mismatch between two patients, the more a DSA production would be susceptible to happend after a transplantation. (See more at https://pubmed.ncbi.nlm.nih.gov/24769079/ or from René J. Duquesnoy)


HLA-EPI aims to establish a "matching" score between a patient and a donor for a transplantation, based on the epitopes carried by their alleles. It is currently available on easy-HLA (https://hla.univ-nantes.fr/index.php#). Its mains advantages are the fact it is pretty easy to use for any user, and the presence of multiples reasearch, which permit to establish mismatching simultaneously many couples of patients.
However, HLA-EPI is currently coded in php, and the it could be optimized in other languages.


The re-edition in bash main to replace the multiple reasearch function wich may be a little bit too slow for too many patients, and sometimes instable on php. It is strongly faster, stable, and can carry as much patients as desired. \n

#Utilisation \n

Currently, HLA.sh respect every aspect of HLA-EPI, this mean it is only possible to evaluate A,B,C,DQ1,DRB1,DQA1,DPB1 alleles. However, the "extended option" is automatic, meaning DRB1 and DQA1 will be calculated (and will return 0 if nothing is given). Imputation hasn't been implemented yet, meaning every alleles fields has to be completed (or the comparaison for this field for this patient and the total count for this patient will be corrupted, there won't have any other consequence). \n

#Input: \n
The aim is to make it easy to use it too. For this are required two files: \n
- A csv file containing all patients haplotype data (as for HLA-EPI initial tool) as : \n group_ID;patient_id;status(R/D);A1;A2;B1;B2;C1;C2;DR1;DR2;DQ1;DQ2;DQA1;DQA2;DP1;DP2 as: \n
group1;patient1;D;0101;0102;5020;0109;0101;0102;0101;5020;0109;0102;5020;0109;5020;0109; \n
group1;patient2;R;0203;0105;5020;0102;5020;0109;0203;0105;0203;0105;0203;0105;0203;0105; \n
- A csv file containing all epitopes-allele relations data as: \n
HLA;allele;eplet;exposition;verification as C;0101;YC,t,f \n



There is no need to have any header, but if there is one in patients file, it has to contain exactly ID_group, or a line full of "0" will appear (without any other consequence actually) \n
However, it is required to respect every field order for both file.

#Output \n

In output is given a csv file containing: \n 
-  The total count of mismatch by allele, then class, then both class together. \n
- The listing of all eplet by allele
- the patient allele (as reminder)
Some eplet are shared by different HLA allele (B and A, etc...) they are compensated, meaning if there is an eplet isn't present in an HLA allele but in an another one it won't be considerated as a mismatch (Meaning it may have some total class I < A + B + C \n

#Note
The epitope data file can be obtained with the script multiplet.sh, which will be soon distributed here

