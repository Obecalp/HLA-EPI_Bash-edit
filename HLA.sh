	#!/bin/bash

# Cherche les allèles correspondant à l'identifiant $i dans le fichier $2, en distinguant les receveur (R) des donneurs (D) et les alleles en fonction de leur position dans le tableau ((A=4, DQ=13), refomante les alleles **:** en ****, puis le stock dans la variable correspondante. X1 et X2 proviennent du donneur, X3 et X4 du receveur.
function attribution(){
	grep "-wE" "$1" "$2" | grep ";$3;" | cut '-d' ';' '-f' "$4" | cut '-d' ':' '-f' '1,2' | sed "s/://g"
	}

#Définition de la variable précision ci dessous, détermine l'argument donné:
function exclusion(){
if [ $(echo "$1" | grep "$2" | wc "-l") -gt 0 ]
then
echo "$3"
fi
}


#établi le mismatch des epitope correspondant aux allèles donnés en $1 et $3 dans la banque de donnée en $3, en tenant compte de la précision demandée en $4	
function score(){

	grep "-wE" "$1" "$2" | grep "-vwE" "HelloWorld$4" | cut '-d' ';' '-f' '3' | grep "-vw" "$(grep "-wE" "$3" "$2" | cut "-d" ';' '-f' '3')" | uniq

}

#définie avec exclusion()
precision=$(exclusion "$3" "E" "|f;[f-t]")$(exclusion "$3" "C" "|t;[f-t]")$(exclusion "$3" "V" "|[f-t];f")$(exclusion "$3" "P" "|[f-t];t")

#header
echo "groupe;HLA A;HLA B;HLA C;CLASSE I;HLA DR;HLA DQ;HLA DQA1; HLA DP;CLASSE II;TOTAL;mm_A;mm_B;mm_C;mm_DR;mm_DQ;mm_DQA1;mm_DR;mm_allI;mm_allII;mm_epiA;mm_epiB;mm_epiC;mm_epiDR;mm_epiDQ;mm_epiDQA1;mm_epiDP;A1_1;A2_1;B1_1;B2_1;C1_1;C2_1;DR1_1;DR2_1;DQ1_1;DQ2_1;DQA11_1;DQA2_1;DP1_1;DP2_1;A1_2;A2_2;B1_2;B2_2;C1_2;C2_2;DR1_2;DR2_2;DQ1_2;DQ2_2;DQA11_2;DQA12_2;DP1_2:DP2_2" > "enrocil.csv"


for i in $(grep ";" "$1" | grep "-v" "ID_group" | cut '-d' ';' '-f' '1'| uniq)
do

A1=$(attribution "$i" "$1" "D" "4")
A2=$(attribution "$i" "$1" "D" "5")
B1=$(attribution "$i" "$1" "D" "6")
B2=$(attribution "$i" "$1" "D" "7")
C1=$(attribution "$i" "$1" "D" "8")
C2=$(attribution "$i" "$1" "D" "9")
DR1=$(attribution "$i" "$1" "D" "10")
DR2=$(attribution "$i" "$1" "D" "11")
DQ1=$(attribution "$i" "$1" "D" "12")
DQ2=$(attribution "$i" "$1" "D" "13")
DQA11=$(attribution "$i" "$1" "D" "14")
DQA12=$(attribution "$i" "$1" "D" "15")
DP1=$(attribution "$i" "$1" "D" "16")
DP2=$(attribution "$i" "$1" "D" "17")




A3=$(attribution "$i" "$1" "R" "4")
A4=$(attribution "$i" "$1" "R" "5")
B3=$(attribution "$i" "$1" "R" "6")
B4=$(attribution "$i" "$1" "R" "7")
C3=$(attribution "$i" "$1" "R" "8")
C4=$(attribution "$i" "$1" "R" "9")
DR3=$(attribution "$i" "$1" "R" "10")
DR4=$(attribution "$i" "$1" "R" "11")
DQ3=$(attribution "$i" "$1" "R" "12")
DQ4=$(attribution "$i" "$1" "R" "13")
DQA13=$(attribution "$i" "$1" "R" "14")
DQA14=$(attribution "$i" "$1" "R" "15")
DP3=$(attribution "$i" "$1" "R" "16")
DP4=$(attribution "$i" "$1" "R" "17")


#découpement du calcul du mismatch alleles (alleleA, etc...
#AD=$(attribution "$i" "$1" "D" "4,5" | sed "s/;/\n/g") => alleles en D
#AR=$(attribution "$i" "$1" "R" "4,5" | sed "s/;/\n/g") => alleles en R
#AR=$(echo "$AD" | grep "$AR" | wc "-l") => calcul score en fonction de la difference
#allele=$((($AR-2)*-1)) #conversion des valeurs

alleleA=$((($(echo "$(attribution "$i" "$1" "D" "4,5" | sed "s/;/\n/g")" | grep "$(attribution "$i" "$1" "R" "4,5" | sed "s/;/\n/g")" | wc "-l")-2)*-1)) #conversion des valeurs


alleleB=$((($(echo "$(attribution "$i" "$1" "D" "6,7" | sed "s/;/\n/g")" | grep "$(attribution "$i" "$1" "R" "6,7" | sed "s/;/\n/g")" | wc "-l")-2)*-1)) #conversion des valeurs


alleleC=$((($(echo "$(attribution "$i" "$1" "D" "8,9" | sed "s/;/\n/g")" | grep "$(attribution "$i" "$1" "R" "8,9" | sed "s/;/\n/g")" | wc "-l")-2)*-1)) #conversion des valeurs

alleleDR=$((($(echo "$(attribution "$i" "$1" "D" "10,11" | sed "s/;/\n/g")" | grep "$(attribution "$i" "$1" "R" "10,11" | sed "s/;/\n/g")" | wc "-l")-2)*-1)) #conversion des valeurs

alleleDQ=$((($(echo "$(attribution "$i" "$1" "D" "12,13" | sed "s/;/\n/g")" | grep "$(attribution "$i" "$1" "R" "12,13" | sed "s/;/\n/g")" | wc "-l")-2)*-1)) #conversion des valeurs

alleleDQA1=$((($(echo "$(attribution "$i" "$1" "D" "14,15" | sed "s/;/\n/g")" | grep "$(attribution "$i" "$1" "R" "14,15" | sed "s/;/\n/g")" | wc "-l")-2)*-1)) #conversion des valeurs

alleleDP=$((($(echo "$(attribution "$i" "$1" "D" "16,17" | sed "s/;/\n/g")" | grep "$(attribution "$i" "$1" "R" "16,17" | sed "s/;/\n/g")" | wc "-l")-2)*-1)) #conversion des valeurs

alleleI=$(($alleleA+$alleleB+$alleleC))
alleleII=$((alleleDR+$alleleDQ+alleleDQA1+alleleDP))






#pour un même identifiant de groupe:

#BaseX: calcul des mismatch d'epitopes contenus dans la banque de donnée de l'allele X 
#scoreX: contient le nombre d'épitope
#epitopeX: nécessaire pour éviter les sauts de case


#Stockage des mismatchs d'épitopes dans les variables baseX
BaseA=$(score "A;$A1;*|A;$A2;*" "$2" "A;$A3;*|A;$A4;*" "$precision")
BaseB=$(score "B;$B1;*|B;$B2;*" "$2" "B;$B3;*|B;$B4;*" "$precision")
BaseC=$(score "C;$C1;*|C;$C2;*" "$2" "C;$C3;*|C;$C4;*" "$precision")
BaseI=$(score "A;$A1;*|A;$A2;*|B;$B1;*|B;$B2;*|C;$C1;*|C;$C2;*" "$2" "A;$A3;*|A;$A4;*|B;$B3;*|B;$B4;*|C;$C3;*|C;$C4;*" "$precision")
BaseDR=$(score "DR;$DR1;*|DR;$DR2;*" "$2" "DR;$DR3;*|DR;$DR4;*" "$precision")
BaseDQ=$(score "DQ;$DQ1;*|DQ;$DQ2;*" "$2" "DQ;$DQ3;*|DQ;$DQ4;*" "$precision")
BaseDQA1=$(score "DQA1;$DQA11;*|DQA1;$DQA12;*" "$2" "DR;$DA13;*|DR;$DA14;*" "$precision")
BaseDP=$(score "DP;$DP1;*|DP;$DP2;*" "$2" "DP;$DP3;*|DP;$DP4;*" "$precision")
BaseII=$(score "DR;$DR1;*|DR;$DR2;*|DQ;$DQ1;*|DQ;$DQ2;*|DQA1;$DQA11;*|DQA1;$DQA12;*|DP;$DP1;*|DP;$DP2;*" "$2" "DR;$DR3;*|DR;$DR4;*|DQ;$DQ3;*|DQ;$DQ4;*|DQA1;$DQA13;*|DQA1;$DQA14;*|DP;$DP3;*|DP;$DP3;*" "$precision")


#Stockage des décomptes dans les variables scoreX
scoreA=$(echo $BaseA | tr ' ' "\n" | grep "-E" "[a-zA-Z0-9]+" | wc "-l")
scoreB=$(echo $BaseB | tr ' ' "\n" | grep "-E" "[a-zA-Z0-9]+" | wc "-l")
scoreC=$(echo $BaseC | tr ' ' "\n" | grep "-E" "[a-zA-Z0-9]+" | wc "-l")
scoreI=$(echo $BaseI | tr ' ' "\n" | grep "-E" "[a-zA-Z0-9]+" | wc "-l")

scoreDR=$(echo $BaseDR | tr ' ' "\n" | grep "-E" "[a-zA-Z0-9]+" | wc "-l")
scoreDQ=$(echo $BaseDQ | tr ' ' "\n" | grep "-E" "[a-zA-Z0-9]+" | wc "-l")
scoreDQA1=$(echo $BaseDQA1 | tr ' ' "\n" | grep "-E" "[a-zA-Z0-9]+" | wc "-l")
scoreDP=$(echo $BaseDP | tr '? ' "\n" | grep "-E" "[a-zA-Z0-9]+" | wc "-l")


scoreII=$(echo $BaseII | tr ' ' "\n" | grep "-E" "[a-zA-Z0-9]+" | wc "-l")
Total=$(($scoreI+$scoreII))

#Listage des epitopes dans les variables epitopeX
epitopeA=$(echo $BaseA)
epitopeB=$(echo $BaseB)
epitopeC=$(echo $BaseC)
epitopeDR=$(echo $BaseDR)
epitopeDQB=$(echo $BaseDQ)
epitopeDQA1=$(echo $BaseDQA1)
epitopeDP=$(echo $BaseDP)


#Rend les résultats sous forme de csv pour chaque id
echo "$i;$scoreA;$scoreB;$scoreC;$scoreI;$scoreDR;$scoreDQ;$scoreDQA1;$scoreDP;$scoreII;$Total;$alleleA;$alleleB;$alleleC;$alleleDR;$alleleDQ;$alleleDQA1;$alleleDP;$alleleI;$alleleII;$epitopeA;$epitopeB;$epitopeC;$epitopeDR;$epitopeDQB;$epitopeDQA1;$epitopeDQP;$A1;$A2;$B1;$B2;$C1;$C2;$DR1;$DR2;$DQ1;$DQ2;$DQA11;$DQA12;$DP1;$DP2;$A3;$A4;$B3;$B4;$C3;$C4;$DR3;$DR4;$DQ3;$DQ4;$DQA13;$DQA14;$DP3;$DP4" >> "enrocil.csv"
done



