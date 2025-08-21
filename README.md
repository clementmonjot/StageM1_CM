# StageM1_CM
This repository contains script to reproduce figures present in M1 internship report

### 1. General use 

The same script was used to obtained all figures present in the internship report. The script need to be adjusted by arguments before being used to create figures. Numerous arguments are used:

a. input argument to read differents rawdata from repository **rawdata/**.

b. control_sample argument to statistically determine which sample will be compared to others.

c. target_gene argument to analyse specific gene.

d. pattern_i argument to analyse specific samples by selecting parts of their names. This argument can be empty if no specific samples want to be analyse.

Figures will be generated in the repository **result/**.


### 2. Figures parameters

## Figure 5

##### Figure 5A

a. Infection.csv 

b. IAV_Control

c. PARP6

```
Rscript Figure_Generator.R Infection.csv IAV_Control PARP6 
```

##### Figure 5B

a. Infection.csv 

b. IAV_Control

c. NP
```
Rscript Figure_Generator.R Infection.csv IAV_Control NP 
```
## Figure 6

##### Figure 6A 

a. Stimulation_HT-DNA

b. Control

c. PARP6
```
Rscript Figure_Generator.R Stimulation_HT-DNA.csv Control PARP6
```
##### Figure 6B 

a. Stimulation_HT-DNA

b. Control

c. PARP8
```
Rscript Figure_Generator.R Stimulation_HT-DNA.csv Control PARP8
```
##### Figure 6C

a. Stimulation_HT-DNA

b. Control

c. IFNB
```
Rscript Figure_Generator.R Stimulation_HT-DNA.csv Control IFNB
```
##### Figure 6D

a. Stimulation_HT-DNA

b. Control

c. TNFa
```
Rscript Figure_Generator.R Stimulation_HT-DNA.csv Control TNFa
```
## Figure 7


##### Figure 7A

a. Stimulation_3pRNA

b. Control

c. PARP6
```
Rscript Figure_Generator.R Stimulation_3pRNA.csv Control PARP6
```
##### Figure 7B

a. Stimulation_3pRNA

b. Control

c. PARP8
```
Rscript Figure_Generator.R Stimulation_3pRNA.csv Control PARP8
```
##### Figure 7C

a. Stimulation_3pRNA

b. Control

c. IFNB
```
Rscript Figure_Generator.R Stimulation_3pRNA.csv Control IFNB 
```
##### Figure 7D

a. Stimulation_3pRNA

b. Control

c. TNFa
```
Rscript Figure_Generator.R Stimulation_3pRNA.csv Control TNFa
```
##### Figure 7E

a. Stimulation_3pRNA

b. 5h Con

c. PARP6

d. 5h
```
Rscript Figure_Generator.R Time_dependant.csv 5h\ Control IFNB 5h
```
  Only 5h specific samples are analysed here. To show other time condition sample analysis,   replace 5h by wanted time.

## Figure 8

##### Figure 8A

a. Overexpression

b. Control

c. RIG1
```
Rscript Figure_Generator.R Overexpression.csv Negative\ Control RIG1
```
##### Figure 8B

a. Overexpression

b. Control

c. cGAS
```
Rscript Figure_Generator.R Overexpression.csv Negative\ Control cGAS
```
##### Figure 8C

a. Overexpression

b. Control

c. PARP8
```
Rscript Figure_Generator.R Overexpression.csv Negative\ Control PARP8
```

#### Figure 8D

a. Overexpression

b. Control

c. IFNB
```
Rscript Figure_Generator.R Overexpression.csv Negative\ Control IFNB
```

