# Trimming Raw Reads Worksheet

<!--- Write name below --->
## Name: 

<!--- For this worksheet, answer the following questions--->

## Q1: What does "cleaning" your reads mean?
Answer:
Filtering the sequence data by removing/trimming low quality reads, adapters, and by other 

### Q2: Open the script called "trim_raw_reads.sh". For every line that says ```insert description here```, replace that text with a description of what the line will do

### Q4: Attach an image of the plot from fastqc showing the average base quality from your reads files
![before-trimming-PHRED-plot-1](./images/beforetrimming1.png)
![before-trimming-PHRED-plot-1](./images/beforetrimming2.png)
---

## The following questions pertain to your first fastp run (without altering the parameters)
### Q5: How many reads were there in the R1 file before filtering?
Answer:
33967298

### Q6: How many bases were there in the R1 file before filtering?
Answer:
101 bases per sequence

### Q7: What proportion of bases were above PHRED score of 20 before filtering?
Answer: 
95% (52/55)

### Q8: What proportion of bases were above PHRED score of 30 before filtering?
Answer: 
72% (40/55)

### Q13: How many reads were there in the R2 file before filtering?
Answer:
33967298

### Q14: How many bases were there in the R2 file before filtering?
Answer:
101 bases per sequence

### Q15: What proportion of bases were above PHRED score of 20 before filtering?
Answer: 
93% (51/55)

### Q16: What proportion of bases were above PHRED score of 30 before filtering?
Answer: 
53% (29/55)


### Q17: How many reads were there in the merged file after filtering?
Answer:
12295196

### Q18: How many bases were there in the merged file after filtering?
Answer:
3500000+

### Q19: What proportion of bases in the merged file were above PHRED score of 20 after filtering?
Answer: 
100%

### Q20: What proportion of bases in the merged file were above PHRED score of 30 after filtering?
Answer: 
90% (38/42)


### Q21: What is the difference between the merged and unmerged files (in principle, not quantitatively)?
Answer:
The merged file is not only a combination of both unmerged files, but is also filtered according to the conditions set by running fastp.

### Q22: Why are the unmerged files for R1 and R2 different lengths?
Answer:


---

## The following questions pertain to Remix 1 (the first time you change fastp parameters)
### Q23: What parameters did you change?
Answer: 
Changed -q from 15 to 30

### Q24: How did you expect this to change the filtering results (be specific)?
Answer: 
Limit the number of sequences left over after cleaning by removing any position with a PHRED score under 30

### Q25: Explain the results. Did the change cause an effect that matched your expectations? Use information from the fastp output to explain.
Answer: 
The number of reads was about cut in half, with every position in the read having a PHRED score of 30 or above (everything located in the green section of the Per base sequence quality map), matching expectations. 
---

## The following questions pertain to Remix 2 (the first time you change fastp parameters)
### Q26: What parameters did you change?
Answer: 
Changed -l from 15 to 20

### Q27: How did you expect this to change the filtering results (be specific)?
Answer: 
Remove any read with a length less than 20 base pairs long

### Q28: Explain the results. Did the change cause an effect that matched your expectations? Use information from the fastp output to explain.
Answer: 
The number of reads kept was again cut in half, and the smallest read was 20 base pairs long, meeting that expectation. Little else changed in terms of the data provided.
