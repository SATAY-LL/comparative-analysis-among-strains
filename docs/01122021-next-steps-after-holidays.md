# Next steps on the analysis of the data

## To compute essentiality

### To ensure that region with few or zero hits are due to biological reasons: 

- [ ] Compute the number of transposons in the upstream and downstream 3KB regions of every gene. Then if the number of transposons is less than the average transposon density of the library * 3KB , then discard the gene for analysis of essentiality. If not then keep the gene for analysis of essentiality.

- [ ] Find a way to get rid of genes that have a low alignment score (gene duplications and gene with repeats), thereby avoiding the situation of having less transposons due to a misalignment. (Look into BWA alignment score)

### Make a list of essential genes from different backgrounds.

- [ ] Have a list of essential genes per genetic background with the adhoc criteria of having few transposons in that region , after filtering out genes with low alignment score and by the neghboring intergenic regions. 

- [ ] Make a heatmap where the rows are all essential genes from all backgrounds and the columns the evolutionary backgrounds. (Like Evelyn plot, but with genes instead of organisms) https://stackoverflow.com/questions/27854243/is-it-possible-to-plot-a-checkerboard-type-plot-in-python

## Fitness calculation

- [ ] Compute the fitness of each gene in the genome from two contributions:
    - contributions of domains
    - contributions of non domains regions
    - excluding the 10% and the last 90% of the gene. (gene Truncations)

- [ ] Add the calculation of the maximum resolution we can achieve by comparing two identical replicates fitnesses difference, fit it to the a function and compute the Full Width at Half maximum (FWHM) of the function. (Floor code)

- Thereby if the difference in fitness between two genes from 2 different replicates is 30% and  the maximum resolution we can get is also 30% , we cant say that the two genes are different. So the fitnes differences between two genes from different libraries have to be greater than the fitness differences between two genes from the same library.