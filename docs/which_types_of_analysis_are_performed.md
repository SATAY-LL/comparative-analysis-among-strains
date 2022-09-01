# Analyses of datasets from SATAY sequencing 

## Computation of the fitness from SATAY using a coarse grained model.
- src/analysis_evolutionary_trajectory.ipynb
  - *getting_r function* to compute  fitness values using a coarse model , where it is applied the intergenic model , and every gene gets a value given by the total number of reads over the total number of insertions in the gene.
    - [x] add test
    
  - normalize the fitness to the wild type values using the formula:

    - normalized_fitness=fitness[gene][background]/fitness[gene][wild type]

 - Visualizations of how fitness values are distributed in the wild type normalized to the HO gene.
 - Visualizations of how fitness values are distributed in essential and non essential genes.
 - Heat map to visualize how fitness values for some polarity genes( normalized to wt merged genotype) change across backgrounds.
 - ROC curve to assess the performance of the low fitness values with essentiality.