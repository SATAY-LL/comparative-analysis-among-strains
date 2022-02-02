# comparative-analysis-among-strains

## Install environment

```bash
git clone https://github.com/SATAY-LL/comparative-analysis-among-strains.git 
cd Data-analysis-multiple-strains
conda env create --file environment.yml
conda activate satay-analysis
```

## Convert all pergene files to tab separated instead of comma separated

`cat oldfile.txt | tr '[,]' '[\t]' > newfile.txt` 

