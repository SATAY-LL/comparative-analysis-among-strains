# comparative-analysis-among-strains

## Install environment

```bash
git clone https://github.com/SATAY-LL/comparative-analysis-among-strains.git 
cd Data-analysis-multiple-strains
conda env create --file environment.yml
conda activate satay-analysis
```

## The pergene files must be tab separated instead of comma separated.

This is a bash code to convert comma separated files to tab separated files, if needed. 

`cat oldfile.txt | tr '[,]' '[\t]' > newfile.txt` 


## To version control the notebooks for the analysis 

1. Convert the modified notebooks to .py files with `jupytext`
    - `jupytext --to py <your-notebook>.ipynb`
2. Add and commit the .py files to the repository
    - `git add <your-notebook>.py`
    - `git commit -m "<commit -message>"`
3. Push the changes to the repository
    
This way of versioning the notebooks are taken from this blog [Version Control With Jupyter Notebook](https://towardsdatascience.com/version-control-with-jupyter-notebook-b9630bc5996e).