# Refactoring Plan: comparative-analysis-among-strains

This document outlines the steps to refactor this codebase into a reusable Python package with testing and modular design.

---

## ✅ Step 1: Assess and Structure the Codebase

1. **Create a new development branch**
   ```bash
   git checkout -b refactor-package
   ```

2. **Organize files**
   Identify:
   - Notebooks that produce results
   - Code blocks that are reused
   - Scripts that load/process/analyze/plot data

3. **Create this folder structure**
   ```
   comparative_analysis/
   ├── comparative_analysis/         # Python package
   │   ├── __init__.py
   │   ├── io.py                     # File reading/writing
   │   ├── analysis.py               # Data processing, stats
   │   ├── plotting.py               # Visualization functions
   ├── notebooks/                    # Research notebooks
   ├── tests/                        # Unit tests
   │   └── test_analysis.py
   ├── setup.py                      # To make it pip-installable
   ├── requirements.txt              # Dependencies
   └── README.md
   ```

---

## ✅ Step 2: Refactor into Functions and Modules

1. **Move repeated logic into functions**
   Example:
   ```python
   def compute_fold_change(control, treated):
       return [c / t for c, t in zip(control, treated)]
   ```

2. **Move functions into modules**
   - File reading → `io.py`
   - Analysis logic → `analysis.py`
   - Plotting → `plotting.py`

3. **Use these functions in your notebooks**
   ```python
   from comparative_analysis.analysis import compute_fold_change
   ```

4. **Avoid hardcoded paths or magic numbers**  
   Use parameters or a config file (`config.py` or `config.yaml`)

---

## ✅ Step 3: Add Unit Tests

1. **Install `pytest`**
   ```bash
   pip install pytest
   ```

2. **Create test files in `tests/` folder**
   Example: `tests/test_analysis.py`
   ```python
   from comparative_analysis.analysis import compute_fold_change

   def test_compute_fold_change():
       result = compute_fold_change([1, 2], [2, 4])
       assert result == [0.5, 0.5]
   ```

3. **Run tests**
   ```bash
   pytest
   ```

---

## ✅ Step 4: Make It a Python Package

1. **Create `setup.py`**
   ```python
   from setuptools import setup, find_packages

   setup(
       name='comparative_analysis',
       version='0.1.0',
       packages=find_packages(),
       install_requires=['numpy', 'pandas', 'matplotlib'],
   )
   ```

2. **Install the package locally**
   ```bash
   pip install -e .
   ```

---

## ✅ Step 5: Document the Codebase

1. **Write docstrings using NumPy or Google style**
2. **Add usage examples in `README.md`**
3. **Optionally create docs using Sphinx or MkDocs later**

---

## ✅ Step 6: Make the Project FAIR and Reusable

1. **Add a license**
   Example: `LICENSE` file with MIT License

2. **Add a `CITATION.cff`**
   Example:
   ```yaml
   cff-version: 1.2.0
   message: "If you use this software, please cite it as below."
   authors:
     - family-names: YourLastName
       given-names: YourFirstName
   title: "Comparative Analysis Among Strains"
   version: 0.1.0
   doi: "10.xxxx/zenodo.1234567"
   date-released: YYYY-MM-DD
   ```

3. **Tag your releases on GitHub**
   ```bash
   git tag -a v0.1.0 -m "Initial refactored release"
   git push origin v0.1.0
   ```

---

## 🛠 Optional: Add GitHub Actions for Continuous Integration

1. Create `.github/workflows/python-app.yml`:
   ```yaml
   name: Python application

   on: [push]

   jobs:
     build:

       runs-on: ubuntu-latest

       steps:
         - uses: actions/checkout@v2
         - name: Set up Python
           uses: actions/setup-python@v2
           with:
             python-version: '3.9'
         - name: Install dependencies
           run: |
             python -m pip install --upgrade pip
             pip install -e .[dev]
         - name: Run tests
           run: |
             pytest
   ```

---

## You're Ready!

Refactor one step at a time. Commit often and track progress in issues or a TODO list.
