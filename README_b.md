# README

This is a README file for the code associated with the paper:

**“Discrete-Time Competing Risks Regression with or without Penalization,” Tomer Meir and Malka Gorfine.**

Check out the paper GitHub repository as well:
[https://github.com/tomer1812/DiscreteTimeSurvivalPenalization](https://github.com/tomer1812/DiscreteTimeSurvivalPenalization)

Check out PyDTS GitHub repository which implements the proposed approach:
[https://github.com/tomer1812/pydts](https://github.com/tomer1812/pydts)


## Instructions

### Environment Setup

1. Set up an environment with Python ≥ 3.9.17 and R ≥ 4.3.1.
   
   Example using Anaconda:
   ```bash
   conda create --name dtsp-env python=3.9.17 r-base=4.3.1
   ```

2. Activate your new Anaconda environment:
   ```bash
   conda activate dtsp-env
   ```

3. Install package dependencies:
   ```bash
   pip install -r requirements.txt
   conda install r-tidyverse r-survival r-FSA r-foreach r-doParallel r-dplyr
   ```

4. Ensure the dependencies are installed in the correct Anaconda environment.


### Running a `.ipynb` File

1. Add your Anaconda environment to the list of interpreters:
   ```bash
   python -m ipykernel install --user --name=dtsp-env --display-name "dtsp-env"
   ```

2. Open Jupyter Notebook.

3. Ensure parameters are set as desired, particularly the output directory or output filename.

4. Set the IPython kernel interpreter to dtsp-env.

5. Run the cells in the notebook.

### Running a `.py` File

1. Open the file.

2. Ensure parameters are set as desired, particularly the output directory or output filename.

3. Save and exit.

4. Run the Python file using the command line:
   ```bash
   python FILENAME.py
   ```

### Running an `.R` File

1. Open the file.

2. Ensure parameters are set as desired, particularly the output directory or output filename.

3. Save and exit.

4. Run the R file using the command line:
   ```bash
   Rscript FILENAME.R
   ```

## List of Files

- **`Implementation-example.R`**: R implementation example of the proposed approach.
- **`Simulations - different n.ipynb`**: Simulations for large sample sizes with two competing events (see Settings 3, 4, 5, 6 from Table S2).
- **`Simulations - different n - 3 competing events.ipynb`**: Simulations for large sample sizes with three competing events (see Settings 7, 8, 9, 10 from Table S2).
- **`Simulations - different d - timing.ipynb`**: Timing comparison analysis (see Web Appendix H).
- **`Simulations - Regularization.ipynb`**: Large sample size LASSO simulations with independent covariates (see Setting 11 from Table S2).
- **`Simulations - Regularization - corr.ipynb`**: Large sample size LASSO simulations with correlated covariates (see Settings 12 from Table S2).
- **`Simulations - SIS-SIS-L.ipynb`**: Sure Independence Screening simulations (see Settings 17, 18, 19 from Table S2).
- **`Small-sample-size.R`**: R example comparing methods for small sample sizes with the ‘exact’ method (see Settings 1, 2 from Table S2).
- **`small_sample_results_reader.py`**: results analysis of the R example comparing methods for small sample sizes with the ‘exact’ method (see Settings 1, 2 from Table S2).
- **`regularization-small-sample-size.py`**: Regularization analysis for small sample sizes with ‘exact’ in Python using PyDTS (see Settings 14, 15).
- **`TP-FP-small-sample-size.py`**: LASSO grid search true positive-false positive analysis for small sample sizes (see Setting 16).
- **`Simulations - TP-FP small sample size results reader.ipynb`**: Results reader for analysing LASSO grid search true positive-false positive analysis for large sample sizes (see Setting 16).
- **`Simulations - Regularization-TP-FP.ipynb`**: LASSO grid search true positive-false positive analysis for large sample sizes (see Setting 13).
- **`Simulations - Regularization-TP-FP results reader.ipynb`**: Results reader for analysing LASSO grid search true positive-false positive analysis for large sample sizes (see Setting 13).
- **`continuous_bias.R`**: Comparison of continuous and discrete analyses (see Web Appendix A).
- **`continuous_bias_large_d.R`**: Comparison of continuous and discrete analyses with a large number of time points (see Web Appendix A).
- **`mimiciv.ipynb`**: MIMIC-IV length-of-stay use-case analysis. Requires specifying the MIMIC-IV 2.0 data directory, output directory, and the `constants.py` file.
- **`constants.py`**: Constants definition for MIMIC-IV analysis.
- **`requirements.R`**: R packages requirements.
- **`requirements.txt`**: Python packages requirements.
- 
