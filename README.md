# Discrete-Time Competing-Risks Regression with or without Penalization

A Python code for the paper: 

[Meir, Tomer and Gorfine, Malka, "Discrete-time Competing-Risks Regression with or without Penalization", 2023.](https://arxiv.org/abs/2303.01186)

The simulations and MIMIC-IV analysis can be replicated by cloning the repository, running the code in the notebooks directory using the docker image built from the Dockerfile of this project. 

An example for an R implementation of the proposed approach is available in src/Implementation-Example.R

The MIMIC-IV (2.0) dataset is accessible at [PhysioNet](https://physionet.org/content/mimiciv/2.0/) and subjected to PhysioNet credentials.

This work is based on PyDTS Python Package:

[Documentation](https://tomer1812.github.io/pydts/)  

[Github](https://github.com/tomer1812/pydts)


## Citations
If you found this work or PyDTS software useful to your research, please cite the papers:

```bibtex
@article{Meir_Gorfine_DTSP_2023,
    author = {Meir, Tomer and Gorfine, Malka},
    doi = {10.48550/arXiv.2303.01186},
    title = {{Discrete-time Competing-Risks Regression with or without Penalization}},
    url = {https://arxiv.org/abs/2303.01186},
    year = {2023}
}

@article{Meir_PyDTS_2022,
    author = {Meir, Tomer and Gutman, Rom, and Gorfine, Malka},
    doi = {10.48550/arXiv.2204.05731},
    title = {{PyDTS: A Python Package for Discrete Time Survival Analysis with Competing Risks}},
    url = {https://arxiv.org/abs/2204.05731},
    year = {2022}
}
```

