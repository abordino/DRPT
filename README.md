# DRPT: Density Ratio Permutation Test

This repository contains the implementation of the **Density Ratio Permutation Test (DRPT)** along with code to reproduce all simulations and experiments presented in the paper [Density Ratio Permutation Tests with connections to distributional shifts and conditional two-sample testing](https://arxiv.org/abs/2505.24529).

Code to reproduce the experiments in the paper [General Frameworks for Conditional Two-Sample Testing](https://arxiv.org/abs/2410.16636) can be found at [this link](https://github.com/suman-cha/Cond2ST?tab=readme-ov-file).

## Repository Structure

The repository is organised into folders corresponding to different experimental settings. Each folder is self-contained and has no dependencies on the others; you can download and run them independently.

- **`diamonds/`**  
  Applies DRPT in a real-world data setting using the diamonds dataset. The 'diamonds' dataset is available in the R-package 'ggplot2'.

- **`stroop/`**  
  Demonstrates DRPT on data from the Stroop effect experiment. The dataset is available at [this link](https://github.com/Lakens/Stroop?tab=readme-ov-file).

- **`frisk/`**  
  Evaluates DRPT performance on the New York frisk dataset. The dataset is available at [this link](https://www.nyc.gov/site/nypd/stats/reports-analysis/stopfrisk.page).

- **`synthetic/`**  
  Tests DRPT on synthetically generated datasets.

---

Feel free to explore each folder based on your use case or interest.

