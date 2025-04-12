ReMFPCA 2.0.0
===========

Updates
-------

- Enhanced documentation for the class `remfpca` by providing more detailed descriptions of parameters.

- Renamed the argument `alpha` to `smooth_tuning` in the remfpca R6 class to better reflect its role in controlling smoothness.


New Additions
-------------
 The new following functions are replaced to GitHub load data:

- Introduced a `joint_power()` function as an alternative method to address smoothness issues in MFPCA, and a `sequential_power()` to tackle both sparsity and smoothness issues.

- Added `electrical_power_data` and `motion_sense_data` as example datasets for multivariate functional data, with variables observed over a one-dimensional domain.

- Added a detailed example for the class `remfpca`, demonstrating the regularized MFPCA approaches.

- Added `plotScores()` function for visualizing functional principal component scores.

- Introduced `reconstructCurve()` function to reconstruct original curves using estimated FPCs.

- Enhanced `plotMFPCA()` with options for better customization and additional plotting styles.

- Added examples to all exported functions for improved usability.

- Extended `estimateMFPCA()` with new argument `scale = TRUE/FALSE` to control scaling before FPCA.

- Improved internal documentation and inline comments for better maintainability.

- Added more extensive unit tests for core functions.

In these updated functions, upon downloading the data files from GitHub into a temporary directory (not the global environment), the target objects are now returned within the function. This modification allows users to save the data into an arbitrary variable of their choice.


Minor improvements and bug fixes
--------------------------------

- Improved error handling for non-strict parameter inputs: the function now infers the intended input format, applies the correction, and issues a warning to inform the user of the changes made.

