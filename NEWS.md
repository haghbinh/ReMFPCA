ReMFPCA 1.1.0
===========

Updates
-------

- ???

New Additions
-------------
 The new following functions are replaced to GitHub load data:
 
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

- ?????

