## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
 
## Resubmission

This is a resubmission. In this version I have:

 * Rewritten the output for R/SPRT.R. The output is now an object from which users may extract the information they are interested in.
 * Removed all cat()/print() that generated not easily suppressed information messages to the console.
 * Updated the corresponding comments, examples and vignettes.
 * Increased version number from 0.1.0 to 0.1.1
