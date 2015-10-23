# snps2assays

This is an R-package to help prepare Fluidigm SNPtype assays (or similar assays) from 
reduced-representation library loci (for example ddRAD or RAD loci).  You pass in a
a data frame of variable sites, a data frame of desired SNPs, and a data frame of the 
consensus sequences of each ddRAD or RAD locus, and it returns the information needed
to order the assays.

## Installing

You can install this directly from GitHub using Hadley Wickham's `devtools` package:
```r
devtools::install_github("eriqande/snps2assays", build_vignettes = TRUE)
```

Once you have it installed, try the usual suspects to learn about the package.  Some useful
entry points will be:
```r
vignette("snps2assays")

package?snps2assays

help(package = "snps2assays")

?assayize
```

## Terms 

As a work of the United States Government, this package is in the
public domain within the United States. Additionally, we waive
copyright and related rights in the work worldwide through the CC0 1.0
Universal public domain dedication.

See TERMS.md for more information.

