# Convert an `"efa_loadings"` table to matrix or a matrix to `"efa_loadings"`

The loadings tables returned by
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
are of class `c("efa_loadings", "LOADINGS")`, which prevents applying
functions on them. This function changes their class to `"matrix"`, and
changes it back to `"efa_loadings"` when done.

## Usage

``` r
.change_class(x, cl = "matrix")
```

## Arguments

- x:

  A table of class `"matrix"` or `"efa_loadings"`.

- cl:

  A character vector with the class to change the table to. Should be
  `c("efa_loadings", "LOADINGS")` or `"matrix"`.

## Value

A table with the loadings, of class either `"efa_loadings"` or
`"matrix"`.

## Author

Andreas Soteriades
