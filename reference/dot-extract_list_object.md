# Extract a list object by its name

Consider a list of named sub-lists. This function extracts, for each
sub-list, the sub-list element that is specified by the user. This
function is useful for extracting results from
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
for each imputation in
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md).

## Usage

``` r
.extract_list_object(alist, object)
```

## Arguments

- alist:

  A list of sub-lists, typically a list of \\m\\ objects of class
  `"efa"`, where \\m\\ is the number of imputations passed to
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md).

- object:

  String of length 1. The name of the object to extract e.g. `"h2"` or
  `"vars_accounted"`.

## Value

A list of length \\m\\, with each element containing the extracted
`object` for the \\k\\th element (\\k = 1,..., m\\).

## Author

Andreas Soteriades
