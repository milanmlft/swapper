# Simulate DE by swapping features

A simple DE simulator that takes a dataset with a grouping factor and
induces DE by swapping some features in one of the groups but not the
other. This should ensure keeping the same data structure as the
original data without having to estimate distribution parameters.

## Usage

``` r
simulateDE(x, which_cols, ...)

# S4 method for class 'ANY'
simulateDE(x, which_cols, prop_DE = 0.01)

# S4 method for class 'SummarizedExperiment'
simulateDE(x, which_cols, ..., use_assay = "counts")

# S4 method for class 'SingleCellExperiment'
simulateDE(x, which_cols, ...)
```

## Arguments

- x:

  A numeric matrix containing features in rows and cells in columns.
  Alternatively, a SummarizedExperiment or SingleCellExperiment object.

- which_cols:

  Character, numeric or logical vector specifying for which columns
  (cells, samples, ...) of `x` the feature swapping should occur.
  Usually, this will be a random subset of columns belonging to a mock
  group.

- ...:

  For the generic, arguments to be passed to specific methods.

  For the `SummarizedExperiment` method, further arguments to be passed
  to the `ANY` method.

  For the `SingleCellExperiment` method, further arguments to be passed
  to the `SummarizedExperiment` method.

- prop_DE:

  Numeric scalar specifying the proportion of features that will be
  simulated as DE. Default: `0.01`.

- use_assay:

  A string or integer scalar specifying the assay of `x` to be used as
  input for the simulation. The default is to use `"counts"`.

## Value

A SummarizedExperiment object with DE induced between the specified
`groups`. The simulated counts are in the `"counts"` assay of the
returned object. The `"sim_group"` column of the `colData` indicates
whether swapping was performed in that column (specified by the
`which_cols` argument). The `rowData` contains the following columns:

- `"is_DE"`: logical vector indicating the ground truth status of each
  feature

- `"swapped_feature"`: character vector indicating the original feature
  with which the feature was swapped

If `x` was a *SummarizedExperiment* object, the original `colData` and
`rowData` are combined with these new columns.

If `x` was a *SingleCellExperiment* object, the output will also be a
*SingleCellExperiment*.

## Details

Note that for *SummarizedExperiment* objects, any additional assays
(containing for example normalized values) are not retained. For simple
library size-factor normalization, size factors for each cell will be
identical before and after simulation. However, for more sophisticated
methods, such as normalization by deconvolution, the size factors might
differ due to the feature swapping. The expected number of DE features
will be equal to `prop_DE * nrow(x)`. Note however that it's possible
that the actual number might be 1 lower than this. This can happen when
a single feature remains un-"swapped" without any other features left to
swap it with.

## Author

Milan Malfait

## Examples

``` r
example_sce <- scuttle::mockSCE()

## Swap features in the "treat1" group
cells_to_swap <- example_sce$Treatment == "treat1"

sim <- simulateDE(example_sce, which_cols = cells_to_swap, prop_DE = 0.1)
sim
#> class: SingleCellExperiment 
#> dim: 2000 200 
#> metadata(0):
#> assays(1): counts
#> rownames(2000): Gene_0001 Gene_0002 ... Gene_1999 Gene_2000
#> rowData names(2): is_DE swapped_feature
#> colnames(200): Cell_001 Cell_002 ... Cell_199 Cell_200
#> colData names(4): Mutation_Status Cell_Cycle Treatment sim_group
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```
