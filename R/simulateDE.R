#' Simulate DE by swapping genes
#'
#' A simple DE simulator that takes a dataset with a grouping factor and induces
#' DE by swapping some genes in one of the groups but not the other. This should
#' ensure keeping the same data structure as the original data without having to
#' estimate distribution parameters.
#'
#' @details
#' Note that for *SummarizedExperiment* objects, any additional assays
#' (containing for example normalized values) are not retained. For simple
#' library size-factor normalization, size factors for each cell will be
#' identical before and after simulation. However, for more sophisticated
#' methods, such as normalization by deconvolution, the size factors might
#' differ due to the gene swapping.
#'
#' @param x A numeric matrix containing features in rows and cells in columns.
#'   Alternatively, a \linkS4class{SummarizedExperiment} or
#'   \linkS4class{SingleCellExperiment} object.
#' @param groups A vector of length equal to `ncol(x)` specifying the group to
#'   which each cell is assigned. DE will be induced between these groups.
#' @param prop_DE Numeric scalar specifying the proportion of genes that will be
#'   simulated as DE. Default: `0.01`.
#' @param ... For the generic, arguments to be passed to specific methods.
#'
#'   For the `SummarizedExperiment` method, further arguments to be passed to
#'   the `ANY` method.
#'
#'   For the `SingleCellExperiment` method, further arguments to be passed to
#'   the `SummarizedExperiment` method.
#'
#' @param use_assay A string or integer scalar specifying the assay of `x` to be
#'   used as input for the simulation. The default is to use `"counts"`.
#'
#' @return
#' A \linkS4class{SummarizedExperiment} object with DE induced between the
#' specified `groups`. The simulated counts are in the `"counts"` assay of the
#' returned object. The `"groups"` used as input for the simulation are recorded
#' in the `"sim_group"` column of the `colData`. The `rowData` contains the
#' following columns:
#'
#' * `"is_DE"`: logical vector indicating the ground truth status of each gene
#' * `"swapped_gene"`: character vector indicating the original gene with which
#' the gene was swapped
#'
#' The expected number of DE genes will be equal to `prop_DE * nrow(x)`. Note
#' however that it's possible that the actual number might be 1 lower than this.
#' This can happen when a single feature remains un-"swapped" without any other
#' features left to swap it with.
#'
#' If `x` was a *SummarizedExperiment* object, the original `colData` and
#' `rowData` are combined with these new columns.
#'
#' If `x` was a *SingleCellExperiment* object, the output will also be a
#' *SingleCellExperiment*.
#'
#' @examples
#' example_sce <- scuttle::mockSCE()
#' sim <- simulateDE(example_sce, groups = example_sce$Treatment, prop_DE = 0.1)
#' sim
#'
#' @author Milan Malfait
#' @name simulateDE
NULL

.simulate_DE <- function(x, groups, prop_DE = 0.01) {
    stopifnot(length(groups) == ncol(x))

    if (prop_DE < 0 || prop_DE > 1) {
        stop("`prop_DE` should be between 0 and 1.", call. = FALSE)
    }

    # TODO: implement method for multiple groups
    if (length(unique(groups)) != 2) {
        stop("Currently only supports simulating DE between 2 groups.",
            call. = FALSE)
    }

    sel_group <- sample(unique(groups), 1)
    sel_group_idx <- which(groups == sel_group)

    n_rows <- nrow(x)
    n_DE <- round(prop_DE * n_rows)
    if (n_DE < 2L) {
        stop("`prop_DE=", prop_DE, "` resulted in less than 2 DE genes.",
            "\nNeed at least 2 DE genes to do the swapping.",
            " Consider increasing `prop_DE`.",
            call. = FALSE)
    }
    de_genes <- sample(n_rows, size = n_DE)

    ## Scramble genes only within the selected group
    scrambling <- .scramble_rows(de_genes)
    sim_cnts <- x
    sim_cnts[de_genes, sel_group_idx] <- x[scrambling$rows, sel_group_idx]

    ## It's possible that a single gene remains unswapped
    row_data <- DataFrame(is_DE = logical(n_rows))
    row_data[de_genes, "is_DE"] <- scrambling$swapped
    row_data[["swapped_gene"]] <- rep(NA_character_, n_rows)
    row_data[de_genes, "swapped_gene"] <- rownames(x)[scrambling$rows]
    rownames(row_data) <- rownames(sim_cnts)

    SummarizedExperiment(
        assays = list(counts = sim_cnts),
        colData = DataFrame(sim_group = groups, row.names = colnames(sim_cnts)),
        rowData = row_data
    )
}

## Helper to scramble rows of a matrix
.scramble_rows <- function(rows) {
    ## Randomly permute rows
    scrambled <- sample(rows)
    non_swapped <- rows == scrambled

    ## Keep scrambling until at most 1 is non-swapped
    while (sum(non_swapped) > 1L) {
        if (sum(non_swapped) == 2L) {
            ## Just swap if there are 2 left
            scrambled[non_swapped] <- rev(scrambled[non_swapped])
        } else {
            scrambled[non_swapped] <- sample(scrambled[non_swapped])
        }
        non_swapped <- rows == scrambled
    }
    list(rows = scrambled, swapped = !non_swapped)
}


#' @export
#' @rdname simulateDE
setGeneric("simulateDE", function(x, groups, ...) standardGeneric("simulateDE"))

#' @export
#' @rdname simulateDE
setMethod("simulateDE", "ANY", .simulate_DE)

#' @export
#' @rdname simulateDE
setMethod("simulateDE", "SummarizedExperiment",
    function(x, groups, ..., use_assay = "counts") {
        out <- .simulate_DE(assay(x, use_assay), groups = groups, ...)

        ## Combine row and col data
        rowData(out) <- cbind(rowData(x), rowData(out))
        colData(out) <- cbind(colData(x), colData(out))
        out
    }
)

#' @export
#' @rdname simulateDE
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("simulateDE", "SingleCellExperiment",
    function(x, groups, ...) {
        out <- callNextMethod(x = x, groups = groups, ...)
        as(out, "SingleCellExperiment")
    }
)
