library(scuttle)
ref <- mockSCE()

test_that("simulateDE works", {
    set.seed(123)
    sim <- simulateDE(ref, groups = ref$Treatment, prop_DE = 0.1)
    expect_s4_class(sim, "SummarizedExperiment")
    expect_equal(dim(sim), dim(ref))
    expect_true(all(c("is_DE", "swapped_gene") %in% colnames(rowData(sim))))
    expect_true("sim_group" %in% colnames(colData(sim)))
})

test_that("simulateDE returns correct number of DE genes", {
    set.seed(42)
    prop <- 0.5
    sim <- simulateDE(ref, groups = ref$Treatment, prop_DE = prop)
    is_de <- rowData(sim)$is_DE
    expected <- round(prop * nrow(ref))
    ## Actual number of DE genes can be 1 lower than 'expected'
    expect_true(any(sum(is_de) == c(expected - 1, expected)))
})

test_that("simulateDE fails for invalid `prop_DE`", {
    expect_error(simulateDE(ref, groups = ref$Treatment, prop_DE = 2),
        "`prop_DE` should be between 0 and 1.")
    expect_error(simulateDE(ref, groups = ref$Treatment, prop_DE = -0.5),
        "`prop_DE` should be between 0 and 1.")
    expect_error(simulateDE(ref, groups = ref$Treatment, prop_DE = 0),
        "resulted in less than 2 DE genes.")
})

test_that("simulateDE behaves consistently with random seed", {
    if (!requireNamespace("withr", quietly = TRUE)) {
        skip("'withr' package not available.")
    } else {
        seed <- 1010101
        x1 <- withr::with_seed(seed, simulateDE(ref, ref$Treatment))
        x2 <- withr::with_seed(seed, simulateDE(ref, ref$Treatment))
        expect_identical(x1, x2)
    }
})
