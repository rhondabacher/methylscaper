library(testthat)

# Mock data for testing
dataIn <- list(
    gch = data.frame(matrix(c(1, 0, 1, 2, 0, 1, 1, 0, 1, 1), nrow = 5)),
    hcg = data.frame(matrix(c(1, 0, 2, 3, 0, 1, 4, 0, 3, 2), nrow = 5))
)
coordinatesObject <- list(weight_start = 1, weight_stop = 5, weight_color = "red", refine_start = 1, refine_stop = 5)


# Test initialOrder function
test_that("initialOrder creates a valid orderObject", {
    orderObject <- initialOrder(dataIn, Method = "PCA")
    expect_type(orderObject, "list")
    expect_true("toClust" %in% names(orderObject))
    expect_true("order1" %in% names(orderObject))
})

# Test buildOrderObjectShiny function
test_that("buildOrderObjectShiny returns initialOrder with weights", {
    coordinatesObject$weight_start <- 1
    coordinatesObject$weight_stop <- 5
    orderObject <- buildOrderObjectShiny(dataIn, "PCA", coordinatesObject, NULL)
    expect_type(orderObject, "list")
    expect_true("toClust" %in% names(orderObject))
    expect_true("order1" %in% names(orderObject))
})

# Test refineOrderShiny function
test_that("refineOrderShiny applies refinement correctly", {
    orderObject <- buildOrderObjectShiny(dataIn, "PCA", coordinatesObject, NULL)
    refinedOrder <- refineOrderShiny(orderObject, "PCA", coordinatesObject)
    
    expect_type(refinedOrder, "integer")  # Expecting an integer vector for order
    expect_equal(length(refinedOrder), length(orderObject$order1))  # Same length as original order
})

# Test refineFunction directly
test_that("refineFunction reorders correctly", {
    orderObject <- initialOrder(dataIn, Method = "PCA")
    refinedOrder <- refineFunction(orderObject, refineStart = 1, refineEnd = 3, Method = "PCA")
    
    expect_type(refinedOrder, "integer")  # Expecting an integer vector for order
    expect_equal(length(refinedOrder), length(orderObject$order1))  # Same length as original order
})

# Test forceReverse function
test_that("forceReverse correctly reverses a subset of the order", {
    orderObject <- initialOrder(dataIn, Method = "PCA")
    reversedOrder <- forceReverse(orderObject, reverseStart = 1, reverseEnd = 3)

    expect_type(reversedOrder, "integer")  # Expecting an integer vector for order
    expect_equal(length(reversedOrder), length(orderObject$order1))  # Same length as original order
    expect_equal(reversedOrder[1:3], rev(orderObject$order1[1:3]))  # Check if the subset is reversed
})


# Test handleBrushCoordinates
test_that("handleBrushCoordinates returns expected list structure", {
    plot_brush <- list(xmin = 0, xmax = 1, ymin = 0.1, ymax = 0.5)
    result <- handleBrushCoordinates(plot_brush, n = 10, m = 5)
    expect_named(result, c("first_row", "last_row", "first_col", "last_col", "weight_color"))
    expect_true(result$weight_color %in% c("red", "yellow"))
})