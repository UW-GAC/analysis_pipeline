context("analysis report tests")

test_that("modelString", {
    outcome <- "x"
    covars <- c("a", "b")
    random <- c("r", "s")
    group_var <- "g"
    expect_equal("x ~ a + b + (1|r) + (1|s) + var(g)", 
                 modelString(outcome, covars, random, group_var))
    expect_equal("x ~ a + (1|r) + var(g)", 
                 modelString(outcome, covars[1], random[1], group_var))
    expect_equal("x ~ (1|r) + (1|s) + var(g)", 
                 modelString(outcome, NULL, random, group_var))
    expect_equal("x ~ a + b + var(g)", 
                 modelString(outcome, covars, NULL, group_var))
    expect_equal("x ~ ", 
                 modelString(outcome, NULL, NULL, NULL))
})
