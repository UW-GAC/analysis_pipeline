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

    # With inverse_normal = TRUE
    expect_equal("rankInvNorm(resid(x)) ~ a + b + (1|r) + (1|s) + var(g)",
                 modelString(outcome, covars, random, group_var, inverse_normal = TRUE))
    expect_equal("rankInvNorm(resid(x)) ~ a + (1|r) + var(g)",
                 modelString(outcome, covars[1], random[1], group_var, inverse_normal = TRUE))
    expect_equal("rankInvNorm(resid(x)) ~ (1|r) + (1|s) + var(g)",
                 modelString(outcome, NULL, random, group_var, inverse_normal = TRUE))
    expect_equal("rankInvNorm(resid(x)) ~ a + b + var(g)",
                 modelString(outcome, covars, NULL, group_var, inverse_normal = TRUE))
    expect_equal("rankInvNorm(resid(x)) ~ ",
                 modelString(outcome, NULL, NULL, NULL, inverse_normal = TRUE))
})

test_that("modelOutcomeString", {
  outcome <- "x"
  expect_equal(modelOutcomeString(outcome), "x")
  expect_equal(modelOutcomeString(outcome, inverse_normal = TRUE),  "rankInvNorm(resid(x))")
})
