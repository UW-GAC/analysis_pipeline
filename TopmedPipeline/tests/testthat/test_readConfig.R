context("config tests")

test_that("setConfigDefaults", {
    config <- setNames(as.character(1:10), letters[1:10])
    req <- c("a","b")
    opt <- c(c="x", d="y")
    expect_equal(suppressMessages(setConfigDefaults(config, req, opt)),
                 c(a="1", b="2", c="3", d="4"))
    expect_equal(suppressMessages(setConfigDefaults(config[1:2], req, opt)),
                 c(a="1", b="2", c="x", d="y"))
    expect_error(suppressMessages(setConfigDefaults(config["a"], req, opt)))
    expect_message(setConfigDefaults(config, req, opt), "unused parameters")
})
