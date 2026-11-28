# Drives the server of mispitools_app() with shiny::testServer, so the
# displayed values are checked without opening a browser.

test_that("individual tab shows the LR from the helper functions", {
  skip_if_not_installed("shiny")
  shiny::testServer(mispitools_app(), {
    session$setInputs(indiv_type = "sex", indiv_MPs = "F",
                      indiv_obsSex = "F", indiv_eps_sex = 0.05,
                      indiv_propF = 0.5)
    expect_equal(output$indiv_lr_value, "LR = 1.900")
    expect_equal(output$indiv_lr_log, "log10(LR) = 0.279")

    session$setInputs(indiv_type = "age", indiv_MPa = 40, indiv_MPr = 5,
                      indiv_obsAge = 42, indiv_epa = 0.05)
    expect_equal(output$indiv_lr_value, sprintf("LR = %.3f", 0.95 / (10 / 80)))

    session$setInputs(indiv_type = "hair", indiv_MPc = "2", indiv_obsC = "4",
                      indiv_epc = 0.1)
    expect_equal(output$indiv_lr_value, sprintf("LR = %.3f", (0.1 / 4) / 0.07))
  })
})

test_that("combine tab multiplies the selected LRs", {
  skip_if_not_installed("shiny")
  shiny::testServer(mispitools_app(), {
    session$setInputs(comb_use_sex = TRUE, comb_use_age = TRUE,
                      comb_use_hair = TRUE,
                      comb_MPs = "F", comb_obsSex = "F", comb_eps = 0.05,
                      comb_MPa = 35, comb_MPr = 5, comb_obsAge = 30,
                      comb_epa = 0.05,
                      comb_MPc = "2", comb_obsC = "2", comb_epc = 0.1)
    # the combine tab uses 51% female, not the 50% of the individual tab
    lr_sex <- 0.95 / 0.51
    lr_age <- 0.95 / (10 / 80)
    lr_hair <- 0.9 / 0.45
    expect_equal(output$comb_lr_sex, sprintf("LR = %.2f", lr_sex))
    expect_equal(output$comb_lr_total,
                 sprintf("LR = %.2f", lr_sex * lr_age * lr_hair))

    # leaving age out drops it from the product
    session$setInputs(comb_use_age = FALSE)
    expect_equal(output$comb_lr_age, "Not included")
    expect_equal(output$comb_lr_total, sprintf("LR = %.2f", lr_sex * lr_hair))
  })
})
