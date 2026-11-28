# Records what the Shiny app displays over a grid of inputs, so a refactor
# of R/app_unified.R can be checked against the version before it.
#
#   Rscript Toys/Ibrahim/app_snapshot.R before.rds   # on main
#   Rscript Toys/Ibrahim/app_snapshot.R after.rds    # on the branch
#   Rscript -e 'all.equal(readRDS("before.rds"), readRDS("after.rds"))'
#
# Uses shiny::testServer, so no browser is needed.

out_file <- commandArgs(TRUE)[1]
if (is.na(out_file)) out_file <- "app_snapshot.rds"

suppressMessages(devtools::load_all(".", quiet = TRUE))
app <- mispitools_app()
snap <- list()

shiny::testServer(app, {
  # individual tab, sex
  for (mps in c("F", "M")) for (obs in c("F", "M"))
    for (eps in c(0, 0.01, 0.05, 0.3)) for (pf in c(0.3, 0.5, 0.51)) {
      session$setInputs(indiv_type = "sex", indiv_MPs = mps, indiv_obsSex = obs,
                        indiv_eps_sex = eps, indiv_propF = pf)
      key <- paste("sex", mps, obs, eps, pf)
      snap[[key]] <<- c(output$indiv_lr_value, output$indiv_lr_log,
                        output$indiv_formula)
    }

  # individual tab, age (includes ranges that hit the 1 and 80 limits)
  for (mpa in c(3, 25, 40, 78)) for (mpr in c(1, 5, 10))
    for (obs in c(1, 20, 35, 40, 50, 79, 85)) for (epa in c(0.01, 0.05)) {
      session$setInputs(indiv_type = "age", indiv_MPa = mpa, indiv_MPr = mpr,
                        indiv_obsAge = obs, indiv_epa = epa)
      key <- paste("age", mpa, mpr, obs, epa)
      snap[[key]] <<- c(output$indiv_lr_value, output$indiv_lr_log,
                        output$indiv_formula)
    }

  # individual tab, hair
  for (mpc in 1:5) for (obs in 1:5) for (epc in c(0.01, 0.1, 0.4)) {
    session$setInputs(indiv_type = "hair", indiv_MPc = as.character(mpc),
                      indiv_obsC = as.character(obs), indiv_epc = epc)
    key <- paste("hair", mpc, obs, epc)
    snap[[key]] <<- c(output$indiv_lr_value, output$indiv_lr_log,
                      output$indiv_formula)
  }

  # combine tab
  for (us in c(TRUE, FALSE)) for (ua in c(TRUE, FALSE)) for (uh in c(TRUE, FALSE))
    for (obs_sex in c("F", "M")) for (obs_age in c(30, 60)) for (obs_c in c("2", "4")) {
      session$setInputs(comb_use_sex = us, comb_use_age = ua, comb_use_hair = uh,
                        comb_MPs = "F", comb_obsSex = obs_sex, comb_eps = 0.05,
                        comb_MPa = 35, comb_MPr = 5, comb_obsAge = obs_age,
                        comb_epa = 0.05, comb_MPc = "2", comb_obsC = obs_c,
                        comb_epc = 0.1)
      key <- paste("comb", us, ua, uh, obs_sex, obs_age, obs_c)
      snap[[key]] <<- c(output$comb_lr_sex, output$comb_lr_age,
                        output$comb_lr_hair, output$comb_lr_total,
                        output$comb_lr_log, output$comb_formula)
    }
})

saveRDS(snap, out_file)
cat(length(snap), "input combinations recorded in", out_file, "\n")
