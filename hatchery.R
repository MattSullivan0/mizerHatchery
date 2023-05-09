# imports + `projectHatchery` and `steadyHatchery` definitions
# remotes::install_github("MattSullivan0/mizerShelf", quiet = TRUE) # small bug fixes to mizerShelf
# remotes::install_github("sizespectrum/mizerExperimental", quiet = TRUE)
library(mizer)
# library(mizerExperimental)
# library(mizerShelf)
library(assertthat) # used in projection function when overloading

# `hatchery_` functions not intended for user calls - only to be called within
# `projectHatchery` and `steadyHatchery`.
hatchery_project_simple <- 
    function(params, 
             n = params@initial_n,
             n_pp = params@initial_n_pp,
             n_other = params@initial_n_other,
             effort = params@initial_effort,
             t = 0, dt = 0.1, steps,
             resource_dynamics_fn = get(params@resource_dynamics),
             other_dynamics_fns = lapply(params@other_dynamics, get),
             rates_fns = lapply(params@rates_funcs, get), ...) {    
    # Handy things ----
    no_sp <- nrow(params@species_params) # number of species
    no_w <- length(params@w) # number of fish size bins
    idx <- 2:no_w
    # Hacky shortcut to access the correct element of a 2D array using 1D 
    # notation
    # This references the egg size bracket for all species, so for example
    # n[w_min_idx_array_ref] = n[,w_min_idx]
    w_min_idx_array_ref <- (params@w_min_idx - 1) * no_sp + (1:no_sp)
    # Matrices for solver
    a <- matrix(0, nrow = no_sp, ncol = no_w)
    b <- matrix(0, nrow = no_sp, ncol = no_w)
    S <- matrix(0, nrow = no_sp, ncol = no_w)

    # Loop over time steps ----
    for (i_time in 1:steps) {
        r <- rates_fns$Rates(
            params, n = n, n_pp = n_pp, n_other = n_other,
            t = t, effort = effort, rates_fns = rates_fns, ...)
        
        # * Update other components ----
        n_other_current <- n_other  # So that the resource dynamics can still 
        # use the current value
        for (component in names(params@other_dynamics)) {
            n_other[[component]] <-
                other_dynamics_fns[[component]](
                    params,
                    n = n,
                    n_pp = n_pp,
                    n_other = n_other_current,
                    rates = r,
                    t = t,
                    dt = dt,
                    component = component,
                    ...
                )
        }
        
        # * Update resource ----
        n_pp <- resource_dynamics_fn(params, n = n, n_pp = n_pp,
                                     n_other = n_other_current, rates = r,
                                     t = t, dt = dt,
                                     resource_rate = params@rr_pp,
                                     resource_capacity = params@cc_pp, ...)
        
        # * Update species ----
        # a_{ij} = - g_i(w_{j-1}) / dw_j dt
        a[, idx] <- sweep(-r$e_growth[, idx - 1, drop = FALSE] * dt, 2,
                          params@dw[idx], "/")
        # b_{ij} = 1 + g_i(w_j) / dw_j dt + \mu_i(w_j) dt
        b[] <- 1 + sweep(r$e_growth * dt, 2, params@dw, "/") + r$mort * dt
        # S_{ij} <- N_i(w_j)
        S[, idx] <- n[, idx, drop = FALSE]
        # Update first size group of n
        n[w_min_idx_array_ref] <-
            (n[w_min_idx_array_ref] + r$rdd * dt / 
                 params@dw[params@w_min_idx]) /
            b[w_min_idx_array_ref]
        # Update n
        # for (i in 1:no_sp) # number of species assumed small, so no need to 
        #                      vectorize this loop over species
        #     for (j in (params@w_min_idx[i]+1):no_w)
        #         n[i,j] <- (S[i,j] - A[i,j]*n[i,j-1]) / B[i,j]
        # This is implemented via Rcpp
        n <- mizer:::inner_project_loop(no_sp = no_sp, no_w = no_w, n = n,
                                A = a, B = b, S = S,
                                w_min_idx = params@w_min_idx)
        
        # * Update time ----
        t <- t + dt

        # ----------------------------------------------------------------------
        # Add hatchery species
        # ----------------------------------------------------------------------
        n <- hatchery_dynamics(params = params, n = n, t = t, dt = dt)
    }
    
    return(list(n = n, n_pp = n_pp, n_other = n_other, rates = r))
}

hatchery_projectToSteady <- function(params,
                            effort = params@initial_effort,
                            distance_func = distanceSSLogN,
                            t_per = 1.5,
                            t_max = 100,
                            dt = 0.1,
                            tol = 0.1 * t_per,
                            return_sim = FALSE,
                            progress_bar = TRUE, ...) {
    
    # function to check hatchery_params is properly defined
    check_hatchery(params)

    params <- validParams(params)
    effort <- validEffortVector(effort, params = params)
    params@initial_effort <- effort
    assert_that(t_max >= t_per,
                tol > 0)
    if ((t_per < dt) || !isTRUE(all.equal((t_per - round(t_per / dt) * dt), 0))) {
        stop("t_per must be a positive multiple of dt")
    }
    t_dimnames <-  seq(0, t_max, by = t_per)
    
    if (is(progress_bar, "Progress")) {
        # We have been passed a shiny progress object
        progress_bar$set(message = "Finding steady state", value = 0)
        proginc <- 1 / ceiling(t_max/t_per)
    }
    
    if (return_sim) {
        # create MizerSim object
        sim <- MizerSim(params, t_dimnames =  t_dimnames)
        sim@n[1, , ] <- params@initial_n
        sim@n_pp[1, ] <- params@initial_n_pp
        sim@n_other[1, ] <- params@initial_n_other
        sim@effort[1, ] <- params@initial_effort
    }
    
    # get functions
    resource_dynamics_fn <- get(params@resource_dynamics)
    other_dynamics_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)
    r <- rates_fns$Rates(
        params, n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 0, 
        effort = effort, rates_fns = rates_fns, ...)
    
    previous <- list(n = params@initial_n,
                     n_pp = params@initial_n_pp,
                     n_other = params@initial_n_other,
                     rates = r)
    
    for (i in 2:length(t_dimnames)) {
        # advance shiny progress bar
        if (is(progress_bar, "Progress")) {
            progress_bar$inc(amount = proginc)
        }
        current <- hatchery_project_simple(params, n = previous$n, n_pp = previous$n_pp,
                                  n_other = previous$n_other, t = 0,
                                  dt = dt, steps = round(t_per / dt),
                                  effort = params@initial_effort,
                                  resource_dynamics_fn = resource_dynamics_fn,
                                  other_dynamics_fns = other_dynamics_fns,
                                  rates_fns = rates_fns)
        if (return_sim) {
            # Store result
            sim@n[i, , ] <- current$n
            sim@n_pp[i, ] <- current$n_pp
            sim@n_other[i, ] <- current$n_other
            sim@effort[i, ] <- params@initial_effort
        }
        
        # Species with no reproduction are going extinct, so stop.
        extinct <- is.na(current$rates$rdd) | current$rates$rdd <= 1e-20
        if (any(extinct)) {
            warning(paste(params@species_params$species[extinct], collapse = ", "),
                    " are going extinct.")
            success <- FALSE
            distance <- NA
            break
        }
        
        distance <- distance_func(params,
                                  current = current,
                                  previous = previous, ...)
        success <- distance < tol
        if (success == TRUE) {
            break
        }
        previous <- current
    }
    if (!success) {
        message("Simulation run did not converge after ",
                (i - 1) * t_per,
                " years. Value returned by the distance function was: ",
                distance)
    } else {
        message("Convergence was achieved in ", (i - 1) * t_per, " years.")
    }
    
    params@initial_n[] <- current$n
    params@initial_n_pp[] <- current$n_pp
    params@initial_n_other[] <- current$n_other
    
    if (return_sim) {
        sim@params <- params
        sel <- 1:i
        sim@n <- sim@n[sel, , , drop = FALSE]
        sim@n_pp <- sim@n_pp[sel, , drop = FALSE]
        sim@n_other <- sim@n_other[sel, , drop = FALSE]
        sim@effort <- sim@effort[sel, , drop = FALSE]
        return(sim)
    } else {
        params@time_modified <- lubridate::now()
        return(params)
    }
}

projectHatchery <- function(object, effort,
                    t_max = 100, dt = 0.1, t_save = 1, t_start = 0,
                    initial_n, initial_n_pp,
                    append = TRUE,
                    progress_bar = TRUE, ...) {
    
    # Set and check initial values ----
    assert_that(t_max > 0)
    if (is(object, "MizerSim")) {
        validObject(object)
        params <- setInitialValues(object@params, object)
        t_start <- getTimes(object)[idxFinalT(object)]
    } else if (is(object, "MizerParams")) {
        params <- validParams(object)
        if (!missing(initial_n)) params@initial_n[] <- initial_n
        if (!missing(initial_n_pp)) params@initial_n_pp[] <- initial_n_pp
    } else {
        stop("The `object` argument must be either a MizerParams or a MizerSim object.")
    }

    # function to check hatchery_params is properly defined
    check_hatchery(params)

    initial_n <- params@initial_n
    initial_n_pp <- params@initial_n_pp
    initial_n_other <- params@initial_n_other
    
    no_sp <- length(params@w_min_idx)
    assert_that(is.array(initial_n),
                is.numeric(initial_n),
                are_equal(dim(initial_n), c(no_sp, length(params@w))))
    assert_that(is.numeric(initial_n_pp),
                length(initial_n_pp) == length(params@w_full))
    
    assert_that(is.null(initial_n_other) || is.list(initial_n_other))
    other_names <- names(params@other_dynamics)
    if (length(other_names) > 0) {
        if (is.null(names(initial_n_other))) {
            stop("The initial_n_other needs to be a named list")
        }
        if (!setequal(names(initial_n_other), other_names)) {
            stop("The names of the entries in initial_n_other do not match ",
                 "the names of the other components of the model.")
        }
    }
    
    # Set effort array ----
    if (missing(effort)) effort <- params@initial_effort
    if (is.null(dim(effort))) { # effort is a vector or scalar
        # Set up the effort array transposed so we can use the recycling rules
        # no point running a simulation with no saved results
        if (t_max < t_save) {
            t_save <- t_max
        }
        times <- seq(t_start, t_start + t_max, by = t_save)
        effort <- validEffortVector(effort, params)
        effort <- t(array(effort, 
                          dim = c(length(effort), length(times)), 
                          dimnames = list(gear = names(effort), 
                                          time = times)))
    } else {
        effort <- validEffortArray(effort, params)
    }
    
    times <- as.numeric(dimnames(effort)[[1]])
    
    # Make the MizerSim object with the right size ----
    # We only save every t_save years
    sim <- MizerSim(params, t_dimnames = times)
    # Set initial population and effort
    sim@n[1, , ] <- initial_n 
    sim@n_pp[1, ] <- initial_n_pp
    sim@n_other[1, ] <- initial_n_other
    sim@effort <- effort
    
    ## Initialise ----
    # get functions
    resource_dynamics_fn <- get(params@resource_dynamics)
    other_dynamics_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)
    
    # Set up progress bar
    if (is(progress_bar, "Progress")) {
        # We have been passed a shiny progress object
        progress_bar$set(message = "Running simulation", value = 0)
        proginc <- 1 / length(times)
    } else if (progress_bar == TRUE) {
        pb <- progress::progress_bar$new(
            format = "[:bar] :percent ETA: :eta",
            total = length(times), width = 60)
        pb$tick(0)
    }
    
    n_list <- list(n = initial_n, n_pp = initial_n_pp,
                   n_other = initial_n_other)
    t <- times[[1]]
    
    ## Loop over time ----
    for (i in 2:length(times)) {
        # number of time steps between saved times
        steps <- round((times[[i]] - t) / dt)
        # advance to next saved time
        n_list <- hatchery_project_simple(
            params, n = n_list$n, n_pp = n_list$n_pp, n_other = n_list$n_other,
            t = t, dt = dt, steps = steps, 
            effort = effort[i - 1, ],
            resource_dynamics_fn = resource_dynamics_fn,
            other_dynamics_fns = other_dynamics_fns,
            rates_fns = rates_fns, ...)
        # Calculate start time for next iteration
        # The reason we don't simply use the next entry in `times` is that
        # those entries may not be separated by exact multiples of dt.
        t <- t + steps * dt
        # Advance progress bar
        if (is(progress_bar, "Progress")) {
            progress_bar$inc(amount = proginc)
        } else if (progress_bar == TRUE) {
            pb$tick()
        }
        
        # Store result
        sim@n[i, , ] <- n_list$n
        sim@n_pp[i, ] <- n_list$n_pp
        sim@n_other[i, ] <- n_list$n_other
    }
    
    # append to previous simulation ----
    if (is(object, "MizerSim") && append) {
        no_t_old <- dim(object@n)[1]
        no_t <- length(times)
        new_t_dimnames <- c(as.numeric(dimnames(object@n)[[1]]),
                            times[2:no_t])
        new_sim <- MizerSim(params, t_dimnames = new_t_dimnames)
        old_indices <- 1:no_t_old
        new_indices <- seq(from = no_t_old + 1, length.out = no_t - 1)
        new_sim@n[old_indices, , ]  <- object@n
        new_sim@n[new_indices, , ]  <- sim@n[2:no_t, , ]
        new_sim@n_pp[old_indices, ] <- object@n_pp
        new_sim@n_pp[new_indices, ] <- sim@n_pp[2:no_t, ]
        new_sim@n_other[old_indices, ]  <- object@n_other
        new_sim@n_other[new_indices, ]  <- sim@n_other[2:no_t, ]
        new_sim@effort[old_indices, ] <- object@effort
        new_sim@effort[new_indices, ] <- sim@effort[2:no_t, ]
        return(new_sim)
    }
    return(sim)
}

steadyHatchery <- function(params, t_max = 100, t_per = 1.5, dt = 0.1,
                   tol = 0.1 * dt, return_sim = FALSE, 
                   preserve = c("reproduction_level", "erepro", "R_max"),
                   progress_bar = TRUE) {
    params <- validParams(params)
    
    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        preserve <- match.arg(preserve)
        old_reproduction_level <- getReproductionLevel(params)
        old_R_max <- params@species_params$R_max
        old_erepro <- params@species_params$erepro
    }
    
    # Force the reproduction to stay at the current level
    params@species_params$constant_reproduction <- getRDD(params)
    old_rdd_fun <- params@rates_funcs$RDD
    params@rates_funcs$RDD <- "constantRDD"
    
    # Force resource to stay at current level
    old_resource_dynamics <- params@resource_dynamics
    params@resource_dynamics <- "resource_constant"
    
    # Force other components to stay at current level
    old_other_dynamics <- params@other_dynamics
    for (res in names(params@other_dynamics)) {
        params@other_dynamics[[res]] <- "constant_other"
    }
    
    object <- hatchery_projectToSteady(params,
                              distance_func = distanceMaxRelRDI,
                              t_per = t_per,
                              t_max = t_max,
                              dt = dt,
                              tol = tol,
                              return_sim = return_sim,
                              progress_bar = progress_bar)
    if (return_sim) {
        params <- object@params
    } else {
        params <- object
    }
    # Restore original RDD and dynamics
    params@rates_funcs$RDD <- old_rdd_fun
    params@other_dynamics <- old_other_dynamics
    params@species_params$constant_reproduction <- NULL
    
    # Set resource dynamics
    params <- setResource(params, resource_dynamics = old_resource_dynamics)
    
    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        if (preserve == "reproduction_level") {
            params <- setBevertonHolt(params, 
                                      reproduction_level = old_reproduction_level)
        } else if (preserve == "R_max") {
            params <- setBevertonHolt(params, 
                                      R_max = old_R_max)
        } else {
            params <- setBevertonHolt(params, erepro = old_erepro)
        }
    }
    
    if (return_sim) {
        object@params <- params
        return(object)
    } else {
        params@time_modified <- lubridate::now()
        return(params)
    }
}

# check hatchery_params is properly defined
check_hatchery <- function(params) {
    # warn if hatchery_params undefined in global scope
    if (!exists("hatchery_params")) {
        message("Warning: `hatchery_params` is not defined, using standard projection.")
    } else {
        # stop if hatchery contains species not present in the model
        for (species in hatchery_params$species) {
            if (!(species %in% species_params(params)$species)) {
                stop("Hatchery species \"", species, "\" is not in the model")
            }
        }
    }
    return(NULL)
}

# specify hatchery dynamics function - called after all other changes to `n` in
# `my_project_simple` function.

hatchery_dynamics <- function(params, n, t, dt) {
    if (exists("hatchery_params")) {

        # support multiple species in `hatchery_params`
        for (row in nrow(hatchery_params)) {
            idx = which(params@species_params$species == 
                        hatchery_params$species[row])

            # size-dependent component of hatchery (w units = grams) -----------
            # compute normal distribution of lobsters to add
            mu <- hatchery_params$mu[row]
            sigma <- hatchery_params$sigma[row]
            # normal distribution NOT in log weight => looks skewed on log scale
            w_dist <- dnorm(params@w, mean=mu, sd=sigma, log=FALSE)

            # ensure the distribution integral is 1 for number density scaling
            w_dist <- w_dist / sum(w_dist * params@dw)

            # time-dependent component of hatchery (t units = years) -----------
            t_dependence <- hatchery_T(t)
            
            # set scaled rates using h(w, t) = W(w) * T(t) ---------------------
            n[idx, ] <- n[idx, ] +
                (hatchery_params$annual_N * w_dist) * t_dependence * dt
        }
    }
    return(n)
}

# hatchery_T defines the time dependence of the hatchery input - this expects
# the integral over 1 year to equal 1 (or average to 1 if multi-year hatching 
# cycle so that the `annual_N` value follows intuition of amount added per year)
hatchery_T <- function(t) {return(1)} # no time dependence

getChange <- function(sim) {
    biomass <- getBiomass(sim)

    return(biomass[length(biomass)] / biomass[1])
}
