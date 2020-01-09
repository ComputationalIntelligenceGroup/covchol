#' Execution of a Gaussian graphical model experiment
#'
#' @param p Vector of dimensions to test.
#' @param d Vector of densities to test.
#' @param r Number of replications for the experiment.
#' @param experiment Function containing the atomic experiment to be executed.
#' As such, it must require at least two arguments, `p` and `d`.
#' @param ... Additional parameters for `experiment`.
#'
#' @details This function executes `r` times the experiment  `experiment` for
#' each value of `p` and `d`, and stores its result as an `rds` object in a
#' newly created directory for each repetition.
#'
#' @return Nothing
#' 
#' @export
execute <- function(p, d, r, experiment, ...) {
	wd <- getwd()
	case_name <- outer(p, d, paste, sep = "_")
	case_fname <- matrix(paste0(case_name, ".rds"), ncol = length(d))
	n_cores <- parallel::detectCores() - 2
	
	cl <- parallel::makeCluster(n_cores, outfile = "")
	doParallel::registerDoParallel(cl)
	
	exp_name <- substitute(experiment)
	repetition <- 0 # For avoiding R CMD check NOTE
	iter <- foreach::foreach(repetition = 1:r, .combine = rbind)
	foreach::"%dopar%"(iter, {
		dir_name <- paste0(exp_name, "_r", repetition)
  	dir.create(paste0(wd, "/", dir_name), showWarnings = FALSE)
  	
		for (i in seq_along(p)) {
			for (j in seq_along(d)) {
				result <- experiment(p = p[i], d = d[j], ...)
				saveRDS(result, file = paste0(dir_name, "/", exp_name, "_", case_fname[i, j])) 
			}
		}
	})
	
	parallel::stopCluster(cl)
}


#' Plot the results of a Gaussian graphical model experiment
#'
#' @param p Vector of dimensions to test.
#' @param d Vector of densities to test.
#' @param r Number of replications for the experiment.
#' @param N Number of samples obtained in each atomic experiment.
#' @param map Function that will map the values obtained from each atomic
#' experiment for each sample.
#' @param reduce Function that will reduce the previously mapped results from
#' each sample, yielding a single value for each atomic experiment.
#' @param show_sd Whether to show the standard deviation for `map` function.
#' @param exp_name Name of the atomic experiment that has been executed.
#' @param ... Additional arguments for `map`.
#' 
#' @details This function plots the results from a previously executed
#' experiment related to a Gaussian graphical model. It expects the results
#' stored in the same format as `excute()`.
#'
#' @return The generated plot.
#' @export
plot_map_reduce <- function(p, d, r, N, map = function(x) {
                              return(x)
                            }, reduce, show_sd = FALSE, exp_name, ...) {
  data <- matrix(
    nrow = length(p), ncol = length(d),
    dimnames = list(p = p, d = d)
  )
  data_sd <- matrix(
    nrow = length(p), ncol = length(d),
    dimnames = list(p = p, d = d)
  )


  for (i in 1:length(p)) {
    exp_res <- array(dim = c(p[i], p[i], N * r))
    for (j in 1:length(d)) {
      for (rep in 1:r) {
        exp_res[, , ((rep - 1) * N + 1):(rep * N)] <-
          readRDS(file = paste0(exp_name, "_r", rep, "/", exp_name, "_", p[i],
		  "_", d[j], ".rds"))
      }
      mapd_mat <- apply(X = exp_res, MARGIN = 3, FUN = map, ...)
      data[i, j] <- reduce(mapd_mat)
      data_sd[i, j] <- stats::sd(mapd_mat)
    }
  }

  wd <- getwd()
  dir.create(paste0(wd, "/plot_", r), showWarnings = FALSE)

  palette <- grDevices::colorRampPalette(colors = c("black", "red"))
  colors <- palette(length(d))

  df <- data %>% as.tbl_cube(met_name = "data") %>% as_tibble()
  df$d <- as.factor(df$d)
  df$data_sd <- data_sd %>% as.tbl_cube(met_name = "data_sd") %>% as_tibble()
  
  pl <- ggplot(df, aes(x = p, y = data, group = d, color = d)) +
    geom_line() +
    geom_point() +
    theme(text = element_text(size = 20), legend.position = "bottom") +
    scale_color_manual(values = colors) +
    xlab("Number of nodes")

  if (show_sd == TRUE) {
    pl <- pl + geom_ribbon(aes(ymin = data - data_sd, ymax = data + data_sd))
  }
  
  return(pl)
}

