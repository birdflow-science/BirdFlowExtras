#' Calc migratory connectivity from BirdFlowR routes
#'
#' **WARNING** this function is experimental and new and likely to change
#' a bunch.
#'
#' Currently, it assumes that all the routes start and end at the same
#' timestep and so uses the first and last timestep of the first route as the
#' starting and ending timesteps and then calculates the migratory connectivity
#' between those timesteps across all routes, dropping routes that don't
#' have both timesteps.
#'
#' It treats the routes as "Telemetry" information when calling
#' `MigConnectivity` functions.
#'
#' Current process is to call `MigConnectivity::estPsi` internally to calculate
#'  the transition probabilities needed to calculate `MC`.
#'  This step is really slow and then pass the result to `calcMC`.
#'
#'  It might make sense to instead call  `estMC` which I think would perform
#'  both steps at once and output uncertainty.
#'
#'  Alternatively, we could see if we can call a more efficient function to
#'  calculate the transitions probabilities as `estPsi` is doing a lot of
#'  work to calculate uncertainties that aren't used with the current process.
#'
#'
#'
#' @param rts Output from [BirdFlowR::route()]
#' @param bf  The BirdFlow model used to make `rts`
#' @return migratory connectivity estimated from the `rts` object.
#' @export
calc_route_mc <- function(rts, bf) {
  # Using MigConnectivity::estPSI treating the routes as tracking data


  ts <- rts$timestep[rts$route_id == rts$route_id[1]]
  origin_t <- ts[1]
  target_t <- ts[length(ts)]
  ids <- unique(rts$route_id)

  origin <- rts[rts$timestep %in% origin_t, ]
  target <- rts[rts$timestep %in% target_t,]

  # Make sure origin and target rows correspond
  # Not necessary for synthetic routes because likely needed for real routes
  retained_ids <- intersect(origin$route_id, target$route_id)
  origin <- origin[match(origin$route_id, retained_ids), ]  |> as.data.frame()
  target <- target[match(target$route_id, retained_ids), ]  |> as.data.frame()
  stopifnot(isTRUE(all.equal(origin$route_id, retained_ids)))
  stopifnot(isTRUE(all.equal(target$route_id, retained_ids)))
  make_pts <- function(x){
    x$geom <- vector(mode = "list", length = nrow(x))
    for(i in seq_len(nrow(x))){
      x$geom[[i]] <- sf::st_point(as.numeric(x[i, c("x","y")]))
    }
    res <- sf::st_as_sf(x, sf_column_name = "geom")
  }

  origin_pts <- make_pts(origin)
  target_pts <- make_pts(target)

  sf::st_crs(origin_pts) <- sf::st_crs(bf)
  sf::st_crs(target_pts) <- sf::st_crs(bf)


  # It seems that MigConnectivity needs origin and target assignment to be 1:n
  # Here I make a crosswalk between i and
  # a new id that ranges 1:n.  I also include a name that is "i" combined with
  # the value of i for each cell.
  # Origin
  origin_sites <- data.frame(i = sort(unique(origin$i)))
  origin_sites$name <- paste0("i", origin_sites$i)
  origin_sites$id <- seq_len(nrow(origin_sites))
  # Add new ID and name to origin
  mv <- match(origin_pts$i, origin_sites$i)
  origin_pts$id <- origin_sites$id[mv]
  origin_pts$name <- origin_sites$name[mv]
  # Target
  target_sites <- data.frame(i = sort(unique(target$i)))
  target_sites$name <- paste0("i", target_sites$i)
  target_sites$id <- seq_len(nrow(target_sites))
  # Add new ID and name to target
  mv <- match(target_pts$i, target_sites$i)
  target_pts$id <- target_sites$id[mv]
  target_pts$name <- target_sites$name[mv]

  # Calculate PSI --- SLOW!
  bf_msg("Calculating psi (transition probabilities) from routes.",
         " This is slow!\n")
  psi <- MigConnectivity::estPsi(originPoints = origin_pts,
                                 targetPoints = target_pts,
                                 originAssignment = origin_pts$id,
                                 targetAssignment = target_pts$id,
                                originNames = origin_sites$name,
                                 targetNames = target_sites$name,
                                 isTelemetry = rep(TRUE, nrow(origin_pts)))


  # Make distance matrices that correspond with origin_sites and target_sites
  bf_dist <- great_circle_distances(bf)
  origin_dist <- bf_dist[origin_sites$i, origin_sites$i]
  target_dist <- bf_dist[target_sites$i, target_sites$i]
  origin_abund <- get_distr(bf, origin_t)[origin_sites$i]

  # Re-standarize to sum to 1 - necessary but maybe problematic
  origin_abund <- origin_abund / sum(origin_abund)

  # Calculate migratory connectivity
  bf_msg("Calculating migratory connectivity.\n")
  mc <- MigConnectivity::calcMC(originDist = origin_dist,
                                targetDist = target_dist,
                                psi = psi$psi$mean,
                                originRelAbund =  origin_abund)


  return(mc)
}

