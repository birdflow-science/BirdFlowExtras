#' Calculate migratory connectivity from BirdFlowR routes
#'
#' **WARNING** this function is experimental and new and likely to change
#' a bunch.
#'
#' @param rts Output from [BirdFlowR::route()]
#' @param bf The BirdFlow model used to make `rts`
#' @param exact logical. Whether to match route time steps exactly to period requested with `...` (TRUE)
#'  or use the route's closest available time steps (FALSE).
#' @inheritDotParams BirdFlowR::lookup_timestep_sequence -x
#' @return migratory connectivity estimated from the `rts` object.
#' @export
#' @examples
#' bf <- BirdFlowModels::amewoo
#' # generate 100 synthetic routes
#' rts <- route(bf, 100, season = "prebreeding")
#' # calculate MC across prebreeding season
#' calc_route_mc(rts, bf, season="prebreeding")
#' # calculate MC across a subset of weeks:
#' calc_route_mc(rts, bf, start=10, end=20)
#' # set exact to false to not enforce exact matches of route timestamps and
#' # requested start and end weeks (in this example end week 30 is after the last
#' # timestamp of the input routes):
#' calc_route_mc_test(rts, bf, start=10, end=30, exact=FALSE)
calc_route_mc <- function(rts, bf, exact=TRUE, ...) {
  # Using MigConnectivity::estPSI treating the routes as tracking data

  ts <- lookup_timestep_sequence(bf, ...)
  origin_t <- ts[1]
  target_t <- ts[length(ts)]
  ids <- unique(rts$data$route_id)

  nearest_timestep <- function(rts, tstep){
    rts$data |>
      dplyr::group_by(route_id) |>
      dplyr::mutate(dt = abs(timestep-tstep)) |>
      dplyr::filter(dt == min(dt)) |>
      dplyr::select(-dt) |>
      dplyr::ungroup() |>
      as.data.frame()
  }

  if(exact){
    origin <- rts$data[rts$data$timestep %in% origin_t, ]
    target <- rts$data[rts$data$timestep %in% target_t,]
  }
  else{
    origin <- nearest_timestep(rts, origin_t)
    target <- nearest_timestep(rts, target_t)
  }

  # Make sure origin and target rows correspond
  # Not necessary for synthetic routes because likely needed for real routes
  retained_ids <- intersect(origin$route_id, target$route_id)
  dropped_ids <- ids[!ids %in% retained_ids]
  if(length(dropped_ids)>0) warning(paste("dropping", length(dropped_ids), "tracks not connecting start and end period ..."))
  origin <- origin[match(origin$route_id, retained_ids), ]  |> as.data.frame()
  target <- target[match(target$route_id, retained_ids), ]  |> as.data.frame()
  stopifnot(isTRUE(all.equal(origin$route_id, retained_ids)))
  stopifnot(isTRUE(all.equal(target$route_id, retained_ids)))

  # Initialize origin and target distributions
  origin_abun <- target_abun <- rep(0, n_active(bf))
  # Tally the origin and target cells
  target_table <- table(target$i)
  origin_table <- table(origin$i)
  # Populate origin and target distributions
  target_abun[as.numeric(names(target_table))]=target_table/sum(target_table)
  origin_abun[as.numeric(names(origin_table))]=origin_table/sum(origin_table)

  # Calculate distance matrices
  dist <- great_circle_distances(bf)  # all active cells

  # matrix product to have all compbinations for two samples for origin_abun
  # note that this outer product remains normalized, because origin_abun is normalized
  origin_abun_prod <- origin_abun %*% t(origin_abun)
  # mu_D effectively equals weighted.mean(origin_dist, origin_abun_prod)
  mu_D <- sum(dist * origin_abun_prod)
  # SD in origin distance between any two given pixels
  sd_D <- sqrt(sum((dist - mu_D)^2 * origin_abun_prod))

  # same calculation, now for target_abun
  target_abun_prod <- target_abun %*% t(target_abun)
  mu_V <- sum(dist * target_abun_prod)
  sd_V <- sqrt(sum((dist - mu_V)^2 * target_abun_prod))

  # construct transition matrix
  psi <- as.matrix(table(origin$i,target$i))
  # normalize the matrix
  psi <- psi/rowSums(psi)

  origin_abun <- origin_abun[as.numeric(rownames(psi))]
  target_abun <- target_abun[as.numeric(colnames(psi))]

  # multiply transition matrix and relative abundance
  psi_abun <- psi * origin_abun

  # standardizing origin distance matrix
  origin_std <- (dist[as.numeric(rownames(psi)), as.numeric(rownames(psi))] - mu_D) / sd_D
  target_std <- (dist[as.numeric(colnames(psi)), as.numeric(colnames(psi))] - mu_V) / sd_V

  # reduce the distance matrix and abundance vectors
  dist <- dist[as.numeric(colnames(psi)), as.numeric(rownames(psi))]

  dim(psi_abun)
  dim(origin_std)
  dim(target_std)

  # calculate MC
  #MC=sum(t(psi_abun) %*% origin_std %*% psi_abun * target_std)
  MC=sum(t(psi_abun) %*% origin_std %*% psi_abun * target_std)

  return(MC)
}

