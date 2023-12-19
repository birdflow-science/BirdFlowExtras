#' Function to calculate migratory connectivity (MC) on a BirdFlow object
#'
#' This function calculates migratory connectivity based on the transitions
#' in a BirdFlowR model using the raster cells as regions.
#' It relies on `MigConnectivity::calcMC()` to do the bulk of the work.
#'
#' @param bf A BirdFlow model
#' @inheritDotParams BirdFlowR::lookup_timestep_sequence -x
#' @return The migratory connectivity of the BirdFlow model over the time period
#' indicated by `...`
#' @export
#' @import BirdFlowR
#' @examples
#' bf <- BirdFlowModels::amewoo
#' mc <- calc_birdflow_mc(bf, season = "prebreeding")
#' print(mc)
#'
calc_birdflow_mc <- function(bf, ...) {

  # Figure out time
  ts <- lookup_timestep_sequence(bf, ...)
  origin_t <- ts[1]
  target_t <- ts[length(ts)]

  # Dynamic masks
  origin_dm <- get_dynamic_mask(bf, origin_t) # origin dynamic cells
  target_dm <- get_dynamic_mask(bf, target_t) # target dynamic cells

  # distance matricies
  dist <- great_circle_distances(bf)  # all active cells
  origin_dist <- dist[origin_dm, origin_dm] # origin dynamic
  target_dist <- dist[target_dm, target_dm] # target dynamic

  # Origin relative abundance
  origin_abun <- get_distr(bf, origin_t)[origin_dm]

  # Transition probabilities
  psi <- t(combine_transitions(bf, ...))

  # Double check dimensions
  stopifnot(isTRUE(all.equal(nrow(psi), sum(origin_dm))))
  stopifnot(isTRUE(all.equal(ncol(psi), sum(target_dm))))

  # Calculate MC
  mc <- MigConnectivity::calcMC(originDist = origin_dist,
                                targetDist = target_dist,
                                originRelAbund = origin_abun,
                                psi = psi)

  return(mc)
}
