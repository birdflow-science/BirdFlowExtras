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

  # Origin and destination relative abundance
  origin_abun <- get_distr(bf, origin_t)[origin_dm]
  target_abun <- get_distr(bf, target_t)[target_dm]

  print("Combining transitions...")
  # Transition probabilities
  psi <- t(combine_transitions(bf, ...))
  print("Done")

  # Double check dimensions
  stopifnot(isTRUE(all.equal(nrow(psi), sum(origin_dm))))
  stopifnot(isTRUE(all.equal(ncol(psi), sum(target_dm))))

  # matrix product to have all compbinations for two samples for origin_abun
  # note that this outer product remains normalized, because origin_abun is normalized
  origin_abun_prod <- origin_abun %*% t(origin_abun)
  # mu_D effectively equals weighted.mean(origin_dist, origin_abun_prod)
  mu_D <- sum(origin_dist * origin_abun_prod)
  # SD in origin distance between any two given pixels
  sd_D <- sqrt(sum((origin_dist - mu_D)^2 * origin_abun_prod))

  # same calculation, now for target_abun
  target_abun_prod <- target_abun %*% t(target_abun)
  mu_V <- sum(target_dist * target_abun_prod)
  sd_V <- sqrt(sum((target_dist - mu_V)^2 * target_abun_prod))

  # multiply transition matrix and relative abundance
  psi_abun <- apply(psi, 2, "*", origin_abun)

  # calculate MC
  MC=0
  print("Starting MC calculation...")
  for (i in 1:dim(origin_dist)[1]){
    print(paste(i,"/",dim(origin_dist)[1]))
    for(j in 1:dim(target_dist)[1]){
      Delta_MC=psi_abun[i,j]*((origin_dist[,i] - mu_D) / sd_D) %*% psi_abun %*% ((target_dist[j,] - mu_V) / sd_V)
      MC=MC+Delta_MC
    }
  }

  return(MC)
}
