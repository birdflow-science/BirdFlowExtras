# BirdFlowExtras 0.0.0.9003
2023-5-5

`calc_route_mc()` and `calc_birdflow_mc()` functions updated. 
- Much faster and more memory efficient
- No longer depend on \pkg{MigConnectivity}  
- Added flexibility as to whether to use marginal or S&T distributions
- `calc_route_mc()` now works with the new  `BirdFlowRoutes` class facilitating
working with both synthetic and real bird movement data.

# BirdFlowExtras 0.0.0.9002
2023-12-19

New `calc_route_mc()` calculates migratory connectivity from `BirdFlowRoutes()`
object.  It doesn't work very well yet though.

# BirdFlowExtras 0.0.0.9001
2023-12-18

First function: `calc_birdflow_mc()` calculates the migratory connectivity 
of a BirdFlow object.
