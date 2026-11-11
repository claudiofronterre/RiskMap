# get_worldpop reports the product unless messages are off

    Code
      adults <- get_worldpop("LBR", 2011, ages = c(20, Inf), dest_dir = dest_dir)
    Message
      WorldPop unconstrained 2000-2020 estimates at 100m, LBR 2011, aged 20 and over: 26 files.

# get_worldpop errors name get_worldpop and give hints

    Code
      get_worldpop("LBR", 2011, ages = c(18, Inf), messages = FALSE)
    Condition
      Error in `get_worldpop()`:
      ! `ages` must start and end on WorldPop age groups.
      x You supplied `c(18, Inf)`.
      i Groups start at 0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, and 80, and the oldest is open-ended.
      i Use, for example, `c(20, Inf)` or `c(15, 49)`.
    Code
      get_worldpop("LBR", 2011, constrained = "yes")
    Condition
      Error in `get_worldpop()`:
      ! `constrained` must be `TRUE` or `FALSE`, not "yes".
    Code
      get_worldpop("LBR", 2035, constrained = TRUE)
    Condition
      Error in `get_worldpop()`:
      ! Constrained estimates can't be requested for 2035.
      i Constrained estimates cover 2015 to 2030 and unconstrained estimates (`constrained = FALSE`) cover 2000 to 2020.

