# Set Control Parameters for Simulation

This function sets control parameters for running simulations,
particularly for MCMC methods. It allows users to specify the number of
simulations, burn-in period, thinning interval, and various other
parameters necessary for the simulation.

## Usage

``` r
set_control_sim(
  n_sim = 12000,
  burnin = 2000,
  thin = 10,
  h = NULL,
  c1.h = 0.01,
  c2.h = 1e-04,
  linear_model = FALSE
)
```

## Arguments

- n_sim:

  Integer. The total number of simulations to run. Default is 12000.

- burnin:

  Integer. The number of initial simulations to discard (burn-in period,
  used for the MCMC algorithm). Default is 2000.

- thin:

  Integer. The interval at which simulations are recorded (thinning
  interval, used for the MCMC algorithm). Default is 10.

- h:

  Numeric. An optional parameter. Must be non-negative if specified.

- c1.h:

  Numeric. A control parameter for the simulation. Must be positive.
  Default is 0.01.

- c2.h:

  Numeric. Another control parameter for the simulation. Must be between
  0 and 1. Default is 1e-04.

- linear_model:

  Logical. If TRUE, the function sets up parameters for a linear model
  and only returns `n_sim`. Default is FALSE.

## Value

A list of control parameters for the simulation with class attribute
"mcmc.RiskMap".

## Details

The function validates the input parameters and ensures they are
appropriate for the simulation that is used in the
[`glgpm`](https://claudiofronterre.github.io/RiskMap/reference/glgpm.md)
fitting function. For non-linear models, it checks that `n_sim` is
greater than `burnin`, that `thin` is positive and a divisor of
`(n_sim - burnin)`, and that `h`, `c1.h`, and `c2.h` are within their
respective valid ranges.

If `linear_model` is TRUE, only `n_sim` and `linear_model` are required,
and the function returns a list containing these parameters.

If `linear_model` is FALSE, the function returns a list containing
`n_sim`, `burnin`, `thin`, `h`, `c1.h`, `c2.h`, and `linear_model`.

## See also

[`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html),
[`forceSymmetric`](https://rdrr.io/pkg/Matrix/man/forceSymmetric-methods.html)

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>

## Examples

``` r
# Example with default parameters
control_params <- set_control_sim()

# Example with custom parameters
control_params <- set_control_sim(n_sim = 15000, burnin = 3000, thin = 20)
```
