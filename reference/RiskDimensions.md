# RiskDimensions

A list containing the bivariate correlations (cormat) of the 9
dimensions on which participants in Fischhoff et al. (1978) rated
different activities and technologies as well as the sample size (N).
This was then analyzed together with ratings of the risks and benefits
of these activities and technologies.

## Usage

``` r
RiskDimensions
```

## Format

A list of 2 with elements "cormat" (9 x 9 matrix of bivariate
correlations) and "N" (scalar). The correlation matrix contains the
following risk dimensions:

- Voluntariness:

  (numeric) - Voluntariness of exposure to the risk.

- Immediacy:

  (numeric) - Immediacy of the risk's effect.

- Known to exposed:

  (numeric) - How well the risk is known to those exposed to it.

- Known to science:

  (numeric) - How well the risk is known to science.

- Controllability:

  (numeric) - Controllability of the risk.

- Newness:

  (numeric) - Newness of the risk.

- Chronic:

  (numeric) - Whether the risk is chronic rather than catastrophic.

- Common:

  (numeric) - Whether the risk is common rather than dreaded.

- Severity of consequences:

  (numeric) - Severity of the consequences.

## Source

Fischhoff, B, Slovic, P, Lichtenstein, S, Read, S, and Combs, B. (1978).
How safe is safe enough? A psychometric study of attitudes towards
technological risks and benefits. Policy Sciences, 9, 127-152. doi:
10.1007/BF00143739
