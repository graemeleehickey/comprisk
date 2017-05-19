# Comparison of joint models for competing risks and longitudinal data

This repository contains the code for each joint model referenced and implemented in the article:

> Hickey GL, Philipson P, Jorgensen A, Kolamunnage-Dona R. A comparison of different joint models for longitudinal and competing risks data: with application to an epilepsy drug randomised control trial. *In submission*.

Each model script is labelled by the first author of the original study. The example study data is available by running

library('joineR') # version 1.1.0 and earlier
    data(epileptic)

and is also available as a file (`epileptic.txt`) in this repository. **Note** that if using `joineR` version 1.2.0 or above, the `epileptic` data has a slightly different column ordering, which the code in this repository relies on. To avoid any issues, please load the data from the `epileptic.txt` file.

## Implementation

To allow for fair comparison, all models were fitted in the same computing environment, namely a HP Desktop PC with an Intel® Core™ i5 3.30GHz processor with 8 GB of RAM running Windows 7 Enterprise SP1 (Microsoft Corporation ©, Redmond, WA). Models 1, 3, and 5 were fitted using 64-bit R version 3.2.3 (R Foundation for Statistical Computing, Vienna, Austria).

The R packages `JM` (version 1.4-2) and `lcmm` (version 1.7.4) were used respectively for Models 3 and 4. Models 1 and 3 rely on R packages `lme` (version 3.1-122) and `survival` (version 1.3-17). Model 2 was fitted using a program compiled from C code in `Cygwin` version 2.4.1 (Red Hat, Inc., Raleigh, NC) that depends on access to the GNU Scientific Library (version 2.1, Free Software Foundation, Inc., Boston, MA).

For all models, the computation time was recorded. In all cases, default settings were used unless stated otherwise in the manuscript (Supplementary Material).

## Structure of repository

- **Model 0**: `separate.R`
- **Model 1**: `williamson2008.R` (Williamson et al., 2008)
- **Model 2**: `elashoff2008.R` and all files in `.\elashoff\` (Elashoff et al., 2008)
- **Model 3**: `rizopoulos2012.R` (Rizopoulos, 2012)
- **Model 4**: `proust-lima2015.R` (Proust-Lima et al., 2015)

An R script (`andrinopoulou2014.R`) and `BUGS` script (`andrinopoulou2014_model`) is also available for fitting a model proposed by Andrinopoulou et al. (2014).

## References

1. D. Rizopoulos, *Joint Models for Longitudinal and Time-to-Event Data: with Applications in R*. Boca Raton, FL: Chapman & Hall/CRC, 2012.

2. C. Proust-Lima, J.-F. Dartigues, and H. Jacqmin-Gadda, “Joint modelling of repeated multivariate cognitive measures and competing risks of dementia and death: a latent process and latent class approach,” *Stat. Med.*, In press., 2015.

3. P. R. Williamson, R. Kolamunnage-Dona, P. Philipson, and A. G. Marson, “Joint modelling of longitudinal and competing risks data,” *Stat. Med.*, vol. **27**, pp. 6426–6438, 2008.

4. R. M. Elashoff, G. Li, and N. Li, “A joint model for longitudinal measurements and survival data in the presence of multiple failure types,” *Biometrics*, vol. **64**, no. 3, pp. 762–771, 2008.

5. E.-R. Andrinopoulou, D. Rizopoulos, J. J. M. Takkenberg, and E. Lesaffre, “Joint modeling of two longitudinal outcomes and competing risk data,” *Stat. Med.*, vol. **33**, no. 18, pp. 3167–3178, 2014.



