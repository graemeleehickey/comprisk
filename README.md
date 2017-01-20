# Comparison of joint models for competing risks and longitudinal data

This repository contains the code for each joint model referenced and implemented in the article:

> Hickey GL, Philipson P, Jorgensen A, Kolamunnage-Dona R. A comparison of different joint models for longitudinal and competing risks data: with application to an epilepsy drug randomised control trial. *In submission*.

Each model script is labelled by the first author of the original study. The example study data is available by running

    library('joineR')
    data(epileptic)

and is also available as a file in this repository.

## Implementation

To allow for fair comparison, all models were fitted in the same computing environment, namely a HP Desktop PC with an Intel® Core™ i5 3.30GHz processor with 8 GB of RAM running Windows 7 Enterprise SP1 (Microsoft Corporation ©, Redmond, WA). Models 1, 3, and 5 were fitted using 64-bit R version 3.2.3 (R Foundation for Statistical Computing, Vienna, Austria).

The R packages JM (version 1.4-2) and lcmm (version 1.7.4) were used respectively for Models 3 and 5. Models 1 and 3 rely on R packages lme (version 3.1-122) and survival (version 1.3-17). Model 2 was fitted using a program compiled from C code in Cygwin version 2.4.1 (Red Hat, Inc., Raleigh, NC) that depends on access to the GNU Scientific Library (version 2.1, Free Software Foundation, Inc., Boston, MA).

Model 4 was fitted using WinBUGS version 1.4.3 (Medical Research Council Biostatistics Unit, Cambridge, UK), which is called from within R using the package R2WinBUGS version 2.1-21. For all models, the computation time was recorded. In all cases, default settings were used unless stated otherwise below.
