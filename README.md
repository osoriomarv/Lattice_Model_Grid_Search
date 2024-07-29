# 3-D Grid Search Lattice Strain Model

This MATLAB script implements a 3-D grid search to fit a lattice strain model for rare earth elements (REE) partitioning between clinopyroxene and melt. The model fits for Young's Modulus, Partition Coefficient, and Ionic Radius using experimental data for univariant +3 REE cations.

## Features

- Implements the lattice strain model equation from Blundy and Wood (1994), Nature
- Uses experimental data from Hauri et al. (1994) for long-hour experiments (69 Hr @ 1430°C and 2.5 GPa)
- Performs a 3-D grid search over Young's Modulus, Partition Coefficient (D0), and Ionic Radius
- Calculates chi-squared values to find the best fit
- Generates plots for chi-squared values and the fitted model

## Key Equations

The script implements the lattice strain model based on the following equations:

1. Lattice Strain Model (Blundy and Wood, 1994):

   $D_i = D_0 \exp\left(\frac{-4\pi E N_A}{RT} \left[\frac{r_0}{2}(r_0 - r_i)^2 - \frac{1}{3}(r_0 - r_i)^3\right]\right)$

   Where:
   - $D_i$: Partition coefficient for element i
   - $D_0$: Strain-free partition coefficient
   - $E$: Young's Modulus
   - $N_A$: Avogadro's number
   - $R$: Gas constant
   - $T$: Temperature in Kelvin
   - $r_0$: Optimum ionic radius
   - $r_i$: Ionic radius of element i

2. Chi-squared calculation:

   $\chi^2 = \sum \frac{(D_{obs} - D_i)^2}{D_{obs}}$

   Where:
   - $D_{obs}$: Observed partition coefficient from experimental data
   - $D_i$: Calculated partition coefficient from the model

The script performs a grid search to find the optimal values of $E$, $D_0$, and $r_0$ that minimize the chi-squared value.

## Data Sources

- REE ionic radii: Shannon and Prewitt
- Partition coefficients: Hauri et al. (1994)
- 
## Outputs

1. Figure 1: Chi-squared values for iterations
2. Figure 2: 3-D Lattice Strain Model for REE^3+

## Parameters

The script performs a grid search over the following ranges:

- Young's Modulus: 350-400 GPa (60 steps)
- Partition Coefficient (D0): 0.5-0.9 (200 steps)
- Ionic Radius: 0.8-1.3 Å (100 steps)

## Notes

- The script uses experimental data for La, Ce, Nd, Sm, Eu, Dy, Er, Yb, and Lu.
- Temperature is set to 1430°C (1703.15 K) based on the experimental conditions.

## References

1. Blundy, J., & Wood, B. (1994). Prediction of crystal-melt partition coefficients from elastic moduli. Nature, 372(6505), 452-454.
2. Hauri, E. H., Wagner, T. P., & Grove, T. L. (1994). Experimental and natural partitioning of Th, U, Pb and other trace elements between garnet, clinopyroxene and basaltic melts. Chemical Geology, 117(1-4), 149-166.
3. Shannon, R. D., & Prewitt, C. T. (1969). Effective ionic radii in oxides and fluorides. Acta Crystallographica Section B: Structural Crystallography and Crystal Chemistry, 25(5), 925-946.
