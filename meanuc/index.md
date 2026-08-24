# Mean Unit Cells

Each zone-axis dataset was reduced to a mean unit cell: the lattice was fit to the drift-corrected image, and all complete unit cells inside a hand-drawn region were averaged. Averaging suppresses noise and the aperiodic parts of the defect contrast, and the result is the experimental projected structure for that orientation. This module compares each mean unit cell against a simulation from the $\tau_4$ model, one page per zone axis.

## Conventions

All lengths on these pages are calibrated with the factor 0.9563 derived from the parent-zone indexing (see [evidence](../structure/evidence.md)). Unit cell overlays use red for the $u$ vector and blue for the $v$ vector, and all intensity images use the perceptually uniform inferno colormap. The atomic model, the simulation, and the experiment on each page share the same field of view and the same cell overlay.

## Simulation

The simulation is a projected-potential approximation. Atom columns are computed by projecting the $\tau_4$ structure along the assigned zone axis, each column is rendered as a Gaussian, and Fe columns are weighted to 3.0 times the peak intensity of Al and Si columns, following the $Z^{1.7}$ scaling of HAADF contrast ($(26/13.5)^{1.7} = 3.0$). A global scan confirmed 3.0 as the best single weight across the series. The Gaussian width $\sigma$ is fit per dataset by maximizing the correlation coefficient with the experimental mean cell over $\sigma$ and the cell origin. Fitted widths are 0.55 to 0.80 Å across the series. Because the registered images have arbitrary frame handedness, the fit also tests a mirrored experimental cell and keeps the better of the two.

## Results

| Dataset | Domain | Zone | Measured cell (Å, deg) | $\sigma$ (Å) | Correlation |
| --- | --- | --- | --- | --- | --- |
| [ZA1](za1.md) | A | $[110]$ | 4.281 × 4.767, 89.9 | 0.55 | 0.95 |
| [ZA2](za2.md) | A | $[130]$ | 1.911 × 4.780, 90.1 | 0.65 | 0.99 |
| [ZA3](za3.md) | A | $[010]$ | 3.031 × 4.767, 90.2 | 0.65 | 0.98 |
| [ZA4](za4.md) | A | $[210]$ | 2.705 × 9.540, 89.7 | 0.55 | 0.97 |
| [ZA5](za5.md) | A | $[010]$ | 3.028 × 4.769, 89.9 | 0.75 | 0.98 |
| [ZA7](za7.md) | B | $[33\bar{1}]$ | 2.949 × 8.671, 89.9 | 0.70 | 1.00 |
| [ZA8](za8.md) | B | $[221]$ | 2.077 × 4.333, 92.6 | 0.80 | 0.97 |
| [ZA9](za9.md) | B | $[331]$ | 8.649 × 2.963, 90.0 | 0.75 | 0.99 |

The two independent $[010]$ datasets (ZA3 and ZA5) and the two mirror-related domain B zones (ZA7 and ZA9) agree with each other, which is a useful internal consistency check on the whole pipeline.
