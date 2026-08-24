# Images

This module documents the source data behind each mean unit cell: the registered image, the fitted lattice, the region used for averaging, and the tilt conditions. One page per zone axis.

## Acquisition

All images were recorded on the TitanX at NCEM at 300 kV in HAADF STEM mode, at 1024 × 1024 pixels. Each orientation was recorded as image pairs with the scan rotated by 90 degrees between members of a pair, and the pairs were registered and averaged to remove scan drift. The registered composites shown here are the input to the lattice fits.

Two practical notes on the metadata. The goniometer tilts quoted on each page come from the stage readings embedded in the microscope files, which are more precise than the values in the folder names, and which recovered the tilt angles for ZA7 where the folder name lists them as unknown. The scan rotation also varies between magnification presets, and the registered images carry arbitrary frame rotations and handedness from the registration step, so no common image frame is assumed anywhere in the analysis; all indexing uses cell metrics and beam directions only.

## Overlay conventions

On each source image the fitted lattice vectors are drawn from the fit origin, $u$ in red and $v$ in blue, lengthened for visibility. The dashed white polygon is the region within which complete unit cells were averaged. The inset shows the resulting mean unit cell tiled 2 × 2.

| Page | Dataset | $\alpha$, $\beta$ (deg) | Assignment |
| --- | --- | --- | --- |
| [ZA1](za1.md) | ZA1_3 | $-7.08$, $-2.91$ | parent $[110]$ |
| [ZA2](za2.md) | ZA2_1 | $18.98$, $-1.26$ | parent $[130]$ |
| [ZA3](za3.md) | ZA3_1 | $37.15$, $-0.22$ | parent $[010]$ |
| [ZA4](za4.md) | ZA4_1 | $-25.66$, $-4.87$ | parent $[210]$ |
| [ZA5](za5.md) | ZA5_2 | $37.15$, $-0.22$ | parent $[010]$ |
| [ZA7](za7.md) | ZA7_1 | $24.70$, $-7.97$ | domain B $[33\bar{1}]$ |
| [ZA8](za8.md) | ZA8_1 | $-24.04$, $-11.99$ | domain B $[221]$ |
| [ZA9](za9.md) | ZA9_1 | $-15.33$, $-11.16$ | domain B $[331]$ |
