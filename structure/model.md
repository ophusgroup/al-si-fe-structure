# The τ₄ Model

The particle is the ternary τ₄ phase, Al₃FeSi₂, in the tetragonal structure determined by Guéneau et al. [](#gueneau1995). It occurs in the particle as two orientation domains of the same structure, described below.

## The phase

| | |
| --- | --- |
| Space group | $I4/mcm$ (No. 140) |
| Lattice parameters | $a = 6.061$ Å, $c = 9.525$ Å |
| Pearson symbol | $tI24$, PdGa$_5$ type, $Z = 4$ |
| Fe | $4a$ $(0, 0, \tfrac{1}{4})$ |
| Al | $4c$ $(0, 0, 0)$ |
| Al/Si mixed | $16l$ $(x, x+\tfrac{1}{2}, z)$, $x = 0.152$, $z = 0.145$ |

The Fe atoms sit in columns along $c$, each Fe at the center of a square antiprism of Al and Si. The $16l$ site carries a mixed Al/Si occupancy in the average structure; in the ordered ground state the Al and Si separate onto distinct sites, lowering the symmetry to $Pbcn$ with an essentially unchanged cell [](#gueneau1995), [](#fang2025).

:::{figure} ../assets/figures/tau4_cell.png
:alt: The tau4 unit cell projected along [100], [110], and [001]
:width: 100%

The τ₄ conventional cell projected along three directions. Fe in red, the Al $4c$ site in dark gray, the mixed Al/Si $16l$ site in light gray.
:::

Use the interactive model below to rotate the structure and snap to the zone axes recorded in the experiment. The projected potential mode reproduces the mean unit cell patterns live for any orientation, which makes clear how the higher index zones like $[331]$ produce their characteristic images: the beam runs nearly along close-packed Fe rows, and the projected columns arrange into the wavy centred pattern seen in the ZA7 and ZA9 data.

:::{anywidget} ../widgets/structure-3d.js
:::

## Composition

EDS quantification on the TitanX, for the same particle as the imaging data:

| Region | Al (at.%) | Si (at.%) | Fe (at.%) |
| --- | --- | --- | --- |
| Thin edge | 52.49 | 29.44 | 18.37 |
| Opposite side | 52.68 | 29.46 | 18.24 |
| Whole particle | 52.70 | 29.31 | 18.41 |

Writing the phase as (Al$_{1-x}$Si$_x$)$_5$Fe following Fang et al. [](#fang2025), the measured Si fraction is $x = 0.36$, inside their computed stability window of $x = 0.33$ to $0.45$. The measured Fe content of 18.3 percent is about 2 percent above the reported homogeneity range. The three measurements agree between thin and thick regions, so absorption is not the cause; a standardless k-factor error of this size is plausible, and we do not read the excess as significant.

## Two orientation domains

The nine zone-axis datasets do not all come from one crystal orientation. Five index on one τ₄ orientation and three on a second τ₄ orientation, and the two groups are internally consistent: within each group the goniometer tilts, the projected cells, and the simulated images all agree with a single rigid lattice. All nine datasets were recorded in one microscope session without remounting the sample, so the particle contains two orientation domains, and each image shows one of them. The domains are volumes of the crystal, not alternating layers, and no image in this dataset shows both at once.

To be concrete: datasets ZA1 through ZA5 image domain A, along its $[110]$, $[130]$, $[010]$, $[210]$, and $[010]$ directions. Datasets ZA7, ZA8, and ZA9 image domain B, along its $[33\bar{1}]$, $[221]$, and $[331]$ directions. Both domains are the same τ₄ structure; only the lattice orientation differs.

:::{figure} ../assets/figures/domain_schematic.png
:alt: Schematic of the two orientation domains sharing the pseudo-cubic Fe subcell
:width: 92%

Cross-section schematic of the two domains. The Fe sublattice of τ₄ is a nearly cubic subcell (4.29 × 4.29 × 4.76 Å), and the tetragonal $c$ axis of domain B points along what is approximately a subcell $a$ axis of domain A. The boundary between the domains has not been imaged.
:::

## What relates the two domains

The measured misorientation between the domains is a rotation of 91.4 degrees about an axis close to the $[110]$ direction of domain A. This is near the classic ferroelastic twin relationship of a tetragonal phase built on a pseudo-cubic subcell: the subcell is cubic to within 11 percent, so the crystal can switch its unique axis to another subcell direction at low cost, the same physics as domain formation in martensites and ferroelectrics.

We are deliberately careful with the word twin. Comparing the measured misorientation against every exact low-index twin operation of τ₄, modulo the full point group, the closest match ($\{112\}$ reflection, equivalently 90 degrees about $\langle 100 \rangle$ of the subcell) still differs by 10 to 15 degrees, mostly in the rotation axis. The expected twin obliquity from the 11 percent subcell tetragonality accounts for about 6 degrees of this; the remainder may reflect boundary relaxation, an intermediate twin chain, or lattice rotation in this heavily faulted particle. The safe statement is that the two domains are related by an approximately 90 degree switch of the tetragonal axis on the shared subcell, and the exact boundary crystallography is not yet determined. Imaging the boundary directly would settle it.

## Stacking disorder and antiphase boundaries

Within the domains, the images show weak diffuse intensity at half-order positions along $c^*$, at reciprocal lattice points with $h+k+l$ odd. These are forbidden in $I4/mcm$ and are exactly the superstructure positions of the ordered $Pbcn$ arrangement. The diffuse peaks are two to three times broader than the fundamentals and streaked along $c^*$: the Al/Si ordering is short range, organized as thin (001) slabs separated by antiphase boundaries with displacement $\tfrac{1}{2}[111]$, the translation lost in the ordering transition. No sharp superstructure reflection appears anywhere in the dataset. This is consistent with the computed antisite cost of 0.07 to 0.10 eV and ordering temperature near 660 K [](#fang2025), and with the profuse planar disorder of the related $\beta$ phase [](#becker2019).

The defect hierarchy of the particle is therefore: orientation domains at the largest scale, (001) antiphase boundaries within domains, and Al/Si site disorder at the finest scale. The domain boundaries and the antiphase boundaries are the targets of the planned first-principles calculations.

## Structure files

CIF files for the τ₄ phase and the other Al-Fe-Si phases considered during the analysis, from the Crystallography Open Database:

- [COD 2010454, τ₄ average structure, I4/mcm](../assets/cif/COD_2010454_delta_FeAl3Si2_I4mcm.cif) (the model used on this site)
- [COD 2010455, τ₄ ordered, Pbcn](../assets/cif/COD_2010455_delta_FeAl3Si2_Pbcn_ordered.cif)
- [COD 2101145, β-Al₄.₅FeSi, Rømming 1994](../assets/cif/COD_2101145_beta_Al4.5FeSi_Romming1994.cif)
- [COD 2002891, β-Al₄.₅FeSi, Hansen 1998](../assets/cif/COD_2002891_beta_Al4.5FeSi_Hansen1998.cif)
- [COD 2106287, α-AlFeSi hexagonal](../assets/cif/COD_2106287_alpha_AlFeSi_P63mmc.cif)
- [COD 2004239, τ₇ Al₃Fe₂Si₃](../assets/cif/COD_2004239_gamma_Fe2Al3Si3_P21n.cif)
- [COD 1533819, λ-AlFeSi](../assets/cif/COD_1533819_lambda_AlFeSi_R-3.cif)
- [COD 2005761, τ₁ Al₂Fe₃Si₃](../assets/cif/COD_2005761_Fe3Al2Si3_triclinic.cif)
- [COD 2005762, τ₈ Al₂Fe₃Si₄](../assets/cif/COD_2005762_Fe3Al2Si4_Cmcm.cif)
- [COD 2107329, β alternative setting](../assets/cif/COD_2107329_beta_Al9Fe2Si2_C2c_split.cif)
- [COD 4330511, Fe₈Al₁₇.₄Si₇.₆](../assets/cif/COD_4330511_Fe8Al17.4Si7.6_P21c.cif)
- [Al, fcc reference](../assets/cif/Al_fcc_Fm-3m.cif)
- [Si, diamond reference](../assets/cif/Si_diamond_Fd-3m.cif)
