# Change log

## v0.2.6
- Support numpy>=2.0.0

## v0.2.5
- Fix symmetry-adapted displacements for noisy and lower-symmetry inputs [[#27]](https://github.com/phonopy/spgrep-modulation/pull/27)
    - Require `spgrep>=0.3.3` for robust treatment for selecting symmetry-adapted bases

## v0.2.4
- Fix {func}`get_translational_subgroup` for Gamma point

## v0.2.3

Initial release to PyPI
- Calculate representation matrices and irreps formed by phonon eigenmodes
- Calculate isotropy subgroups of irreps of space groups on the fly
- Generate modulated structures in selected order-parameter directions of isotropy subgroups
