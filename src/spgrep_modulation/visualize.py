"""Utilities for visualization in notebook."""

from __future__ import annotations

import math
from itertools import product

import numpy as np
from nglview import NGLWidget, show_pymatgen
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element, Species
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.phonopy import get_pmg_structure

COLOR_SCHEMES = {
    "jmol": {
        "H": [255, 255, 255],
        "He": [217, 255, 255],
        "Li": [204, 128, 255],
        "Be": [194, 255, 0],
        "B": [255, 181, 181],
        "C": [144, 144, 144],
        "N": [48, 80, 248],
        "O": [255, 13, 13],
        "F": [144, 224, 80],
        "Ne": [179, 227, 245],
        "Na": [171, 92, 242],
        "Mg": [138, 255, 0],
        "Al": [191, 166, 166],
        "Si": [240, 200, 160],
        "P": [255, 128, 0],
        "S": [255, 255, 48],
        "Cl": [31, 240, 31],
        "Ar": [128, 209, 227],
        "K": [143, 64, 212],
        "Ca": [61, 255, 0],
        "Sc": [230, 230, 230],
        "Ti": [191, 194, 199],
        "V": [166, 166, 171],
        "Cr": [138, 153, 199],
        "Mn": [156, 122, 199],
        "Fe": [224, 102, 51],
        "Co": [240, 144, 160],
        "Ni": [80, 208, 80],
        "Cu": [200, 128, 51],
        "Zn": [125, 128, 176],
        "Ga": [194, 143, 143],
        "Ge": [102, 143, 143],
        "As": [189, 128, 227],
        "Se": [255, 161, 0],
        "Br": [166, 41, 41],
        "Kr": [92, 184, 209],
        "Rb": [112, 46, 176],
        "Sr": [0, 255, 0],
        "Y": [148, 255, 255],
        "Zr": [148, 224, 224],
        "Nb": [115, 194, 201],
        "Mo": [84, 181, 181],
        "Tc": [59, 158, 158],
        "Ru": [36, 143, 143],
        "Rh": [10, 125, 140],
        "Pd": [0, 105, 133],
        "Ag": [192, 192, 192],
        "Cd": [255, 217, 143],
        "In": [166, 117, 115],
        "Sn": [102, 128, 128],
        "Sb": [158, 99, 181],
        "Te": [212, 122, 0],
        "I": [148, 0, 148],
        "Xe": [66, 158, 176],
        "Cs": [87, 23, 143],
        "Ba": [0, 201, 0],
        "La": [112, 212, 255],
        "Ce": [255, 255, 199],
        "Pr": [217, 255, 199],
        "Nd": [199, 255, 199],
        "Pm": [163, 255, 199],
        "Sm": [143, 255, 199],
        "Eu": [97, 255, 199],
        "Gd": [69, 255, 199],
        "Tb": [48, 255, 199],
        "Dy": [31, 255, 199],
        "Ho": [0, 255, 156],
        "Er": [0, 230, 117],
        "Tm": [0, 212, 82],
        "Yb": [0, 191, 56],
        "Lu": [0, 171, 36],
        "Hf": [77, 194, 255],
        "Ta": [77, 166, 255],
        "W": [33, 148, 214],
        "Re": [38, 125, 171],
        "Os": [38, 102, 150],
        "Ir": [23, 84, 135],
        "Pt": [208, 208, 224],
        "Au": [255, 209, 35],
        "Hg": [184, 184, 208],
        "Tl": [166, 84, 77],
        "Pb": [87, 89, 97],
        "Bi": [158, 79, 181],
        "Po": [171, 92, 0],
        "At": [117, 79, 69],
        "Rn": [66, 130, 150],
        "Fr": [66, 0, 102],
        "Ra": [0, 125, 0],
        "Ac": [112, 171, 250],
        "Th": [0, 186, 255],
        "Pa": [0, 161, 255],
        "U": [0, 143, 255],
        "Np": [0, 128, 255],
        "Pu": [0, 107, 255],
        "Am": [84, 92, 242],
        "Cm": [120, 92, 227],
        "Bk": [138, 79, 227],
        "Cf": [161, 54, 212],
        "Es": [179, 31, 212],
        "Fm": [179, 31, 186],
        "Md": [179, 13, 166],
        "No": [189, 13, 135],
        "Lr": [199, 0, 102],
        "Rf": [204, 0, 89],
        "Db": [209, 0, 79],
        "Sg": [217, 0, 69],
        "Bh": [224, 0, 56],
        "Hs": [230, 0, 46],
        "Mt": [235, 0, 38],
    }
}


class ColorScheme:
    """Color scheme for atoms."""

    def __init__(self, scheme="jmol") -> None:
        """Initialize color scheme."""
        # TODO: Add more color_scheme
        self._color_scheme = COLOR_SCHEMES[scheme]

    def get_color(self, site: PeriodicSite) -> list[float]:
        """Return color based on species."""
        key = self._get_key(site.specie)
        color = [c / 256 for c in self._color_scheme[key]]
        return color

    def get_hex_color(self, site: PeriodicSite) -> str:
        """Return color hex based on species."""
        key = self._get_key(site.specie)
        color = self._color_scheme[key]
        hex = f"#{color[0]:02x}{color[1]:02x}{color[2]:02x}"
        return hex

    def _get_key(self, specie: Element | Species) -> str:
        if isinstance(specie, Element):
            return str(specie)
        else:
            # with oxidation state
            return str(specie.element)


def viewer(
    cell: PhonopyAtoms,
    show_unitcell: bool = True,
    show_axes: bool = True,
    width: int | None = None,
    height: int | None = None,
) -> NGLWidget:
    """Return NGLView's widget for structure.

    Parameters
    ----------
    cell: PhonopyAtoms
    show_unitcell: Iff true, show frame of unit cell
    show_axes: Iff true, show a, b, and c axes
    width: in pixel
    height: in pixel

    Returns
    -------
    view: NGLWidget
    """
    # Show image atoms near unit cell
    structure = get_pmg_structure(cell)
    locals_and_ghosts = get_local_and_ghost_sites(structure)
    structure_display = Structure.from_sites(locals_and_ghosts)

    view = show_pymatgen(structure_display)
    view.clear()
    view.center()

    cc = ColorScheme(scheme="jmol")

    view = _add_sites(view, cc, [site for site in structure_display])

    if show_unitcell:
        view.add_unitcell()

    if show_axes:
        view = _add_axes(view, structure.lattice.matrix)

    # ref: https://github.com/nglviewer/nglview/issues/900
    view.control.spin([1, 0, 0], -math.pi / 2)
    view.control.spin([0, 0, 1], math.pi * 0.45)
    view.control.spin([0, 1, 0], math.pi * 0.1)
    view.camera = "perspective"

    if (width is not None) and (height is not None):
        view._set_size(w=f"{width}px", h=f"{height}px")

    return view


def get_local_and_ghost_sites(structure: Structure, eps: float = 1e-8) -> list[PeriodicSite]:
    """Return structure with ghost atoms near unit cell."""
    # Wrap frac_coords in [0, 1)
    wrapped_sites = []
    for site in structure:
        frac_coords = np.remainder(site.frac_coords, 1)
        wrapped_sites.append(
            PeriodicSite(
                species=site.species,
                coords=frac_coords,
                lattice=structure.lattice,
                properties=site.properties,
            )
        )

    locals_and_ghosts = []
    for site in wrapped_sites:
        for jimage in product([0, 1 - eps], repeat=3):
            # Skip original site
            if np.allclose(jimage, 0):
                continue

            new_frac_coords = site.frac_coords + np.array(jimage)
            if np.all(new_frac_coords < 1 + eps):
                new_site = PeriodicSite(
                    species=site.species,
                    coords=new_frac_coords,
                    lattice=structure.lattice,
                    properties=site.properties,
                )
                locals_and_ghosts.append(new_site)

    for site in wrapped_sites:
        locals_and_ghosts.append(site)

    return locals_and_ghosts


def _add_sites(
    view: NGLWidget,
    cc: ColorScheme,
    sites: list[PeriodicSite],
) -> NGLWidget:
    for i, si in enumerate(sites):
        # ref: https://github.com/nglviewer/nglview/issues/913
        # selection=[i] is equivalent to f"@{i}" in "selection language".
        # See https://nglviewer.org/ngl/api/manual/usage/selection-language.html
        hex_color = cc.get_hex_color(si.site)
        view.add_spacefill(
            radius=0.5,
            selection=[i],
            color=hex_color,
        )

    return view


def _add_axes(view: NGLWidget, matrix: np.ndarray) -> NGLWidget:
    # Ref: https://github.com/pyiron/pyiron_atomistics/blob/c5df5e87745d7b575463f7b2a0b588e18007dc40/pyiron_atomistics/atomistics/structure/_visualize.py#L388-L403
    axes_start = matrix[0] + np.array([0, -2, 0])
    arrow_radius = 0.1
    text_size = 1
    text_color = [0, 0, 0]  # black
    arrow_colors = [
        [1.0, 0.0, 0.0],  # Red
        [0.0, 1.0, 0.0],  # Green
        [0.0, 0.0, 1.0],  # Blue
    ]
    arrow_names = ["a", "b", "c"]

    for i, (arrow_color, arrow_name) in enumerate(zip(arrow_colors, arrow_names)):
        start = list(axes_start)
        basis_i = matrix[i]
        shift = basis_i / np.linalg.norm(basis_i)
        end = list(axes_start + shift)

        view.shape.add_arrow(start, end, arrow_color, arrow_radius, f"{arrow_name}-axis")
        view.shape.add_text(end, text_color, text_size, arrow_name)

    return view
