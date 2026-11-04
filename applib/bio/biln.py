"""
BILN <-> HELM conversion and 2D structure rendering.

BILN (Boehringer Ingelheim Line Notation) conversion is provided by the
vendored MIT-licensed script `biln_converter.py` (see BILN_LICENSE.txt),
originally from https://github.com/rochoa85/BILN-converter.

Structure rendering converts BILN -> HELM -> RDKit molecule using the
`helmkit` package (https://github.com/adaliaramon/helmkit), which must be
installed separately (`pip install helmkit`).
"""
import logging

from applib.bio import biln_converter

logger = logging.getLogger(__name__)

# The vendored script references a module-level `logger` only when invoked as
# `__main__`; define it here so conversion errors are logged instead of
# raising NameError.
biln_converter.logger = logger


def biln_to_helm(biln):
    """Convert a BILN string to HELM notation.

    Raises ValueError if the BILN string cannot be parsed.
    """
    if not biln or not biln.strip():
        raise ValueError('BILN sequence is empty')
    try:
        return biln_converter.BILN(biln=biln.strip()).getHELM()
    except SystemExit as err:
        raise ValueError(f'Invalid BILN sequence') from err
    except Exception as err:
        raise ValueError(f'Invalid BILN sequence: {err}') from err


def helm_to_biln(helm):
    """Convert a HELM string to BILN notation.

    Raises ValueError if the HELM string cannot be parsed.
    """
    if not helm or not helm.strip():
        raise ValueError('HELM sequence is empty')
    try:
        return biln_converter.BILN(helm=helm.strip()).getBILN()
    except SystemExit as err:
        raise ValueError(f'Invalid HELM sequence') from err
    except Exception as err:
        raise ValueError(f'Invalid HELM sequence: {err}') from err


def _mol_to_svg(mol, width=500, height=400):
    """Render an RDKit molecule to an SVG string."""
    from rdkit.Chem import rdMolDraw2D
    d2d = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(d2d, mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()


def _helm_to_mol(helm):
    """Convert HELM to an RDKit molecule using helmkit."""
    try:
        from helmkit import Molecule
    except ImportError as err:
        raise RuntimeError('helmkit is not installed (pip install helmkit)') from err
    molecule = Molecule(helm)
    return molecule.mol


def render_biln_svg(biln, width=500, height=400):
    """Render a BILN peptide to an SVG string (BILN -> HELM -> RDKit)."""
    helm = biln_to_helm(biln)
    return _mol_to_svg(_helm_to_mol(helm), width=width, height=height)


def render_helm_svg(helm, width=500, height=400):
    """Render a HELM peptide to an SVG string."""
    return _mol_to_svg(_helm_to_mol(helm), width=width, height=height)
