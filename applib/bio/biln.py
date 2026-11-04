"""
BILN / HELM / FASTA conversion and 2D structure rendering, via pyPept.

pyPept (https://github.com/Boehringer-Ingelheim/pyPept) uses BILN as its
master notation and converts between BILN, HELM and FASTA. It requires
BioPython and RDKit and is installed with:

    pip install git+https://github.com/Boehringer-Ingelheim/pyPept.git

pyPept is imported lazily so that the rest of the application can import this
module without the package being installed.
"""
import logging
import re

logger = logging.getLogger(__name__)

_STANDARD_AA = set('ACDEFGHIKLMNPQRSTVWY')


def _converter():
    from pyPept.converter import Converter
    return Converter


def _parser():
    from pyPept.biln.parser import BILNParser
    return BILNParser


def biln_to_helm(biln):
    """Convert a BILN string to HELM notation. Raises ValueError on failure."""
    biln = (biln or '').strip()
    if not biln:
        raise ValueError('BILN sequence is empty')
    try:
        return _converter()(biln=biln).get_helm()
    except Exception as err:
        raise ValueError(f'Invalid BILN sequence: {err}') from err


def helm_to_biln(helm):
    """Convert a HELM string to BILN notation. Raises ValueError on failure."""
    helm = (helm or '').strip()
    if not helm:
        raise ValueError('HELM sequence is empty')
    try:
        return _converter()(helm=helm).get_biln()
    except Exception as err:
        raise ValueError(f'Invalid HELM sequence: {err}') from err


def fasta_to_biln(fasta):
    """Convert a FASTA one-letter sequence to BILN (natural amino acids only).

    Raises ValueError if the sequence contains non-standard residues.
    """
    fasta = ''.join(str(fasta or '').split()).upper()
    if not fasta:
        return ''
    invalid = [ch for ch in fasta if ch not in _STANDARD_AA]
    if invalid:
        raise ValueError(f'FASTA contains non-standard residues: {"".join(invalid)}')
    return '-'.join(fasta)


def biln_to_fasta(biln):
    """Convert a BILN sequence to a FASTA one-letter sequence.

    Only the main peptide chain is considered (branched chains after '.' are
    ignored). Only the 20 natural amino acids can be represented in FASTA, so
    a modified/non-natural residue raises ValueError (the caller clears FASTA).
    """
    biln = (biln or '').strip()
    if not biln:
        return ''

    main_chain = biln.split('.')[0]
    residues = []
    for token in main_chain.split('-'):
        res = re.sub(r'\([^)]*\)', '', token).strip()
        if res in _STANDARD_AA:
            residues.append(res)
        else:
            raise ValueError(f'Residue "{res}" cannot be represented in FASTA')
    return ''.join(residues)


def is_valid_biln(biln):
    """Return True if the BILN string is valid according to pyPept's parser."""
    biln = (biln or '').strip()
    if not biln:
        return False
    try:
        return bool(_parser().IsValidSequence(biln))
    except Exception:
        return False


def sync_biln_fasta(biln, fasta):
    """Sync BILN (master) and FASTA on save.

    BILN is the master notation: when present it overwrites FASTA. When only
    FASTA is provided (import), BILN is derived from it. If either conversion
    fails, the derived field is cleared (empty).
    """
    biln = (biln or '').strip()
    fasta = (fasta or '').strip()

    if biln:
        try:
            fasta = biln_to_fasta(biln)
        except Exception:
            fasta = ''
        return biln, fasta

    if fasta:
        try:
            biln = fasta_to_biln(fasta)
        except Exception:
            biln = ''
        return biln, fasta

    return '', ''


def biln_to_mol(biln):
    """Convert a BILN string to an RDKit molecule using pyPept."""
    from pyPept.sequence import Sequence
    from pyPept.molecule import Molecule
    seq = Sequence(biln)
    return Molecule(seq).get_molecule(fmt='ROMol')


def write_2d_svg(pep_mol, width=500, height=400):
    """Render an RDKit molecule (e.g. from pyPept) to a 2D SVG string."""
    from rdkit.Chem import rdMolDraw2D
    d2d = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(d2d, pep_mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()


def render_biln_svg(biln, width=500, height=400):
    """Render a BILN peptide to an SVG string."""
    return write_2d_svg(biln_to_mol((biln or '').strip()), width=width, height=height)


def render_helm_svg(helm, width=500, height=400):
    """Render a HELM peptide to an SVG string."""
    return render_biln_svg(helm_to_biln(helm), width=width, height=height)
