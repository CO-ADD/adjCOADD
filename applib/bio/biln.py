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


def sync_sequences(biln, helm, fasta, changed=None):
    """Synchronise BILN / HELM / FASTA from the edited field.

    Only the fields derived from the source field are regenerated. The source
    is the highest-priority non-empty field in ``changed``; priority is
    BILN > HELM > FASTA. If no edited field is non-empty, the values are
    returned unchanged (no sync). On conversion failure the derived field is
    cleared (empty).

    :param changed: iterable of field names edited in this save
        (subset of {'biln', 'helm', 'fasta'}).
    :return: (biln, helm, fasta)
    """
    biln = (biln or '').strip()
    helm = (helm or '').strip()
    fasta = (fasta or '').strip()
    changed = set(changed or [])
    values = {'biln': biln, 'helm': helm, 'fasta': fasta}

    # Pick the highest-priority edited field that is non-empty.
    source = None
    for field in ('biln', 'helm', 'fasta'):
        if field in changed and values[field]:
            source = field
            break

    if source is None:
        # Nothing was edited (or the edited field was cleared) -> no sync.
        return biln, helm, fasta

    if source == 'biln':
        try:
            helm = biln_to_helm(biln)
        except Exception:
            helm = ''
        try:
            fasta = biln_to_fasta(biln)
        except Exception:
            fasta = ''
    elif source == 'helm':
        try:
            biln = helm_to_biln(helm)
        except Exception:
            biln = ''
        if biln:
            try:
                fasta = biln_to_fasta(biln)
            except Exception:
                fasta = ''
        else:
            fasta = ''
    elif source == 'fasta':
        try:
            biln = fasta_to_biln(fasta)
        except Exception:
            biln = ''
        if biln:
            try:
                helm = biln_to_helm(biln)
            except Exception:
                helm = ''
        else:
            helm = ''

    return biln, helm, fasta


def biln_to_mol(biln):
    """Convert a BILN string to an RDKit molecule using pyPept.

    The monomer library is read from the database (dpeptide.Monomer) when
    available, falling back to pyPept's bundled SDF library otherwise.
    """
    from pyPept.sequence import Sequence
    from pyPept.molecule import Molecule

    monomer_df = _monomer_info_from_db()
    if monomer_df is not None:
        import pyPept.sequence as pseq
        # Sequence.__init__ looks up get_monomer_info in the pyPept.sequence
        # module namespace, so patching it makes pyPept read monomers from the
        # database instead of the bundled SDF file.
        pseq.get_monomer_info = lambda path: monomer_df

    seq = Sequence(biln)
    return Molecule(seq).get_molecule(fmt='ROMol')


def _monomer_info_from_db():
    """Return the pyPept monomer DataFrame from the database, or None."""
    try:
        from dpeptide.models import Monomer
        df = Monomer.spec_from_database()
        if df is not None and not df.empty:
            return df
    except Exception:
        logger.exception('Failed to load monomer library from database')
    return None


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
