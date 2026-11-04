"""Import pyPept monomer SDF into dpeptide.Monomer.

Usage:
    python manage.py import_pypept_monomers [--path dpeptide/monomers/pyPept_Monomers.sdf] [--upload]

Without --upload the command only validates/parses the SDF and reports what
would be imported.
"""
import os

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from rdkit import Chem

from apputil.models import ApplicationUser
from dpeptide.models import Monomer


def _parse_int_list(value, length=4):
    items = [v.strip() for v in (value or '').split(',')]
    out = [None if v in ('', 'None') else int(v) for v in items]
    while len(out) < length:
        out.append(None)
    return out[:length]


def _parse_str_list(value, length=4):
    items = [v.strip() for v in (value or '').split(',')]
    out = [None if v in ('', 'None') else v for v in items]
    while len(out) < length:
        out.append(None)
    return out[:length]


class Command(BaseCommand):
    help = "Import pyPept monomer SDF into dpeptide.Monomer"

    def add_arguments(self, parser):
        parser.add_argument(
            '--path',
            default=os.path.join('dpeptide', 'monomers', 'pyPept_Monomers.sdf'),
            help='Path to the pyPept monomers SDF',
        )
        parser.add_argument('--upload', action='store_true', default=False)
        parser.add_argument('--user', default='J.Zuegg', help='User name for created/updated records')

    def handle(self, *args, **options):
        path = options['path']
        if not os.path.isabs(path):
            path = os.path.join(os.getcwd(), path)
        if not os.path.isfile(path):
            raise CommandError(f'Monomer SDF not found: {path}')

        owner = ApplicationUser.objects.get(name= options['user'])
        created = updated = 0

        supplier = Chem.SDMolSupplier(path)
        for mol in supplier:
            if mol is None:
                continue
            symbol = mol.GetProp('symbol').strip()

            entry = Monomer.objects.filter(symbol=symbol).first()
            if entry is None:
                entry = Monomer(symbol=symbol)
                created += 1
            else:
                updated += 1

            entry.monomer_name = mol.GetProp('m_name').strip()
            entry.monomer_abbr = mol.GetProp('m_abbr').strip()
            entry.monomer_type = mol.GetProp('m_type').strip()
            entry.monomer_subtype = mol.GetProp('m_subtype').strip()
            entry.nat_analog = mol.GetProp('natAnalog').strip()
            entry.pdb_name = mol.GetProp('pdbName').strip()
            entry.rgroups = _parse_str_list(mol.GetProp('m_Rgroups'))
            entry.rgroup_idx = _parse_int_list(mol.GetProp('m_RgroupIdx'))
            entry.attachment_idx = _parse_int_list(mol.GetProp('m_attachmentPointIdx'))
            entry.mol = mol
            entry.molfile = Chem.MolToMolBlock(mol)

            if options['upload']:
                with transaction.atomic(using='dpeptide'):
                    entry.save(user=owner)

        self.stdout.write(self.style.SUCCESS(
            f'Parsed {created + updated} monomers '
            f'({created} new, {updated} existing). '
            f'Upload: {options["upload"]}'
        ))
        if not options['upload']:
            self.stdout.write('Run with --upload to save to the database.')
