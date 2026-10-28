import csv
from pathlib import Path

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from apputil.models import Menu, Menu_Item


def _to_bool(value):
    return str(value).strip().lower() in ('true', '1', 'yes', 'y', 't')


def _to_int(value):
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return 0


class Command(BaseCommand):
    help = "Load navigation menus and menu items from the provided CSV files"

    def add_arguments(self, parser):
        parser.add_argument(
            '--menu-file', default='data/menu/app_menu.csv',
            help="Path to the Menu CSV file (default: data/menu/app_menu.csv)",
        )
        parser.add_argument(
            '--item-file', default='data/menu/app_menu_item.csv',
            help="Path to the Menu_Item CSV file (default: data/menu/app_menu_item.csv)",
        )
        parser.add_argument(
            '--base-dir', default=None,
            help="Base directory for relative file paths (default: settings.BASE_DIR)",
        )

    def _resolve(self, path, base_dir):
        p = Path(path)
        if not p.is_absolute():
            p = Path(base_dir) / p
        return p

    def handle(self, *args, **options):
        base_dir = options['base_dir'] or settings.BASE_DIR

        menu_file = self._resolve(options['menu_file'], base_dir)
        item_file = self._resolve(options['item_file'], base_dir)

        # ---- Load Menus ---------------------------------------------------
        if not menu_file.exists():
            raise CommandError(f"Menu CSV file not found: {menu_file}")

        n_menu = 0
        with open(menu_file, newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                menu_id = (row.get('menu_id') or '').strip()
                if not menu_id:
                    continue
                menu, created = Menu.objects.update_or_create(
                    menu_id=menu_id,
                    defaults={
                        'menu_name': (row.get('menu_name') or menu_id).strip(),
                        'menu_label': (row.get('menu_label') or row.get('menu_name') or menu_id).strip(),
                        'icon': (row.get('icon') or '').strip(),
                        'order': _to_int(row.get('order')),
                        'permission': (row.get('permission') or '').strip(),
                    },
                )
                n_menu += 1
                self.stdout.write(f"  Menu {menu_id}: {'created' if created else 'updated'}")

        # ---- Load Menu Items ---------------------------------------------
        if not item_file.exists():
            raise CommandError(f"Menu Item CSV file not found: {item_file}")

        n_item = 0
        with open(item_file, newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                item_id = (row.get('item_id') or '').strip()
                if not item_id:
                    continue

                menu_id = (row.get('menu_id') or '').strip()
                try:
                    menu = Menu.objects.get(menu_id=menu_id)
                except Menu.DoesNotExist:
                    raise CommandError(f"Menu {menu_id} not found for item {item_id} - load menus first")

                parent = None
                parent_id = (row.get('parent') or '').strip()
                if parent_id:
                    try:
                        parent = Menu_Item.objects.get(item_id=parent_id)
                    except Menu_Item.DoesNotExist:
                        raise CommandError(f"Parent item {parent_id} not found for item {item_id}")

                item, created = Menu_Item.objects.update_or_create(
                    item_id=item_id,
                    defaults={
                        'menu_id': menu,
                        'parent': parent,
                        'label': (row.get('label') or '').strip(),
                        'url_name': (row.get('url_name') or '').strip(),
                        'url_path': (row.get('url_path') or '').strip(),
                        'icon': (row.get('icon') or '').strip(),
                        'order': _to_int(row.get('order')),
                        'separator': _to_bool(row.get('separator')),
                        'disabled': _to_bool(row.get('disabled')),
                        'modal': _to_bool(row.get('modal')),
                        'permission': (row.get('permission') or '').strip(),
                    },
                )
                n_item += 1
                self.stdout.write(f"  Item {item_id}: {'created' if created else 'updated'}")

        self.stdout.write(self.style.SUCCESS(f"Loaded {n_menu} menus and {n_item} menu items"))
