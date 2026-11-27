from django.conf import settings
from django.db.utils import DatabaseError
from django.urls import reverse, NoReverseMatch

from apputil.models import Menu, Menu_Item


def global_settings(request):
    return {
        'DEVELOPMENT': settings.DEVELOPMENT,
        'VERSION': settings.VERSION,
        'DBHOST_NAME': settings.DBHOST_NAME,
        'DJANGO_VER': settings.DJANGO_VER,
        'PYTHON_VER': settings.PYTHON_VER,
    }


def navbar_menus(request):
    """
    Build the navigation menus for the top navbar from the Menu/Menu_Item
    models, filtered by user permission and soft-delete status.
    """
    user = getattr(request, 'user', None)
    is_authenticated = bool(getattr(user, 'is_authenticated', False))
    is_super = is_authenticated and getattr(user, 'is_superuser', False)
    is_admin = is_authenticated and (is_super or (user.has_permission('Admin') if hasattr(user, 'has_permission') else False))

    def allowed(permission):
        if permission == 'Superuser':
            return is_super
        if permission == 'Admin':
            return is_admin
        return True

    if not is_authenticated:
        return {'nav_menus': []}

    try:
        menus = list(Menu.objects.filter(astatus__gte=0).order_by('order', 'menu_name'))
        items = list(Menu_Item.objects.filter(astatus__gte=0).order_by('order', 'label'))
    except DatabaseError:
        # Tables not migrated yet - render an empty menu.
        return {'nav_menus': []}

    by_menu = {}
    for item in items:
        by_menu.setdefault(item.menu_id_id, []).append(item)

    nav_menus = []
    for menu in menus:
        if not allowed(menu.permission):
            continue
        menu_items = by_menu.get(menu.menu_id, [])
        children_by_parent = {}
        for item in menu_items:
            if item.parent_id is not None and allowed(item.permission):
                children_by_parent.setdefault(item.parent_id, []).append(item)

        def build(item):
            href = '#'
            if item.url_name:
                try:
                    href = reverse(item.url_name)
                except NoReverseMatch:
                    href = item.url_path or '#'
            elif item.url_path:
                href = item.url_path
            return {
                'label': item.label,
                'href': href,
                'separator': item.separator,
                'disabled': item.disabled,
                'modal': item.modal,
                'children': [build(child) for child in children_by_parent.get(item.item_id, [])],
            }

        nav_menus.append({
            'label': menu.menu_label,
            'items': [build(item) for item in menu_items if item.parent_id is None and allowed(item.permission)],
        })

    return {'nav_menus': nav_menus}