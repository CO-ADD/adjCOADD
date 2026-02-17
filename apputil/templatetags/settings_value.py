from django import template
from django.conf import settings

register = template.Library()

class AssignNode(template.Node):
    def __init__(self, name, value):
        self.name = name
        self.value = value

    def render(self, context):
        context[self.name] = getattr(settings, self.value.resolve(context, True), "")
        return ''

# settings value
@register.simple_tag
def settings_value(name):
    if name in ['DEVELOPMENT','DBHOST_NAME','VERSION','DJANGO_VER','PYTHON_VER']:
        return getattr(settings, name, "")
    else:
        return ""