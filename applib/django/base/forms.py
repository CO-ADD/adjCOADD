"""
Shared base ModelForm providing consistent, full-width, modern field rendering
for all Create/Update forms across the application.

- Adds Bootstrap form-control / form-select classes to standard widgets so
  text inputs, selects and textareas fill the available width.
- Renders long CharFields (max_length >= TEXTAREA_MIN_LENGTH) as multi-line
  Textareas.
"""
from django import forms

# Widgets that keep their native appearance (no form-control class)
_PLAIN_WIDGETS = (
    forms.CheckboxInput,
    forms.CheckboxSelectMultiple,
    forms.RadioSelect,
    forms.HiddenInput,
)

# CharField variants that must stay single-line regardless of max_length
_SINGLE_LINE_CHARFIELDS = (
    forms.EmailField,
    forms.URLField,
    forms.SlugField,
    forms.RegexField,
    forms.UUIDField,
    forms.GenericIPAddressField,
)


class Base_ModelForm(forms.ModelForm):
    """
    ModelForm base class for consistent form field rendering.
    """

    # CharFields whose max_length is at least this value are shown as textareas.
    TEXTAREA_MIN_LENGTH = 200

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._style_fields()

    def _style_fields(self):
        for field in self.fields.values():
            widget = field.widget

            if isinstance(widget, forms.HiddenInput):
                continue

            # Large text fields -> multi-line textarea
            if (
                isinstance(field, forms.CharField)
                and not isinstance(field, _SINGLE_LINE_CHARFIELDS)
                and not isinstance(widget, (forms.Textarea, forms.PasswordInput))
                and field.max_length is not None
                and field.max_length >= self.TEXTAREA_MIN_LENGTH
            ):
                attrs = dict(widget.attrs)
                attrs.setdefault('rows', '3')
                field.widget = forms.Textarea(attrs=attrs)
                widget = field.widget

            # Checkboxes / radios keep their native appearance
            if isinstance(widget, _PLAIN_WIDGETS):
                continue

            # Bootstrap full-width styling
            attrs = dict(widget.attrs)
            classes = attrs.get('class', '')
            if isinstance(widget, forms.Select):
                if 'form-select' not in classes and 'form-control' not in classes:
                    attrs['class'] = f'{classes} form-select'.strip()
            else:
                if 'form-control' not in classes:
                    attrs['class'] = f'{classes} form-control'.strip()
            widget.attrs = attrs
