from django.apps import apps
from django.db.models import ForeignKey

# Get the Author model
Author = apps.get_model('myapp', 'Author')  # Replace 'myapp' with your app name

# Get all models in the app
all_models = apps.get_models()

# Find models with Author as a foreign key
models_with_author_fk = []
for model in all_models:
    for field in model._meta.fields:
        if isinstance(field, ForeignKey) and field.related_model == Author:
            models_with_author_fk.append(model.__name__)


# Get the Author model
Author = apps.get_model('myapp', 'Author')  # Replace 'myapp' with your app name

# Replace 'old_author_name' and 'new_author_name' with the actual names
old_author_name = 'Old Author Name'
new_author_name = 'New Author Name'

# Find the old and new authors
old_author = Author.objects.get(name=old_author_name)
new_author = Author.objects.get(name=new_author_name)

# Get all models in the app
all_models = apps.get_models()

# Update author_id in models with Author as a foreign key
for model in all_models:
    for instance in model.objects.filter(author=old_author):
        instance.author = new_author
        instance.save()