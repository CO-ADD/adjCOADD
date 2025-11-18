from dsample.models import Project
from dsample.utils.summary import update_project_summary

qryPrj = Project.objects.all()
for djPrj in qryPrj:
    print(djPrj)
    update_project_summary(djPrj)
djPrj.save()
