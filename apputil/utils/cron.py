#from django_cron import CronJobBase, Schedule
from django.core.management import call_command
from apputil.models import ApplicationLog

# #-----------------------------------------------------------------
# class Backup_adjCOADD(CronJobBase):
# #-----------------------------------------------------------------
#     RUN_EVERY_MINS = 1

#     schedule = Schedule(run_every_mins=RUN_EVERY_MINS )
#     code = 'Backup_adjCOADD_code'

#     def do(self):
#         pass

#-----------------------------------------------------------------
def Backup_adjCOADD():
#-----------------------------------------------------------------
    print("Backup DB and Media")
    try:
        call_command('dbbackup','-z')
        call_command('mediabackup','-z')
        ApplicationLog.add('Backup', 'Backup_adjCOADD','type', None,"object","dec","finish")
    except:
        ApplicationLog.add('Backup', 'Backup_adjCOADD','type', None,"object","dec","failed")
