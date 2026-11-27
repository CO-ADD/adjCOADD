"""
Background tasks for dsample (CastDB).

Uses Django's `django.tasks` framework (Django 6.1+).
Requires a task backend and a worker, see adjcoadd/settings.py (TASKS).
"""
from django.tasks import task

from apputil.models import ApplicationLog
from applib.mol.abase import get_ABase_Batches, upload_ABase_Batches


@task
def task_update_ABase_Batches(test=0, upload=True, overwrite=False):
    """
    Fetch ABase compound batches from the Oracle ABase database and
    upload them into the CastDB database.

    Parameters
    ----------
    test : int
        Fetch First N rows only (0 = all rows).
    upload : bool
        Persist the records (False = dry run).
    overwrite : bool
        Update existing records as well as new ones.
    """
    try:
        dfABase = get_ABase_Batches(test=test)
        result = upload_ABase_Batches(dfABase, upload=upload, overwrite=overwrite)
        ApplicationLog.add(
            'Task', 'Update ABase', 'Info', None,
            str(test), f"ABase upload completed {result}", 'Completed',
        )
        return result
    except Exception as err:
        ApplicationLog.add(
            'Task', 'Update ABase', 'Error', None,
            str(test), f"ABase upload failed: {err}", 'Failed',
        )
        raise
