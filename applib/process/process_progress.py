from django.core.cache import cache

from applib.logging.validation_log import Validation_Log


#=================================================================================================
# Progress-reporting process base class for Process
#=================================================================================================
class Process_Progress:
    """
    Base class for long-running ``*_Process`` steps that report progress.

    Subclass this instead of writing a module-level function when a process
    should expose progress (e.g. for a polling endpoint or a progress modal).

    Common plumbing provided here:

      * ``report_progress(done, total, item, message)`` - emits a progress dict
        to an optional callback and writes it to Django's cache under
        ``progress_key`` so a JSON endpoint can poll it immediately.
      * ``finish()`` - convenience call marking completion.
      * ``execute()`` - classmethod: instantiate + run in one call, mirroring
        the old function-call style.

    Subclasses must:

      * set ``process_name`` (used for the Validation_Log label),
      * implement ``run()`` returning the ``Validation_Log``.
    """

    process_name = 'Process'
    progress_timeout = 300  # seconds a progress snapshot stays in cache

    def __init__(self, request, dirname=None, file_list=None, appuser=None,
                 progress_cb=None, progress_key=None, progress_timeout=None):
        self.request = request
        self.dirname = dirname
        self.file_list = file_list
        self.appuser = appuser
        self.progress_cb = progress_cb
        self.progress_key = progress_key or self._default_progress_key()
        if progress_timeout is not None:
            self.progress_timeout = progress_timeout
        self.valLog = Validation_Log(self.process_name)
        self.progress = {'done': 0, 'total': 0, 'percent': 0, 'item': '', 'message': ''}

    # ----------------------------------------------------------------
    @staticmethod
    def _default_progress_key():
        import uuid
        return f"process_progress:{uuid.uuid4().hex}"

    # ----------------------------------------------------------------
    @property
    def n_files(self):
        return len(self.file_list) if self.file_list else 0

    # ----------------------------------------------------------------
    def report_progress(self, done, total, item='', message=''):
        """Emit a progress snapshot to the cache and/or the callback."""
        data = {
            'done': done,
            'total': total,
            'percent': int(done * 100 / total) if total else 0,
            'item': item or '',
            'message': message or '',
        }
        self.progress = data
        if self.progress_cb is not None:
            try:
                self.progress_cb(data)
            except Exception as err:
                print(f" [{self.process_name}] progress_cb error: {err}")
        if self.progress_key:
            try:
                cache.set(self.progress_key, data, timeout=self.progress_timeout)
            except Exception as err:
                print(f" [{self.process_name}] cache error: {err}")
        return data

    # ----------------------------------------------------------------
    def finish(self):
        """Mark the process as complete."""
        return self.report_progress(1, 1, '', 'Finished')

    # ----------------------------------------------------------------
    @classmethod
    def execute(cls, *args, **kwargs):
        """Instantiate and run in one call (function-style compatibility)."""
        return cls(*args, **kwargs).run()

    # ----------------------------------------------------------------
    def run(self):
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement run()")