from datetime import datetime
from pathlib import Path


class ProgressLogger:
    def __init__(self, config, output_dir=None):
        self.path = self._progress_log_path(config, output_dir)
        self._file = None

        if self.path is not None:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._file = self.path.open("w", encoding="utf-8")

    @staticmethod
    def _progress_log_path(config, output_dir=None):
        if getattr(config, "debug_path", None):
            return Path(config.debug_path) / "pyTimSTEM_progress.txt"

        if output_dir:
            return Path(output_dir) / "pyTimSTEM_progress.txt"

        return None

    def write(self, message):
        if self._file is None:
            return

        timestamp = datetime.now().isoformat(timespec="seconds")
        self._file.write(f"[{timestamp}] {message}\n")
        self._file.flush()

    def close(self):
        if self._file is not None:
            self._file.close()
            self._file = None


# add plots for vz, supercell-vz, aberration function, aperature function, probe function 