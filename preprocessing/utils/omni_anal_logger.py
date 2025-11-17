"""
Module: omni_anal_logger.py

Purpose:
    Lightweight utility for structured logging with elapsed time tracking,
    nested operations, and automatic inclusion of log level, file, and function names.
"""

import time
from datetime import datetime
import inspect
import os

class OmniAnalLogger:
    def __init__(self):
        self.indent_level = 0
        self.indent_str = "    "

    def _timestamp(self) -> str:
        now = datetime.now()
        ms = int(now.microsecond / 1000)
        return f"[{now.hour:02}:{now.minute:02}:{now.second:02}.{ms:03}]"

    def _indent(self) -> str:
        return self.indent_str * self.indent_level

    def _caller_context(self, stack_depth: int = 3) -> str:
        """
        Returns a string like 'module.py.function_name' for the original caller.
        """
        frame = inspect.stack()[stack_depth]
        module = inspect.getmodule(frame[0])
        filename = os.path.basename(module.__file__) if module and hasattr(module, '__file__') else "unknown"
        function = frame.function
        return f"{filename}.{function}"

    def _log(self, level: str, message: str) -> None:
        """
        Internal method to format and print a log message.

        Parameters:
            level (str): Log level (INFO, WARNING, ERROR)
            message (str): The message to log
        """
        context = self._caller_context()
        print(f"{self._timestamp()} [{level.upper()}] [{self._indent()}{context}] {message}")

    def info(self, message: str) -> None:
        """
        Log an informational message.

        Parameters:
            message (str): The message to print.
        """
        self._log("INFO", message)

    def warning(self, message: str) -> None:
        """
        Log a warning message.

        Parameters:
            message (str): The message to print.
        """
        self._log("WARNING", message)

    def debug(self, message: str) -> None:
        """
        Log an dubug message.

        Parameters:
            message (str): The message to print.
        """
        self._log("DEBUG", message)

    def start(self, message: str) -> float:
        """
        Log the start of a long-running operation.

        Parameters:
            message (str): Description of the operation.

        Returns:
            float: Start time (use with .end())
        """
        self.info(f"Starting: {message}")
        self.indent_level += 1
        return time.time()

    def end(self, start_time: float) -> None:
        """
        Log the end of an operation started with .start()

        Parameters:
            start_time (float): The start time returned by .start()
        """
        elapsed = time.time() - start_time
        seconds = int(elapsed)
        milliseconds = int((elapsed - seconds) * 1000)
        self.indent_level = max(0, self.indent_level - 1)
        self.info(f"Finished in {seconds:02}.{milliseconds:03}s")

    def time_block(self, message: str):
        """
        Context manager for timing a code block.

        Parameters:
            message (str): Description of the block being timed.

        Usage:
            with omni_anal_logger.time_block("Step name"):
                do_work()
        """
        return _TimerContext(self, message)


class _TimerContext:
    def __init__(self, omni_anal_logger: OmniAnalLogger, message: str):
        self.logger = omni_anal_logger
        self.message = message
        self.start_time = None

    def __enter__(self):
        self.start_time = self.logger.start(self.message)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.logger.end(self.start_time)


# Instantiate a default global logger for convenience
omni_anal_logger = OmniAnalLogger()
