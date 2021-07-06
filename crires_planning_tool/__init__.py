import logging
import tqdm

try:
    import colorlog
except ImportError:
    # install colorlog for colored log output
    colorlog = None

class TqdmLoggingHandler(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

tqdm.tqdm.get_lock()

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

console = TqdmLoggingHandler()
console.setLevel(logging.DEBUG)
if colorlog is not None:
    console.setFormatter(
        colorlog.ColoredFormatter("%(log_color)s%(levelname)s - %(message)s")
    )
logger.addHandler(console)

del logging
del colorlog

__all__ = ["etc_form", "interactive_graph", "nasa_exoplanets_archive", "transit_planner"]
