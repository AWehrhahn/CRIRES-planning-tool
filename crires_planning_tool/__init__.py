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
logging.captureWarnings(True)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

__console__ = TqdmLoggingHandler()
__console__.setLevel(logging.WARNING)
if colorlog is not None:
    __console__.setFormatter(
        colorlog.ColoredFormatter("%(log_color)s%(levelname)s - %(name)s - %(message)s")
    )
logger.addHandler(__console__)


del logging
del colorlog

__all__ = ["etc_form", "interactive_graph", "nasa_exoplanets_archive", "transit_planner"]
