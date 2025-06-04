import logging


class LogConfig:
    BasicFormat = "%(message)s"
    Format = "[ %(asctime)s | %(levelname)-3s] %(message)s"
    FullFormat = "[ %(asctime)s | %(levelname)-3s | %(filename)s] %(message)s"
    Level = logging.INFO
    # Level = logging.DEBUG
    # Level = logging.WARNING

    @staticmethod
    def get_package_name(name):
        return ".".join(name.split(".")[:-1])


logger = logging.getLogger(__name__)
logger.setLevel(LogConfig.Level)
handler = logging.StreamHandler()
formatter = logging.Formatter(LogConfig.FullFormat)
handler.setFormatter(formatter)
logger.addHandler(handler)
