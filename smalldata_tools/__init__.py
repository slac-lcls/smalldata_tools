import logging

class LogConfig:
    BasicFormat = '%(message)s'
    Format = '[ %(asctime)s | %(levelname)-3s] %(message)s'
    FullFormat = '[ %(asctime)s | %(name)-13s | %(levelname)-8s] %(message)s'
    Level = logging.INFO

    @staticmethod
    def get_package_name(name):
        return '.'.join(name.split('.')[:-1])

logger = logging.getLogger(__name__)
logger.setLevel(LogConfig.Level)
