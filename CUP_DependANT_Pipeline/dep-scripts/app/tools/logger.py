import logging
from logging.handlers import RotatingFileHandler
import sys
from app.tools.config_parser import config


if config.logging_level != "STDOUT""":
    log_dir = config.this_repo / "logs"
    if not log_dir.is_dir():
        log_dir.mkdir(parents=True, exist_ok=True)


def show_only_error(record):
    return record.levelname == "ERROR" or record.levelname == "CRITICAL"


def create_logger() -> logging.Logger:
    logger = logging.getLogger("pearlai")
    formatter = logging.Formatter("{asctime} - {levelname} - {message}",
                                  style="{",
                                  datefmt="%Y-%m-%d %H:%M")
    if config.logging_level != "STDOUT":
        file_handler = RotatingFileHandler(str(log_dir / "info.log"), mode="a", encoding="utf-8", maxBytes=1000000)
        error_handler = RotatingFileHandler(str(log_dir / "error.log"), mode="a", encoding="utf-8", maxBytes=1000000)
        error_handler.addFilter(show_only_error)

        logger.addHandler(file_handler)
        logger.addHandler(error_handler)
        logger.setLevel(config.logging_level)

        file_handler.setFormatter(formatter)
        error_handler.setFormatter(formatter)
    else:
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)
    return logger


def clinical_trials_logger() -> logging.Logger:
    logger = logging.getLogger("clinical trials")
    formatter = logging.Formatter("{asctime} - {levelname} - {message}",
                                  style="{",
                                  datefmt="%Y-%m-%d %H:%M")
    if config.logging_level != "STDOUT":
        file_handler = RotatingFileHandler(str(log_dir / "clinical_trials_info.log"),
                                           mode="a", encoding="utf-8", maxBytes=1000000)
        error_handler = RotatingFileHandler(str(log_dir / "clinical_trials_error.log"),
                                            mode="a", encoding="utf-8", maxBytes=1000000)
        error_handler.addFilter(show_only_error)

        logger.addHandler(file_handler)
        logger.addHandler(error_handler)
        logger.setLevel(config.logging_level)

        file_handler.setFormatter(formatter)
        error_handler.setFormatter(formatter)
    else:
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)
    return logger
