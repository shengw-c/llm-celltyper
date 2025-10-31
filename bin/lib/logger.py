#!/usr/bin/env python3
"""
Logging configuration module for the cell type annotation pipeline.
Provides centralized logging setup with file and console handlers.
"""

import logging
import os
import sys
from pathlib import Path
from typing import Optional
from datetime import datetime


class PipelineLogger:
    """Manages logging configuration for the annotation pipeline."""
    
    _loggers = {}
    
    @classmethod
    def get_logger(
        cls,
        name: str,
        log_dir: str = "logs",
        level: int = logging.INFO,
        log_file: Optional[str] = None,
        cell_type: Optional[str] = None
    ) -> logging.Logger:
        """
        Get or create a logger with file and console handlers.
        
        Args:
            name: Logger name (typically __name__ of the module)
            log_dir: Directory to store log files
            level: Logging level (default: INFO)
            log_file: Custom log file name (optional, will auto-generate if None)
            cell_type: Cell type being processed (will be included in log filename)
        
        Returns:
            Configured logger instance
        """
        # Return existing logger if already configured
        if name in cls._loggers:
            return cls._loggers[name]
        
        # Create logger
        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.propagate = False
        
        # Remove existing handlers to avoid duplicates
        logger.handlers.clear()
        
        # Create log directory if it doesn't exist
        os.makedirs(log_dir, exist_ok=True)
        
        # Generate log file name if not provided
        if log_file is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            module_name = name.split('.')[-1]
            
            # Include cell_type in filename if provided
            if cell_type:
                # Replace spaces and special chars with underscores
                cell_type_safe = cell_type.replace(' ', '_').replace('/', '_')
                log_file = f"{cell_type_safe}_{timestamp}.log"
            else:
                log_file = f"{module_name}_{timestamp}.log"
        
        log_path = os.path.join(log_dir, log_file)
        
        # Create formatters
        file_formatter = logging.Formatter(
            '%(asctime)s | %(name)s | %(levelname)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        console_formatter = logging.Formatter(
            '%(asctime)s | %(levelname)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # File handler - detailed logging
        file_handler = logging.FileHandler(log_path, mode='a')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
        
        # Console handler - cleaner output
        # Use stderr instead of stdout to avoid polluting stdout in scripts
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(level)
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
        
        # Store logger
        cls._loggers[name] = logger
        
        logger.info(f"Logger initialized. Log file: {log_path}")
        
        return logger
    
    @classmethod
    def set_level(cls, name: str, level: int):
        """Set logging level for a specific logger."""
        if name in cls._loggers:
            cls._loggers[name].setLevel(level)
    
    @classmethod
    def get_all_log_files(cls, log_dir: str = "logs") -> list:
        """Get list of all log files in the log directory."""
        if not os.path.exists(log_dir):
            return []
        return [os.path.join(log_dir, f) for f in os.listdir(log_dir) if f.endswith('.log')]
