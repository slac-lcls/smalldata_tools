"""
SharedAzavCache
===============

Draft helper for storing azimuthalBinning precomputed arrays in MPI shared memory.

This is a sketch intended to mirror the SharedGeoCache interface used in psana.
Integration points should call into this cache to get or build shared arrays.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import logging
import re
from typing import Any, Dict, Iterable, Optional, Tuple

import numpy as np


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SharedAzavKey:
    """Immutable key identifying one azimuthal-binning cache group."""
    det_name: str
    azav_id: str


class SharedAzavCache:
    """Thin wrapper around a MPISharedMemory-like object for azimuthalBinning."""

    def __init__(self, shared_mem: Any = None, logger_obj: Optional[logging.Logger] = None):
        self.shared_mem = shared_mem
        self.logger = logger_obj or logger
        self._local_meta: Dict[str, Dict[str, Any]] = {}

    @property
    def enabled(self) -> bool:
        return self.shared_mem is not None

    @staticmethod
    def _safe_name(name: str) -> str:
        return re.sub(r"[^0-9A-Za-z_]+", "_", name)

    @staticmethod
    def _stable_json(obj: Any) -> str:
        return json.dumps(obj, sort_keys=True, separators=(",", ":"))

    @classmethod
    def build_azav_id(
        cls,
        geom_id: Optional[str],
        params: Optional[Dict[str, Any]] = None,
    ) -> str:
        """Build a short hash for azimuthal-binning cache invalidation."""
        payload = {
            "geom_id": geom_id or "",
            "params": params or {},
        }
        digest = hashlib.sha1(cls._stable_json(payload).encode("utf-8")).hexdigest()
        return digest[:12]

    @classmethod
    def make_key(cls, det_name: str, azav_id: str) -> SharedAzavKey:
        return SharedAzavKey(
            det_name=cls._safe_name(det_name),
            azav_id=cls._safe_name(azav_id),
        )

    @staticmethod
    def _prefix(key: SharedAzavKey) -> str:
        return f"azav_{key.det_name}_{key.azav_id}"

    def record_meta(self, key: SharedAzavKey, meta: Dict[str, Any]) -> None:
        """Record small metadata locally (not shared)."""
        self._local_meta[self._prefix(key)] = meta

    def get_meta(self, key: SharedAzavKey) -> Optional[Dict[str, Any]]:
        return self._local_meta.get(self._prefix(key))

    def get_or_allocate(
        self,
        key: SharedAzavKey,
        name: str,
        shape: Iterable[int],
        dtype: Any,
        zero_init: bool = False,
    ) -> Tuple[np.ndarray, bool]:
        """Return a shared array, allocating if needed.

        Returns:
            (array, created)
        """
        if not self.enabled:
            raise RuntimeError("SharedAzavCache is not enabled (no shared_mem).")

        dtype = np.dtype(dtype)
        shape_tuple = tuple(int(dim) for dim in shape)
        full_name = f"{self._prefix(key)}_{name}"

        if self.shared_mem.has_array(full_name):
            handle = (
                self.shared_mem.get_handle(full_name)
                if hasattr(self.shared_mem, "get_handle")
                else None
            )
            if handle and (handle.shape != shape_tuple or handle.dtype != dtype):
                raise ValueError(
                    f"Shared array {full_name} shape/dtype mismatch: "
                    f"{handle.shape}/{handle.dtype} vs {shape_tuple}/{dtype}"
                )
            return self.shared_mem.get_array(full_name), False

        array = self.shared_mem.allocate_array(
            full_name, shape_tuple, dtype, zero_init=zero_init
        )
        return array, True

    def get_if_present(self, key: SharedAzavKey, name: str) -> Optional[np.ndarray]:
        """Return a shared array if it exists, otherwise None."""
        if not self.enabled:
            return None
        full_name = f"{self._prefix(key)}_{name}"
        if not self.shared_mem.has_array(full_name):
            return None
        return self.shared_mem.get_array(full_name)

    def barrier(self) -> None:
        if self.enabled and hasattr(self.shared_mem, "barrier"):
            self.shared_mem.barrier()
