"""Shared numerical checks for Python input adapters; no pandas dependency."""

from numbers import Integral

import numpy as np


MAX_INDEX = np.iinfo(np.int32).max


def reject(invalid, label, requirement, *, rows=None):
    positions = np.argwhere(invalid)
    if positions.size:
        index = tuple(int(i) for i in positions[0])
        if rows is not None:
            row = rows[index[0]]
            if isinstance(row, np.generic):
                row = row.item()
            location = f"row {row!r}"
        else:
            location = f"index {index[0] if len(index) == 1 else index}"
        raise ValueError(f"{label} must {requirement}; invalid value at {location}")


def numeric(values, label, *, rows=None):
    values = np.asarray(values)
    if values.dtype.kind not in "iuf":
        raise ValueError(f"{label} must have a real numeric dtype")
    with np.errstate(over="ignore", invalid="ignore"):
        values = values.astype(np.float64, copy=False)
    # At least one dimension also lets scalar inputs share error reporting.
    reject(~np.isfinite(np.atleast_1d(values)), label,
           "contain only finite, non-missing values", rows=rows)
    return values


def integers(values, label, minimum=0, *, rows=None):
    values = numeric(values, label, rows=rows)
    checked = np.atleast_1d(values)
    reject(checked != np.floor(checked), label, "contain integer values", rows=rows)
    reject((checked < minimum) | (checked > MAX_INDEX), label,
           f"contain integers in [{minimum}, {MAX_INDEX}]", rows=rows)
    return values.astype(np.int64)


def nonnegative(values, label, *, positive=False, rows=None):
    values = numeric(values, label, rows=rows)
    checked = np.atleast_1d(values)
    reject(checked <= 0 if positive else checked < 0, label,
           "be strictly positive" if positive else "be nonnegative", rows=rows)
    return values


def integer_scalar(value, label, minimum=0, maximum=MAX_INDEX):
    if (isinstance(value, bool) or not isinstance(value, Integral)
            or not minimum <= value <= maximum):
        raise ValueError(f"{label} must be an integer in [{minimum}, {maximum}]")
    return int(value)


def index_base(value):
    if isinstance(value, bool) or not isinstance(value, Integral) or value not in (0, 1):
        raise ValueError("node_index_base must be 0 or 1")
    return int(value)
