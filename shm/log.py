"""Checks of a sensor log: metadata consistency, saturation, limit exceedances and robust scores."""
from __future__ import annotations

import numpy as np
import pandas as pd

COLUMNS = ["sensor_id", "value", "timestamp", "sensor_type", "location"]
LIMITS = {"strain": 50.0, "temperature": 100.0, "vibration": 30.0}


def load_log(path):
    """Read the sensor log and parse its timestamps."""
    log = pd.read_csv(path)
    missing = [column for column in COLUMNS if column not in log.columns]
    if missing:
        raise ValueError(f"missing columns: {missing}")
    log["timestamp"] = pd.to_datetime(log["timestamp"])
    return log


def exceedances(log, limits=LIMITS):
    """Readings above the limit of their own sensor type, listed type by type in the order of `limits`."""
    unknown = sorted(set(log["sensor_type"]) - set(limits))
    if unknown:
        raise ValueError(f"no limit defined for sensor types: {unknown}")
    flagged = log[log["value"] > log["sensor_type"].map(limits)]
    order = flagged["sensor_type"].map({name: position for position, name in enumerate(limits)})
    return flagged.iloc[np.argsort(order.to_numpy(), kind="stable")]


def robust_scores(log):
    """Modified z-score of each reading within its sensor type (Iglewicz and Hoaglin, 1993)."""
    def score(values):
        deviation = values - values.median()
        mad = deviation.abs().median()
        if mad > 0:
            return 0.6745 * deviation / mad
        # more than half of the values are equal, so fall back to the mean absolute deviation
        return deviation / (1.253314 * deviation.abs().mean())

    return log.groupby("sensor_type")["value"].transform(score)


def check_log(log, limits=LIMITS):
    """Facts needed before the readings can be interpreted."""
    types = log.groupby("sensor_id")["sensor_type"].nunique()
    places = log.groupby("sensor_id")["location"].nunique()
    by_type = {}
    for sensor_type, values in log.groupby("sensor_type")["value"]:
        top = values.max()
        by_type[sensor_type] = {
            "readings": int(len(values)), "minimum": float(values.min()), "maximum": float(top),
            "readings_at_maximum": int((values == top).sum()),
            "above_limit": int((values > limits[sensor_type]).sum()), "limit": float(limits[sensor_type]),
        }
    steps = log["timestamp"].sort_values().diff().dt.total_seconds().dropna()
    return {
        "readings": int(len(log)),
        "sensors": int(log["sensor_id"].nunique()),
        "start": log["timestamp"].min().isoformat(),
        "end": log["timestamp"].max().isoformat(),
        "interval_s": sorted(float(step) for step in steps.unique()),
        "sensors_reporting_several_types": sorted(types[types > 1].index),
        "sensors_reporting_several_locations": sorted(places[places > 1].index),
        "by_type": by_type,
        "largest_robust_score": float(np.abs(robust_scores(log)).max()),
    }
