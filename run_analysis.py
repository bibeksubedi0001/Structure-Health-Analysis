"""Checks of the sensor log, and a damage detection study on a simulated five-storey shear frame."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from shm import (LIMITS, NoveltyDetector, ShearFrame, ambient_record, check_log, exceedances, identify_modes,
                 load_log, locate_damage, mac)

# Liberation Serif has the metrics of Times New Roman and stands in where that font is not installed
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"]})

ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "reports" / "visualizations"
RESULTS = ROOT / "results"

STOREYS = 5
FLOOR_MASS = 50_000.0            # kg
STOREY_STIFFNESS = 80e6          # N/m
DAMPING = 0.02                   # ratio, the same in every mode
SAMPLING = 50.0                  # Hz
DURATION = 600.0                 # s per record
NOISE = 0.05                     # measurement noise as a share of the signal RMS
TEMPERATURE_SPREAD = 0.015       # storey stiffness differs by up to this share between records
FALSE_ALARM = 0.01
LOSSES = (0.02, 0.03, 0.05, 0.10, 0.20)


def finish(fig, name):
    fig.tight_layout()
    fig.savefig(FIGURES / name, dpi=150, bbox_inches="tight")
    plt.close(fig)


def percent(value, digits=0):
    return f"{100 * value:.{digits}f} %"


def annotate(ax, image, matrix, fmt):
    """Write each value in its cell, in white where the cell colour is dark."""
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            red, green, blue, _ = image.cmap(image.norm(matrix[i, j]))
            dark = 0.299 * red + 0.587 * green + 0.114 * blue < 0.5
            ax.text(j, i, format(matrix[i, j], fmt).replace("-", "\u2212"), ha="center", va="center", fontsize=9,
                    color="white" if dark else "black")


def analyse_log(path):
    """Check the sensor log, write the readings above their limit and plot the log."""
    log = load_log(path)
    facts = check_log(log)
    flagged = exceedances(log)
    flagged.to_csv(path.parent / "anomalies.csv", index=False, date_format="%Y-%m-%dT%H:%M:%S")
    # the first version's plots applied the limit of each sensor's first reading to all of its readings
    first_type_limit = log.groupby("sensor_id")["sensor_type"].transform("first").map(LIMITS)
    facts["flagged_by_first_version_plots"] = int((log["value"] > first_type_limit).sum())
    facts["flagged"] = int(len(flagged))

    sensors = sorted(log["sensor_id"].unique())
    fig, axes = plt.subplots(len(LIMITS), 1, figsize=(10, 8), sharex=True)
    for ax, (sensor_type, limit) in zip(axes, LIMITS.items()):
        part = log[log["sensor_type"] == sensor_type]
        for number, sensor in enumerate(sensors):
            rows = part[part["sensor_id"] == sensor]
            ax.scatter(rows["timestamp"], rows["value"], s=28, color=f"C{number}", label=sensor)
        above = part[part["value"] > limit]
        ax.scatter(above["timestamp"], above["value"], s=90, facecolors="none", edgecolors="tab:red", linewidths=1.2)
        ax.axhline(limit, color="tab:red", linestyle="--", linewidth=1)
        ax.set_ylabel("value")
        ax.set_title(f"{sensor_type}: {len(above)} of {len(part)} readings above the limit of {limit:g}", fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].legend(ncol=len(sensors), fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.45))
    axes[-1].set_xlabel("time")
    finish(fig, "log_readings.png")

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    for ax, column, title in zip(axes, ("sensor_type", "location"), ("sensor type", "location")):
        table = pd.crosstab(log["sensor_id"], log[column])
        image = ax.imshow(table.to_numpy(), cmap="Blues", vmin=0)
        annotate(ax, image, table.to_numpy(), "d")
        ax.set_xticks(range(table.shape[1]), table.columns)
        ax.set_yticks(range(table.shape[0]), table.index)
        ax.set_title(f"Readings by sensor and {title}", fontsize=10)
    finish(fig, "log_metadata.png")
    return facts


def frame_study(rng, n_baseline, n_healthy, n_damaged):
    """Identify the modes of the frame from simulated records and test damage detection and location."""
    frame = ShearFrame([FLOOR_MASS] * STOREYS, [STOREY_STIFFNESS] * STOREYS, DAMPING)
    model_frequencies, model_shapes = frame.modes()
    floors = np.arange(1, STOREYS + 1)

    def identify(structure, temperature=True):
        factor = 1 + TEMPERATURE_SPREAD * rng.uniform(-1, 1) if temperature else 1.0
        record = ambient_record(structure.scaled(factor), DURATION, SAMPLING, rng, NOISE)
        return identify_modes(record, SAMPLING, model_frequencies, model_shapes)

    def frequencies_of(structure, count, temperature=True):
        return np.array([identify(structure, temperature).frequencies for _ in range(count)])

    single = identify(frame, temperature=False)
    shape_mac = [mac(single.shapes[:, mode], model_shapes[:, mode]) for mode in range(STOREYS)]
    constant = frequencies_of(frame, n_baseline, temperature=False)
    baseline = frequencies_of(frame, n_baseline)
    healthy = frequencies_of(frame, n_healthy)
    all_modes = NoveltyDetector(baseline, FALSE_ALARM)
    first_mode = NoveltyDetector(baseline[:, 0], FALSE_ALARM)
    reference = baseline.mean(axis=0)
    fractions = frame.strain_energy_fractions()

    shape = (len(LOSSES), STOREYS)
    detected, detected_first, located = np.zeros(shape), np.zeros(shape), np.zeros(shape)
    profile = np.zeros(shape + (STOREYS,))
    index = {}
    for i, loss in enumerate(LOSSES):
        for j, storey in enumerate(floors):
            frequencies = frequencies_of(frame.with_damage(storey, loss), n_damaged)
            index[(loss, storey)] = all_modes.distance(frequencies)
            detected[i, j] = all_modes.alarms(frequencies).mean()
            detected_first[i, j] = first_mode.alarms(frequencies[:, 0]).mean()
            estimate = np.array([locate_damage(row, reference, fractions) for row in frequencies])
            located[i, j] = (estimate.argmax(axis=1) == j).mean()
            # averaging the records first removes most of the temperature effect and the identification scatter
            profile[i, j] = locate_damage(frequencies.mean(axis=0), reference, fractions)
    exact_change = np.array([frame.with_damage(storey, 0.10).modes()[0] / model_frequencies - 1 for storey in floors])

    # figures
    lines, first_value = single.spectrum
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.semilogy(lines, first_value, color="tab:blue", linewidth=0.9, label="first singular value of the spectral matrix")
    for mode, (model, found) in enumerate(zip(model_frequencies, single.frequencies)):
        ax.axvline(model, color="0.5", linestyle="--", linewidth=0.8, label="model frequency" if mode == 0 else None)
        ax.text(found, 1.01, f"{found:.3f} Hz", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8)
    ax.set(xlim=(0, 15), xlabel="frequency (Hz)", ylabel="(m/s2)2 / Hz")
    ax.set_title("Frequency-domain decomposition of one 10-minute record (labels: identified frequencies)", pad=18)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="lower right")
    finish(fig, "fdd_spectrum.png")

    fig, axes = plt.subplots(1, STOREYS, figsize=(11, 4.2), sharey=True)
    for mode, ax in enumerate(axes):
        model = np.append(0.0, model_shapes[:, mode] / np.abs(model_shapes[:, mode]).max())
        found = single.shapes[:, mode] / np.abs(single.shapes[:, mode]).max()
        ax.plot(model, np.arange(STOREYS + 1), color="0.4", label="model")
        ax.plot(found, floors, "o", color="tab:red", label="identified")
        ax.axvline(0, color="0.8", linewidth=0.8)
        ax.set_title(f"mode {mode + 1}\n{single.frequencies[mode]:.2f} Hz, MAC {shape_mac[mode]:.4f}", fontsize=9)
        ax.set_xlim(-1.15, 1.15)
    axes[0].set_ylabel("floor")
    axes[0].legend(fontsize=8, loc="lower right")
    finish(fig, "mode_shapes.png")

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    image = ax.imshow(100 * exact_change, cmap="Blues_r", vmax=0)
    annotate(ax, image, 100 * exact_change, ".2f")
    ax.set_xticks(range(STOREYS), [f"mode {mode}" for mode in floors])
    ax.set_yticks(range(STOREYS), [f"storey {storey}" for storey in floors])
    ax.set_title("Change of each natural frequency (%)\nfor a 10 % stiffness loss in one storey", fontsize=10)
    fig.colorbar(image, ax=ax, fraction=0.046)
    finish(fig, "frequency_changes.png")

    groups = [("baseline", all_modes.distance(baseline)), ("healthy", all_modes.distance(healthy))]
    groups += [(f"{percent(loss)} loss", index[(loss, 1)]) for loss in LOSSES]
    fig, ax = plt.subplots(figsize=(11, 4.6))
    start = 0
    for label, values in groups:
        colour = "0.45" if label == "baseline" else "tab:blue" if label == "healthy" else "tab:orange"
        ax.semilogy(np.arange(start, start + len(values)), values, ".", color=colour, markersize=4)
        ax.text(start + len(values) / 2, 1.02, label, ha="center", va="bottom", fontsize=8,
                transform=ax.get_xaxis_transform())
        start += len(values)
        ax.axvline(start - 0.5, color="0.85", linewidth=0.8)
    ax.axhline(all_modes.limit, color="tab:red", linestyle="--", label=f"limit for a {percent(FALSE_ALARM)} false-alarm rate")
    ax.set(xlabel="record", ylabel="squared Mahalanobis distance", xlim=(0, start))
    ax.set_title("Novelty index from five frequencies; damage in storey 1", fontsize=10, pad=18)
    ax.legend(fontsize=8, loc="lower right")
    finish(fig, "novelty_index.png")

    labels = [percent(loss) for loss in LOSSES]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, rates, title in zip(axes, (detected, detected_first),
                                ("Five frequencies, Mahalanobis distance", "First frequency only")):
        for j, storey in enumerate(floors):
            ax.plot(range(len(LOSSES)), 100 * rates[:, j], "o-", label=f"storey {storey}")
        ax.set_xticks(range(len(LOSSES)), labels)
        ax.set(xlabel="stiffness loss in one storey", title=title)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("records flagged (%)")
    axes[0].legend(fontsize=8)
    finish(fig, "detection_rate.png")

    ten = LOSSES.index(0.10)
    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.2, 1]})
    for j, storey in enumerate(floors):
        left.plot(range(len(LOSSES)), 100 * located[:, j], "o-", label=f"storey {storey}")
    left.set_xticks(range(len(LOSSES)), labels)
    left.set(xlabel="stiffness loss in one storey", ylabel="single records with the largest\nestimate in the right storey (%)",
             ylim=(0, 105), title="Damage location from one record")
    left.grid(alpha=0.3)
    left.legend(fontsize=8)
    image = right.imshow(100 * profile[ten], cmap="Oranges", vmin=0)
    annotate(right, image, 100 * profile[ten], ".1f")
    right.set_xticks(range(STOREYS), floors)
    right.set_yticks(range(STOREYS), floors)
    right.set(xlabel="storey of the estimate", ylabel="damaged storey",
              title=f"Estimated loss (%) from the mean of {n_damaged} records\ntrue loss 10 %")
    fig.colorbar(image, ax=right, fraction=0.046)
    finish(fig, "damage_location.png")

    def table(matrix):
        return {percent(loss): [float(value) for value in row] for loss, row in zip(LOSSES, matrix)}

    return {
        "frame": {"storeys": STOREYS, "floor_mass_kg": FLOOR_MASS, "storey_stiffness_N_per_m": STOREY_STIFFNESS,
                  "damping_ratio": DAMPING},
        "records": {"duration_s": DURATION, "sampling_Hz": SAMPLING, "noise_ratio": NOISE,
                    "temperature_stiffness_spread": TEMPERATURE_SPREAD, "baseline": n_baseline, "healthy": n_healthy,
                    "per_damage_case": n_damaged},
        "model_frequencies_Hz": model_frequencies.tolist(),
        "one_record": {"frequencies_Hz": single.frequencies.tolist(), "damping_ratio": single.damping.tolist(),
                       "error_percent": (100 * (single.frequencies / model_frequencies - 1)).tolist(), "mac": shape_mac},
        "scatter_percent": {"constant_temperature": (100 * constant.std(axis=0, ddof=1) / constant.mean(axis=0)).tolist(),
                            "varying_temperature": (100 * baseline.std(axis=0, ddof=1) / baseline.mean(axis=0)).tolist()},
        "limit": {"five_frequencies": float(all_modes.limit), "first_frequency": float(first_mode.limit)},
        "false_alarm_rate": {"nominal": FALSE_ALARM, "five_frequencies": float(all_modes.alarms(healthy).mean()),
                             "first_frequency": float(first_mode.alarms(healthy[:, 0]).mean())},
        "detection_rate": {"five_frequencies": table(detected), "first_frequency": table(detected_first)},
        "location_rate": table(located),
        "estimated_loss_in_damaged_storey_from_mean_frequencies": table(
            np.array([[profile[i, j, j] for j in range(STOREYS)] for i in range(len(LOSSES))])),
        "largest_estimated_loss_in_another_storey_from_mean_frequencies": table(
            np.array([[np.delete(profile[i, j], j).max() for j in range(STOREYS)] for i in range(len(LOSSES))])),
        "frequency_change_percent_for_10_percent_loss": (100 * exact_change).tolist(),
    }


def report(command, log, study):
    def percent_row(label, values, digits=0):
        return f"| {label} | " + " | ".join(percent(value, digits) for value in values) + " |"

    header = "| Stiffness loss | " + " | ".join(f"Storey {storey}" for storey in range(1, STOREYS + 1)) + " |"
    rule = "|---" * (STOREYS + 1) + "|"
    one = study["one_record"]
    lines = ["# Results", "", f"Command: `{command}`", "", "## Sensor log", "", "| | |", "|---|---|",
             f"| Readings | {log['readings']} from {log['sensors']} sensor ids, one every {log['interval_s'][0]:g} s |",
             f"| Period | {log['start']} to {log['end']} |",
             f"| Sensor ids that report more than one sensor type | {len(log['sensors_reporting_several_types'])} of {log['sensors']} |",
             f"| Sensor ids that report more than one location | {len(log['sensors_reporting_several_locations'])} of {log['sensors']} |",
             f"| Readings above the limit of their own type | {log['flagged']} |",
             f"| Readings flagged in the plots of the first version | {log['flagged_by_first_version_plots']} |",
             f"| Largest robust score of a reading | {log['largest_robust_score']:.2f} (3.5 is the usual outlier limit) |",
             "", "| Sensor type | Readings | Minimum | Maximum | Readings equal to the maximum | Limit | Above the limit |",
             "|---|---|---|---|---|---|---|"]
    lines += [f"| {name} | {row['readings']} | {row['minimum']:g} | {row['maximum']:g} | {row['readings_at_maximum']} | "
              f"{row['limit']:g} | {row['above_limit']} |" for name, row in log["by_type"].items()]
    records = study["records"]
    lines += ["", "## Shear-frame study", "",
              f"Five storeys, {FLOOR_MASS / 1000:g} t per floor, {STOREY_STIFFNESS / 1e6:g} MN/m per storey, "
              f"{percent(DAMPING)} damping. Records of {DURATION / 60:g} minutes at {SAMPLING:g} Hz with {percent(NOISE)} "
              f"measurement noise; the stiffness of the whole frame differs by up to {percent(TEMPERATURE_SPREAD, 1)} "
              f"between records.", "",
              "### One record against the model", "",
              "| Mode | Model frequency (Hz) | Identified (Hz) | Error | Identified damping | MAC of the mode shape |",
              "|---|---|---|---|---|---|"]
    lines += [f"| {mode + 1} | {study['model_frequencies_Hz'][mode]:.3f} | {one['frequencies_Hz'][mode]:.3f} | "
              f"{one['error_percent'][mode]:+.2f} % | {100 * one['damping_ratio'][mode]:.2f} % | {one['mac'][mode]:.4f} |"
              for mode in range(STOREYS)]
    lines += ["", f"### Scatter of the identified frequencies over {records['baseline']} healthy records", "",
              "| | " + " | ".join(f"Mode {mode}" for mode in range(1, STOREYS + 1)) + " |", rule,
              "| Constant temperature | " + " | ".join(f"{value:.2f} %" for value in study["scatter_percent"]["constant_temperature"]) + " |",
              "| Varying temperature | " + " | ".join(f"{value:.2f} %" for value in study["scatter_percent"]["varying_temperature"]) + " |",
              "", f"### False alarms on {records['healthy']} further healthy records", "", "| Detector | Limit | False alarms |",
              "|---|---|---|",
              f"| Five frequencies | {study['limit']['five_frequencies']:.1f} | {percent(study['false_alarm_rate']['five_frequencies'], 1)} |",
              f"| First frequency | {study['limit']['first_frequency']:.1f} | {percent(study['false_alarm_rate']['first_frequency'], 1)} |",
              "", f"### Records flagged, five frequencies ({records['per_damage_case']} records per case)", "", header, rule]
    lines += [percent_row(loss, row) for loss, row in study["detection_rate"]["five_frequencies"].items()]
    lines += ["", "### Records flagged, first frequency only", "", header, rule]
    lines += [percent_row(loss, row) for loss, row in study["detection_rate"]["first_frequency"].items()]
    lines += ["", "### Single records in which the largest estimated loss is in the damaged storey", "", header, rule]
    lines += [percent_row(loss, row) for loss, row in study["location_rate"].items()]
    lines += ["", f"### Estimated loss in the damaged storey, from the mean frequencies of {records['per_damage_case']} records",
              "", header, rule]
    lines += [percent_row(loss, row, 1)
              for loss, row in study["estimated_loss_in_damaged_storey_from_mean_frequencies"].items()]
    lines += ["", "### Largest estimated loss in any other storey, same records", "", header, rule]
    lines += [percent_row(loss, row, 1)
              for loss, row in study["largest_estimated_loss_in_another_storey_from_mean_frequencies"].items()]
    lines += ["", "### Exact change of the model frequencies for a 10 % loss", "",
              "| Damaged storey | " + " | ".join(f"Mode {mode}" for mode in range(1, STOREYS + 1)) + " |", rule]
    lines += [f"| {storey + 1} | " + " | ".join(f"{value:.2f} %" for value in row) + " |"
              for storey, row in enumerate(study["frequency_change_percent_for_10_percent_loss"])]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, default=ROOT / "backend" / "data.csv")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--baseline", type=int, default=60, help="healthy records used to set the limits")
    parser.add_argument("--healthy", type=int, default=200, help="further healthy records for the false-alarm rate")
    parser.add_argument("--damaged", type=int, default=40, help="records per damage case")
    args = parser.parse_args()
    FIGURES.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(exist_ok=True)

    log = analyse_log(args.log)
    study = frame_study(np.random.default_rng(args.seed), args.baseline, args.healthy, args.damaged)
    command = "python run_analysis.py" + "".join(
        f" --{name} {getattr(args, name)}" for name in ("seed", "baseline", "healthy", "damaged")
        if getattr(args, name) != parser.get_default(name))
    text = report(command, log, study)
    (RESULTS / "summary.md").write_text(text, encoding="utf-8")
    (RESULTS / "summary.json").write_text(json.dumps({"command": command, "sensor_log": log, "frame_study": study},
                                                     indent=2), encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
