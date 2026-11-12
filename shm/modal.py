"""Output-only modal identification, damage location and novelty detection."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, signal, stats


def spectral_matrix(acceleration, fs, nperseg=2048):
    """Welch cross-spectral density matrix: frequencies and an array (frequencies, channels, channels)."""
    channels = np.asarray(acceleration, dtype=float).T
    frequencies, spectrum = signal.csd(channels[:, None, :], channels[None, :, :], fs=fs, nperseg=nperseg)
    return frequencies, np.moveaxis(spectrum, -1, 0)


def mac(a, b):
    """Modal assurance criterion between two mode shapes: 1 for the same shape, 0 for orthogonal shapes."""
    return float(abs(np.vdot(a, b)) ** 2 / (np.vdot(a, a).real * np.vdot(b, b).real))


@dataclass
class ModalEstimate:
    frequencies: np.ndarray           # Hz
    damping: np.ndarray               # ratio, NaN where the peak could not be fitted
    shapes: np.ndarray                # one column per mode
    spectrum: tuple                   # spectral lines (Hz) and first singular value


def fit_peak(frequencies, values):
    """Natural frequency and damping ratio from the spectral lines around one resonance peak.

    For one mode under white-noise force the acceleration spectrum S satisfies
    f^4 / S = A [(f^2 - fn^2)^2 + (2 zeta fn f)^2], a quadratic in f^2, so a weighted linear fit is enough.
    Returns None when the fitted curve has no minimum inside the lines.
    """
    y = np.asarray(frequencies, dtype=float) ** 2
    z = y**2 / np.asarray(values, dtype=float)
    centre, half = y.mean(), (y.max() - y.min()) / 2
    # the estimate of each spectral line has a constant relative error, hence the weights 1 / z
    a, b, c = np.polyfit((y - centre) / half, z / z.min(), 2, w=z.min() / z)
    if a <= 0 or abs(b / (2 * a)) > 1:
        return None
    vertex = centre - half * b / (2 * a)
    natural_squared = np.sqrt(vertex**2 + max(half**2 * (c / a - (b / (2 * a)) ** 2), 0.0))
    return float(np.sqrt(natural_squared)), float(np.sqrt(max((1 - vertex / natural_squared) / 2, 0.0)))


def identify_modes(acceleration, fs, reference_frequencies, reference_shapes, nperseg=2048, band=0.12, drop=8.0,
                   mac_limit=0.8):
    """Frequency-domain decomposition (Brincker et al., 2001), tracking each mode near its reference.

    The peak of the first singular value gives the mode shape. The frequency and damping come from a fit to
    the lines around the peak that stay above 1 / drop of its height and keep the same singular vector.
    """
    frequencies, spectrum = spectral_matrix(acceleration, fs, nperseg)
    vectors, values, _ = np.linalg.svd(spectrum)
    first_value, first_vector = values[:, 0], vectors[:, :, 0]
    identified, damping, shapes = [], [], []
    for target, shape in zip(reference_frequencies, np.asarray(reference_shapes).T):
        lines = np.flatnonzero(np.abs(frequencies - target) <= band * target)
        # weighting by the similarity to the reference shape keeps a neighbouring mode from being picked
        similarity = np.abs(first_vector[lines] @ shape) ** 2 / (shape @ shape)
        peak = lines[int(np.argmax(first_value[lines] * similarity))]

        def in_peak(line):
            same_shape = abs(np.vdot(first_vector[line], first_vector[peak])) ** 2 >= mac_limit
            return same_shape and first_value[line] >= first_value[peak] / drop

        low = high = peak
        while low > lines[0] and in_peak(low - 1):
            low -= 1
        while high < lines[-1] and in_peak(high + 1):
            high += 1
        fitted = fit_peak(frequencies[low:high + 1], first_value[low:high + 1]) if high - low >= 3 else None
        identified.append(fitted[0] if fitted else frequencies[peak])
        damping.append(fitted[1] if fitted else np.nan)
        vector = first_vector[peak]
        vector = np.real(vector * np.exp(-1j * np.angle(vector[np.argmax(np.abs(vector))])))
        shapes.append(vector * np.sign(vector @ shape))
    return ModalEstimate(np.array(identified), np.array(damping), np.array(shapes).T, (frequencies, first_value))


def locate_damage(frequencies, baseline_frequencies, energy_fractions):
    """Storey stiffness losses that best explain a change of the squared frequencies.

    To first order the relative change of the squared frequency of a mode is minus the sum over storeys of
    (stiffness loss x share of the modal strain energy in that storey). The losses are found by non-negative
    least squares.
    """
    change = (np.asarray(frequencies) / np.asarray(baseline_frequencies)) ** 2 - 1
    return optimize.nnls(np.asarray(energy_fractions), -change)[0]


class NoveltyDetector:
    """Squared Mahalanobis distance from a healthy population, with a limit for a chosen false-alarm rate."""

    def __init__(self, baseline, false_alarm=0.01):
        baseline = np.asarray(baseline, dtype=float)
        if baseline.ndim == 1:
            baseline = baseline[:, None]
        n, p = baseline.shape
        self.mean = baseline.mean(axis=0)
        self.precision = np.linalg.inv(np.atleast_2d(np.cov(baseline, rowvar=False)))
        # distance of a new healthy observation from n earlier ones follows a scaled F distribution
        self.limit = p * (n - 1) * (n + 1) / (n * (n - p)) * stats.f.ppf(1 - false_alarm, p, n - p)

    def distance(self, features):
        features = np.asarray(features, dtype=float)
        delta = features.reshape(-1, len(self.mean)) - self.mean
        return np.einsum("ij,jk,ik->i", delta, self.precision, delta)

    def alarms(self, features):
        return self.distance(features) > self.limit
