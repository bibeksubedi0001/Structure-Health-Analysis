"""Shear-frame model of a multi-storey building: modes, sensitivities and ambient response."""
from __future__ import annotations

import numpy as np
from scipy import linalg, signal


class ShearFrame:
    """Lumped floor masses on storey springs, fixed at the base, with the same damping ratio in every mode."""

    def __init__(self, masses, stiffnesses, damping_ratio=0.02):
        self.masses = np.asarray(masses, dtype=float)
        self.stiffnesses = np.asarray(stiffnesses, dtype=float)
        self.damping_ratio = float(damping_ratio)
        if self.masses.shape != self.stiffnesses.shape or self.masses.ndim != 1:
            raise ValueError("one mass and one storey stiffness are needed per floor")

    @property
    def n(self):
        return len(self.masses)

    def stiffness_matrix(self):
        k = self.stiffnesses
        return np.diag(k + np.append(k[1:], 0.0)) - np.diag(k[1:], 1) - np.diag(k[1:], -1)

    def modes(self):
        """Natural frequencies in Hz, lowest first, and mass-normalised mode shapes in columns."""
        eigenvalues, shapes = linalg.eigh(self.stiffness_matrix(), np.diag(self.masses))
        return np.sqrt(eigenvalues) / (2 * np.pi), shapes * np.sign(shapes[-1])

    def with_damage(self, storey, loss):
        """Frame with the stiffness of one storey (1 is the lowest) reduced by the fraction `loss`."""
        stiffnesses = self.stiffnesses.copy()
        stiffnesses[storey - 1] *= 1 - loss
        return ShearFrame(self.masses, stiffnesses, self.damping_ratio)

    def scaled(self, factor):
        """Frame with every storey stiffness multiplied by `factor`, as for a change of temperature."""
        return ShearFrame(self.masses, self.stiffnesses * factor, self.damping_ratio)

    def strain_energy_fractions(self):
        """Share of the strain energy of each mode (rows) that is stored in each storey (columns)."""
        shapes = self.modes()[1]
        drift = np.diff(np.vstack([np.zeros(self.n), shapes]), axis=0)
        energy = self.stiffnesses[:, None] * drift**2
        return (energy / energy.sum(axis=0)).T

    def response(self, force, dt):
        """Floor displacements and accelerations for forces held constant over each time step.

        Each mode is a single oscillator that is advanced with its exact discrete-time solution, so the
        result has no time-stepping error for this loading.
        """
        frequencies, shapes = self.modes()
        modal_force = np.asarray(force, dtype=float) @ shapes
        modal_displacement = np.empty_like(modal_force)
        modal_acceleration = np.empty_like(modal_force)
        for mode, frequency in enumerate(frequencies):
            omega, zeta = 2 * np.pi * frequency, self.damping_ratio
            system = (np.array([[0.0, 1.0], [-omega**2, -2 * zeta * omega]]), np.array([[0.0], [1.0]]),
                      np.eye(2), np.zeros((2, 1)))
            a, b, c, d, _ = signal.cont2discrete(system, dt, method="zoh")
            numerator, denominator = signal.ss2tf(a, b, c, d)
            displacement = signal.lfilter(numerator[0], denominator, modal_force[:, mode])
            velocity = signal.lfilter(numerator[1], denominator, modal_force[:, mode])
            modal_displacement[:, mode] = displacement
            modal_acceleration[:, mode] = modal_force[:, mode] - 2 * zeta * omega * velocity - omega**2 * displacement
        return modal_displacement @ shapes.T, modal_acceleration @ shapes.T


def ambient_record(frame, duration, fs, rng, noise_ratio=0.05, force_rms=100.0, settle=30.0):
    """Floor accelerations under independent white-noise floor forces, with measurement noise added."""
    steps, skip = int(round(duration * fs)), int(round(settle * fs))
    force = rng.normal(0.0, force_rms, (steps + skip, frame.n))
    acceleration = frame.response(force, 1 / fs)[1][skip:]
    noise = rng.normal(0.0, 1.0, acceleration.shape) * noise_ratio * acceleration.std(axis=0)
    return acceleration + noise
