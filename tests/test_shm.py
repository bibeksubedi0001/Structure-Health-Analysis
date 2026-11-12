"""Checks of the frame model, the identification and the detectors against theory and synthetic data."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import signal, stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shm import (NoveltyDetector, ShearFrame, ambient_record, check_log, exceedances, fit_peak, identify_modes,
                 load_log, locate_damage, mac, robust_scores, spectral_matrix)


def uniform_frame():
    return ShearFrame([5e4] * 5, [8e7] * 5, 0.02)


class FrameModel(unittest.TestCase):
    def test_uniform_frame_matches_the_closed_form_frequencies(self):
        frequencies, _ = uniform_frame().modes()
        mode = np.arange(1, 6)
        exact = 2 * np.sqrt(8e7 / 5e4) * np.sin((2 * mode - 1) * np.pi / 22) / (2 * np.pi)
        np.testing.assert_allclose(frequencies, exact, rtol=1e-12)

    def test_mode_shapes_are_mass_normalised(self):
        frame = ShearFrame([2.0, 3.0, 1.5], [400.0, 300.0, 500.0])
        frequencies, shapes = frame.modes()
        np.testing.assert_allclose(shapes.T @ np.diag(frame.masses) @ shapes, np.eye(3), atol=1e-12)
        np.testing.assert_allclose(shapes.T @ frame.stiffness_matrix() @ shapes,
                                   np.diag((2 * np.pi * frequencies) ** 2), atol=1e-9)

    def test_damage_and_scaling_change_the_right_stiffness(self):
        frame = uniform_frame()
        np.testing.assert_allclose(frame.with_damage(2, 0.25).stiffnesses / 8e7, [1, 0.75, 1, 1, 1])
        np.testing.assert_allclose(frame.scaled(0.81).modes()[0], 0.9 * frame.modes()[0], rtol=1e-12)
        np.testing.assert_allclose(frame.stiffnesses, 8e7)

    def test_strain_energy_fractions_predict_small_frequency_changes(self):
        frame = uniform_frame()
        fractions = frame.strain_energy_fractions()
        np.testing.assert_allclose(fractions.sum(axis=1), 1.0, rtol=1e-12)
        for storey in range(1, 6):
            change = frame.with_damage(storey, 0.001).modes()[0] / frame.modes()[0] - 1
            np.testing.assert_allclose(change, -0.0005 * fractions[:, storey - 1], atol=5e-7)

    def test_response_matches_a_direct_state_space_solution(self):
        frame = ShearFrame([2.0, 3.0, 1.5], [400.0, 300.0, 500.0], 0.03)
        dt, steps = 0.01, 400
        force = np.random.default_rng(0).normal(size=(steps, 3))
        displacement, acceleration = frame.response(force, dt)
        mass, stiffness = np.diag(frame.masses), frame.stiffness_matrix()
        frequencies, shapes = frame.modes()
        damping = mass @ shapes @ np.diag(2 * 0.03 * 2 * np.pi * frequencies) @ shapes.T @ mass
        a = np.block([[np.zeros((3, 3)), np.eye(3)],
                      [-np.linalg.solve(mass, stiffness), -np.linalg.solve(mass, damping)]])
        b = np.vstack([np.zeros((3, 3)), np.linalg.inv(mass)])
        ad, bd, *_ = signal.cont2discrete((a, b, np.eye(6), np.zeros((6, 3))), dt)
        state = np.zeros(6)
        for step in range(steps):
            np.testing.assert_allclose(displacement[step], state[:3], atol=1e-10)
            expected = np.linalg.solve(mass, force[step] - damping @ state[3:] - stiffness @ state[:3])
            np.testing.assert_allclose(acceleration[step], expected, atol=1e-9)
            state = ad @ state + bd @ force[step]

    def test_constant_force_settles_at_the_static_displacement(self):
        frame = ShearFrame([2.0, 3.0, 1.5], [400.0, 300.0, 500.0], 0.1)
        load = np.array([1.0, -2.0, 3.0])
        displacement, acceleration = frame.response(np.tile(load, (8000, 1)), 0.01)
        np.testing.assert_allclose(displacement[-1], np.linalg.solve(frame.stiffness_matrix(), load), rtol=1e-8)
        np.testing.assert_allclose(acceleration[-1], 0.0, atol=1e-8)

    def test_free_vibration_has_the_damped_period_and_decay(self):
        zeta, damped_frequency, dt = 0.05, 2.0, 0.005
        natural = 2 * np.pi * damped_frequency / np.sqrt(1 - zeta**2)
        frame = ShearFrame([1000.0], [1000.0 * natural**2], zeta)
        force = np.zeros((2000, 1))
        force[0] = 1e4
        displacement = frame.response(force, dt)[0][:, 0]
        period = int(round(1 / (damped_frequency * dt)))
        ratio = np.exp(-2 * np.pi * zeta / np.sqrt(1 - zeta**2))
        np.testing.assert_allclose(displacement[1 + period:1500 + period], ratio * displacement[1:1500],
                                   atol=1e-9 * np.abs(displacement).max())

    def test_white_noise_response_has_the_variance_of_random_vibration_theory(self):
        mass, frequency, zeta, dt, force_rms = 1000.0, 1.0, 0.05, 0.02, 50.0
        omega = 2 * np.pi * frequency
        frame = ShearFrame([mass], [mass * omega**2], zeta)
        force = np.random.default_rng(1).normal(0.0, force_rms, (1_000_000, 1))
        displacement = frame.response(force, dt)[0][5000:, 0]
        # sigma^2 = pi S0 / (c k) with S0 = force_rms^2 dt / (2 pi) for a force held constant over dt
        expected = force_rms**2 * dt / (4 * zeta * omega**3 * mass**2)
        self.assertAlmostEqual(displacement.var() / expected, 1.0, delta=0.06)

    def test_ambient_record_has_the_requested_length_and_noise(self):
        frame = uniform_frame()
        clean = ambient_record(frame, 100.0, 50.0, np.random.default_rng(2), noise_ratio=0.0)
        noisy = ambient_record(frame, 100.0, 50.0, np.random.default_rng(2), noise_ratio=0.2)
        self.assertEqual(clean.shape, (5000, 5))
        np.testing.assert_allclose((noisy - clean).std(axis=0) / clean.std(axis=0), 0.2, rtol=0.05)


class Identification(unittest.TestCase):
    def test_mac(self):
        shape = np.array([1.0, 2.0, 3.0])
        self.assertAlmostEqual(mac(shape, -2.5 * shape), 1.0, places=12)
        self.assertAlmostEqual(mac(shape, 1j * shape), 1.0, places=12)
        self.assertAlmostEqual(mac(shape, np.array([3.0, 0.0, -1.0])), 0.0, places=12)

    def test_peak_fit_is_exact_for_a_single_mode(self):
        lines = np.linspace(4.6, 5.4, 40)
        natural, zeta = 5.0, 0.02
        spectrum = 3.0 * lines**4 / ((natural**2 - lines**2) ** 2 + (2 * zeta * natural * lines) ** 2)
        frequency, damping = fit_peak(lines, spectrum)
        self.assertAlmostEqual(frequency, natural, places=9)
        self.assertAlmostEqual(damping, zeta, places=9)

    def test_peak_fit_rejects_a_curve_without_a_peak(self):
        lines = np.linspace(4.6, 5.4, 40)
        trough = lines**4 / (700.0 - (lines**2 - 25.0) ** 2)
        self.assertIsNone(fit_peak(lines, trough))
        self.assertIsNone(fit_peak(lines[:15], 1 / (30.0 - lines[:15])))

    def test_spectral_matrix_is_hermitian_and_matches_welch(self):
        data = np.random.default_rng(3).normal(size=(4096, 3))
        frequencies, spectrum = spectral_matrix(data, 50.0, nperseg=512)
        np.testing.assert_allclose(spectrum, np.conj(np.swapaxes(spectrum, 1, 2)), atol=1e-12)
        reference = signal.welch(data[:, 1], fs=50.0, nperseg=512)[1]
        np.testing.assert_allclose(spectrum[:, 1, 1].real, reference, rtol=1e-10)
        self.assertEqual(spectrum.shape, (len(frequencies), 3, 3))

    def test_identification_recovers_the_model_modes(self):
        frame = uniform_frame()
        model_frequencies, model_shapes = frame.modes()
        record = ambient_record(frame, 600.0, 50.0, np.random.default_rng(4))
        estimate = identify_modes(record, 50.0, model_frequencies, model_shapes)
        np.testing.assert_allclose(estimate.frequencies, model_frequencies, rtol=0.01)
        for mode in range(5):
            self.assertGreater(mac(estimate.shapes[:, mode], model_shapes[:, mode]), 0.99)
            self.assertGreater(estimate.shapes[:, mode] @ model_shapes[:, mode], 0.0)
        self.assertTrue(np.all((estimate.damping > 0.01) & (estimate.damping < 0.035)))

    def test_modes_are_tracked_when_two_frequencies_move_close_together(self):
        frame = uniform_frame()
        model_frequencies, model_shapes = frame.modes()
        damaged = frame.with_damage(3, 0.3)
        record = ambient_record(damaged, 600.0, 50.0, np.random.default_rng(5))
        estimate = identify_modes(record, 50.0, model_frequencies, model_shapes)
        np.testing.assert_allclose(estimate.frequencies, damaged.modes()[0], rtol=0.012)


class DamageLocation(unittest.TestCase):
    def test_single_damaged_storey_is_located_from_exact_frequencies(self):
        frame = uniform_frame()
        healthy, fractions = frame.modes()[0], frame.strain_energy_fractions()
        np.testing.assert_allclose(locate_damage(healthy, healthy, fractions), 0.0, atol=1e-12)
        for storey in range(1, 6):
            losses = locate_damage(frame.with_damage(storey, 0.05).modes()[0], healthy, fractions)
            self.assertEqual(int(np.argmax(losses)), storey - 1)
            self.assertAlmostEqual(losses[storey - 1], 0.05, delta=0.005)
            self.assertLess(np.delete(losses, storey - 1).max(), 0.005)

    def test_two_damaged_storeys_are_separated(self):
        frame = uniform_frame()
        damaged = frame.with_damage(2, 0.05).with_damage(5, 0.04)
        losses = locate_damage(damaged.modes()[0], frame.modes()[0], frame.strain_energy_fractions())
        np.testing.assert_allclose(losses, [0.0, 0.05, 0.0, 0.0, 0.04], atol=0.006)


class Novelty(unittest.TestCase):
    def test_limit_gives_the_nominal_false_alarm_rate(self):
        rng = np.random.default_rng(6)
        mixing = rng.normal(size=(5, 5))
        rate = np.mean([NoveltyDetector(rng.normal(size=(60, 5)) @ mixing, 0.05).alarms(rng.normal(size=(300, 5)) @ mixing)
                        for _ in range(300)])
        self.assertAlmostEqual(rate, 0.05, delta=0.006)

    def test_single_feature_limit_is_the_t_prediction_limit(self):
        values = np.random.default_rng(7).normal(10.0, 2.0, 25)
        detector = NoveltyDetector(values, 0.01)
        self.assertAlmostEqual(detector.limit, stats.t.ppf(0.995, 24) ** 2 * 26 / 25, places=9)
        self.assertAlmostEqual(detector.distance(values[:1])[0], (values[0] - values.mean()) ** 2 / values.var(ddof=1),
                               places=9)

    def test_a_shifted_observation_is_flagged(self):
        rng = np.random.default_rng(8)
        detector = NoveltyDetector(rng.normal(size=(80, 5)), 0.01)
        self.assertTrue(detector.alarms(np.full(5, 5.0))[0])
        self.assertFalse(detector.alarms(np.zeros(5))[0])


class SensorLog(unittest.TestCase):
    def test_exceedances_reproduce_the_stored_anomalies(self):
        log = load_log(ROOT / "backend" / "data.csv")
        stored = pd.read_csv(ROOT / "backend" / "anomalies.csv", parse_dates=["timestamp"])
        pd.testing.assert_frame_equal(exceedances(log).reset_index(drop=True), stored)

    def test_log_checks(self):
        facts = check_log(load_log(ROOT / "backend" / "data.csv"))
        self.assertEqual(facts["readings"], 100)
        self.assertEqual(facts["interval_s"], [60.0])
        self.assertEqual(len(facts["sensors_reporting_several_types"]), 5)
        self.assertEqual(len(facts["sensors_reporting_several_locations"]), 5)
        self.assertEqual(facts["by_type"]["vibration"]["readings_at_maximum"], 18)
        self.assertEqual({name: row["above_limit"] for name, row in facts["by_type"].items()},
                         {"strain": 16, "temperature": 0, "vibration": 22})
        self.assertLess(facts["largest_robust_score"], 3.5)

    def test_robust_scores_single_out_an_outlier(self):
        log = pd.DataFrame({"sensor_type": ["strain"] * 7 + ["vibration"] * 5,
                            "value": [10, 11, 9, 10, 12, 10, 50, 5, 5, 5, 5, 9]})
        scores = robust_scores(log).to_numpy()
        self.assertGreater(scores[6], 3.5)
        self.assertLess(np.abs(scores[:6]).max(), 3.5)
        self.assertEqual(int(np.argmax(scores[7:])), 4)
        self.assertTrue(np.all(np.isfinite(scores)))

    def test_bad_logs_are_rejected(self):
        log = load_log(ROOT / "backend" / "data.csv")
        with self.assertRaises(ValueError):
            exceedances(log.assign(sensor_type="humidity"))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "log.csv"
            log.drop(columns="location").to_csv(path, index=False)
            with self.assertRaises(ValueError):
                load_log(path)


if __name__ == "__main__":
    unittest.main()
