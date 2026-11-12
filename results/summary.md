# Results

Command: `python run_analysis.py`

## Sensor log

| | |
|---|---|
| Readings | 100 from 5 sensor ids, one every 60 s |
| Period | 2025-01-12T08:45:19 to 2025-01-12T10:24:19 |
| Sensor ids that report more than one sensor type | 5 of 5 |
| Sensor ids that report more than one location | 5 of 5 |
| Readings above the limit of their own type | 38 |
| Readings flagged in the plots of the first version | 29 |
| Largest robust score of a reading | 3.34 (3.5 is the usual outlier limit) |

| Sensor type | Readings | Minimum | Maximum | Readings equal to the maximum | Limit | Above the limit |
|---|---|---|---|---|---|---|
| strain | 40 | 3.14 | 96.28 | 1 | 50 | 16 |
| temperature | 29 | 0.28 | 98.83 | 1 | 100 | 0 |
| vibration | 31 | 0.93 | 50 | 18 | 30 | 22 |

## Shear-frame study

Five storeys, 50 t per floor, 80 MN/m per storey, 2 % damping. Records of 10 minutes at 50 Hz with 5 % measurement noise; the stiffness of the whole frame differs by up to 1.5 % between records.

### One record against the model

| Mode | Model frequency (Hz) | Identified (Hz) | Error | Identified damping | MAC of the mode shape |
|---|---|---|---|---|---|
| 1 | 1.812 | 1.804 | -0.44 % | 2.37 % | 1.0000 |
| 2 | 5.289 | 5.283 | -0.12 % | 2.66 % | 0.9999 |
| 3 | 8.338 | 8.329 | -0.10 % | 2.13 % | 0.9999 |
| 4 | 10.711 | 10.711 | -0.01 % | 2.09 % | 0.9971 |
| 5 | 12.217 | 12.212 | -0.03 % | 1.87 % | 0.9998 |

### Scatter of the identified frequencies over 60 healthy records

| | Mode 1 | Mode 2 | Mode 3 | Mode 4 | Mode 5 |
|---|---|---|---|---|---|
| Constant temperature | 0.22 % | 0.15 % | 0.11 % | 0.11 % | 0.10 % |
| Varying temperature | 0.51 % | 0.49 % | 0.48 % | 0.45 % | 0.47 % |

### False alarms on 200 further healthy records

| Detector | Limit | False alarms |
|---|---|---|
| Five frequencies | 18.4 | 0.5 % |
| First frequency | 7.2 | 0.5 % |

### Records flagged, five frequencies (40 records per case)

| Stiffness loss | Storey 1 | Storey 2 | Storey 3 | Storey 4 | Storey 5 |
|---|---|---|---|---|---|
| 2 % | 8 % | 28 % | 18 % | 15 % | 12 % |
| 3 % | 25 % | 40 % | 72 % | 57 % | 15 % |
| 5 % | 88 % | 100 % | 98 % | 98 % | 85 % |
| 10 % | 100 % | 100 % | 100 % | 100 % | 100 % |
| 20 % | 100 % | 100 % | 100 % | 100 % | 100 % |

### Records flagged, first frequency only

| Stiffness loss | Storey 1 | Storey 2 | Storey 3 | Storey 4 | Storey 5 |
|---|---|---|---|---|---|
| 2 % | 0 % | 0 % | 0 % | 0 % | 0 % |
| 3 % | 0 % | 2 % | 2 % | 0 % | 0 % |
| 5 % | 20 % | 8 % | 5 % | 0 % | 0 % |
| 10 % | 78 % | 72 % | 35 % | 2 % | 0 % |
| 20 % | 100 % | 100 % | 100 % | 42 % | 0 % |

### Single records in which the largest estimated loss is in the damaged storey

| Stiffness loss | Storey 1 | Storey 2 | Storey 3 | Storey 4 | Storey 5 |
|---|---|---|---|---|---|
| 2 % | 88 % | 65 % | 48 % | 52 % | 50 % |
| 3 % | 95 % | 75 % | 78 % | 72 % | 75 % |
| 5 % | 100 % | 98 % | 92 % | 100 % | 98 % |
| 10 % | 100 % | 100 % | 100 % | 100 % | 100 % |
| 20 % | 100 % | 100 % | 100 % | 100 % | 100 % |

### Estimated loss in the damaged storey, from the mean frequencies of 40 records

| Stiffness loss | Storey 1 | Storey 2 | Storey 3 | Storey 4 | Storey 5 |
|---|---|---|---|---|---|
| 2 % | 1.7 % | 2.2 % | 1.8 % | 1.9 % | 2.0 % |
| 3 % | 3.0 % | 2.7 % | 3.0 % | 3.0 % | 2.7 % |
| 5 % | 5.2 % | 4.5 % | 4.8 % | 5.0 % | 5.1 % |
| 10 % | 10.1 % | 10.3 % | 9.9 % | 9.7 % | 10.8 % |
| 20 % | 21.8 % | 20.7 % | 19.0 % | 19.3 % | 20.6 % |

### Largest estimated loss in any other storey, same records

| Stiffness loss | Storey 1 | Storey 2 | Storey 3 | Storey 4 | Storey 5 |
|---|---|---|---|---|---|
| 2 % | 0.2 % | 0.6 % | 0.0 % | 0.2 % | 0.0 % |
| 3 % | 0.1 % | 0.3 % | 0.3 % | 0.0 % | 0.0 % |
| 5 % | 0.1 % | 0.0 % | 0.3 % | 0.6 % | 0.2 % |
| 10 % | 0.0 % | 0.5 % | 0.4 % | 0.4 % | 0.4 % |
| 20 % | 0.0 % | 2.0 % | 2.8 % | 2.2 % | 1.3 % |

### Exact change of the model frequencies for a 10 % loss

| Damaged storey | Mode 1 | Mode 2 | Mode 3 | Mode 4 | Mode 5 |
|---|---|---|---|---|---|
| 1 | -1.93 % | -1.59 % | -1.05 % | -0.52 % | -0.14 % |
| 2 | -1.63 % | -0.16 % | -0.61 % | -1.91 % | -0.92 % |
| 3 | -1.14 % | -0.59 % | -1.62 % | -0.16 % | -1.72 % |
| 4 | -0.59 % | -1.93 % | -0.16 % | -1.17 % | -1.39 % |
| 5 | -0.16 % | -1.17 % | -1.95 % | -1.48 % | -0.47 % |
