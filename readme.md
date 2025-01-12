# **Structural Health Sensor Analysis and Visualization**

## **Overview**
This project is designed to analyze and visualize sensor data for monitoring the health of structures like bridges, towers, and buildings. It provides various features, including anomaly detection, sensor comparisons, and data visualization.

The project is implemented in Python using **Pandas** for data manipulation, **Matplotlib** for visualizations, and file operations for dynamic data handling.


## **Technologies Used**
- **Python**: Core programming language.
- **Pandas**: For data manipulation and processing.
- **Matplotlib**: For creating line graphs, comparison plots, and anomaly visualizations.
- **OS Module**: For file and folder operations.

---

## **Folder Structure**
```
Structural_Health_Sensor_Analysis/
│
├── backend/
│   ├── data.csv                  # Sensor data file
│   ├── data_operation.py   # Python script for data analysis and visualization
│
├── reports/
│   ├── visualizations/           # Folder for saved graphs
│
├── README.md                     # Project documentation
```

---

## **Setup Instructions**
1. Clone the repository:
   ```bash
   git clone "https://github.com/bibeksubedi0001/Structure-Health-Analysis"
   cd Structural_Health_Sensor_Analysis
   ```

2. Install dependencies:
   ```bash
   pip install pandas matplotlib
   ```

3. Place the `data.csv` file in the `backend/` folder. This file should contain the sensor data with the following columns:
   - `sensor_id`
   - `sensor_type`
   - `location`
   - `value`
   - `timestamp`

4. Run the script:
   ```bash
   python backend/visualization_script.py
   ```

---

## **How to Use**
### **1. Load Sensor Data**
The script will load and preview the data from `data.csv`.

### **2. Generate Visualizations**
- **Line Graphs**: Individual graphs for each sensor, saved in the `reports/visualizations/` folder.
- **Sensor Comparisons**: A single graph comparing all sensors.
- **Anomalies**: Highlighted graphs for each sensor, showing values exceeding predefined thresholds.

### **3. Adjust Thresholds**
Modify the `thresholds` dictionary in the script to change anomaly detection limits:
```python
thresholds = {
    "strain": 50,
    "temperature": 100,
    "vibration": 30
}
```

---

## **Outputs**
- **Line Graphs**: `line_graph_<sensor_id>.png`
- **Comparison Graph**: `comparison_all_sensors.png`
- **Anomaly Graphs**: `anomalies_<sensor_id>.png`

All outputs are saved in the `reports/visualizations/` folder.

---

## **Example Data Format**
Here’s an example of the `data.csv` file format:
```csv
sensor_id,sensor_type,location,value,timestamp
sensor_1,strain,Bridge A,45.5,2025-01-12T12:00:00
sensor_2,vibration,Tower C,30.2,2025-01-12T12:01:00
sensor_3,temperature,Building B,102.3,2025-01-12T12:02:00
```
**Author**

Name: Bibek Subedi

Email: 078bce035.bibek@pcampus.edu.np