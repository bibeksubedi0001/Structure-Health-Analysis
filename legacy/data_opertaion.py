import pandas as pd
from datetime import datetime
import os
import matplotlib.pyplot as plt

# Constants for file paths
DATA_FILE = "backend/data.csv"
VISUALIZATION_FOLDER = "reports/visualizations"
os.makedirs(VISUALIZATION_FOLDER, exist_ok=True)

# Utility Functions
def print_separator():
    """Prints a separator line for console output."""
    print("\n" + "-" * 80 + "\n")

# 1. Load Sensor Data
def load_data():
    """Load sensor data from the CSV file and return it as a DataFrame.

    Returns:
        pd.DataFrame: Loaded sensor data.
    """
    try:
        print("Attempting to load data...")
        if not os.path.exists(DATA_FILE):
            print(f"Error: {DATA_FILE} does not exist. Please ensure the file is available.")
            return pd.DataFrame()
        
        data = pd.read_csv(DATA_FILE)
        print("Data loaded successfully.")
        print_separator()
        print("Preview of loaded data:")
        print(data.head())
        print_separator()
        return data
    except pd.errors.EmptyDataError:
        print("Error: The data file is empty.")
        return pd.DataFrame()
    except Exception as e:
        print(f"An unexpected error occurred while loading data: {e}")
        return pd.DataFrame()

# 2. Save Sensor Data
def save_data(data):
    """Save the sensor data to the CSV file.

    Args:
        data (pd.DataFrame): The data to be saved.
    """
    try:
        print("Attempting to save data...")
        data.to_csv(DATA_FILE, index=False)
        print("Data saved successfully to", DATA_FILE)
    except Exception as e:
        print(f"An error occurred while saving data: {e}")

# 3. Plot Line Graph for All Sensors
def plot_all_sensors(data):
    """Plot line graphs of sensor readings over time for all sensors.

    Args:
        data (pd.DataFrame): The sensor data.
    """
    print("Plotting line graphs for all sensors...")
    unique_sensors = data["sensor_id"].unique()

    for sensor_id in unique_sensors:
        filtered_data = data[data["sensor_id"] == sensor_id]
        filtered_data["timestamp"] = pd.to_datetime(filtered_data["timestamp"], errors="coerce")
        filtered_data = filtered_data.dropna(subset=["timestamp"])

        plt.figure(figsize=(10, 6))
        plt.plot(filtered_data["timestamp"], filtered_data["value"], marker="o", linestyle="-", label=f"{sensor_id}")
        plt.title(f"Sensor Readings Over Time for {sensor_id}")
        plt.xlabel("Timestamp")
        plt.ylabel("Sensor Value")
        plt.legend()
        plt.grid(True)

        output_file = os.path.join(VISUALIZATION_FOLDER, f"line_graph_{sensor_id}.png")
        plt.savefig(output_file)
        plt.close()
        print(f"Line graph saved for {sensor_id} to {output_file}")

# 4. Compare Sensors in a Single Graph
def compare_sensors(data):
    """Compare all sensors on a single line graph.

    Args:
        data (pd.DataFrame): The sensor data.
    """
    print("Comparing all sensors on a single graph...")
    unique_sensors = data["sensor_id"].unique()

    plt.figure(figsize=(12, 8))
    for sensor_id in unique_sensors:
        filtered_data = data[data["sensor_id"] == sensor_id]
        filtered_data["timestamp"] = pd.to_datetime(filtered_data["timestamp"], errors="coerce")
        filtered_data = filtered_data.dropna(subset=["timestamp"])

        plt.plot(filtered_data["timestamp"], filtered_data["value"], marker="o", linestyle="-", label=f"{sensor_id}")

    plt.title("Comparison of Sensor Readings Over Time")
    plt.xlabel("Timestamp")
    plt.ylabel("Sensor Value")
    plt.legend()
    plt.grid(True)

    output_file = os.path.join(VISUALIZATION_FOLDER, "comparison_all_sensors.png")
    plt.savefig(output_file)
    plt.show()
    print(f"Comparison graph saved to {output_file}")

# 5. Highlight Anomalies for All Sensors
def highlight_anomalies_all_sensors(data, thresholds):
    """Highlight anomalies for all sensors on individual graphs.

    Args:
        data (pd.DataFrame): The sensor data.
        thresholds (dict): Dictionary of thresholds for each sensor type.
    """
    print("Highlighting anomalies for all sensors...")
    unique_sensors = data["sensor_id"].unique()

    for sensor_id in unique_sensors:
        sensor_data = data[data["sensor_id"] == sensor_id]
        sensor_type = sensor_data["sensor_type"].iloc[0]
        threshold = thresholds.get(sensor_type, None)

        if threshold is None:
            print(f"No threshold defined for sensor_type: {sensor_type}. Skipping {sensor_id}.")
            continue

        sensor_data["timestamp"] = pd.to_datetime(sensor_data["timestamp"], errors="coerce")
        sensor_data = sensor_data.dropna(subset=["timestamp"])

        anomalies = sensor_data[sensor_data["value"] > threshold]

        plt.figure(figsize=(10, 6))
        plt.plot(sensor_data["timestamp"], sensor_data["value"], label=f"{sensor_id} Readings", linestyle="-", marker="o")
        plt.scatter(anomalies["timestamp"], anomalies["value"], color="red", label="Anomalies", zorder=5)
        plt.title(f"Sensor Readings with Anomalies Highlighted ({sensor_id})")
        plt.xlabel("Timestamp")
        plt.ylabel("Sensor Value")
        plt.legend()
        plt.grid(True)

        output_file = os.path.join(VISUALIZATION_FOLDER, f"anomalies_{sensor_id}.png")
        plt.savefig(output_file)
        plt.close()
        print(f"Anomalies graph saved for {sensor_id} to {output_file}")

# Main Script for Visualization
def main():
    print("=== Sensor Data Visualization ===")
    sensor_data = load_data()

    if sensor_data.empty:
        print("No data to visualize. Exiting.")
        return

    # Line graphs for all sensors
    plot_all_sensors(sensor_data)

    # Comparison of all sensors
    compare_sensors(sensor_data)

    # Highlight anomalies for all sensors
    thresholds = {
        "strain": 50,
        "temperature": 100,
        "vibration": 30
    }
    highlight_anomalies_all_sensors(sensor_data, thresholds)

if __name__ == "__main__":
    main()
