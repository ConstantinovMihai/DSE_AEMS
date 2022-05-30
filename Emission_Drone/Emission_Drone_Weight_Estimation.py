"""
Initial Weight Estimation for the emission and air pollution weight estimation
"""

def get_payload_mass():
    sensors_dict = {"CO-CX":8, "NO-B4":13, "SO2-A4":6, "IRC-AT":15, "OPC-R2":30, "NO2-AE":6, "O3-AH":6}
    mass_all_sensors = sum(sensors_dict.values())
    return mass_all_sensors

