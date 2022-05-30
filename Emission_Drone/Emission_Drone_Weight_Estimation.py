"""
Initial Weight Estimation for the emission and air pollution weight estimation
"""


def get_payload_mass():
    sensors_dict = {"CO-CX": 8, "NO-B4": 13, "SO2-A4": 6, "IRC-AT": 15, "OPC-R2": 30, "NO2-AE": 6, "O3-AH": 6}
    support_dict = {"fan": 80, "airpump": 10, "sensor_adapters": 40}  # masses in grams
    mass_support = sum(support_dict.values())
    mass_all_sensors = sum(sensors_dict.values())
    return mass_all_sensors + mass_support

def get_initial_mass_est():
    payload_mass = get_payload_mass()
    weights_dict = {"structure":800, "motors":500, "propellers":60, "flight_controllers":49, \
            "ESC's":120, "Arduino":25, "battery":800, "wiring":100, "payload":payload_mass}
    return sum(weights_dict.values())

class Drone:
    def __init__(self, total_mass=0, payload=get_payload_mass()):
        self.total_mass = total_mass
        self.payload = payload

    def MassIteration(self):
        pass

