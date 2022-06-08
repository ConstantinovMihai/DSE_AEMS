"""
Initial sizing for the emission and air pollution weight estimation
"""


class Drone:
    @staticmethod
    def get_payload_mass():
        sensors_dict = {"CO-CX": 8, "NO-B4": 13, "SO2-A4": 6, "IRC-AT": 15, "OPC-R2": 30, "NO2-AE": 6, "O3-AH": 6}
        support_dict = {"fan": 80, "airpump": 10, "sensor_adapters": 100}  # masses in grams
        return sum(support_dict.values()) + sum(sensors_dict.values())

    @staticmethod
    def get_initial_mass_est():
        payload_mass = Drone.get_payload_mass()
        weights_dict = {"structure": 800, "motors": 352, "propellers": 60, "flight_controllers": 49,\
                        "ESC's": 88, "Arduino": 25, "battery": 792, "wiring": 100, "Obstacle avoidance": 90,\
                       "payload": payload_mass}
        return weights_dict

    def __init__(self):
        self.mass_dict = Drone.get_initial_mass_est()

    def propellers(self):
        pass

drone = Drone()
print(sum(drone.mass_dict.values()))