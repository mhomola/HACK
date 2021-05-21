from Subsystem_design.common_constants import Constants

class FlightEnvelope(Constants):

    def __init__(self):
        super().__init__()

        # Conversions
        self.conv1 = 3.28084  # 1 m = conv1 ft
        self.conv2 = 0.00194032  # 1 kg/m^3 = conv2 slug/ft^3
        self.conv3 = 2.20462  # 1 kg = conv3 punds
        self.conv4 = 1.94384  # 1 m/s = conv4 knots


