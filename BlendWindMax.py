# ----------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# BlendWindMax.py
# ----------------------------------------------------------------------------

ToolType = "numeric"
WeatherElementEdited = "Wind"
ScreenList = ["Wind"]

import numpy as np
import SmartScript
import sys

# Define the models you want to check for the maximum wind speed.
MODELS_TO_CHECK = ["NBM", "HRRR", "NAM12", "GFS", "UWWRF1", "MPASRT", "LRRFS"]

VariableList = [
    ("Blend Weight toward Max:", 0.50, "scale", [0.0, 1.0], 0.10),
]

class Tool (SmartScript.SmartScript):
    def __init__(self, dbss):
        self._dbss = dbss
        SmartScript.SmartScript.__init__(self, dbss)

    def execute(self, Wind, GridTimeRange, varDict):
        weight = varDict["Blend Weight toward Max:"]

        currSpeed = Wind[0]
        currDir = Wind[1]

        # Initialize the composite max grid with the current wind speed
        maxSpeed = np.copy(currSpeed)

        models_used = []

        # Loop through our model list to find the highest speeds
        for modelName in MODELS_TO_CHECK:
            # Check up to 12 runs back for HRRR, otherwise check 3 runs back
            if modelName == "HRRR":
                runOffsets = range(0, -12, -1) # [0, -1, -2, ..., -11]
            else:
                runOffsets = [0, -1, -2]

            for runOffset in runOffsets:
                dbid = self.findDatabase(modelName, runOffset)
                if dbid is None:
                    continue

                # Wrapped in a try/except to bypass the GFE 'exprName' bug
                # if a database doesn't contain the requested parameter.
                try:
                    modWind = self.getGrids(dbid, "Wind", "SFC", GridTimeRange, mode="Max", noDataError=0)
                except Exception:
                    modWind = None

                if modWind is not None:
                    modSpeed = modWind[0]
                    # Compare and keep the highest values at every grid point
                    maxSpeed = np.maximum(maxSpeed, modSpeed)
                    models_used.append(modelName)
                    break # Found the most recent available run for this model, move to next model

        # If no models were found at all, just return the current wind
        if not models_used:
            self.statusBarMsg("No models found. Returning current wind.", "U")
            return Wind

        self.statusBarMsg("Blended Max Wind using: " + ", ".join(models_used), "R")

        # Blend the current wind toward the composite maximum using the slider weight
        blendedSpeed = (maxSpeed * weight) + (currSpeed * (1.0 - weight))

        return (blendedSpeed, currDir)