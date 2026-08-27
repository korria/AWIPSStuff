# ----------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# BlendGustMax.py
# ----------------------------------------------------------------------------

ToolType = "numeric"
WeatherElementEdited = "WindGust"
ScreenList = ["Wind", "WindGust"]

import numpy as np
import SmartScript
import sys

# Define the models you want to check for the maximum gust.
MODELS_TO_CHECK = ["NBM", "NBM90Pct", "HRRR", "NAM12", "GFS", "UWWRF1", "MPASRT", "LRRFS"]

VariableList = [
    ("Blend Weight toward Max:", 0.50, "scale", [0.0, 1.0], 0.10),
    ("Fallback Gust Multiplier:", 1.40, "scale", [1.0, 2.0], 0.05)
]

class Tool (SmartScript.SmartScript):
    def __init__(self, dbss):
        self._dbss = dbss
        SmartScript.SmartScript.__init__(self, dbss)

    def execute(self, WindGust, Wind, GridTimeRange, varDict):
        weight = varDict["Blend Weight toward Max:"]
        multiplier = varDict["Fallback Gust Multiplier:"]

        currSpeed = Wind[0]
        currGust = WindGust

        # Initialize the max grid. Use the current gust, but ensure it's at least
        # the fallback multiplier applied to the current wind speed.
        maxGust = np.maximum(currGust, currSpeed * multiplier)

        models_used = []

        # Loop through our model list to find the highest gusts
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
                    modGust = self.getGrids(dbid, "WindGust", "SFC", GridTimeRange, mode="Max", noDataError=0)
                except Exception:
                    modGust = None

                if modGust is not None:
                    # Keep the highest values at every grid point
                    maxGust = np.maximum(maxGust, modGust)
                    if modelName not in models_used:
                        models_used.append(modelName)
                    break

        if models_used:
            self.statusBarMsg("Blended Max Gust using: " + ", ".join(models_used), "R")
        else:
            self.statusBarMsg("No model gusts found. Using fallback multiplier.", "U")

        # Blend the current gust toward the composite maximum
        blendedGust = (maxGust * weight) + (currGust * (1.0 - weight))

        # Constraint: Ensure gust is at least Wind Speed + 1 kt
        diff = blendedGust - currSpeed
        mask = diff < 1.0
        blendedGust[mask] = (currSpeed + 1.0)[mask]

        return blendedGust