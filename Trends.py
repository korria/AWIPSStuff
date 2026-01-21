#--------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# Trends - version 2.0
#
# Show trends of values over many models/runs with optional observed data overlay.
#
# 1.0 - initial version - Timothy Barker  
# 2.0 - Modernized with PySide6 and observation support
# ----------------------------------------------------------------------------
#
# C O N F I G U R A T I O N   S E C T I O N
#
# ----------------------------------------------------------------------------

ModelConfig = [
    "Fcst           |#000000|X",
    "GFS            |#7f0000|GFS",
    "GFSBC          |#ff0000|GFS",
    "ADJMEX         |#7f3f00|GFS",
    "ADJMEXBC       |#ff7f00|GFS",
    "ADJMAV         |#bfbf00|GFS",
    "ADJMAVBC       |#ffff00|GFS",
    "ECMWF          |#00007f|EC",
    "ECMWFBC        |#0000ff|EC",
    "ADJECE         |#3f007f|EC",
    "ADJECEBC       |#7f007f|EC",
    "ADJECS         |#bf00bf|EC",
    "ADJECSBC       |#ff00ff|EC",
    "ADJECH         |#7f00ff|EC",
    "ADJECM         |#7f7fff|EC",
    "ADJECL         |#ffff7f|EC",
    "CMCnh          |#7fbfff|CAN",
    "CMCnhBC        |#3f7fff|CAN",
    "CMCreg         |#7fbfff|CAN",
    "CMCregBC       |#3f7fff|CAN",
    "NAM12          |#007f00|NAM",
    "NAM12BC        |#00ff00|NAM",
    "ADJMET         |#007f7f|NAM",
    "ADJMETBC       |#00bfbf|NAM",
    "SREF           |#003f3f|NAM",
    "SREFBC         |#007f7f|NAM",
    "HIRESWarw      |#1f7f3f|NAM",
    "HIRESWarwBC    |#3fff7f|NAM",
    "HIRESWnmm      |#3fff1f|NAM",
    "HIRESWnmmBC    |#7fff3f|NAM",
    "NBM            |#3f3fff|CONS",
    "NBMEXP         |#3030f0|CONS",
    "CONSAll        |#3f3f3f|CONS",
    "BCCONSAll      |#7fffff|CONS",
    "CONSRaw        |#7f7f7f|CONS",
    "BCCONSRaw      |#ff7fff|CONS",
    "CONSMOS        |#bfbfbf|CONS",
    "BCCONSMOS      |#ffff7f|CONS",
    "SuperBlend     |#7F007f|CONS",
    "WPCGuide       |#3f7fff|WPC",
    "WPCGuideBC     |#7fbfff|WPC",
]

# Observation database options
ObsSourceConfig = {
    "Obs": {"color": "#00FF00", "dash": [6, 3], "label": "Observed"},
    "URMA": {"color": "#FF6600", "dash": [6, 3], "label": "URMA Analysis"},
    "RTMA": {"color": "#9900FF", "dash": [6, 3], "label": "RTMA Analysis"},
}

LabelSizeStrings = ["8", "10", "12", "14"]

# ----------------------------------------------------------------------------
#
# E N D   C O N F I G U R A T I O N   S E C T I O N
#
# ============================================================================

ToolType = "numeric"
WeatherElementEdited = "None"
ScreenList = ["SCALAR", "VECTOR"]

import numpy as np
import SmartScript
import calendar
import math
import time
import sys
import copy

import BOIVerifyUtility
import TimeRange

# PySide6 imports are deferred to avoid conflicts with AWIPS Jep at module load time
# They will be imported when the dialog is actually created
PySide6_QtWidgets = None
PySide6_QtCore = None
PySide6_QtGui = None

def _import_pyside6():
    """Deferred import of PySide6 to avoid Jep/shiboken conflicts at module load."""
    global PySide6_QtWidgets, PySide6_QtCore, PySide6_QtGui
    if PySide6_QtWidgets is None:
        from PySide6 import QtWidgets as PySide6_QtWidgets_mod
        from PySide6 import QtCore as PySide6_QtCore_mod
        from PySide6 import QtGui as PySide6_QtGui_mod
        PySide6_QtWidgets = PySide6_QtWidgets_mod
        PySide6_QtCore = PySide6_QtCore_mod
        PySide6_QtGui = PySide6_QtGui_mod
    return PySide6_QtWidgets, PySide6_QtCore, PySide6_QtGui


class Tool(SmartScript.SmartScript):
    def __init__(self, dbss):
        SmartScript.SmartScript.__init__(self, dbss)
        self._dbss = dbss

    def preProcessTool(self):
        self.VU = BOIVerifyUtility.BOIVerifyUtility(self._dbss, None)

        # Get the model configuration info into self
        self.MODELS = []
        self.COLORS = {}
        self.primaryGroup = {}
        self.state = {}
        self.lineWidth = {}
        self.pointState = {}
        self.labelState = {}
        self.labelSize = {}
        self.groups = {}
        
        for modelInfo in ModelConfig:
            (rawMod, rawColor, rawGroup) = modelInfo.split("|")
            mod = rawMod.strip()
            color = rawColor.strip()
            self.MODELS.append(mod)
            self.COLORS[mod] = color
            group = rawGroup.strip()
            self.primaryGroup[mod] = group
            if group in self.groups:
                self.groups[group].append(mod)
            else:
                self.groups[group] = [mod]
            self.state[mod] = 1
            self.lineWidth[mod] = 1
            self.pointState[mod] = 1
            self.labelState[mod] = 0
            self.labelSize[mod] = 1

        # Observation source settings
        self.obsSource = "None"  # "None", "Obs", "URMA", "RTMA"
        self.obsData = {}
        self.obsFull = {}
        self.showObs = False

        self.readStatus()
        
        # Setup the graphing settings
        self.outdata = {}
        self.outfull = {}
        self.saveGrids = {}
        self.saveElement = "none"
        self.savePStart = 0
        self.savePEnd = 0
        self.mintime = self._gmtime().unixTime()
        self.maxtime = self._gmtime().unixTime() + 3600
        self.graphStart = 0
        self.graphmin = 0.0
        self.allmax = 1.0
        self.labelText = "inches"
        self.titleText = "Trends"
        self.numpts = 0
        self.absMin = 0.0
        self.absMax = 100.0
        self.perStart = self._gmtime().unixTime()
        self.perEnd = self._gmtime().unixTime() + 3600
        self.minTick = 1
        self.element = "MaxT"
        self.cursorReadout = 0

    def preProcessGrid(self):
        # Import PySide6 at runtime to avoid Jep conflicts
        QtWidgets, QtCore, QtGui = _import_pyside6()
        
        # Create the PySide6 application if needed
        self.app = QtWidgets.QApplication.instance()
        if self.app is None:
            self.app = QtWidgets.QApplication(sys.argv)
        
        self.dlg = TrendDialog(
            title="Model Trends",
            callback_method=self.getStuff,
            popup_callback=self.popupHandler,
            parent_tool=self
        )
        
        self.getStuff("Update")
        self.dlg.exec()
        self.cancel()

    def execute(self):
        "Model/Forecast trend over time"
        return

    def getStuff(self, buttonName):
        QtWidgets, QtCore, QtGui = _import_pyside6()
        
        if buttonName == "Close":
            self.saveStatus()
            return
        
        if buttonName == "Resize":
            self.drawGraph(
                self.outdata, self.outfull,
                max(self.graphStart, self.mintime), self.maxtime,
                self.graphmin, self.allmax,
                label=self.labelText, title=self.titleText,
                colors=self.COLORS, numpts=self.numpts
            )
            return
        
        if buttonName == "Redraw":
            self.drawGraph(
                self.outdata, self.outfull,
                max(self.graphStart, self.mintime), self.maxtime,
                self.graphmin, self.allmax,
                label=self.labelText, title=self.titleText,
                colors=self.COLORS, numpts=self.numpts
            )
            return

        # Disable the Update button while working
        self.dlg.updateButton.setText("WORKING")
        self.dlg.updateButton.setEnabled(False)
        self.dlg.updateButton.setStyleSheet("background-color: red; color: white;")
        QtWidgets.QApplication.processEvents()

        # Get the current editArea
        editArea = self.getActiveEditArea()
        ea = self.encodeEditArea(editArea)
        eaFlat = np.ravel(ea)
        self.numpts = np.sum(ea)
        
        if self.numpts == 0:
            ea = (self.getTopo() * 0) + 1
            eaFlat = np.ravel(ea)
            self.numpts = np.sum(ea)

        # Get selected Parm
        selectedParms = self.selectedParms()
        if selectedParms is None or not selectedParms:
            self.dlg.updateButton.setEnabled(True)
            self.dlg.updateButton.setText("Update")
            self.dlg.updateButton.setStyleSheet("")
            return
        
        (selectedName, selectedLevel, selectModelID) = selectedParms[0]
        parmName = selectedName

        selectTR = TimeRange.TimeRange(self._dbss.getParmOp().getSelectionTimeRange())
        
        # Get the mutable database id
        mutid = self.mutableID()
        mutName = mutid.modelName()

        startTime = selectTR.startTime()
        durHours = int((selectTR.duration() / 3600.0) + 0.5)
        (ntr, gridTimes) = self.getGridTimes(mutid, parmName, "SFC", startTime, durHours)
        start = gridTimes[0].startTime().unixTime()
        end = gridTimes[-1].endTime().unixTime()
        self.perStart = start
        self.perEnd = end
        durHours = int(((end - start) / 3600.0) + 0.5)
        (gyea, gmon, gday, ghou, gmin, gsec, gwda, gyda, gdst) = self._gmtime().timetuple()
        zz = calendar.timegm((gyea, gmon, gday, 0, 0, 0, 0, 0, 0))
        hourstart = int(round((start - zz) / 3600.0))
        hourend = int(round((end - zz) / 3600.0))
        FullTimeRange = self.createTimeRange(hourstart, hourend, "Zulu")

        reqElement = parmName
        accumVal = 0
        if reqElement in ["QPF", "SnowAmt"]:
            accumVal = 1
            parmTitle = "%d-hr %s" % (durHours, parmName)
        else:
            parmTitle = parmName
        self.titleText = "Trend of %s" % parmTitle

        gridInfos = self.getGridInfo(mutid, parmName, "SFC", FullTimeRange)
        if not gridInfos:
            self.dlg.updateButton.setEnabled(True)
            self.dlg.updateButton.setText("Update")
            self.dlg.updateButton.setStyleSheet("")
            return

        gridMax = gridInfos[0].maxLimit()
        if accumVal == 1:
            self.absMax = gridMax * 20.0
        else:
            self.absMax = gridMax
        self.absMin = gridInfos[0].minLimit()
        inUnits = gridInfos[0].units()
        inPrecision = gridInfos[0].precision()

        if inPrecision < 1:
            self.minTick = 1
        else:
            self.minTick = 10.0 ** (-inPrecision)

        self.labelText = inUnits

        if inPrecision == 0:
            self.valueFormat = "%.0f"
        else:
            self.valueFormat = "%." + ("%d" % inPrecision) + "f"

        # Get information about the currently displayed parm
        parmObj = self.getParm(mutid, parmName, "SFC")
        if parmObj is not None:
            from com.raytheon.viz.gfe.rsc import DiscreteDisplayUtil
            ctInfo = DiscreteDisplayUtil.buildColorMapParameters(parmObj)
            if ctInfo is not None:
                colorTable = ctInfo.getColorMapName()
                displayMinval = ctInfo.getColorMapMin()
                displayMaxval = ctInfo.getColorMapMax()

        lev = "SFC"
        self.outdata = {}
        self.outfull = {}
        outkeys = []
        self.allmin = 10000.0
        self.allmax = 0.0
        self.mintime = self.perStart
        self.maxtime = self.perEnd

        # Clear the cache if not the same as before
        if (self.saveElement != reqElement) or (self.savePStart != self.perStart) or (self.savePEnd != self.perEnd):
            del self.saveGrids
            self.saveGrids = {}
            self.saveFull = {}
            self.savePStart = self.perStart
            self.savePEnd = self.perEnd
            self.saveElement = reqElement
            useCache = 0
        else:
            useCache = 1

        modsWithData = []
        for mod in self.MODELS:
            maxvers = -40
            if mod in ("Fcst", "ISC", "NBM", "NBMEXP"):
                maxvers = -1

            element = reqElement
            self.element = reqElement
            tempoff = 0
            if self.state[mod] == 0:
                tempoff = 1
                self.state[mod] = 1

            # Try getting the most recent data by reading grids from GFE
            for vers in range(0, maxvers, -1):
                modID = self.findDatabase(mod, version=vers)
                if modID is None:
                    break
                if not modID.isValid():
                    break
                modelName = modID.modelName()

                parmObj = self.getParm(modID, element, lev)
                if parmObj is None:
                    continue

                if mod not in ("Fcst", "ISC", "NBM", "NBMEXP"):
                    modTimeUnix = modID.modelTime().unixTime()
                else:
                    if self._gmtime().unixTime() > self.perStart:
                        modTimeUnix = self.perStart
                    else:
                        modTimeUnix = self._gmtime().unixTime()

                outkey = "%s,%d" % (mod, modTimeUnix)

                if outkey in self.saveGrids:
                    finalGrid = self.saveGrids[outkey]
                    fullFlag = self.saveFull[outkey]
                else:
                    gridInfos = self.getGridInfo(modID, element, lev, FullTimeRange)
                    numAdded = 0
                    if reqElement in ("QPF", "SnowAmt"):
                        finalGrid = self.empty()
                    elif reqElement == "MaxT":
                        finalGrid = self.newGrid(-500)
                    elif reqElement == "MinT":
                        finalGrid = self.newGrid(+500)
                    else:
                        finalGrid = self.empty()

                    fullDuration = FullTimeRange.duration()
                    totTime = 0
                    for gridInfo in gridInfos:
                        gridTR = gridInfo.gridTime()
                        duration = gridTR.duration()
                        intTR = gridTR.intersection(FullTimeRange)
                        inside = intTR.duration()
                        if (duration < 1) or (inside < 1):
                            continue
                        pctInside = float(inside) / float(duration)
                        
                        grid = self.getGrids(modID, element, lev, gridTR, noDataError=0)
                        if grid is None:
                            continue
                        if not isinstance(grid, np.ndarray) and len(grid) == 2:
                            (mag, direc) = grid
                            grid = mag

                        if reqElement in ("QPF", "SnowAmt"):
                            finalGrid += (grid * pctInside)
                        elif reqElement == "MaxT":
                            finalGrid = np.maximum(finalGrid, grid)
                        elif reqElement == "MinT":
                            finalGrid = np.minimum(finalGrid, grid)
                        else:
                            finalGrid += grid
                        totTime += inside
                        numAdded += 1

                    if numAdded == 0:
                        continue

                    if (reqElement not in ("QPF", "MaxT", "MinT", "SnowAmt")) and (numAdded > 1):
                        finalGrid /= float(numAdded)
                    if reqElement in ("QPF", "SnowAmt"):
                        if totTime >= fullDuration:
                            fullFlag = 1
                        else:
                            fullFlag = 0
                    else:
                        fullFlag = 1

                    self.saveGrids[outkey] = copy.copy(finalGrid)
                    self.saveFull[outkey] = fullFlag

                statsList = self.getStats(finalGrid, eaFlat)
                self.outdata[outkey] = copy.copy(statsList)
                self.outfull[outkey] = fullFlag
                self.allmax = max(self.allmax, statsList[1])
                self.allmin = min(self.allmin, statsList[0])
                self.mintime = min(self.mintime, modTimeUnix)
                self.maxtime = max(self.maxtime, modTimeUnix)
                
                if mod not in modsWithData:
                    modsWithData.append(mod)
                if len(modsWithData) < 10:
                    self.drawGraph(
                        self.outdata, self.outfull,
                        max(self.graphStart, self.mintime), self.maxtime,
                        self.graphmin, self.allmax,
                        label=self.labelText, title=self.titleText,
                        colors=self.COLORS, numpts=self.numpts
                    )

            outkeys = list(self.outdata.keys())

            if useCache == 0:
                archiveMod = mod
                if mod == "Fcst":
                    archiveMod = "Official"
                if reqElement in ("QPF", "SnowAmt"):
                    modeVal = "Sum"
                elif reqElement == "MaxT":
                    modeVal = "Max"
                elif reqElement == "MinT":
                    modeVal = "Min"
                else:
                    modeVal = "TimeWtAverage"

                cases = self.getModelCases(element, archiveMod, self.perStart, self.perEnd, accum=accumVal)
                casekeys = sorted(list(cases.keys()))
                casekeys.reverse()
                
                for casekey in casekeys:
                    (basestr, startstr, endstr) = casekey.split(",")
                    outkey = "%s,%s" % (mod, basestr)
                    base = int(basestr)

                    if outkey not in outkeys:
                        t0 = time.time()
                        if outkey in self.saveGrids:
                            verGrid = self.saveGrids[outkey]
                        else:
                            start = int(startstr)
                            end = int(endstr)
                            recs = cases[casekey]
                            verGrid = self.VU.getVerGrids(
                                archiveMod, base, element, self.perStart, self.perEnd,
                                mode=modeVal, recList=recs
                            )
                            if verGrid is None:
                                continue
                            if not isinstance(verGrid, np.ndarray) and len(verGrid) == 2:
                                (mag, direc) = verGrid
                                verGrid = mag
                            if reqElement in ("QPF", "SnowAmt"):
                                totTime = 0
                                for rec in recs:
                                    st = self.VU.fncStime[rec]
                                    en = self.VU.fncEtime[rec]
                                    totTime += (en - st)
                                if totTime >= FullTimeRange.duration():
                                    fullFlag = 1
                                else:
                                    fullFlag = 0
                            else:
                                fullFlag = 1

                            self.saveGrids[outkey] = copy.copy(verGrid)
                            self.saveFull[outkey] = fullFlag

                        statsList = self.getStats(verGrid, eaFlat)
                        self.outdata[outkey] = copy.copy(statsList)
                        self.outfull[outkey] = fullFlag
                        self.allmax = max(self.allmax, statsList[1])
                        self.allmin = min(self.allmin, statsList[0])
                        self.mintime = min(self.mintime, base)
                        self.maxtime = max(self.maxtime, base)
                        outkeys = list(self.outdata.keys())

                        if accumVal == 1:
                            self.graphmin = 0.0
                        else:
                            self.graphmin = self.allmin
                        if mod not in modsWithData:
                            modsWithData.append(mod)
                        if len(modsWithData) < 10:
                            self.drawGraph(
                                self.outdata, self.outfull,
                                max(self.graphStart, self.mintime), self.maxtime,
                                self.graphmin, self.allmax,
                                label=self.labelText, title=self.titleText,
                                colors=self.COLORS, numpts=self.numpts
                            )

                if len(modsWithData) >= 10:
                    self.drawGraph(
                        self.outdata, self.outfull,
                        max(self.graphStart, self.mintime), self.maxtime,
                        self.graphmin, self.allmax,
                        label=self.labelText, title=self.titleText,
                        colors=self.COLORS, numpts=self.numpts
                    )
            else:
                gridkeys = sorted(list(self.saveGrids.keys()))
                gridkeys.reverse()
                for key in gridkeys:
                    (modstr, basestr) = key.split(",")
                    if modstr == mod:
                        verGrid = self.saveGrids[key]
                        base = int(basestr)
                        statsList = self.getStats(verGrid, eaFlat)
                        self.outdata[key] = copy.copy(statsList)
                        self.outfull[key] = self.saveFull[key]
                        self.allmax = max(self.allmax, statsList[1])
                        self.allmin = min(self.allmin, statsList[0])
                        self.mintime = min(self.mintime, base)
                        self.maxtime = max(self.maxtime, base)
                        outkeys = list(self.outdata.keys())

            if accumVal == 1:
                self.graphmin = 0.0
            else:
                self.graphmin = self.allmin
            self.drawGraph(
                self.outdata, self.outfull,
                max(self.graphStart, self.mintime), self.maxtime,
                self.graphmin, self.allmax,
                label=self.labelText, title=self.titleText,
                colors=self.COLORS, numpts=self.numpts
            )

            if tempoff == 1:
                self.state[mod] = 0
                self.drawGraph(
                    self.outdata, self.outfull,
                    max(self.graphStart, self.mintime), self.maxtime,
                    self.graphmin, self.allmax,
                    label=self.labelText, title=self.titleText,
                    colors=self.COLORS, numpts=self.numpts
                )

            QtWidgets.QApplication.processEvents()

        # Fetch observation data if enabled
        if self.showObs and self.obsSource != "None":
            self.fetchObsData(reqElement, eaFlat, FullTimeRange, accumVal)

        if accumVal == 1:
            self.graphmin = 0.0
        else:
            self.graphmin = self.allmin
        self.drawGraph(
            self.outdata, self.outfull,
            max(self.graphStart, self.mintime), self.maxtime,
            self.graphmin, self.allmax,
            label=self.labelText, title=self.titleText,
            colors=self.COLORS, numpts=self.numpts
        )

        # All done - return to waiting
        self.dlg.updateButton.setEnabled(True)
        self.dlg.updateButton.setText("Update")
        self.dlg.updateButton.setStyleSheet("")
        return

    def fetchObsData(self, element, eaFlat, FullTimeRange, accumVal):
        """Fetch observation data from Obs, URMA, or RTMA database."""
        self.obsData = {}
        self.obsFull = {}
        
        obsSource = self.obsSource
        if obsSource == "None":
            return
        
        # Map element names for observations if needed
        obsElement = element
        if element == "MaxT":
            obsElement = "T"  # Observations use T, we'll take max
        elif element == "MinT":
            obsElement = "T"  # Observations use T, we'll take min
        
        # Determine mode based on element
        if element in ("QPF", "SnowAmt"):
            modeVal = "Sum"
        elif element == "MaxT":
            modeVal = "Max"
        elif element == "MinT":
            modeVal = "Min"
        else:
            modeVal = "TimeWtAverage"
        
        try:
            # Check if the observation database exists and has data
            if not self.VU.checkFile(obsElement, obsSource):
                # Try alternate names
                if obsSource == "Obs":
                    altSources = ["Metar", "METAR", "Obs"]
                elif obsSource == "URMA":
                    altSources = ["URMA25", "URMA", "URMAe"]
                elif obsSource == "RTMA":
                    altSources = ["RTMA25", "RTMA", "RTMAe"]
                else:
                    altSources = []
                
                found = False
                for altSource in altSources:
                    if self.VU.checkFile(obsElement, altSource):
                        obsSource = altSource
                        found = True
                        break
                
                if not found:
                    return
            
            # Get observation cases
            obsCases = self.getModelCases(obsElement, obsSource, self.perStart, self.perEnd, accum=accumVal)
            
            if not obsCases:
                return
            
            casekeys = sorted(list(obsCases.keys()))
            
            for casekey in casekeys:
                (basestr, startstr, endstr) = casekey.split(",")
                base = int(basestr)
                start = int(startstr)
                end = int(endstr)
                
                recs = obsCases[casekey]
                obsGrid = self.VU.getVerGrids(
                    obsSource, base, obsElement, self.perStart, self.perEnd,
                    mode=modeVal, recList=recs
                )
                
                if obsGrid is None:
                    continue
                
                if not isinstance(obsGrid, np.ndarray) and len(obsGrid) == 2:
                    (mag, direc) = obsGrid
                    obsGrid = mag
                
                # Determine full coverage flag
                if element in ("QPF", "SnowAmt"):
                    totTime = 0
                    for rec in recs:
                        st = self.VU.fncStime[rec]
                        en = self.VU.fncEtime[rec]
                        totTime += (en - st)
                    if totTime >= FullTimeRange.duration():
                        fullFlag = 1
                    else:
                        fullFlag = 0
                else:
                    fullFlag = 1
                
                # Store the observation data with a special key
                obsKey = "OBS_%s,%d" % (self.obsSource, base)
                statsList = self.getStats(obsGrid, eaFlat)
                self.obsData[obsKey] = copy.copy(statsList)
                self.obsFull[obsKey] = fullFlag
                
                # Update min/max
                self.allmax = max(self.allmax, statsList[1])
                self.allmin = min(self.allmin, statsList[0])
                
        except Exception as e:
            print("Error fetching observation data: %s" % str(e))
            sys.stdout.flush()

    def getStats(self, grid, areaFlat):
        flatGrid = np.ravel(grid)
        data = np.compress(areaFlat, flatGrid)
        sdata = np.sort(data)
        datapts = data.shape[0]
        minvalue = sdata[0]
        maxvalue = sdata[datapts - 1]
        total = np.add.reduce(sdata)
        avg = float(total) / float(datapts)
        p01 = sdata[self.pointNumber(0.01, datapts)]
        p05 = sdata[self.pointNumber(0.05, datapts)]
        p10 = sdata[self.pointNumber(0.10, datapts)]
        p25 = sdata[self.pointNumber(0.25, datapts)]
        p50 = sdata[self.pointNumber(0.50, datapts)]
        p75 = sdata[self.pointNumber(0.75, datapts)]
        p90 = sdata[self.pointNumber(0.90, datapts)]
        p95 = sdata[self.pointNumber(0.95, datapts)]
        p99 = sdata[self.pointNumber(0.99, datapts)]
        return (minvalue, maxvalue, avg, p01, p05, p10, p25, p50, p75, p90, p95, p99)

    def pointNumber(self, percentage, datapts):
        return min(max(int(percentage * datapts) - 1, 0), datapts)

    def getModelCases(self, parm, model, perStart, perEnd, accum=0):
        cases = {}
        if not self.VU.checkFile(parm, model):
            return cases

        totalTime = perEnd - perStart
        endsin = (self.VU.fncEtime[:] > perStart) & (self.VU.fncEtime[:] <= perEnd)
        startsin = (self.VU.fncStime[:] < perEnd) & (self.VU.fncStime[:] >= perStart)
        crosses = (self.VU.fncStime[:] < perStart) & (self.VU.fncEtime[:] > perEnd)
        recmatch = endsin | startsin | crosses

        if recmatch.any():
            recnumberList = list(np.compress(recmatch, self.VU.fncRecs[:]))
            baselist = []
            for rec in recnumberList:
                recnum = int(rec)
                if self.VU.fncBtime[recnum] not in baselist:
                    baselist.append(self.VU.fncBtime[recnum])

            for base in baselist:
                reclist = []
                for rec in recnumberList:
                    recnum = int(rec)
                    if self.VU.fncBtime[recnum] == base:
                        reclist.append(recnum)
                key = "%d,%d,%d" % (base, perStart, perEnd)
                cases[key] = reclist
        return cases

    def drawGraph(self, graphData, fullData, xmin, xmax, ymin, ymax,
                  label="value", title="Graph", scaleDist=1,
                  windowHeight=100, histcolor="blue", axescolor="black",
                  xticknum=25, yticknum=25, showWidth=1, showTotal=1, showMax=1,
                  showValueLine=0, valueLine=" ", showPercentileLine=0,
                  percentileLine=" ", colors={}, numpts=100):
        
        typeStr = ("Min", "Max", "Average", "1st percentile", "5th percentile",
                   "10th percentile", "25th percentile", "50th percentile",
                   "75th percentile", "90th percentile", "95th percentile",
                   "99th percentile")

        graphType = self.dlg.graphType
        graphTypeStr = typeStr[graphType]

        # Setup the y-axis min/max based on which stat we are showing
        allkeys = list(graphData.keys())
        ymin = 999999.0
        ymax = 0.0
        for key in allkeys:
            value = graphData[key][graphType]
            ymin = min(ymin, value)
            ymax = max(ymax, value)
        
        # Include observation data in y-range calculation
        if self.showObs and self.obsData:
            for key in self.obsData.keys():
                value = self.obsData[key][graphType]
                ymin = min(ymin, value)
                ymax = max(ymax, value)
        
        if ymax == ymin:
            ymax = ymin + 1.0
        yrange = ymax - ymin
        if yrange < 0:
            yrange = 1.0

        gxmin = xmin
        gxmax = xmax

        # Find good number of tick marks in vertical
        ytick = max(self.niceNumDec(yrange / yticknum, 1), self.minTick)

        # Based on tick marks - set max/min
        gymax = min(((int(float(ymax) / float(ytick)) + 1) * ytick), self.absMax)
        gymin = max((int(float(ymin) / float(ytick)) * ytick), self.absMin)

        # Clear and redraw
        self.dlg.graphWidget.clear()
        
        # Setup graph coordinates
        self.dlg.graphWidget.setGraphCoords(gxmin, gxmax, gymin, gymax)

        # Draw axes
        self.dlg.graphWidget.drawAxes(gxmin, gxmax, gymin, gymax, ytick, label, self._gmtime)

        # Title
        self.dlg.graphWidget.drawTitle(title, graphTypeStr, numpts, 
                                        self.getPeriodString(self.perStart, self.perEnd))

        # Put the labels in the right order
        modOrder = self.getModelOrder(list(graphData.keys()))

        # Draw model data
        allkeys = sorted(list(graphData.keys()))
        allkeys.reverse()
        
        displayNum = -1
        for mod in modOrder:
            displayNum += 1
            lastx = -1
            lasty = -1
            if mod in colors:
                histcolor = colors[mod]
            else:
                histcolor = "blue"

            # Draw model label
            self.dlg.graphWidget.drawModelLabel(
                mod, displayNum, histcolor,
                self.state[mod] == 1,
                self.changeModelStatus,
                self.changeGroupStatus,
                self.postModelPopup
            )

            prevkey = "--"
            for key in allkeys:
                (modstr, basestr) = key.split(",")
                if mod != modstr:
                    continue
                base = int(basestr)
                value = graphData[key][graphType]
                (sx, sy) = self.dlg.graphWidget.graphcoord(base, value)

                if self.dlg.graphWidget.pointInGraph(sx, sy):
                    if self.pointState[mod] == 1 and self.state[mod] == 1:
                        dataLabelText = self.valueFormat % value
                        fontSizeStr = LabelSizeStrings[self.labelSize[mod]]
                        self.dlg.graphWidget.drawDataPoint(
                            sx, sy, histcolor, mod,
                            dataLabelText if self.labelState[mod] == 1 else None,
                            int(fontSizeStr)
                        )

                if lastx >= 0:
                    lw = self.lineWidth[mod]
                    full = fullData.get(prevkey, 1) == 1
                    if self.state[mod] == 1:
                        self.dlg.graphWidget.drawLine(
                            lastx, lasty, sx, sy,
                            histcolor, lw, full, mod
                        )

                lastx = sx
                lasty = sy
                prevkey = key

        # Draw observation line if enabled
        if self.showObs and self.obsData:
            self.drawObsLine(graphType)

        self.dlg.graphWidget.update()

    def drawObsLine(self, graphType):
        """Draw the observation data as a dashed line."""
        if not self.obsData:
            return
        
        obsConfig = ObsSourceConfig.get(self.obsSource, {})
        obsColor = obsConfig.get("color", "#00FF00")
        obsDash = obsConfig.get("dash", [6, 3])
        obsLabel = obsConfig.get("label", "Observed")
        
        # Sort observation keys by time
        obsKeys = sorted(list(self.obsData.keys()))
        
        lastx = -1
        lasty = -1
        prevkey = "--"
        
        for key in obsKeys:
            parts = key.split(",")
            if len(parts) >= 2:
                base = int(parts[1])
            else:
                continue
            
            value = self.obsData[key][graphType]
            (sx, sy) = self.dlg.graphWidget.graphcoord(base, value)
            
            if self.dlg.graphWidget.pointInGraph(sx, sy):
                # Draw observation point
                self.dlg.graphWidget.drawDataPoint(sx, sy, obsColor, "OBS", None, 10)
            
            if lastx >= 0:
                full = self.obsFull.get(prevkey, 1) == 1
                self.dlg.graphWidget.drawLine(
                    lastx, lasty, sx, sy,
                    obsColor, 2, False, "OBS",  # Always dashed for obs
                    dash_pattern=obsDash
                )
            
            lastx = sx
            lasty = sy
            prevkey = key
        
        # Add observation source label
        self.dlg.graphWidget.drawObsLabel(obsLabel, obsColor)

    def getPeriodString(self, perStart, perEnd):
        DAYS = ("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
        MONTHS = ("Jan", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
        (syea, smon, sday, shou, smin, ssec, swda, syda, sdst) = self._gmtime(perStart).timetuple()
        (eyea, emon, eday, ehou, emin, esec, ewda, eyda, edst) = self._gmtime(perEnd).timetuple()
        if self.element in ("MaxT", "MinT", "MaxRH", "MinRH", "TdMrn", "TdAft"):
            if (eyda == syda) or ((eyda == (syda + 1)) and (self.element in ("MaxT", "MinRH", "TdAft"))):
                subText = "%s %s %d, %4.4d" % (DAYS[swda], MONTHS[smon], sday, syea)
            else:
                subText = "%s %s %d -- %s %s %d, %4.4d" % (DAYS[swda], MONTHS[smon], sday, DAYS[ewda], MONTHS[emon], eday, syea)
        elif (ewda == swda) or ((eyda == (syda + 1)) and (ehou == 0)):
            subText = "%2.2dZ - %2.2dZ %s %s %d, %4.4d" % (shou, ehou, DAYS[swda], MONTHS[smon], sday, syea)
        else:
            subText = "%4.4d/%2.2d/%2.2d %2.2dZ --- %4.4d/%2.2d/%2.2d %2.2dZ" % (syea, smon, sday, shou, eyea, emon, eday, ehou)
        return subText

    def changeModelStatus(self, modelStr):
        if self.state[modelStr] == 1:
            self.state[modelStr] = 0
        else:
            self.state[modelStr] = 1
        self.drawGraph(
            self.outdata, self.outfull,
            max(self.graphStart, self.mintime), self.maxtime,
            self.graphmin, self.allmax,
            label=self.labelText, title=self.titleText,
            colors=self.COLORS, numpts=self.numpts
        )

    def changeGroupStatus(self, modelStr):
        group = self.primaryGroup[modelStr]
        newState = 0 if self.state[modelStr] == 1 else 1
        for mod in self.groups[group]:
            self.state[mod] = newState
        self.drawGraph(
            self.outdata, self.outfull,
            max(self.graphStart, self.mintime), self.maxtime,
            self.graphmin, self.allmax,
            label=self.labelText, title=self.titleText,
            colors=self.COLORS, numpts=self.numpts
        )

    def postModelPopup(self, modelStr, pos):
        self.dlg.showModelPopup(modelStr, pos)

    def popupHandler(self, changeType, mod=None):
        if mod is None:
            mod = self.dlg.currentPopupModel
        
        if changeType == "LineChange":
            lw = self.dlg.popupLineWidth
            self.lineWidth[mod] = lw
        elif changeType == "ColorChange":
            color = QColorDialog.getColor(QColor(self.COLORS[mod]), self.dlg, "%s Color" % mod)
            if color.isValid():
                self.COLORS[mod] = color.name()
        elif changeType == "DataChange":
            self.pointState[mod] = 1 if self.pointState[mod] == 0 else 0
        elif changeType == "LabelChange":
            self.labelState[mod] = 1 if self.labelState[mod] == 0 else 0
        elif changeType == "LabelSizeChange":
            self.labelSize[mod] = self.dlg.popupLabelSize
        elif changeType == "TimeSplit":
            xtime = self.dlg.timeSplitX
            xr = 6 * (60 * 60)
            xround = (int(xtime / xr)) * xr
            self.graphStart = xround
        elif changeType == "TimeReset":
            self.graphStart = self.mintime
        elif changeType == "ObsSourceChange":
            self.obsSource = self.dlg.obsSourceCombo.currentText()
            self.showObs = (self.obsSource != "None")
            # Trigger data refresh
            self.getStuff("Update")
            return

        self.drawGraph(
            self.outdata, self.outfull,
            max(self.graphStart, self.mintime), self.maxtime,
            self.graphmin, self.allmax,
            label=self.labelText, title=self.titleText,
            colors=self.COLORS, numpts=self.numpts
        )

    def getModelOrder(self, allkeys):
        uniqueModels = []
        for key in allkeys:
            (modstr, basestr) = key.split(",")
            if modstr not in uniqueModels:
                uniqueModels.append(modstr)
        outOrder = []
        if "Fcst" in uniqueModels:
            outOrder.append("Fcst")
            uniqueModels.remove("Fcst")
        for mod in self.MODELS:
            if mod in uniqueModels:
                outOrder.append(mod)
        return outOrder

    def saveStatus(self):
        print("inside saveStatus")
        sys.stdout.flush()
        self.saveObject("state", self.state, "TrendsTool")
        self.saveObject("lineWidth", self.lineWidth, "TrendsTool")
        self.saveObject("pointState", self.pointState, "TrendsTool")
        self.saveObject("labelState", self.labelState, "TrendsTool")
        self.saveObject("labelSize", self.labelSize, "TrendsTool")
        self.saveObject("colors", self.COLORS, "TrendsTool")
        self.saveObject("obsSource", self.obsSource, "TrendsTool")
        return

    def readStatus(self):
        print("inside readStatus")
        
        lf = self.getUserFile("state", "TrendsTool")
        if lf.exists():
            tempstate = self.getObject("state", "TrendsTool")
            for key in tempstate:
                if key in self.state:
                    self.state[key] = tempstate[key]

        lf = self.getUserFile("lineWidth", "TrendsTool")
        if lf.exists():
            templw = self.getObject("lineWidth", "TrendsTool")
            for key in templw:
                if key in self.lineWidth:
                    self.lineWidth[key] = templw[key]

        lf = self.getUserFile("pointState", "TrendsTool")
        if lf.exists():
            tempps = self.getObject("pointState", "TrendsTool")
            for key in tempps:
                if key in self.pointState:
                    self.pointState[key] = tempps[key]

        lf = self.getUserFile("labelState", "TrendsTool")
        if lf.exists():
            templs = self.getObject("labelState", "TrendsTool")
            for key in templs:
                if key in self.labelState:
                    self.labelState[key] = templs[key]

        lf = self.getUserFile("labelSize", "TrendsTool")
        if lf.exists():
            templs = self.getObject("labelSize", "TrendsTool")
            for key in templs:
                if key in self.labelSize:
                    self.labelSize[key] = templs[key]

        lf = self.getUserFile("colors", "TrendsTool")
        if lf.exists():
            tempcolors = self.getObject("colors", "TrendsTool")
            for key in tempcolors:
                if key in self.COLORS:
                    self.COLORS[key] = tempcolors[key]

        lf = self.getUserFile("obsSource", "TrendsTool")
        if lf.exists():
            self.obsSource = self.getObject("obsSource", "TrendsTool")
            self.showObs = (self.obsSource != "None")

        sys.stdout.flush()
        return

    def getUserFile(self, name, category):
        from com.raytheon.uf.common.localization import PathManagerFactory, LocalizationContext
        LocalizationType = LocalizationContext.LocalizationType
        LocalizationLevel = LocalizationContext.LocalizationLevel
        pathMgr = PathManagerFactory.getPathManager()
        path = 'gfe/userPython/' + category + '/' + name
        lc = pathMgr.getContext(LocalizationType.valueOf('CAVE_STATIC'), LocalizationLevel.valueOf('USER'))
        lf = pathMgr.getLocalizationFile(lc, path)
        return lf

    def niceNumDec(self, val, roundit):
        if val == 0:
            return 1
        e = math.floor(math.log10(val))
        a = 10.0 ** e
        f = val / a
        if roundit > 0:
            if f < 1.5:
                nf = 1
            elif f < 3.0:
                nf = 2
            elif f < 7.0:
                nf = 5
            else:
                nf = 10
        else:
            if f <= 1.0:
                nf = 1
            elif f <= 2.0:
                nf = 2.0
            elif f <= 5.0:
                nf = 5.0
            else:
                nf = 10.0
        return nf * a



# Cache for dynamically created UI classes
_GraphWidget = None
_TrendDialog = None


def _create_ui_classes():
    """Create UI classes after PySide6 is imported to avoid Jep conflicts."""
    global _GraphWidget, _TrendDialog
    
    if _GraphWidget is not None and _TrendDialog is not None:
        return _GraphWidget, _TrendDialog
    
    # Import PySide6
    QtWidgets, QtCore, QtGui = _import_pyside6()
    
    class GraphWidget(QtWidgets.QWidget):
        """Custom widget for drawing the trend graph using QPainter."""
        
        def __init__(self, parent=None):
            super().__init__(parent)
            self.setMinimumSize(800, 600)
            self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
            self.setMouseTracking(True)
            
            # Graph coordinate system
            self.xmin = 0
            self.xmax = 1
            self.ymin = 0
            self.ymax = 1
            self.xmult = 1
            self.xoff = 0
            self.ymult = 1
            self.yoff = 0
            
            # Margins
            self.left = 60
            self.right = 165
            self.top = 70
            self.bottom = 70
            
            # Drawing lists
            self.lines = []
            self.points = []
            self.texts = []
            self.labels = []
            self.axes_items = []
            self.obs_items = []
            
            # Callbacks
            self.model_click_callback = None
            self.model_right_click_callback = None
            self.model_middle_click_callback = None
            
            # Model label positions for click detection
            self.model_label_rects = {}
            
            # Cursor readout
            self.show_readout = False
            self.readout_pos = QtCore.QPointF(0, 0)
            self.readout_value = ""
            self.valueFormat = "%.1f"
            
            # Time conversion function
            self._gmtime = None
            
        def setValueFormat(self, fmt):
            self.valueFormat = fmt
            
        def setGraphCoords(self, xmin, xmax, ymin, ymax):
            self.xmin = xmin
            self.xmax = xmax
            self.ymin = ymin
            self.ymax = ymax
            
            width = self.width()
            height = self.height()
            
            sxmin = self.left
            sxmax = width - self.right
            symin = self.bottom
            symax = height - self.top
            
            if xmax != xmin:
                self.xmult = (sxmax - sxmin) / (xmax - xmin)
            else:
                self.xmult = 1
            self.xoff = sxmin
            
            if ymax != ymin:
                self.ymult = (symax - symin) / (ymax - ymin)
            else:
                self.ymult = 1
            self.yoff = symin
            self.windowHeight = height
            
        def graphcoord(self, x, y):
            sx = ((x - self.xmin) * self.xmult) + self.xoff
            sy = self.windowHeight - (((y - self.ymin) * self.ymult) + self.yoff)
            return sx, sy
        
        def valcoord(self, sx, sy):
            x = ((sx - self.xoff) / self.xmult) + self.xmin
            y = ((-sy + self.windowHeight - self.yoff) / self.ymult) + self.ymin
            return x, y
        
        def pointInGraph(self, sx, sy):
            width = self.width()
            height = self.height()
            return (sx >= self.left and sx <= width - self.right and
                    sy >= self.top and sy <= height - self.bottom)
        
        def clear(self):
            self.lines = []
            self.points = []
            self.texts = []
            self.labels = []
            self.axes_items = []
            self.obs_items = []
            self.model_label_rects = {}
            
        def drawAxes(self, xmin, xmax, ymin, ymax, ytick, label, gmtime_func):
            self._gmtime = gmtime_func
            width = self.width()
            height = self.height()
            
            # Draw box around graph
            (nx, ny) = self.graphcoord(xmin, ymin)
            (xx, xy) = self.graphcoord(xmax, ymax)
            self.axes_items.append(('rect', nx, xy, xx - nx, ny - xy, 'black'))
            
            # Y-axis ticks and labels
            numticks = int((ymax - ymin) / ytick) + 1
            labeldigits = max(-math.floor(math.log10(ytick)) if ytick > 0 else 0, 0)
            
            for i in range(numticks):
                y = ymin + (i * ytick)
                (tx, ty) = self.graphcoord(xmin, y)
                (tx2, ty2) = self.graphcoord(xmax, y)
                
                # Ticks on left
                self.axes_items.append(('line', tx - 3, ty, tx + 3, ty, 'black', 1))
                # Ticks on right  
                self.axes_items.append(('line', tx2 - 3, ty2, tx2 + 3, ty2, 'black', 1))
                
                # Labels
                if i % 3 == 0:
                    if labeldigits == 0:
                        labelstring = "%d" % y
                    else:
                        fmt = "%%.%df" % labeldigits
                        labelstring = fmt % y
                    self.axes_items.append(('text', tx - 7, ty, labelstring, 'black', 'right'))
                    self.axes_items.append(('text', tx2 + 7, ty2, labelstring, 'black', 'left'))
            
            # X-axis ticks and labels (time-based)
            xtick = 6 * 60 * 60  # 6 hours
            numticks = int((xmax - xmin) / xtick) + 1
            
            for i in range(numticks):
                x = xmin + (i * xtick)
                (tx, ty) = self.graphcoord(x, ymin)
                (tx2, ty2) = self.graphcoord(x, ymax)
                
                if self._gmtime:
                    (xyea, xmon, xday, xhou, xmin_t, xsec, xwda, xyda, xdst) = self._gmtime(x).timetuple()
                    
                    # Longer tick at midnight
                    ticklen = 5 if xhou == 0 else 3
                    self.axes_items.append(('line', tx, ty - ticklen, tx, ty + ticklen, 'black', 1))
                    self.axes_items.append(('line', tx2, ty2 - ticklen, tx2, ty2 + ticklen, 'black', 1))
                    
                    # Grid line at midnight
                    if xhou == 0 and i > 0 and i < numticks - 1:
                        self.axes_items.append(('line', tx, ty, tx, ty2, '#A0A0A0', 1))
                    
                    # Date label at noon
                    if xhou == 12:
                        labelstring = "%2.2d/%2.2d" % (xmon, xday)
                        self.axes_items.append(('text', tx, ty + 15, labelstring, 'black', 'center'))
            
            # Y-axis label (vertical)
            (xdum, midy) = self.graphcoord(xmin, (ymin + ymax) / 2.0)
            self.axes_items.append(('vtext', nx - 40, midy, label, 'black'))
            self.axes_items.append(('vtext', xx + 40, midy, label, 'black'))
            
            # X-axis label
            (midx, ydum) = self.graphcoord((xmax + xmin) / 2.0, 0)
            self.axes_items.append(('text', midx, ny + 35, "Model Run", 'black', 'center'))
            
        def drawTitle(self, title, graphTypeStr, numpts, periodString):
            width = self.width()
            xmid = width / 2
            
            self.texts.append(('title', xmid, 10, title, 'black', 24, True))
            
            if numpts > 1:
                pointText = "(%s of %d points)" % (graphTypeStr, numpts)
                self.texts.append(('subtitle', xmid, 36, pointText, 'black', 10, False))
                yoff = 12
            else:
                yoff = 0
            
            self.texts.append(('period', xmid, 36 + yoff, periodString, 'black', 12, False))
            
        def drawModelLabel(self, mod, displayNum, color, visible, click_cb, middle_cb, right_cb):
            width = self.width()
            x = width - 115
            y = 50 + (displayNum * 16)
            
            textColor = color if visible else 'grey'
            self.labels.append((mod, x, y, textColor))
            
            # Store click area
            self.model_label_rects[mod] = QtCore.QRectF(x, y - 8, 100, 16)
            
            # Store callbacks
            self.model_click_callback = click_cb
            self.model_middle_click_callback = middle_cb
            self.model_right_click_callback = right_cb
            
        def drawDataPoint(self, sx, sy, color, tag, label_text=None, font_size=10):
            self.points.append((sx, sy, color, tag, label_text, font_size))
            
        def drawLine(self, x1, y1, x2, y2, color, width, solid, tag, dash_pattern=None):
            # Clip line to graph area
            w = self.width()
            h = self.height()
            left = self.left
            right = w - self.right
            top = self.top
            bottom = h - self.bottom
            
            clipped = self.clipLine(left, top, right, bottom, x1, y1, x2, y2)
            if clipped[0] is not None:
                self.lines.append((clipped[0], clipped[1], clipped[2], clipped[3], 
                                  color, width, solid, tag, dash_pattern))
        
        def drawObsLabel(self, label, color):
            """Draw observation source label in the legend area."""
            width = self.width()
            x = width - 115
            # Find the next available y position after model labels
            y = 50 + (len(self.labels) * 16) + 20
            self.obs_items.append(('label', x, y, label, color))
        
        def clipLine(self, left, top, right, bottom, x1, y1, x2, y2):
            """Liang-Barsky line clipping algorithm."""
            dx = x2 - x1
            dy = y2 - y1
            dt0, dt1 = 0, 1

            checks = ((-dx, x1 - left),
                      (dx, right - x1),
                      (-dy, y1 - top),
                      (dy, bottom - y1))

            for p, q in checks:
                if p == 0:
                    if q < 0:
                        return None, None, None, None
                    continue
                dt = q / (p * 1.0)
                if p < 0:
                    if dt > dt1:
                        return None, None, None, None
                    dt0 = max(dt0, dt)
                else:
                    if dt < dt0:
                        return None, None, None, None
                    dt1 = min(dt1, dt)
            if dt0 > 0:
                x1 += dt0 * dx
                y1 += dt0 * dy
            if dt1 < 1:
                x2 = x1 + dt1 * dx
                y2 = y1 + dt1 * dy
            return x1, y1, x2, y2
        
        def paintEvent(self, event):
            painter = QtGui.QPainter(self)
            painter.setRenderHint(QtGui.QPainter.Antialiasing)
            
            # White background
            painter.fillRect(self.rect(), QtCore.Qt.white)
            
            # Draw axes items
            for item in self.axes_items:
                if item[0] == 'rect':
                    _, x, y, w, h, color = item
                    painter.setPen(QtGui.QPen(QtGui.QColor(color), 1))
                    painter.drawRect(int(x), int(y), int(w), int(h))
                elif item[0] == 'line':
                    _, x1, y1, x2, y2, color, width = item
                    painter.setPen(QtGui.QPen(QtGui.QColor(color), width))
                    painter.drawLine(int(x1), int(y1), int(x2), int(y2))
                elif item[0] == 'text':
                    _, x, y, text, color, align = item
                    painter.setPen(QtGui.QColor(color))
                    font = QtGui.QFont("Helvetica", 10)
                    painter.setFont(font)
                    fm = painter.fontMetrics()
                    tw = fm.horizontalAdvance(text)
                    if align == 'right':
                        painter.drawText(int(x - tw), int(y + 4), text)
                    elif align == 'center':
                        painter.drawText(int(x - tw / 2), int(y + 4), text)
                    else:
                        painter.drawText(int(x), int(y + 4), text)
                elif item[0] == 'vtext':
                    _, x, y, text, color = item
                    painter.save()
                    painter.translate(x, y)
                    painter.rotate(-90)
                    painter.setPen(QtGui.QColor(color))
                    font = QtGui.QFont("Helvetica", 10)
                    painter.setFont(font)
                    fm = painter.fontMetrics()
                    tw = fm.horizontalAdvance(text)
                    painter.drawText(int(-tw / 2), 0, text)
                    painter.restore()
            
            # Draw lines
            for line in self.lines:
                x1, y1, x2, y2, color, width, solid, tag, dash_pattern = line
                pen = QtGui.QPen(QtGui.QColor(color), width)
                if not solid:
                    if dash_pattern:
                        pen.setDashPattern(dash_pattern)
                    else:
                        pen.setStyle(QtCore.Qt.DashLine)
                painter.setPen(pen)
                painter.drawLine(int(x1), int(y1), int(x2), int(y2))
            
            # Draw points
            for point in self.points:
                sx, sy, color, tag, label_text, font_size = point
                painter.setPen(QtGui.QPen(QtGui.QColor(color), 1))
                painter.setBrush(QtGui.QBrush(QtGui.QColor(color)))
                painter.drawRect(int(sx - 2), int(sy - 2), 4, 4)
                
                if label_text:
                    font = QtGui.QFont("Helvetica", font_size)
                    painter.setFont(font)
                    painter.drawText(int(sx), int(sy - 5), label_text)
            
            # Draw titles
            for text in self.texts:
                ttype, x, y, txt, color, size, bold = text
                font = QtGui.QFont("Helvetica", size)
                if bold:
                    font.setBold(True)
                painter.setFont(font)
                painter.setPen(QtGui.QColor(color))
                fm = painter.fontMetrics()
                tw = fm.horizontalAdvance(txt)
                painter.drawText(int(x - tw / 2), int(y + size), txt)
            
            # Draw model labels
            for label in self.labels:
                mod, x, y, color = label
                font = QtGui.QFont("Helvetica", 11)
                painter.setFont(font)
                painter.setPen(QtGui.QColor(color))
                painter.drawText(int(x), int(y), mod)
            
            # Draw observation items
            for item in self.obs_items:
                if item[0] == 'label':
                    _, x, y, text, color = item
                    font = QtGui.QFont("Helvetica", 11)
                    font.setItalic(True)
                    painter.setFont(font)
                    painter.setPen(QtGui.QColor(color))
                    # Draw a dashed line sample
                    pen = QtGui.QPen(QtGui.QColor(color), 2)
                    pen.setStyle(QtCore.Qt.DashLine)
                    painter.setPen(pen)
                    painter.drawLine(int(x - 30), int(y - 4), int(x - 5), int(y - 4))
                    # Draw label
                    painter.setPen(QtGui.QColor(color))
                    painter.drawText(int(x), int(y), text)
            
            # Draw cursor readout
            if self.show_readout:
                font = QtGui.QFont("Helvetica", 12)
                font.setBold(True)
                painter.setFont(font)
                painter.setPen(QtGui.QColor('black'))
                painter.drawText(int(self.readout_pos.x() + 5), int(self.readout_pos.y() - 5), self.readout_value)
            
        def mousePressEvent(self, event):
            pos = event.position()
            
            # Check if click is on a model label
            for mod, rect in self.model_label_rects.items():
                if rect.contains(pos):
                    if event.button() == QtCore.Qt.LeftButton and self.model_click_callback:
                        self.model_click_callback(mod)
                    elif event.button() == QtCore.Qt.MiddleButton and self.model_middle_click_callback:
                        self.model_middle_click_callback(mod)
                    elif event.button() == QtCore.Qt.RightButton and self.model_right_click_callback:
                        self.model_right_click_callback(mod, event.globalPosition().toPoint())
                    return
            
            # Check if click is in graph area
            if self.pointInGraph(pos.x(), pos.y()):
                if event.button() == QtCore.Qt.LeftButton:
                    self.show_readout = True
                    self.readout_pos = pos
                    (x, y) = self.valcoord(pos.x(), pos.y())
                    self.readout_value = self.valueFormat % y
                    self.update()
        
        def mouseMoveEvent(self, event):
            pos = event.position()
            
            # Update cursor based on position
            in_label = False
            for mod, rect in self.model_label_rects.items():
                if rect.contains(pos):
                    self.setCursor(QtCore.Qt.PointingHandCursor)
                    in_label = True
                    break
            
            if not in_label:
                if self.pointInGraph(pos.x(), pos.y()):
                    self.setCursor(QtCore.Qt.CrossCursor)
                else:
                    self.setCursor(QtCore.Qt.ArrowCursor)
            
            # Update readout if showing
            if self.show_readout:
                if self.pointInGraph(pos.x(), pos.y()):
                    self.readout_pos = pos
                    (x, y) = self.valcoord(pos.x(), pos.y())
                    self.readout_value = self.valueFormat % y
                else:
                    self.show_readout = False
                self.update()
        
        def mouseReleaseEvent(self, event):
            if event.button() == QtCore.Qt.LeftButton:
                self.show_readout = False
                self.update()
        
        def resizeEvent(self, event):
            super().resizeEvent(event)
            # Recalculate graph coordinates on resize
            self.setGraphCoords(self.xmin, self.xmax, self.ymin, self.ymax)


    class TrendDialog(QtWidgets.QDialog):
        """Main dialog window for the Trends tool using PySide6."""
        
        def __init__(self, title="Model Trends", callback_method=None, popup_callback=None, parent_tool=None):
            super().__init__()
            
            self.setWindowTitle(title)
            self.setMinimumSize(1000, 700)
            self.resize(1200, 800)
            
            self.callback_method = callback_method
            self.popup_callback = popup_callback
            self.parent_tool = parent_tool
            
            # Graph type (0=Min, 1=Max, 2=Avg, etc.)
            self.graphType = 2
            
            # Popup menu state
            self.currentPopupModel = ""
            self.popupLineWidth = 1
            self.popupLabelSize = 1
            self.timeSplitX = 0
            
            self.setupUI()
            
        def setupUI(self):
            layout = QtWidgets.QVBoxLayout(self)
            
            # Menu bar
            menuBar = QtWidgets.QMenuBar(self)
            layout.setMenuBar(menuBar)
            
            # Graph menu
            graphMenu = menuBar.addMenu("Graph")
            
            graphTypeGroup = QtGui.QActionGroup(self)
            graphTypeGroup.setExclusive(True)
            
            graphTypes = [
                ("Max", 1), ("99th Percentile", 11), ("95th Percentile", 10),
                ("90th Percentile", 9), ("75th Percentile", 8), ("50th Percentile", 7),
                ("Avg", 2), ("25th Percentile", 6), ("10th Percentile", 5),
                ("5th Percentile", 4), ("1st Percentile", 3), ("Min", 0)
            ]
            
            for name, value in graphTypes:
                action = QtGui.QAction(name, self)
                action.setCheckable(True)
                action.setData(value)
                if value == 2:  # Default to Avg
                    action.setChecked(True)
                action.triggered.connect(lambda checked, v=value: self.setGraphType(v))
                graphTypeGroup.addAction(action)
                graphMenu.addAction(action)
            
            # Observation source selection in menu bar
            obsMenu = menuBar.addMenu("Observations")
            
            obsGroup = QtGui.QActionGroup(self)
            obsGroup.setExclusive(True)
            
            obsSources = ["None", "Obs", "URMA", "RTMA"]
            for source in obsSources:
                action = QtGui.QAction(source, self)
                action.setCheckable(True)
                action.setData(source)
                if source == "None":
                    action.setChecked(True)
                action.triggered.connect(lambda checked, s=source: self.setObsSource(s))
                obsGroup.addAction(action)
                obsMenu.addAction(action)
            
            # Also add a combo box for quick access
            obsLayout = QtWidgets.QHBoxLayout()
            obsLayout.addWidget(QtWidgets.QLabel("Observation Source:"))
            self.obsSourceCombo = QtWidgets.QComboBox()
            self.obsSourceCombo.addItems(["None", "Obs", "URMA", "RTMA"])
            self.obsSourceCombo.currentTextChanged.connect(self.onObsSourceChanged)
            obsLayout.addWidget(self.obsSourceCombo)
            obsLayout.addStretch()
            
            obsWidget = QtWidgets.QWidget()
            obsWidget.setLayout(obsLayout)
            layout.addWidget(obsWidget)
            
            # Graph widget
            self.graphWidget = GraphWidget(self)
            layout.addWidget(self.graphWidget, 1)
            
            # Button bar
            buttonLayout = QtWidgets.QHBoxLayout()
            buttonLayout.addStretch()
            
            self.updateButton = QtWidgets.QPushButton("Update")
            self.updateButton.setFixedWidth(100)
            self.updateButton.clicked.connect(self.onUpdate)
            buttonLayout.addWidget(self.updateButton)
            
            closeButton = QtWidgets.QPushButton("Close")
            closeButton.setFixedWidth(100)
            closeButton.clicked.connect(self.onClose)
            buttonLayout.addWidget(closeButton)
            
            buttonLayout.addStretch()
            layout.addLayout(buttonLayout)
            
            # Context menus
            self.createContextMenus()
            
        def createContextMenus(self):
            # Model popup menu
            self.modelPopup = QtWidgets.QMenu(self)
            self.modelNameAction = self.modelPopup.addAction("Model Name")
            self.modelNameAction.setEnabled(False)
            self.modelPopup.addSeparator()
            
            # Line width submenu
            lineWidthMenu = self.modelPopup.addMenu("Line Width")
            self.lineWidthGroup = QtGui.QActionGroup(self)
            for w in [1, 2, 3]:
                action = QtGui.QAction("%d px" % w, self)
                action.setCheckable(True)
                action.setData(w)
                if w == 1:
                    action.setChecked(True)
                action.triggered.connect(lambda checked, width=w: self.onLineWidthChange(width))
                self.lineWidthGroup.addAction(action)
                lineWidthMenu.addAction(action)
            
            # Color action
            colorAction = self.modelPopup.addAction("Color...")
            colorAction.triggered.connect(self.onColorChange)
            
            # Data points toggle
            self.dataPointsAction = self.modelPopup.addAction("Show Data Points")
            self.dataPointsAction.setCheckable(True)
            self.dataPointsAction.setChecked(True)
            self.dataPointsAction.triggered.connect(self.onDataPointsChange)
            
            # Data labels toggle
            self.dataLabelsAction = self.modelPopup.addAction("Show Data Labels")
            self.dataLabelsAction.setCheckable(True)
            self.dataLabelsAction.setChecked(False)
            self.dataLabelsAction.triggered.connect(self.onDataLabelsChange)
            
            # Label size submenu
            labelSizeMenu = self.modelPopup.addMenu("Label Size")
            self.labelSizeGroup = QtGui.QActionGroup(self)
            for i, size in enumerate(LabelSizeStrings):
                action = QtGui.QAction("%s pt" % size, self)
                action.setCheckable(True)
                action.setData(i)
                if i == 1:
                    action.setChecked(True)
                action.triggered.connect(lambda checked, idx=i: self.onLabelSizeChange(idx))
                self.labelSizeGroup.addAction(action)
                labelSizeMenu.addAction(action)
            
            # X-axis popup menu
            self.xAxisPopup = QtWidgets.QMenu(self)
            timeSplitAction = self.xAxisPopup.addAction("Start Time Axis Here")
            timeSplitAction.triggered.connect(self.onTimeSplit)
            timeResetAction = self.xAxisPopup.addAction("Reset Time Axis to all times")
            timeResetAction.triggered.connect(self.onTimeReset)
            
        def setGraphType(self, value):
            self.graphType = value
            if self.callback_method:
                self.callback_method("Redraw")
                
        def setObsSource(self, source):
            self.obsSourceCombo.setCurrentText(source)
            
        def onObsSourceChanged(self, source):
            if self.popup_callback:
                self.popup_callback("ObsSourceChange")
        
        def showModelPopup(self, modelStr, pos):
            self.currentPopupModel = modelStr
            
            # Update menu state
            if self.parent_tool:
                color = self.parent_tool.COLORS.get(modelStr, "#000000")
                self.modelNameAction.setText(modelStr)
                
                # Update line width selection
                lw = self.parent_tool.lineWidth.get(modelStr, 1)
                for action in self.lineWidthGroup.actions():
                    if action.data() == lw:
                        action.setChecked(True)
                
                # Update data points state
                self.dataPointsAction.setChecked(self.parent_tool.pointState.get(modelStr, 1) == 1)
                
                # Update data labels state
                self.dataLabelsAction.setChecked(self.parent_tool.labelState.get(modelStr, 0) == 1)
                
                # Update label size
                ls = self.parent_tool.labelSize.get(modelStr, 1)
                for action in self.labelSizeGroup.actions():
                    if action.data() == ls:
                        action.setChecked(True)
            
            self.modelPopup.popup(pos)
            
        def onLineWidthChange(self, width):
            self.popupLineWidth = width
            if self.popup_callback:
                self.popup_callback("LineChange", self.currentPopupModel)
                
        def onColorChange(self):
            if self.popup_callback:
                self.popup_callback("ColorChange", self.currentPopupModel)
                
        def onDataPointsChange(self):
            if self.popup_callback:
                self.popup_callback("DataChange", self.currentPopupModel)
                
        def onDataLabelsChange(self):
            if self.popup_callback:
                self.popup_callback("LabelChange", self.currentPopupModel)
                
        def onLabelSizeChange(self, idx):
            self.popupLabelSize = idx
            if self.popup_callback:
                self.popup_callback("LabelSizeChange", self.currentPopupModel)
                
        def onTimeSplit(self):
            if self.popup_callback:
                self.popup_callback("TimeSplit")
                
        def onTimeReset(self):
            if self.popup_callback:
                self.popup_callback("TimeReset")
        
        def onUpdate(self):
            if self.callback_method:
                self.callback_method("Update")
                
        def onClose(self):
            if self.callback_method:
                self.callback_method("Close")
            self.accept()
            
        def resizeEvent(self, event):
            super().resizeEvent(event)
            if self.callback_method:
                self.callback_method("Resize")

    # Cache the classes
    _GraphWidget = GraphWidget
    _TrendDialog = TrendDialog
    
    return GraphWidget, TrendDialog


def GraphWidget(*args, **kwargs):
    """Factory function for GraphWidget - ensures PySide6 is imported first."""
    GraphWidgetClass, _ = _create_ui_classes()
    return GraphWidgetClass(*args, **kwargs)


def TrendDialog(*args, **kwargs):
    """Factory function for TrendDialog - ensures PySide6 is imported first."""
    _, TrendDialogClass = _create_ui_classes()
    return TrendDialogClass(*args, **kwargs)
