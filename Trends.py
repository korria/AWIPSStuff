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
# 2.0 - Added observation overlay support (Obs, URMA, RTMA) with dashed lines
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

# Observation database options - color and dash pattern for each source
ObsSourceConfig = {
    "None": {"color": None, "dash": None, "label": None},
    "Obs": {"color": "#00DD00", "dash": (6, 4), "label": "Observed"},
    "URMA": {"color": "#FF6600", "dash": (6, 4), "label": "URMA Analysis"},
    "RTMA": {"color": "#9900FF", "dash": (6, 4), "label": "RTMA Analysis"},
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
import AppDialog

try:
    import tkinter as tk
except ImportError:
    import Tkinter as tk
try:
    import tkinter.colorchooser as tkColor
except ImportError:
    import tkColorChooser as tkColor

import TimeRange


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
        self.obsSource = "None"
        self.obsData = {}
        self.obsFull = {}

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
        self.dlg = ToolDialog(
            "Model Trends",
            callbackMethod=self.getStuff,
            popupCallback=self.popupHandler
        )
        self.dlg.histCanv.bind("<Configure>", self.resizecanvas)
        self.dlg.obsSource.set(self.obsSource)
        
        self.getStuff("Update")
        self.dlg.mainloop()
        self.cancel()

    def execute(self):
        "Model/Forecast trend over time"
        return

    def getStuff(self, buttonName):
        if buttonName == "Close":
            self.saveStatus()
            return
        
        if buttonName == "Resize":
            self.drawGraph(self.outdata, self.outfull,
                max(self.graphStart, self.mintime), self.maxtime,
                self.graphmin, self.allmax,
                label=self.labelText, title=self.titleText,
                colors=self.COLORS, numpts=self.numpts)
            return
        
        if buttonName == "Redraw":
            self.drawGraph(self.outdata, self.outfull,
                max(self.graphStart, self.mintime), self.maxtime,
                self.graphmin, self.allmax,
                label=self.labelText, title=self.titleText,
                colors=self.COLORS, numpts=self.numpts)
            return

        # Disable the Update button while working
        self.dlg.updateButton.configure(text="WORKING")
        normalBG = self.dlg.updateButton.cget("bg")
        self.dlg.updateButton.configure(bg="red")
        normalFG = self.dlg.updateButton.cget("fg")
        self.dlg.updateButton.configure(fg="white")
        self.dlg.updateButton.configure(state=tk.DISABLED)
        self.dlg.update_idletasks()

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
            self.dlg.updateButton.configure(state=tk.NORMAL)
            self.dlg.updateButton.configure(text="Update")
            self.dlg.updateButton.configure(bg=normalBG)
            self.dlg.updateButton.configure(fg=normalFG)
            return
        
        (selectedName, selectedLevel, selectModelID) = selectedParms[0]
        parmName = selectedName

        selectTR = TimeRange.TimeRange(self._dbss.getParmOp().getSelectionTimeRange())
        mutid = self.mutableID()

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
            self.dlg.updateButton.configure(state=tk.NORMAL)
            self.dlg.updateButton.configure(text="Update")
            self.dlg.updateButton.configure(bg=normalBG)
            self.dlg.updateButton.configure(fg=normalFG)
            return

        gridMax = gridInfos[0].maxLimit()
        self.absMax = gridMax * 20.0 if accumVal == 1 else gridMax
        self.absMin = gridInfos[0].minLimit()
        inUnits = gridInfos[0].units()
        inPrecision = gridInfos[0].precision()

        self.minTick = 1 if inPrecision < 1 else 10.0 ** (-inPrecision)
        self.labelText = inUnits
        self.valueFormat = "%.0f" if inPrecision == 0 else "%." + ("%d" % inPrecision) + "f"

        parmObj = self.getParm(mutid, parmName, "SFC")
        if parmObj is not None:
            from com.raytheon.viz.gfe.rsc import DiscreteDisplayUtil
            ctInfo = DiscreteDisplayUtil.buildColorMapParameters(parmObj)

        lev = "SFC"
        self.outdata = {}
        self.outfull = {}
        self.allmin = 10000.0
        self.allmax = 0.0
        self.mintime = self.perStart
        self.maxtime = self.perEnd

        if (self.saveElement != reqElement) or (self.savePStart != self.perStart) or (self.savePEnd != self.perEnd):
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
            maxvers = -1 if mod in ("Fcst", "ISC", "NBM", "NBMEXP") else -40
            element = reqElement
            self.element = reqElement
            tempoff = 0
            if self.state[mod] == 0:
                tempoff = 1
                self.state[mod] = 1

            for vers in range(0, maxvers, -1):
                modID = self.findDatabase(mod, version=vers)
                if modID is None or not modID.isValid():
                    break

                parmObj = self.getParm(modID, element, lev)
                if parmObj is None:
                    continue

                if mod not in ("Fcst", "ISC", "NBM", "NBMEXP"):
                    modTimeUnix = modID.modelTime().unixTime()
                else:
                    modTimeUnix = self.perStart if self._gmtime().unixTime() > self.perStart else self._gmtime().unixTime()

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
                    fullFlag = 1 if totTime >= fullDuration or reqElement not in ("QPF", "SnowAmt") else 0

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
                    self.drawGraph(self.outdata, self.outfull,
                        max(self.graphStart, self.mintime), self.maxtime,
                        self.graphmin, self.allmax,
                        label=self.labelText, title=self.titleText,
                        colors=self.COLORS, numpts=self.numpts)

            outkeys = list(self.outdata.keys())

            if useCache == 0:
                archiveMod = "Official" if mod == "Fcst" else mod
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
                        if outkey in self.saveGrids:
                            verGrid = self.saveGrids[outkey]
                        else:
                            recs = cases[casekey]
                            verGrid = self.VU.getVerGrids(archiveMod, base, element, self.perStart, self.perEnd, mode=modeVal, recList=recs)
                            if verGrid is None:
                                continue
                            if not isinstance(verGrid, np.ndarray) and len(verGrid) == 2:
                                (mag, direc) = verGrid
                                verGrid = mag
                            if reqElement in ("QPF", "SnowAmt"):
                                totTime = sum(self.VU.fncEtime[rec] - self.VU.fncStime[rec] for rec in recs)
                                fullFlag = 1 if totTime >= FullTimeRange.duration() else 0
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

                        self.graphmin = 0.0 if accumVal == 1 else self.allmin
                        if mod not in modsWithData:
                            modsWithData.append(mod)
                        if len(modsWithData) < 10:
                            self.drawGraph(self.outdata, self.outfull,
                                max(self.graphStart, self.mintime), self.maxtime,
                                self.graphmin, self.allmax,
                                label=self.labelText, title=self.titleText,
                                colors=self.COLORS, numpts=self.numpts)

                if len(modsWithData) >= 10:
                    self.drawGraph(self.outdata, self.outfull,
                        max(self.graphStart, self.mintime), self.maxtime,
                        self.graphmin, self.allmax,
                        label=self.labelText, title=self.titleText,
                        colors=self.COLORS, numpts=self.numpts)
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

            self.graphmin = 0.0 if accumVal == 1 else self.allmin
            self.drawGraph(self.outdata, self.outfull,
                max(self.graphStart, self.mintime), self.maxtime,
                self.graphmin, self.allmax,
                label=self.labelText, title=self.titleText,
                colors=self.COLORS, numpts=self.numpts)

            if tempoff == 1:
                self.state[mod] = 0
                self.drawGraph(self.outdata, self.outfull,
                    max(self.graphStart, self.mintime), self.maxtime,
                    self.graphmin, self.allmax,
                    label=self.labelText, title=self.titleText,
                    colors=self.COLORS, numpts=self.numpts)

            self.dlg.update_idletasks()

        # Fetch observation data if enabled
        if self.obsSource != "None":
            self.fetchObsData(reqElement, eaFlat, FullTimeRange, accumVal)

        self.graphmin = 0.0 if accumVal == 1 else self.allmin
        self.drawGraph(self.outdata, self.outfull,
            max(self.graphStart, self.mintime), self.maxtime,
            self.graphmin, self.allmax,
            label=self.labelText, title=self.titleText,
            colors=self.COLORS, numpts=self.numpts)

        self.dlg.updateButton.configure(state=tk.NORMAL)
        self.dlg.updateButton.configure(text="Update")
        self.dlg.updateButton.configure(bg=normalBG)
        self.dlg.updateButton.configure(fg=normalFG)
        return

    def fetchObsData(self, element, eaFlat, FullTimeRange, accumVal):
        """Fetch observation data from Obs, URMA, or RTMA database."""
        self.obsData = {}
        self.obsFull = {}
        
        if self.obsSource == "None":
            return
        
        modeVal = {"QPF": "Sum", "SnowAmt": "Sum", "MaxT": "Max", "MinT": "Min"}.get(element, "TimeWtAverage")
        
        try:
            foundSource = None
            altNames = {"Obs": ["Obs", "Metar", "METAR", "Observed"],
                       "URMA": ["URMA", "URMA25", "URMAe"],
                       "RTMA": ["RTMA", "RTMA25", "RTMAe"]}.get(self.obsSource, [self.obsSource])
            
            for altSource in altNames:
                if self.VU.checkFile(element, altSource):
                    foundSource = altSource
                    break
            
            if foundSource is None:
                print("Could not find observation data for %s in %s" % (element, self.obsSource))
                return
            
            obsCases = self.getModelCases(element, foundSource, self.perStart, self.perEnd, accum=accumVal)
            if not obsCases:
                return
            
            for casekey in sorted(obsCases.keys()):
                (basestr, startstr, endstr) = casekey.split(",")
                base = int(basestr)
                recs = obsCases[casekey]
                
                obsGrid = self.VU.getVerGrids(foundSource, base, element, self.perStart, self.perEnd, mode=modeVal, recList=recs)
                if obsGrid is None:
                    continue
                if not isinstance(obsGrid, np.ndarray) and len(obsGrid) == 2:
                    obsGrid = obsGrid[0]
                
                if element in ("QPF", "SnowAmt"):
                    totTime = sum(self.VU.fncEtime[rec] - self.VU.fncStime[rec] for rec in recs)
                    fullFlag = 1 if totTime >= FullTimeRange.duration() else 0
                else:
                    fullFlag = 1
                
                obsKey = "OBS,%d" % base
                statsList = self.getStats(obsGrid, eaFlat)
                self.obsData[obsKey] = copy.copy(statsList)
                self.obsFull[obsKey] = fullFlag
                self.allmax = max(self.allmax, statsList[1])
                self.allmin = min(self.allmin, statsList[0])
                
        except Exception as e:
            print("Error fetching observation data: %s" % str(e))
            import traceback
            traceback.print_exc()

    def getStats(self, grid, areaFlat):
        flatGrid = np.ravel(grid)
        data = np.compress(areaFlat, flatGrid)
        sdata = np.sort(data)
        datapts = data.shape[0]
        minvalue = sdata[0]
        maxvalue = sdata[datapts - 1]
        avg = float(np.add.reduce(sdata)) / float(datapts)
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
                reclist = [int(rec) for rec in recnumberList if self.VU.fncBtime[int(rec)] == base]
                cases["%d,%d,%d" % (base, perStart, perEnd)] = reclist
        return cases

    def drawGraph(self, graphData, fullData, xmin, xmax, ymin, ymax,
                  label="value", title="Graph", scaleDist=1, windowHeight=100,
                  histcolor="blue", axescolor="black", xticknum=25, yticknum=25,
                  showWidth=1, showTotal=1, showMax=1, showValueLine=0, valueLine=" ",
                  showPercentileLine=0, percentileLine=" ", colors={}, numpts=100):
        
        typeStr = ("Min", "Max", "Average", "1st percentile", "5th percentile",
                   "10th percentile", "25th percentile", "50th percentile",
                   "75th percentile", "90th percentile", "95th percentile", "99th percentile")

        graphType = self.dlg.graphType.get()
        graphTypeStr = typeStr[graphType]

        allkeys = list(graphData.keys())
        ymin, ymax = 999999.0, 0.0
        for key in allkeys:
            value = graphData[key][graphType]
            ymin, ymax = min(ymin, value), max(ymax, value)
        
        if self.obsData:
            for key in self.obsData.keys():
                value = self.obsData[key][graphType]
                ymin, ymax = min(ymin, value), max(ymax, value)
        
        if ymax == ymin:
            ymax = ymin + 1.0
        yrange = max(ymax - ymin, 1.0)

        ytick = max(self.niceNumDec(yrange / yticknum, 1), self.minTick)
        gymax = min(((int(float(ymax) / float(ytick)) + 1) * ytick), self.absMax)
        gymin = max((int(float(ymin) / float(ytick)) * ytick), self.absMin)

        self.dlg.histCanv.delete(tk.ALL)
        self.dlg.histCanv.create_text(0, 0, text="start", font="Helvetica 14", anchor=tk.SE, justify=tk.RIGHT, fill="black", tags="cursorReadout", state=tk.HIDDEN)

        self.left = 50.0
        self.right = self.dlg.curwidth - 165.0
        self.bot = 60.0
        self.top = self.dlg.curheight - 54.0
        self.fbot = 54.0
        self.ftop = self.dlg.curheight - 60.0
        self.setgraph(xmin, xmax, gymin, gymax, self.left, self.right, self.bot, self.top, self.dlg.curheight)

        graphArea = self.dlg.histCanv.create_polygon(self.left, self.fbot, self.left, self.ftop, self.right, self.ftop, self.right, self.fbot, self.left, self.fbot, fill="", outline="", tag="graphArea")
        self.dlg.histCanv.tag_bind(graphArea, '<Enter>', self.targetCursor)
        self.dlg.histCanv.tag_bind(graphArea, '<Leave>', self.normCursor)
        self.dlg.histCanv.tag_bind(graphArea, "<Button-1>", self.showReadOut)
        self.dlg.histCanv.tag_bind(graphArea, "<Button1-Motion>", self.moveReadOut)
        self.dlg.histCanv.tag_bind(graphArea, "<ButtonRelease-1>", self.hideReadOut)

        self.valhaxes(self.dlg.histCanv, xmin, xmax, gymin, gymax, ytick, label, axescolor)

        xmid = (self.left + self.right) / 2.0
        self.dlg.histCanv.create_text(xmid, 2, anchor=tk.N, text=title, justify=tk.CENTER, font="helvetica 24 bold", fill=axescolor)

        if numpts > 1:
            self.dlg.histCanv.create_text(xmid, 26, anchor=tk.N, text="(%s of %d points)" % (graphTypeStr, numpts), justify=tk.CENTER, font="helvetica 10", fill=axescolor)
            yoffset = 12
        else:
            yoffset = 0

        self.dlg.histCanv.create_text(xmid, 26 + yoffset, anchor=tk.N, text=self.getPeriodString(self.perStart, self.perEnd), justify=tk.CENTER, font="helvetica 12", fill=axescolor)

        modOrder = self.getModelOrder(list(graphData.keys()))
        allkeys = sorted(list(graphData.keys()))
        allkeys.reverse()
        displayNum = -1
        lastBoxY = 0
        
        for mod in modOrder:
            displayNum += 1
            lastx, lasty = -1, -1
            histcolor = colors.get(mod, "blue")

            modtextId = self.dlg.histCanv.create_text(self.dlg.curwidth - 115, 40 + (displayNum * 13), anchor=tk.SW, text=mod, justify=tk.LEFT, font="helvetica 12", fill=histcolor if self.state[mod] == 1 else "grey", tags="%sButton" % mod)

            bbox = self.dlg.histCanv.bbox(modtextId)
            (x1, y1, x2, y2) = bbox
            yt = y1 + 3 if lastBoxY == 0 else lastBoxY
            buttonArea = self.dlg.histCanv.create_polygon(x1, yt, x1, y2 - 1, x2, y2 - 1, x2, yt, x1, yt, fill="", outline="", tag="buttonArea")
            lastBoxY = y2 + 1

            def handler(event, self=self, modelStr=mod):
                return self.changeModelStatus(event, modelStr)
            self.dlg.histCanv.tag_bind(buttonArea, '<Button-1>', handler)
            def handler(event, self=self, modelStr=mod):
                return self.changeGroupStatus(event, modelStr)
            self.dlg.histCanv.tag_bind(buttonArea, '<Button-2>', handler)
            def handler(event, self=self, modelStr=mod):
                return self.postModelPopup(event, modelStr)
            self.dlg.histCanv.tag_bind(buttonArea, '<Button-3>', handler)
            self.dlg.histCanv.tag_bind(buttonArea, '<Enter>', self.pickCursor)
            self.dlg.histCanv.tag_bind(buttonArea, '<Leave>', self.normCursor)

            prevkey = "--"
            for key in allkeys:
                (modstr, basestr) = key.split(",")
                if mod != modstr:
                    continue
                base = int(basestr)
                value = graphData[key][graphType]
                (sx, sy) = self.graphcoord(base, value)
                
                if (sx >= self.left) and (sx <= self.right) and (sy >= self.fbot) and (sy <= self.ftop):
                    self.dlg.histCanv.create_text(sx, sy, anchor=tk.S, text=self.valueFormat % value, justify=tk.CENTER, font="helvetica %s" % LabelSizeStrings[self.labelSize[mod]], fill=histcolor, tags="%sLabel" % mod)
                    self.dlg.histCanv.create_line(sx - 2, sy - 2, sx - 2, sy + 2, sx + 2, sy + 2, sx + 2, sy - 2, sx - 2, sy - 2, fill=histcolor, tags="%sPoint" % mod)
                
                if lastx >= 0:
                    lw = self.lineWidth[mod]
                    (x1, y1, x2, y2) = self.clipLine(self.left, self.fbot, self.right, self.ftop, lastx, lasty, sx, sy)
                    if x1 is not None and prevkey != "--":
                        dash = None if fullData.get(prevkey, 1) == 1 else (lw, lw)
                        self.dlg.histCanv.create_line(x1, y1, x2, y2, fill=histcolor, width=lw, dash=dash, tags=mod)

                lastx, lasty, prevkey = sx, sy, key

            if (self.pointState[mod] == 0) or (self.state[mod] == 0):
                self.dlg.histCanv.itemconfigure("%sPoint" % mod, state=tk.HIDDEN)
            if (self.labelState[mod] == 0) or (self.pointState[mod] == 0) or (self.state[mod] == 0):
                self.dlg.histCanv.itemconfigure("%sLabel" % mod, state=tk.HIDDEN)
            if self.state[mod] == 0:
                self.dlg.histCanv.itemconfigure(mod, state=tk.HIDDEN)

        # Draw observation line
        if self.obsSource != "None" and self.obsData:
            self.drawObsLine(graphType, displayNum + 1)

        self.dlg.histCanv.lift("buttonArea")
        modOrder.reverse()
        for mod in modOrder:
            self.dlg.histCanv.lift("%sLabel" % mod)
        self.dlg.histCanv.lift("graphArea")
        self.dlg.update_idletasks()

    def drawObsLine(self, graphType, labelPosition):
        """Draw the observation data as a dashed line."""
        obsConfig = ObsSourceConfig.get(self.obsSource, {})
        obsColor = obsConfig.get("color", "#00DD00")
        obsDash = obsConfig.get("dash", (6, 4))
        obsLabel = obsConfig.get("label", "Observed")
        
        self.dlg.histCanv.create_text(self.dlg.curwidth - 115, 40 + (labelPosition * 13), anchor=tk.SW, text=obsLabel, justify=tk.LEFT, font="helvetica 12 italic", fill=obsColor, tags="ObsButton")
        
        obsKeys = sorted(list(self.obsData.keys()))
        lastx, lasty, prevkey = -1, -1, "--"
        
        for key in obsKeys:
            parts = key.split(",")
            if len(parts) < 2:
                continue
            base = int(parts[1])
            value = self.obsData[key][graphType]
            (sx, sy) = self.graphcoord(base, value)
            
            if (sx >= self.left) and (sx <= self.right) and (sy >= self.fbot) and (sy <= self.ftop):
                self.dlg.histCanv.create_oval(sx - 3, sy - 3, sx + 3, sy + 3, fill=obsColor, outline=obsColor, tags="ObsPoint")
            
            if lastx >= 0:
                (x1, y1, x2, y2) = self.clipLine(self.left, self.fbot, self.right, self.ftop, lastx, lasty, sx, sy)
                if x1 is not None:
                    self.dlg.histCanv.create_line(x1, y1, x2, y2, fill=obsColor, width=2, dash=obsDash, tags="ObsLine")
            
            lastx, lasty, prevkey = sx, sy, key

    def getPeriodString(self, perStart, perEnd):
        DAYS = ("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
        MONTHS = ("Jan", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
        (syea, smon, sday, shou, smin, ssec, swda, syda, sdst) = self._gmtime(perStart).timetuple()
        (eyea, emon, eday, ehou, emin, esec, ewda, eyda, edst) = self._gmtime(perEnd).timetuple()
        if self.element in ("MaxT", "MinT", "MaxRH", "MinRH", "TdMrn", "TdAft"):
            if (eyda == syda) or ((eyda == (syda + 1)) and (self.element in ("MaxT", "MinRH", "TdAft"))):
                return "%s %s %d, %4.4d" % (DAYS[swda], MONTHS[smon], sday, syea)
            return "%s %s %d -- %s %s %d, %4.4d" % (DAYS[swda], MONTHS[smon], sday, DAYS[ewda], MONTHS[emon], eday, syea)
        elif (ewda == swda) or ((eyda == (syda + 1)) and (ehou == 0)):
            return "%2.2dZ - %2.2dZ %s %s %d, %4.4d" % (shou, ehou, DAYS[swda], MONTHS[smon], sday, syea)
        return "%4.4d/%2.2d/%2.2d %2.2dZ --- %4.4d/%2.2d/%2.2d %2.2dZ" % (syea, smon, sday, shou, eyea, emon, eday, ehou)

    def pickCursor(self, event):
        self.dlg.histCanv.config(cursor="hand2")

    def normCursor(self, event):
        self.dlg.histCanv.config(cursor="")

    def targetCursor(self, event):
        self.dlg.histCanv.config(cursor="tcross")

    def showReadOut(self, event):
        (valx, valy) = self.valcoord(event.x, event.y)
        self.cursorReadout = 1
        self.dlg.histCanv.itemconfigure("cursorReadout", text=self.valueFormat % valy)
        self.dlg.histCanv.lift("cursorReadout")
        self.dlg.histCanv.lift("buttonArea")
        self.dlg.histCanv.coords("cursorReadout", event.x, event.y)
        self.dlg.histCanv.itemconfigure("cursorReadout", state=tk.NORMAL)

    def moveReadOut(self, event):
        if (event.x < self.left) or (event.x > self.right) or (event.y < self.fbot) or (event.y > self.ftop):
            self.cursorReadout = 0
            self.dlg.histCanv.itemconfigure("cursorReadout", state=tk.HIDDEN)
        else:
            (valx, valy) = self.valcoord(event.x, event.y)
            self.dlg.histCanv.itemconfigure("cursorReadout", text=self.valueFormat % valy)
            self.dlg.histCanv.coords("cursorReadout", event.x, event.y)
            if self.cursorReadout == 0:
                self.dlg.histCanv.lift("cursorReadout")
                self.dlg.histCanv.lift("buttonArea")
                self.dlg.histCanv.itemconfigure("cursorReadout", state=tk.NORMAL)
                self.cursorReadout = 1

    def hideReadOut(self, event):
        self.dlg.histCanv.itemconfigure("cursorReadout", state=tk.HIDDEN)
        self.cursorReadout = 0

    def clipLine(self, left, top, right, bottom, x1, y1, x2, y2):
        dx, dy = x2 - x1, y2 - y1
        dt0, dt1 = 0, 1
        for p, q in ((-dx, x1 - left), (dx, right - x1), (-dy, y1 - top), (dy, bottom - y1)):
            if p == 0:
                if q < 0:
                    return None, None, None, None
                continue
            dt = q / p
            if p < 0:
                if dt > dt1:
                    return None, None, None, None
                dt0 = max(dt0, dt)
            else:
                if dt < dt0:
                    return None, None, None, None
                dt1 = min(dt1, dt)
        if dt0 > 0:
            x1, y1 = x1 + dt0 * dx, y1 + dt0 * dy
        if dt1 < 1:
            x2, y2 = x1 + dt1 * dx, y1 + dt1 * dy
        return x1, y1, x2, y2

    def changeModelStatus(self, event, modelStr):
        if self.state[modelStr] == 1:
            self.dlg.histCanv.itemconfigure(modelStr, state=tk.HIDDEN)
            self.dlg.histCanv.itemconfigure("%sPoint" % modelStr, state=tk.HIDDEN)
            self.dlg.histCanv.itemconfigure("%sLabel" % modelStr, state=tk.HIDDEN)
            self.dlg.histCanv.itemconfigure("%sButton" % modelStr, fill='grey')
            self.state[modelStr] = 0
        else:
            self.dlg.histCanv.itemconfigure(modelStr, state=tk.NORMAL)
            if self.pointState[modelStr] == 1:
                self.dlg.histCanv.itemconfigure("%sPoint" % modelStr, state=tk.NORMAL)
            if self.labelState[modelStr] == 1:
                self.dlg.histCanv.itemconfigure("%sLabel" % modelStr, state=tk.NORMAL)
            self.dlg.histCanv.itemconfigure("%sButton" % modelStr, fill=self.COLORS[modelStr])
            self.state[modelStr] = 1

    def changeGroupStatus(self, event, modelStr):
        group = self.primaryGroup[modelStr]
        newState = 0 if self.state[modelStr] == 1 else 1
        for mod in self.groups[group]:
            if newState == 0:
                self.dlg.histCanv.itemconfigure(mod, state=tk.HIDDEN)
                self.dlg.histCanv.itemconfigure("%sPoint" % mod, state=tk.HIDDEN)
                self.dlg.histCanv.itemconfigure("%sLabel" % mod, state=tk.HIDDEN)
                self.dlg.histCanv.itemconfigure("%sButton" % mod, fill='grey')
            else:
                self.dlg.histCanv.itemconfigure(mod, state=tk.NORMAL)
                if self.pointState[mod] == 1:
                    self.dlg.histCanv.itemconfigure("%sPoint" % mod, state=tk.NORMAL)
                if self.labelState[mod] == 1:
                    self.dlg.histCanv.itemconfigure("%sLabel" % mod, state=tk.NORMAL)
                self.dlg.histCanv.itemconfigure("%sButton" % mod, fill=self.COLORS[mod])
            self.state[mod] = newState

    def postModelPopup(self, event, modelStr):
        self.dlg.ModelPopup.entryconfigure(0, label=modelStr, foreground=self.COLORS[modelStr], activeforeground=self.COLORS[modelStr])
        self.dlg.lineWidth.set(self.lineWidth[modelStr])
        self.dlg.dataPoints.set(self.pointState[modelStr])
        self.dlg.dataLabels.set(self.labelState[modelStr])
        self.dlg.dataLabelSize.set(self.labelSize[modelStr])
        self.dlg.ModelPopup.post(event.x_root, event.y_root)
        self.dlg.ModelPopup.grab_set()

    def postAxisPopup(self, event):
        self.timeSplitX = event.x
        self.dlg.XAxisPopup.post(event.x_root, event.y_root)
        self.dlg.XAxisPopup.grab_set()

    def popupHandler(self, changeType):
        mod = self.dlg.ModelPopup.entrycget(0, "label")
        if changeType == "LineChange":
            self.lineWidth[mod] = self.dlg.lineWidth.get()
        elif changeType == "ColorChange":
            self.COLORS[mod] = self.dlg.lineColor.get()
        elif changeType == "ColorEdit":
            (o1, o2) = tkColor.askcolor(self.COLORS[mod], title="%s Color" % mod)
            if o1 and o2:
                self.COLORS[mod] = "#%02x%02x%02x" % (int(o1[0]), int(o1[1]), int(o1[2]))
                self.dlg.lineColor.set(self.COLORS[mod])
        elif changeType == "DataChange":
            self.pointState[mod] = self.dlg.dataPoints.get()
        elif changeType == "LabelChange":
            self.labelState[mod] = self.dlg.dataLabels.get()
        elif changeType == "LabelSizeChange":
            self.labelSize[mod] = self.dlg.dataLabelSize.get()
        elif changeType == "TimeSplit":
            (xtime, ytime) = self.valcoord(self.timeSplitX, 0)
            self.graphStart = (int(xtime / (6 * 3600))) * (6 * 3600)
        elif changeType == "TimeReset":
            self.graphStart = self.mintime
        elif changeType == "ObsSourceChange":
            self.obsSource = self.dlg.obsSource.get()
            self.getStuff("Update")
            return

        self.drawGraph(self.outdata, self.outfull, max(self.graphStart, self.mintime), self.maxtime, self.graphmin, self.allmax, label=self.labelText, title=self.titleText, colors=self.COLORS, numpts=self.numpts)

    def getModelOrder(self, allkeys):
        uniqueModels = []
        for key in allkeys:
            modstr = key.split(",")[0]
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
        for name, obj in [("state", self.state), ("lineWidth", self.lineWidth), ("pointState", self.pointState), ("labelState", self.labelState), ("labelSize", self.labelSize), ("colors", self.COLORS), ("obsSource", self.obsSource)]:
            self.saveObject(name, obj, "TrendsTool")

    def readStatus(self):
        for name, target in [("state", self.state), ("lineWidth", self.lineWidth), ("pointState", self.pointState), ("labelState", self.labelState), ("labelSize", self.labelSize), ("colors", self.COLORS)]:
            lf = self.getUserFile(name, "TrendsTool")
            if lf.exists():
                temp = self.getObject(name, "TrendsTool")
                for key in temp:
                    if key in target:
                        target[key] = temp[key]
        lf = self.getUserFile("obsSource", "TrendsTool")
        if lf.exists():
            self.obsSource = self.getObject("obsSource", "TrendsTool")

    def getUserFile(self, name, category):
        from com.raytheon.uf.common.localization import PathManagerFactory, LocalizationContext
        pathMgr = PathManagerFactory.getPathManager()
        lc = pathMgr.getContext(LocalizationContext.LocalizationType.valueOf('CAVE_STATIC'), LocalizationContext.LocalizationLevel.valueOf('USER'))
        return pathMgr.getLocalizationFile(lc, 'gfe/userPython/' + category + '/' + name)

    def niceNumDec(self, val, roundit):
        if val == 0:
            return 1
        e = math.floor(math.log10(val))
        a = 10.0 ** e
        f = val / a
        if roundit > 0:
            nf = 1 if f < 1.5 else 2 if f < 3.0 else 5 if f < 7.0 else 10
        else:
            nf = 1 if f <= 1.0 else 2.0 if f <= 2.0 else 5.0 if f <= 5.0 else 10.0
        return nf * a

    def valhaxes(self, canvas, xmin, xmax, ymin, ymax, ytick, label, axescolor="black"):
        (nx, ny), (xx, xy) = self.graphcoord(xmin, ymin), self.graphcoord(xmax, ymax)
        canvas.create_line(nx, ny, xx, ny, xx, xy, nx, xy, nx, ny, fill=axescolor)
        self.vtick(canvas, xmin, 3, 3, ymin, ymax, ytick, label=1, labelinterval=3, labeloffset=-7.0, labelanchor=tk.E, color=axescolor)
        self.vtick(canvas, xmax, 3, 3, ymin, ymax, ytick, label=1, labelinterval=3, labeloffset=7.0, labelanchor=tk.W, color=axescolor)
        xtick = 6 * 3600
        self.hgridDate(canvas, ymin, ymax, xmin, xmax, xtick)
        self.htick(canvas, ymin, 3, 3, xmin, xmax, xtick, label=1, labeltype="date", labelinterval=12, labeloffset=-5.0, labelanchor=tk.N, color=axescolor)
        self.htick(canvas, ymax, 3, 3, xmin, xmax, xtick, label=0, color=axescolor)
        (midx, _) = self.graphcoord((xmax + xmin) / 2.0, 0)
        canvas.create_text(midx, ny + 35, anchor=tk.N, text="Model Run", fill=axescolor)
        canvas.create_polygon(nx, ny, xx, ny, xx, ny + 35, nx, ny + 35, fill="", outline="", tags="axisarea")
        self.dlg.histCanv.tag_bind("axisarea", '<Button-3>', lambda e: self.postAxisPopup(e))
        (_, midy) = self.graphcoord(xmin, (ymin + ymax) / 2.0)
        labtext = "\n".join(label)
        canvas.create_text(nx - 35, midy, anchor=tk.E, text=labtext, fill=axescolor)
        canvas.create_text(xx + 35, midy, anchor=tk.W, text=labtext, fill=axescolor)

    def htick(self, canvas, yval, widabv, widblo, minx, maxx, xint, label=0, labeloffset=5, labelwidabv=2, labelwidblo=2, labelinterval=1, labeltype="value", skipfirst=0, skiplast=0, labelanchor=tk.N, negative=0, color="black"):
        numticks = int((maxx - minx) / xint) + 1
        for i in range(numticks):
            x = minx + (i * xint)
            (tx, ty) = self.graphcoord(x, yval)
            (_, _, _, xhou, _, _, _, _, _) = self._gmtime(x).timetuple()
            extra = labelwidabv if xhou == 0 else 0
            canvas.create_line(tx, ty - widabv - extra, tx, ty + widblo + extra, fill=color)
            if label and xhou == labelinterval and not ((skipfirst and i == 0) or (skiplast and i == numticks - 1)):
                (_, gmon, gday, _, _, _, _, _, _) = self._gmtime(x).timetuple()
                canvas.create_text(tx, ty - labeloffset, anchor=labelanchor, text="%2.2d/%2.2d" % (gmon, gday) if labeltype == "date" else str(int(x)), fill=color)

    def hgridDate(self, canvas, ymin, ymax, minx, maxx, xint, gridinterval=0, color="#A0A0A0"):
        numticks = int((maxx - minx) / xint) + 1
        for i in range(1, numticks - 1):
            x = minx + (i * xint)
            (_, _, _, xhou, _, _, _, _, _) = self._gmtime(x).timetuple()
            if xhou == gridinterval:
                (tx, ny), (_, xy) = self.graphcoord(x, ymin), self.graphcoord(x, ymax)
                canvas.create_line(tx, ny, tx, xy, fill=color)

    def vtick(self, canvas, xval, widlft, widrgt, miny, maxy, yint, label=0, labeloffset=5, labelwidlft=2, labelwidrgt=2, labelinterval=1, skipfirst=0, skiplast=0, labelanchor=tk.W, negative=0, color="black"):
        numticks = int((maxy - miny) / yint) + 1
        labeldigits = max(-math.floor(math.log10(yint)), 0)
        for i in range(numticks):
            y = miny + (i * yint)
            (tx, ty) = self.graphcoord(xval, y)
            canvas.create_line(tx - widlft, ty, tx + widrgt, ty, fill=color)
            if label and (i % labelinterval == 0) and not ((skipfirst and i == 0) or (skiplast and i == numticks - 1)):
                labelstring = ("%d" % y) if labeldigits == 0 else ("%%.%df" % labeldigits) % y
                canvas.create_text(tx + labeloffset, ty, anchor=labelanchor, text=labelstring, fill=color)
                canvas.create_line(tx - widlft - labelwidlft, ty, tx + widrgt + labelwidrgt, ty, fill=color)

    def setgraph(self, xmin, xmax, ymin, ymax, sxmin, sxmax, symin, symax, windowHeight):
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.xmult = (sxmax - sxmin) / (xmax - xmin)
        self.xoff = sxmin
        self.ymult = (symax - symin) / (ymax - ymin)
        self.yoff = symin
        self.yflip = windowHeight

    def graphcoord(self, x, y):
        return ((x - self.xmin) * self.xmult) + self.xoff, self.yflip - (((y - self.ymin) * self.ymult) + self.yoff)

    def valcoord(self, sx, sy):
        return ((sx - self.xoff) / self.xmult) + self.xmin, ((-sy + self.yflip - self.yoff) / self.ymult) + self.ymin

    def resizecanvas(self, event):
        self.dlg.curwidth, self.dlg.curheight = float(event.width), float(event.height)
        self.drawGraph(self.outdata, self.outfull, max(self.graphStart, self.mintime), self.maxtime, self.graphmin, self.allmax, label=self.labelText, title=self.titleText, colors=self.COLORS, numpts=self.numpts)


class ToolDialog(AppDialog.AppDialog):
    def __init__(self, title="Tk", callbackMethod=None, popupCallback=None, **kwargs):
        self.__callbackMethod = callbackMethod
        self.__popupCallback = popupCallback
        AppDialog.AppDialog.__init__(self, **kwargs)
        self.title(title)
        self.wm_attributes("-topmost", 1)
        self.protocol("WM_DELETE_WINDOW", self.closeCB)
        self.update_idletasks()
        self.minsize(self.winfo_width(), self.winfo_height() + 30)

    def buttonbox(self):
        buttonFrame = tk.Frame(self)
        buttonCenter = tk.Frame(buttonFrame)
        self.updateButton = tk.Button(buttonCenter, text="Update", command=self.updateCB, width=10)
        self.updateButton.pack(side=tk.LEFT, pady=5, padx=10)
        tk.Button(buttonCenter, text="Close", width=10, command=self.closeCB).pack(side=tk.LEFT, pady=5, padx=10)
        buttonCenter.pack(side=tk.TOP, fill=tk.X, expand=tk.TRUE)
        buttonFrame.pack(side=tk.BOTTOM)

    def body(self, master):
        self.createMenu()
        bodyFrame = tk.Frame(master)
        self.buildGraph(bodyFrame)
        bodyFrame.pack(side=tk.TOP, fill=tk.BOTH, expand=tk.TRUE)
        return bodyFrame

    def createMenu(self):
        self.menuBar = tk.Menu(self)
        self.winfo_toplevel()["menu"] = self.menuBar
        
        self.graphMenu = tk.Menu(self.menuBar)
        self.menuBar.add_cascade(label="Graph", menu=self.graphMenu)
        self.graphType = tk.IntVar(value=2)
        for lbl, val in [("Max", 1), ("99th Percentile", 11), ("95th Percentile", 10), ("90th Percentile", 9), ("75th Percentile", 8), ("50th Percentile", 7), ("Avg", 2), ("25th Percentile", 6), ("10th Percentile", 5), ("5th Percentile", 4), ("1st Percentile", 3), ("Min", 0)]:
            self.graphMenu.add_radiobutton(label=lbl, value=val, variable=self.graphType, command=self.setGraphType)

        self.obsMenu = tk.Menu(self.menuBar)
        self.menuBar.add_cascade(label="Observations", menu=self.obsMenu)
        self.obsSource = tk.StringVar(value="None")
        for lbl, val in [("None", "None"), ("Obs (METAR)", "Obs"), ("URMA Analysis", "URMA"), ("RTMA Analysis", "RTMA")]:
            self.obsMenu.add_radiobutton(label=lbl, value=val, variable=self.obsSource, command=self.setObsSource)

        self.ModelPopup = tk.Menu(self, tearoff=0)
        self.LineWidthMenu = tk.Menu(self.ModelPopup, tearoff=0)
        self.ColorMenu = tk.Menu(self.ModelPopup, tearoff=0)
        self.LabelSizeMenu = tk.Menu(self.ModelPopup, tearoff=0)

        self.dataPoints, self.dataLabels, self.dataLabelSize = tk.IntVar(value=1), tk.IntVar(value=0), tk.IntVar(value=1)
        self.ModelPopup.add_command(label="Model Name")
        self.ModelPopup.add_separator()
        self.ModelPopup.add_cascade(label="Line Width", menu=self.LineWidthMenu)
        self.ModelPopup.add_cascade(label="Color", menu=self.ColorMenu)
        self.ModelPopup.add_checkbutton(label="Show Data Points", variable=self.dataPoints, command=lambda: self.__popupCallback("DataChange"))
        self.ModelPopup.add_checkbutton(label="Show Data Labels", variable=self.dataLabels, command=lambda: self.__popupCallback("LabelChange"))
        self.ModelPopup.add_cascade(label="Data Label Size", menu=self.LabelSizeMenu)
        
        self.lineWidth = tk.IntVar(value=1)
        for w in [1, 2, 3]:
            self.LineWidthMenu.add_radiobutton(label="%d px" % w, value=w, variable=self.lineWidth, command=lambda: self.__popupCallback("LineChange"))
        
        self.lineColor = tk.StringVar(value="black")
        boxchar = "\u2588"
        for c in ["#000000", "#ffffff", "#ff0000", "#ffa500", "#ffff00", "#00ff00", "#228b22", "#00ffff", "#0000ff", "#a020f0", "#ffc0cb", "#ff00ff"]:
            self.ColorMenu.add_radiobutton(value=c, label=boxchar, foreground=c, activeforeground=c, hidemargin=tk.TRUE, variable=self.lineColor, command=lambda: self.__popupCallback("ColorChange"))
        self.ColorMenu.add_command(label="Custom...", command=lambda: self.__popupCallback("ColorEdit"))

        for i, size in enumerate(LabelSizeStrings):
            self.LabelSizeMenu.add_radiobutton(label="%s pt" % size, value=i, variable=self.dataLabelSize, command=lambda: self.__popupCallback("LabelSizeChange"))

        self.XAxisPopup = tk.Menu(self, tearoff=0)
        self.XAxisPopup.add_command(label="Start Time Axis Here", command=lambda: self.__popupCallback("TimeSplit"))
        self.XAxisPopup.add_command(label="Reset Time Axis to all times", command=lambda: self.__popupCallback("TimeReset"))

    def setGraphType(self):
        self.__callbackMethod("Redraw")

    def setObsSource(self):
        self.__popupCallback("ObsSourceChange")

    def updateCB(self):
        self.__callbackMethod("Update")

    def closeCB(self):
        self.__callbackMethod("Close")
        self.cancel()

    def buildGraph(self, master, windowWidth=1000, windowHeight=800):
        self.histCanv = tk.Canvas(master=master, width=windowWidth, height=windowHeight, bd=2, relief=tk.RIDGE, bg="white")
        self.histCanv.pack(side=tk.TOP, expand=tk.YES, anchor=tk.N, fill=tk.BOTH)
        self.curwidth, self.curheight = windowWidth, windowHeight
