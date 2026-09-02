# ----------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# Diurnal
#
# Original Author: Tom LeFebvre 2.7
# 
# ADAPTED FOR BOI
# Version: 2.8 - Aug 11, 2022 - Added AptT and Max/MinAptT - KA
#          2.9 - Oct 30, 2023 - Added WBGT - KA
#
# ----------------------------------------------------------------------------

# The MenuItems list defines the GFE menu item(s) under which the
# Procedure is to appear.
# Possible items are: Populate, Edit, Consistency, Verify, Hazards
MenuItems = ["Populate"]

from math import *
import AbsTime
import SmartScript
import TimeRange, DatabaseID
import ProcessVariableList
import time,calendar,datetime,dateutil,pytz,copy
import sys,LogStream
import DiurnalConfig
import numpy as np
import datetime as dt
from datetime import timedelta, date

## For documentation on the available commands,
##   see the SmartScript Utility, which can be viewed from
##   the Edit Actions Dialog Utilities window

class Procedure (SmartScript.SmartScript):
    def __init__(self, dbss):
        SmartScript.SmartScript.__init__(self, dbss)

    def makeTimeRange(self, start, end):
        """Creates a time range.
        Args:
            double start - start of time range in seconds since the epoch began
            double end - end of time range in seconds since the epoch began
        Returns:
            Time range appropriate for AWIPS version
        """

        startTime = AbsTime.AbsTime(start)
        endTime = AbsTime.AbsTime(end)

        return TimeRange.TimeRange(startTime, endTime)

    def makeMaxTimeRange(self):
        
        import TimeRange, AbsTime
        tr = TimeRange.allTimes()

        return tr

    def getParmMinMaxLimits(self, modelName, weName):

        parm = self.getParm(modelName, weName, "SFC")
        return parm.getGridInfo().getMinValue(), parm.getGridInfo().getMaxValue()
    
    def getParmTimeConstraints(self, parm):

        parmStart = parm.getGridInfo().getTimeConstraints().getStartTime()
        parmDuration = parm.getGridInfo().getTimeConstraints().getDuration()

        return parmStart, parmDuration
    
    def timeRangeSort(self, a, b):
        if a.startTime().unixTime() < b.startTime().unixTime():
            return -1
        else:
            return 1

    def getWEInventory(self, modelName, WEName, timeRange=None):

        if timeRange is None:
            timeRange = TimeRange.allTimes()

        trList = []
        try:
            gridInfo = self.getGridInfo(modelName, WEName, "SFC", timeRange)
        except:
            return trList

        for g in gridInfo:
            if timeRange.overlaps(g.gridTime()):
                trList.append(g.gridTime())

        return trList

    def _cubicSpline(self, grids, times, interpTimes):
        """
        Main cubic spline method that accepts a list of grids and int time
        along with a list of times for which grids are to be calculated.
        This method returns the corresponding list of grids that matches the
        interpTimes list.
        """

        # STEP 1: Create coefficients for cubic spline curve
        # zCoefs : List of cubic spline coefficient grids computed to fit the
        #          curve defined by grids and times
        # n      : length of grids - 1.
        # Determine coefficients
        if grids == [] or times == []:
            self.statusBarMsg("No grids or times to interpolate spline.", "S")
            return

        gridShape = grids[0].shape

        timeGrids = []
        for t in times:
            tGrid = np.zeros(gridShape) + t
            timeGrids.append(tGrid)

        n = len(grids) - 1
        zCoefs = self._spline3_coef(n, timeGrids, grids)

        # Create interpolated grids using coefficients
        # interpTimes : List of times for which we want interpolated grids
        # gridList : List of interpolated Grids

        # Create interpolated grids
        gridList = []
        for interpTime in interpTimes:
            x = np.zeros(gridShape) + interpTime  # make a grid of times
            xGrid = self._spline3_eval(n, timeGrids, grids, zCoefs, x)
            gridList.append(xGrid)

        return gridList

    def _spline3_coef(self, n, t, y):
        """
        This method calculates the spline coefficients that are later used to
        calculate grids at the interpolation times.  This method is just a helper
        method to _cubicSpline and should not be called directly.
        """

        gridShape = y[0].shape
        # These will get filled in later with grids as values
        # They are just place holders
        h = [0] * n
        b = [0] * n
        u = [0] * n
        v = [0] * n
        z = [0] * (n + 1)
        # Calculate h and b
        #   range 0 thru n-1
        for i in range(n):
            h[i] = t[i + 1] - t[i]
            b[i] = (y[i + 1] - y[i]) / h[i]
        # Calculate u and v as functions of h and b
        #   range 1 thru n-1
        u[1] = 2 * (h[0] + h[1])
        v[1] = 6 * (b[1] - b[0])
        for i in range(2, n):
            u[i] = 2.0 * (h[i] + h[i - 1]) - h[i - 1] ** 2.0 / u[i - 1]
            v[i] = 6.0 * (b[i] - b[i - 1]) - h[i - 1] * v[i - 1] / u[i - 1]
        # Calculate z
        #   range 0 thru n
        z[n] = np.zeros(gridShape)
        for i in range(n - 1, 0, -1):
            z[i] = (v[i] - h[i] * z[i + 1]) / u[i]
        z[0] = np.zeros(gridShape)
        return z

    def _spline3_eval(self, n, t, y, z, x):
        """
        This method accepts the spline coefficients and calculates a grid.
        This method is a help method to _cubicSpline and should not be
        called directly
        """

        for i in range(n - 1, 0, -1):
            if x[0][0] - t[i][0][0] >= 0:
                break
        h = t[i + 1] - t[i]
        tmp = (z[i] / 2) + (x - t[i]) * (z[i + 1] - z[i]) / (6 * h)
        tmp = (
            -(h / 6) * (z[i + 1] + 2 * z[i]) + (y[i + 1] - y[i]) / h + (x - t[i]) * tmp
        )

        return y[i] + (x - t[i]) * tmp

    def removeExcludedTRs(self, trList):
        """
        Returns a new list ot timeRanges with the excluded periods
        removed. Uses self._excludeTRList defined else where.
        """

        newTRList = []

        for tr in trList:
            found = False
            for exTR in self._excludeTRList:
                if tr.overlaps(exTR):
                    found = True
            if not found:
                newTRList.append(tr)

        return newTRList

    def calcObsHourlyAvg(self, modelName, weName):
        """
        This method fetches all of the grid in the obs database,
        calculates the average at each hour and returns a list of
        grids representing the average of each hour.
        """

        trList = self.getWEInventory(modelName, weName)

        # Remove some timeRanges if they to be excluded
        if self._excludeTRList:
            trList = self.removeExcludedTRs(trList)

        sumList = [self.empty()] * 24  # one for each hour
        gridCount = [0] * 24  # number of grids for each hour
        for tr in trList:
            grid = self.getGrids(modelName, weName, "SFC", tr)
            hour = time.gmtime(tr.startTime().unixTime()).tm_hour

            sumList[hour] += grid
            gridCount[hour] += 1

        avgList = [self.empty()] * 24

        for i in range(len(avgList)):
            if gridCount[i] > 0:
                avgList[i] = sumList[i] / gridCount[i]

        return avgList

    def getAvailableDBIDs(self, modelName):
        """
        Returns a list of IFP database IDs that match the specified
        modelName.  The list is returned in latest-first time order.
        """

        allDBs = self.availableDatabases()
        dbList = []
        # Find all the dbIDs that match this model name and are
        # not D2D databases.
        for dbID in allDBs:
            dbIDModelName = dbID.modelName()
            if dbIDModelName == modelName and dbID.modelIdentifier().find("D2D") == -1:
                dbList.append(dbID)

        return dbList

    def getWEFrequency(self, modelName, weName):
        """
        Returns the frequency (in hours) of the specified modelName and weather element
        """

        try:
            parm = self.getParm(modelName, weName, "SFC")
            parmRepeat = parm.getGridInfo().getTimeConstraints().getRepeatInterval()
            return parmRepeat / 3600
        except:
            self.statusBarMsg(
                modelName
                + " "
                + weName
                + " not found when looking up forecast interval.",
                "U",
            )
            return 0

    def expandTimeRange(self, timeRange, timeStep):
        start = timeRange.startTime().unixTime()
        end = timeRange.endTime().unixTime()

        inc = int(timeStep / 2) * 3600
        start -= inc
        end += inc
        if timeStep == 6:
            end -= 3600
        tr = self.makeTimeRange(start, end)

        return tr
    
    def highlightTRs(self, weName, highlightDict, executeTR):
        trList = list(highlightDict.keys())
        trList.sort(key=lambda a: a.startTime().unixTime())
    
        for tr in trList:
            color, modelName = highlightDict[tr]
            if not executeTR.overlaps(tr):
                continue
            
            start = tr.startTime().unixTime()
            end = tr.endTime().unixTime()
            if start < executeTR.startTime().unixTime():
                start = executeTR.startTime().unixTime()
            if end > executeTR.endTime().unixTime():
                end = executeTR.endTime().unixTime()

            tr = self.makeTimeRange(start, end)
            self.highlightGrids("Fcst", weName, "SFC", tr, color)
        return

    def extractMessages(self, highlightDict):

        trList = list(highlightDict.keys())

        msgDict = {}

        if len(trList) == 0:
            return msgDict

        trList.sort(key=lambda a: a.startTime().unixTime())
        runningTR = trList[0]
        runningColor, runningModel = highlightDict[runningTR]
        for tr in trList:
            color, model = highlightDict[tr]
            adjacent = (tr.startTime() == runningTR.endTime() or
                        tr.overlaps(runningTR))
            sameMsg = (color == runningColor and model == runningModel)
            if (adjacent and sameMsg) or trList.index(tr) == 0:
                runningTR = runningTR.combineWith(tr)
            else:
                msgDict[runningTR] = runningModel
                runningTR = tr
                runningColor, runningModel = highlightDict[tr]
     
        msgDict[runningTR] = runningModel

        return msgDict

    def timeRangeString(self, timeRange):
        start = timeRange.startTime().unixTime()
        end = timeRange.endTime().unixTime()

        startDay = time.strftime("%a", time.gmtime(start))
        startHour = str(time.gmtime(start).tm_hour).zfill(2)
        endDay = time.strftime("%a", time.gmtime(end))
        endHour = str(time.gmtime(end).tm_hour).zfill(2)

        timeStr = startDay + " " + startHour + "Z-" + \
                  endDay + " " + endHour + "Z" 
    
        return timeStr
    
    
    def reportMissingGrids(self, weName, highlightDict, executeTR):

        msgDict = self.extractMessages(highlightDict)

        trList = list(msgDict.keys())
        trList.sort(key=lambda a: a.startTime().unixTime())
        for tr in trList:
            if not executeTR.overlaps(tr):
                continue
            timeStr = self.timeRangeString(tr)
            modelName = msgDict[tr]
            msg = "Diurnal - " + modelName + " " + weName + " grids used at: " + timeStr
            self.statusBarMsg(msg, "R")

        return

    def modelInfoStr(self, modelName, modelTime):
        
        day = str(time.gmtime(modelTime).tm_mday).zfill(2)
        hour = str(time.gmtime(modelTime).tm_hour).zfill(2)
        modelStr = day + "/" + hour + "Z " + modelName

        return modelStr


    def getModelInventory(self, modelList, weName, timeRange):
        """
        This method returns a dictionary containing the databaseID
        for each inventory timeRange that overlaps the specified
        timeRange.  If a particular grid is missing, a substitute
        will be found with the same valid time.
        """

        primaryModelName = modelList[0]

        # get the list of database IDs that matches the primary model
        dbIDList = self.getAvailableDBIDs(primaryModelName)
        if not dbIDList:
            print("ERROR! Primary model " + primaryModelName + " not found.")
            return []

        primaryDB = dbIDList[0]
        modelTime = int(primaryDB.modelTime().unixTime())

        # if the modelTime is zero, we're probably using a singleton database. Just use
        # the start time of the timeRange
        if modelTime == 0:
            modelTime = timeRange.startTime().unixTime()
        elif modelTime > timeRange.startTime().unixTime():
            modelTime -= 6 * 3600

        stopTime = int(timeRange.endTime().unixTime())

        trList = []

        # make a list of timeRanges when we expect grids.
        gridTimeStep = int(self.getWEFrequency(primaryModelName, weName))
        for t in range(modelTime, stopTime, gridTimeStep * 3600):
            tr = self.makeTimeRange(t, t + 3600)
            trList.append(tr)

        # make the list of dbIDs that we may need
        dbList = []
        dbInvDict = {}
        for modelName in modelList:
            dbIDList = self.getAvailableDBIDs(modelName)

            if dbIDList:
                dbID = dbIDList[0]
                dbList.append(dbID)
                dbInv = self.getWEInventory(dbID, weName)
                dbIDname = repr(dbID)
                dbInvDict[dbIDname] = dbInv  # save the inventory for this db

                # Secretly add the previous version to the list as a
                # fallback choice but only after the first model
                if len(dbIDList) > 1:  # has more than one version
                    dbID = dbIDList[1]
                    dbList.append(dbID)
                    dbInv = self.getWEInventory(dbID, weName)

                    dbIDname = repr(dbID)
                    dbInvDict[dbIDname] = dbInv

        invList = []
        self._highlightDict = {}
        # Fetch a dbID for each timeRange checking the inventory to ensure
        # that we'll get a grid when we ask for it
        for tr in trList:
            for dbID in dbList:

                modelName = dbID.modelName()

                dbName = repr(dbID)

                if modelName == "":
                    print("Error in determining modelName.")

                gridTimeStep = self.getWEFrequency(modelName, weName)
                expandedTR = self.expandTimeRange(tr, gridTimeStep)
                if tr in dbInvDict[dbName]:  # found the tr in this inventory
                    invList.append((tr, dbID))  # save the dbID for this timeRange
                    # if we used any substitute model source record the choice
                    # in the highlightDict for later highlighting
                    if dbList.index(dbID) > 0:  # not the primary choice
                        day = str(time.gmtime(dbID.modelTime().unixTime()).tm_mday)
                        hour = str(time.gmtime(dbID.modelTime().unixTime()).tm_hour)
                        modelStr = self.modelInfoStr(
                            modelName, dbID.modelTime().unixTime()
                        )
                        self._highlightDict[expandedTR] = (
                            self.configDict["highlightColors"][dbList.index(dbID) - 1],
                            modelStr,
                        )
                    break

        return invList
    
    def makeHourlyGrids(self, modelList, weName, timeRange):
        """
        This method returns hourly grids from the specified model.
        A cublic spline technique is used to interpolation temporally
        from the grids that are available.  The specified timeRange is
        ignored when a (model is specified (not obs).
        """

        hourlyGridDict = {}
        primaryModelName = modelList[0]

        # If the model is "Obs" fetch Obs data and average each hour.
        if (
            primaryModelName in self.configDict["obsModelList"]
        ):  # Obs grids are averaged per hour
            avgList = self.calcObsHourlyAvg(primaryModelName, weName)
            # Duplicate the hourly Obs over the entire master list
            # of timeRanges.
            start = timeRange.startTime().unixTime()
            end = timeRange.endTime().unixTime() + 3600  # inclusive
            for t in range(start, end, 3600):
                hour = time.gmtime(t).tm_hour
                tr = self.makeTimeRange(t, t + 3600)
                hourlyGridDict[tr] = avgList[hour]

            return hourlyGridDict

        # If the source is not Obs, fetch the inventory and make hourly grids
        gridInv = self.getModelInventory(modelList, weName, timeRange)

        startTimeList = []
        gridList = []
        trList = []
        for tr, dbID in gridInv:
            startTimeList.append(tr.startTime().unixTime())
            grid = self.getGrids(dbID, weName, "SFC", tr)
            gridList.append(grid)
            trList.append(tr)

        if not trList:
            self.statusBarMsg("No model grids found for this timeRange.", "S")
            return hourlyGridDict

        # make a list of times one hour apart
        interpTimes = list(
            range(
                trList[0].startTime().unixTime(),
                trList[-1].startTime().unixTime(),
                3600,
            )
        )

        # interpolate the grids to one hour and save in the global variable
        hourlyGridList = self._cubicSpline(gridList, startTimeList, interpTimes)

        # create and save the timeRanges for each of the interpolated grids
        # and save them in the dictionary
        trList = []
        for i in range(len(interpTimes)):
            start = interpTimes[i]
            tr = self.makeTimeRange(start, start + 3600)
            hourlyGridDict[tr] = hourlyGridList[i]

        # replace the master list with the timeRanges calculated above

        return hourlyGridDict
    

    # Fetches and returns the names of the time periods for the
    # GUI display.
    def getTimePeriods(self):
        #definedPeriods = self.getConfigItem("TimeScalePeriods",[])
        definedPeriods = self.getConfigItem("PublishTimes",[])
        trNames = []
        for trName in definedPeriods:
            trNames.append(trName)

        return trNames

    def makeExecuteTR(self, startName, endName, selectedTimeRange):
        """
        Using the specified timeRange info, determines the timeRange
        over which the tool should run.
        """

        # check for DST config flag and adjust if set
        dstOffset = 0
        if self.configDict.get("adjustDST", False):
            localTime = time.localtime()
            dst = localTime.tm_isdst
            dstOffset = dst * 3600
            print("adjusted time one hour")

        # Fetch the timeRanges and apply the offset only to named TRs
        if startName == "Selected TimeRange":
            start = selectedTimeRange.startTime().unixTime()
        else:
            start = self.getTimeRange(startName).startTime().unixTime() + dstOffset
        if endName == "Selected TimeRange":
            end = selectedTimeRange.endTime().unixTime()
        else:
            end = self.getTimeRange(endName).endTime().unixTime() + dstOffset

        executeTR = self.makeTimeRange(start, end)

        return executeTR

    def getMinMaxTRLists(self, weName, executeTR):
        """
        Get the inventory (timeRanges) for the min and max elements.
        The grids in the Fcst database are all that really matter
        so just fetch those and return them.
        """

        minWeName = "Min" + weName
        maxWeName = "Max" + weName

        # get the Fcst min and max inventory
        fcstMinTRList = self.getWEInventory("Fcst", minWeName)
        fcstMaxTRList = self.getWEInventory("Fcst", maxWeName)

        return fcstMinTRList, fcstMaxTRList

    def checkFcstMinMaxGrids(self, weName, executeTR, minDict, maxDict):
        """
        Check the Fcst set of min and max grids to see any are missing
        and a crude check on integrity for each grid that does exist.
        Returns True if all checks were successful, and False otherwise.
        """

        minTRList = sorted(list(minDict.keys()), key=lambda a: a.startTime().unixTime())
        maxTRList = sorted(list(maxDict.keys()), key=lambda a: a.startTime().unixTime())
        minWE = "Min" + weName
        maxWE = "Max" + weName

        minParm = self.getParm("Fcst", minWE, "SFC")
        maxParm = self.getParm("Fcst", maxWE, "SFC")

        oneDay = 24 * 3600
        execStart = int(executeTR.startTime().unixTime() / oneDay) * oneDay
        execEnd = int(executeTR.endTime().unixTime() / oneDay) * oneDay

        completeInv = True

        # Check min parm inventory
        parmStart, parmDuration = self.getParmTimeConstraints(minParm)

        for t in range(execStart, execEnd, oneDay):
            start = t + parmStart
            end = start + parmDuration
            tr = self.makeTimeRange(start, end)
            if tr not in minTRList:
                completeInv = False
                msg = "Missing Fcst {} grid at {}. Aborting.".format(minWE, str(tr))
                self.statusBarMsg(msg, "S")

        # Check max parm inventory
        parmStart, parmDuration = self.getParmTimeConstraints(maxParm)

        for t in range(execStart, execEnd, oneDay):
            start = t + parmStart
            end = start + parmDuration
            tr = self.makeTimeRange(start, end)
            if tr not in maxTRList:
                completeInv = False
                msg = "Missing Fcst {} grid at {}. Aborting.".format(maxWE, str(tr))
                self.statusBarMsg(msg, "S")

        # return now if something is wrong
        if not completeInv:
            return completeInv

        # get the min and max allowed values for MinT and MaxT
        minTMin, minTmax = self.getParmMinMaxLimits("Fcst", minWE)
        maxTMin, maxTmax = self.getParmMinMaxLimits("Fcst", maxWE)

        # Check grid values for scratch grids
        for tr in minTRList:
            # make sure we're in range
            if (
                tr.startTime().unixTime() > execEnd
                or tr.endTime().unixTime() < execStart
            ):
                continue

            grid = minDict[tr]
            # calculate the average value over the entire grid
            average = np.sum(np.sum(grid)) / np.product(list(grid.shape))
            if (
                average == minTMin
            ):  # if the average is the min value it's a scratch grid
                msg = "Possible scratch detected grid at {}. Aborting!".format(str(tr))
                self.statusBarMsg(msg, "S")
                completeInv = False

        for tr in maxTRList:
            # make sure we're in range
            if (
                tr.startTime().unixTime() > execEnd
                or tr.endTime().unixTime() < execStart
            ):
                continue

            grid = maxDict[tr]
            minVal = min(grid.flat)
            # calculate the average value over the entire grid
            average = np.sum(np.sum(grid)) / np.product(list(grid.shape))
            if (
                average == maxTMin
            ):  # if the average is the min value it's a scratch grid
                msg = "Possible scratch detected grid at {}. Aborting!".format(str(tr))
                self.statusBarMsg(msg, "S")
                completeInv = False

        return completeInv

    # generates a min grid from the hourly grids
    def calcMinMaxFromHourlies(self, timeRange, minOrMax):
        """
        generates a min or max grid from the hourly grids
        """

        trList = list(self._modelHrGridDict.keys())
        trList.sort(key=lambda a: a.startTime().unixTime())

        # See if we have any overlapping hourly grids. If not, return None
        if not trList:
            return None

        hourlyStart = trList[0].startTime().unixTime()
        hourlyEnd = trList[-1].endTime().unixTime()
        if (
            timeRange.endTime().unixTime() < hourlyStart
            or timeRange.startTime().unixTime() > hourlyEnd
        ):
            return None

        if minOrMax == "min":
            grid = self.newGrid(9999.0)
        else:
            grid = self.newGrid(-9999.0)

        for tr in trList:
            if tr.overlaps(timeRange):
                if minOrMax == "min":
                    mask = self._modelHrGridDict[tr] < grid
                else:
                    mask = self._modelHrGridDict[tr] > grid
                grid[mask] = (self._modelHrGridDict[tr])[mask]

        return grid    

    def getMinMaxGrids(self, modelName, weName, minTRList, maxTRList):
        """
        Fetches the min and max grids for the specified modelName and weather
        element.  Returns a dictionary of min grids and max grids, respectively.
        Both the min/max grids in the database and the hourly min/max grids are
        retrieved, if only one exists use that.  If both exist then optionally
        combine them so that we get the lowest min and highest max.
        """

        minWE = "Min" + weName
        maxWE = "Max" + weName
        # get the min grids
        minInv = self.getWEInventory(modelName, minWE)
        minGridDict = {}
        for tr in minTRList:
            # make a min grid from the hourlies
            grid = self.calcMinMaxFromHourlies(tr, "min")

            # Fetch the model min grid if we can
            modelGrid = None
            if tr in minInv:
                modelGrid = self.getGrids(modelName, minWE, "SFC", tr, mode="First")
                # Special exception for the Fcst database only
                if modelName == "Fcst":
                    minGridDict[tr] = modelGrid
                    continue

            if grid is None and modelGrid is None:  # no grids - keeping going
                continue

            if grid is None and modelGrid is not None:
                minGridDict[tr] = modelGrid  # use model grid if hourlies missing
                continue

            if modelGrid is None and grid is not None:
                minGridDict[tr] = grid  # use hourly min if model min missing
                continue

            # Both grids exist
            # Pick the min of the modelGrid and the hourly min
            if not self.configDict["useHourlyMinMax"]:
                mask = modelGrid < grid
                grid[mask] = (modelGrid)[mask]

            minGridDict[tr] = grid

        # get the max grids
        maxInv = self.getWEInventory(modelName, maxWE)
        maxGridDict = {}
        for tr in maxTRList:
            # make a max grid from the hourlies
            grid = self.calcMinMaxFromHourlies(tr, "max")

            # Fetch the model max grid if we can
            modelGrid = None
            if tr in maxInv:
                modelGrid = self.getGrids(modelName, maxWE, "SFC", tr, mode="First")
                # Special exception for the Fcst database only
                if modelName == "Fcst":
                    maxGridDict[tr] = modelGrid
                    continue

            if grid is None and modelGrid is None:  # no grids - keeping going
                continue

            if grid is None and modelGrid is not None:
                maxGridDict[tr] = modelGrid  # use model grid if hourlies missing
                continue

            if modelGrid is None and grid is not None:
                maxGridDict[tr] = grid  # use hourly max if model max missing
                continue

            # Both grids exist
            # Pick the max of the modelGrid and the hourly max
            if not self.configDict["useHourlyMinMax"]:
                mask = modelGrid > grid
                grid[mask] = (modelGrid)[mask]

            maxGridDict[tr] = grid

        return minGridDict, maxGridDict
        

    def calcFTrendGrids(self):
        """
        Using the specified grid list and timerange list, this method
        figures out for each timeRange the time-series trend of the
        grids (before/after min/max) and returns a grid list of the
        same size where 0 is trending down and 1 is trending up.
        """

        # make a dict where we will store the fTrends grids
        self._fTrendDict = {}

        # Process all of the Min timeranges

        trList = list(self._modelHrGridDict.keys())
        trList.sort(key=lambda a: a.startTime().unixTime())

        minTRList = list(self._modelMinGridDict.keys())
        minTRList.sort(key=lambda a: a.startTime().unixTime())

        # calculate the min value over each min timeRange
        for minTR in minTRList:
            minGrid = self.newGrid(9999.0)  # init running min

            # a grid that represents the time of the min value
            timeOfMin = self.empty()

            # Loop over each hourly timeRange to find the time of the min
            for tr in trList:
                if minTR.overlaps(tr):

                    # Make a mask at points where hourly min is less than running min
                    mask = self._modelHrGridDict[tr] <= minGrid

                    # update the running min
                    minGrid[mask] = self._modelHrGridDict[tr][mask]

                    # update the time of the min
                    timeOfMin[mask] = tr.startTime().unixTime()

            for tr in trList:
                if minTR.overlaps(tr):
                    start = tr.startTime().unixTime()
                    self._fTrendDict[tr] = np.where(start <= timeOfMin, -1, 1)

                    #################################################################
                    if self.configDict["DEBUG"]:
                        self.createGrid(
                            "Fcst",
                            "fTrend",
                            "SCALAR",
                            self._fTrendDict[tr].astype(np.float32),
                            tr,
                        )

        # Process all of the Max timeranges
        maxTRList = list(self._modelMaxGridDict.keys())
        maxTRList.sort(key=lambda a: a.startTime().unixTime())

        # calculate the max value over each max timeRange
        for maxTR in maxTRList:
            maxGrid = self.newGrid(-9999.0)  # init running max

            # a grid that represents the time of the max value
            timeOfMax = self.empty()

            # Loop over each hourly timeRange to find the time of the max
            for tr in trList:
                if maxTR.overlaps(tr):

                    # Make a mask at points where hourly max is np.less than running max
                    mask = self._modelHrGridDict[tr] >= maxGrid

                    # update the running max
                    maxGrid[mask] = self._modelHrGridDict[tr][mask]

                    # update the time of the max
                    timeOfMax[mask] = tr.startTime().unixTime()

            for tr in trList:
                if maxTR.overlaps(tr):
                    start = tr.startTime().unixTime()
                    self._fTrendDict[tr] = np.where(start >= timeOfMax, -1, 1)

                    #########################################################################
                    if self.configDict["DEBUG"]:
                        self.createGrid(
                            "Fcst",
                            "fTrend",
                            "SCALAR",
                            self._fTrendDict[tr].astype(np.float32),
                            tr,
                        )

        return

    def calcAdjacentTimeRanges(self, minOrMax, timeRange):
        """
        This method returns either the index of the overlapping min or max
        timeRange list or returns the indexes on either side of the specified
        timeRange.  "minOrMax must have the values "min" or "max".  Other
        values will return no indexes.  All indexes are returned as a list
        containing either one index or two.  If the specified timeRange is
        before the first or after the last in the min or max timeRange list,
        and empty is returned.
        """

        # get the right list
        if minOrMax == "min":
            trList = list(self._modelMinGridDict.keys())
        elif minOrMax == "max":
            trList = list(self._modelMaxGridDict.keys())
        else:
            return []

        trList.sort(key=lambda a: a.startTime().unixTime())
        # Check for before start and after end cases
        if timeRange.endTime().unixTime() < trList[0].startTime().unixTime() - (
            6 * 3600
        ):
            return []
        elif timeRange.startTime().unixTime() >= trList[-1].endTime().unixTime():
            return [trList[-1], trList[-1]]

        for i in range(len(trList)):
            # if there's an overlap, return the current index
            if trList[i].overlaps(timeRange):
                return [trList[i], trList[i]]
            # Since the trList is in order, as soon as we're past timeRange we stop
            elif trList[i].startTime().unixTime() >= timeRange.endTime().unixTime():
                return [trList[i - 1], trList[i]]

        # we should never get here
        print("Oh NO!!!!  we got here anyway!!!!!!")
        return []
        

    def calcDiurnalGrids(self, weName, executeTR, maskArea, fullMaskArea):
        """
        This method does the calculation of the hourly grids based on the model min/max
        model trends (fTrends) and the Fcst min/max.
        """

        trList = list(self._modelHrGridDict.keys())
        trList.sort(key=lambda a: a.startTime().unixTime())

        # if no bits set, set mask to entire area
        if not maskArea.any():
            maskArea = self.newGrid(True, np.bool)

        minLimit, maxLimit = self.getParmMinMaxLimits("Fcst", weName)
        for tr in trList:
            # Only generate grid within the specified timeRange
            if not executeTR.overlaps(tr):
                continue

            # find the indexes that overlap or are on either side of this tr
            minTRList = self.calcAdjacentTimeRanges("min", tr)
            maxTRList = self.calcAdjacentTimeRanges("max", tr)

            # an empty list means we're too far away from a min or max
            # grid to use it so skip this timeRange
            if not minTRList or not maxTRList:
                continue

            # Poke in model min/max values based on the fTrends grid
            # for this hour.
            minModelGrid = np.where(
                self._fTrendDict[tr] == 1,
                self._modelMinGridDict[minTRList[0]],
                self._modelMinGridDict[minTRList[1]],
            )

            maxModelGrid = np.where(
                self._fTrendDict[tr] == 1,
                self._modelMaxGridDict[maxTRList[1]],
                self._modelMaxGridDict[maxTRList[0]],
            )

            #####################################################################
            if self.configDict["DEBUG"]:
                self.createGrid(
                    "Fcst",
                    weName + "modelMinUsed",
                    "SCALAR",
                    minModelGrid.astype(np.float32),
                    tr,
                    minAllowedValue=-100.0,
                    maxAllowedValue=150.0,
                )
                self.createGrid(
                    "Fcst",
                    weName + "modelMaxUsed",
                    "SCALAR",
                    maxModelGrid.astype(np.float32),
                    tr,
                    minAllowedValue=-100.0,
                    maxAllowedValue=150.0,
                )

            # Make sure that minGrid < maxGrid everywhere
            mask = maxModelGrid <= minModelGrid
            minModelGrid[mask] = (maxModelGrid - 0.5)[mask]

            # calculate the f-value
            f = (self._modelHrGridDict[tr] - minModelGrid) / (
                maxModelGrid - minModelGrid
            )

            minFcstGrid = np.where(
                self._fTrendDict[tr] == 1,
                self._fcstMinGridDict[minTRList[0]],
                self._fcstMinGridDict[minTRList[1]],
            )

            maxFcstGrid = np.where(
                self._fTrendDict[tr] == 1,
                self._fcstMaxGridDict[maxTRList[1]],
                self._fcstMaxGridDict[maxTRList[0]],
            )

            ########################################################################
            # Display the Fcst Min and Max used for each tr in the final computation
            if self.configDict["DEBUG"]:
                self.createGrid(
                    "Fcst",
                    weName + "FcstMinUsed",
                    "SCALAR",
                    minFcstGrid.astype(np.float32),
                    tr,
                )
                self.createGrid(
                    "Fcst",
                    weName + "FcstMaxUsed",
                    "SCALAR",
                    maxFcstGrid.astype(np.float32),
                    tr,
                )

            # Linear function
            mask = f > 1.0
            f[mask] = (self._modelHrGridDict[tr] - maxModelGrid + 1.0)[mask]
            mask = f < 0.0
            f[mask] = (self._modelHrGridDict[tr] - minModelGrid)[mask]

            # Calculate and display the hourly grid
            grid = minFcstGrid + (maxFcstGrid - minFcstGrid) * f

            # For linear method  *************************************
            mask = f > 1.0
            grid[mask] = (maxFcstGrid + (f - 1.0))[mask]
            mask = f < 0.0
            grid[mask] = (minFcstGrid + f)[mask]

            if self.configDict["DEBUG"]:
                self.createGrid(
                    "Fcst", "fValue", "SCALAR", f.astype(np.float32), tr, precision=2
                )

            # Clip the grid - just in case
            grid.clip(minLimit, maxLimit, grid)

            # A bug in the GFE requires (re)setting all grids to gray
            self.highlightGrids("Fcst", weName, "SFC", tr, "gray50")

            try:
                oldGrid = self.getGrids("Fcst", weName, "SFC", tr, mode="First")
                oldArea = ~maskArea
                grid[oldArea] = (oldGrid)[oldArea]
            except:
                pass

            # Check for custom smooth area(s)
            if self._smoothAreas:
                for ea in self._smoothAreas:
                    weight = self._smoothWtDict.get(ea, 0)
                    if weight:
                        smMask = self.encodeEditArea(ea)
                        thisMask = smMask & fullMaskArea
                        if not thisMask.any():
                            print("{} does not intersect current edit area", format(ea))
                            continue

                        grid = self.smoothGrid(grid, weight, thisMask)

            self.createGrid("Fcst", weName, "SCALAR", grid.astype(np.float32), tr)

        return

    # This method replace any configuration items with those found in the DiurnalConfig.
    def replaceConfigValues(self):

        try:
            DiurnalConfig.configDict
        except:
            msg = "No config dictionary found in DiurnalConfig. Using base defaults."
            self.statusBarMsg(msg, "S")
            return

        # Replace the values for each item found
        for k in DiurnalConfig.configDict.keys():
            if k in self.configDict:
                self.configDict[k] = DiurnalConfig.configDict[k]

        return

    def overwriteWithClimo(self, primaryModel, executeTR, editArea):
        """
        In cases where the primary model selected runs out of grids before the
        selected timeRange and the user selects climatology as the secondary
        model, we substitute with the climatology algorithm.
        """

        trList = list(self._highlightDict.keys())
        trList.sort(key=lambda a: a.startTime().unixTime())

        start = -1
        end = -1

        for tr in trList:
            if not executeTR.overlaps(tr):
                continue

            color, modelName = self._highlightDict[tr]
            # extract the model name
            parts = modelName.split()
            model = parts[-1]
            if model != primaryModel:
                if start == -1:
                    start = tr.startTime().unixTime()  # save the startTime
                end = tr.endTime().unixTime()
                # reset the highlightDict entry for this timeRange
                self._highlightDict[tr] = ("blue", "Climatology")

                # Adjust start or end times based on execute time
                if start < executeTR.startTime().unixTime():
                    start = executeTR.startTime().unixTime()
                if end > executeTR.endTime().unixTime():
                    end = executeTR.endTime().unixTime()

        if start == -1:  # All primary models found so don't do anything
            return

        climoTR = self.makeTimeRange(start, end)

        # construct a varDict and call DiurnalFromClimo over this timeRange
        varDict = {}
        varDict["Beginning"] = "Selected TimeRange"
        varDict["Ending"] = "Selected TimeRange"

        self.callProcedure(
            "DiurnalTFromClimo", varDict=varDict, timeRange=climoTR, editArea=editArea
        )

        return
    
    def fTrend(self, timeRange):
        """
        Uses the fTable to figure out the trend at the current month/hour
        Returns 1 if the trend is increasing, -1 if the trend is decreasing
        """

        # get the hour and month as integers
        start = timeRange.startTime().unixTime()
        hour = time.gmtime(start).tm_hour
        month = time.gmtime(start).tm_mon - 1

        # fetch the month "column" using array slicing
        hourList = np.array(self.configDict.get("fTable"))[:, month]

        # Figure out the last and next hour.  Make sure they're all 0->23
        lastHour = hour - 1
        if lastHour < 0:
            lastHour = lastHour + 24

        nextHour = hour + 1
        if nextHour >= 24:
            nextHour = nextHour - 24

        # Calculate the trend
        trend = cmp(hourList[nextHour] - hourList[hour], 0)

        # If the trend is flat, use last hour to decide
        if trend == 0:
            trend = cmp(hourList[lastHour] - hourList[hour], 0)

        return trend

    def cmp(self, x, y):
        """
        Replacement for built-in function cmp that was removed in Python 3.
        """

        if x < y:
            return -1
        elif x > y:
            return 1
        else:
            return 0

    def makeDateStr(self, dayTime):
        
        gTime = time.gmtime(dayTime)
        year = gTime.tm_year
        month = gTime.tm_mon
        day = gTime.tm_mday
        timeStr = str(year).zfill(2) + "/" + str(month).zfill(2) + "/" + str(day).zfill(2)

        return timeStr

    def destroyDialog(self):  
        if AWIPS_ENVIRON == 'AWIPS1':
            self._master.destroy()
        else:
            self._master.destroy()
        return

    def cancelCommand(self):
        self._continue = False
        self.destroyDialog()
        return

    def continueCommand(self):
        self._continue = True
        self.destroyDialog()
        return
    
    # Pop a dialog to see if user wishes to continue
    def userWantsToContinue(self):

        variableList = []
        variableList.append(("Small Edit Area Detected. Continue?", "No",
                             "radio", ["Yes", "No"]))
        myVarDict = {}
        processVarList = ProcessVariableList.ProcessVariableList(
            "Continue?", variableList, myVarDict) 
        status = processVarList.status()
        if status.upper() != "OK": 
            self.cancel()
            return False

        answer = myVarDict["Small Edit Area Detected. Continue?"]
        if answer == "Yes":
            return True

        return False

    #===========================================================================
    #
    #  loadIfNot - check to see if the specified parm (from the mutable model)
    #    is loaded in GFE.  If not - go ahead and load it.
    #
    def loadIfNot(self, parmName, parmLevel="SFC"):
        mutID=self.mutableID()
        loadedInfo=self.loadedParms()
        for loaded in loadedInfo:
            (name,level,id)=loaded
            if ((name==parmName)and(level==parmLevel)and(id==mutID)):
                return
        #
        self.loadParm(mutID,parmName,parmLevel)
        return
    #=========================================================================
    #
    #  makeTd - make all the Td grids that where a T and RH grid can be found
    #           in the fcsttr timerange.
    #
    def makeTd(self,dbid,fcsttr):
        #self.logtime("Making Td grids")
        allTInfo=self.getGridInfo(dbid,"T","SFC",fcsttr)
        #self._SmartScript__pythonGrids=[]
        #self.unCacheElements(["T","RH"])
        if len(allTInfo)>0:
            cnt=0
            for info in allTInfo:
                tr=info.gridTime()
                cnt=cnt+1
                #begin=tr.startTime().unixTime()
                #ending=tr.endTime().unixTime()
                #cachekey="%s-%d-%d"%("T",begin,ending)
                #if cachekey in self.mycache:
                #   T=self.mycache[cachekey]
                #   #self.logtime("Got T for %s from cache"%tr)
                #else:
                T=self.getGrids(dbid,"T","SFC",tr,noDataError=0)
                if T is None:
                    continue
                #cachekey="%s-%d-%d"%("RH",begin,ending)
                #if cachekey in self.mycache:
                #   RH=self.mycache[cachekey]
                #   #self.logtime("Got RH for %s from cache"%tr)
                #else:
                RH=self.getGrids(dbid,"RH","SFC",tr,noDataError=0)
                if RH is None:
                    continue

                Tc = .556 * (T - 32.0)
                rh = np.clip(RH, 0.001, 99.999) / 100.0
                x = (np.log(rh) / 17.67) + (Tc / (Tc + 243.5))
                tdc = (243.5 * x) / (1 - x)
                Td = (1.8 * tdc) + 32.0

                Td=np.clip(Td,-80.0,120.0)
                self.createGrid(dbid,"Td","SCALAR",Td,tr)
                #cachekey="Td-%d-%d"%(begin,ending)
                #self.mycache[cachekey]=Td
                #self.callSmartTool("DoNothing","Td",None,tr) 
        return
    #=========================================================================
    #
    #  makeHI - make all the HeatIndex grids where a T and Td grid can be found
    #           in the fcsttr timerange.
    #
    def makeHI(self,dbid,fcsttr):
        #self.logtime("Making HeatIndex grids")
        allTInfo=self.getGridInfo(dbid,"T","SFC",fcsttr)
        #self._SmartScript__pythonGrids=[]
        #self.unCacheElements(["T","Td"])
        if len(allTInfo)>0:
            cnt=0
            for info in allTInfo:
                tr=info.gridTime()
                #hou=tr.startTime().hour
                #if ((hou==0)or(cnt==0)):
                #    day=tr.startTime().day
                #    mon=tr.startTime().month
                #    #self.logtime("   day %d/%d"%(mon,day))
                cnt=cnt+1
                #begin=tr.startTime().unixTime()
                #ending=tr.endTime().unixTime()
                #cachekey="%s-%d-%d"%("T",begin,ending)
                #if cachekey in self.mycache:
                #   T=self.mycache[cachekey]
                #   #self.logtime("Got T for %s from cache"%tr)
                #else:
                T=self.getGrids(dbid,"T","SFC",tr,noDataError=0)
                if T is None:
                    continue
                #
                #  HeatIndex only valid when T greater than 80
                #  so check for entire grid less than 80 degrees
                #  saving time of doing the calculations
                #
                maxT=maximum.reduce(maximum.reduce(T))
                if maxT<80:
                    HeatIndex=T
                    self.createGrid(dbid,"HeatIndex","SCALAR",HeatIndex,tr)
                    continue
                
                #cachekey="%s-%d-%d"%("Td",begin,ending)
                #if cachekey in self.mycache:
                #   Td=self.mycache[cachekey]
                #   #self.logtime("Got Td for %s from cache"%tr)
                #else:
                Td=self.getGrids(dbid,"Td","SFC",tr,noDataError=0)
                if Td is None:
                    continue
                Tc = .556 * (T - 32.0)
                Tdc = .556 * (Td - 32.0)
                #print "shape of T is %s"%str(T.shape)
                #print "shape of Tc is %s"%str(Tc.shape)
                #print "shape of Tdc is %s"%str(Tdc.shape)
                Vt = 6.11 * power(10,(Tc * 7.5 / (Tc + 237.3)))
                Vd = 6.11 * power(10,(Tdc * 7.5 / (Tdc + 237.3)))   
                RH = (Vd / Vt) * 100.0
    
                A = -42.379
                B =  2.04901523 * T
                C = 10.14333127 * RH
                D = -0.22475541 * T * RH
                E = -0.00683783 * power(T, 2)
                F = -0.05481717 * power(RH, 2)
                G =  0.00122874 * power(T, 2) * RH
                H =  0.00085282 * T * power(RH, 2)
                I = -0.00000199 * power(T, 2) * power(RH, 2)
    
                HeatIndex = A + B + C + D + E + F + G + H + I

                # make the adjustments for low humidity
                rhLessThan13 = less(RH, 13.0)
                T80to112 = logical_and(greater_equal(T, 80), less_equal(T, 112))
                downMask = logical_and(rhLessThan13, T80to112)

                # make array that is T where conditions are true and 100, otherwise
                adjustT = where(downMask, T, 80.0)
                HeatIndex = where(downMask, HeatIndex - (((13.0 - RH) / 4.0) * \
                              sqrt((17.0 - abs(adjustT - 95.0)) / 17)),
                              HeatIndex)
 
                # make the adjustments for high humidity
                rhGreater85 = greater(RH, 85.0)
                T80to87 = logical_and(greater_equal(T, 80.0), less_equal(T, 87.0))
                HeatIndex = where(logical_and(rhGreater85, T80to87),
                          HeatIndex + (((RH - 85.0) / 10.0) * ((87 - T) / 5.0)),
                          HeatIndex)
                                   
                HeatIndex = where(less(T, 80), T, HeatIndex)
                self.createGrid(dbid,"HeatIndex","SCALAR",HeatIndex,tr)
        return
    #=========================================================================
    #
    #  makeWC - make all the WindChill grids where a T and Wind grid can be found
    #           in the fcsttr timerange.
    #
    def makeWC(self,dbid,fcsttr):
        #self.logtime("Making WindChill grids")
        allTInfo=self.getGridInfo(dbid,"T","SFC",fcsttr)
        #self._SmartScript__pythonGrids=[]
        #self.unCacheElements(["T","Wind"])
        if len(allTInfo)>0:
            cnt=0
            for info in allTInfo:
                tr=info.gridTime()
                #hou=tr.startTime().hour
                #if ((hou==0)or(cnt==0)):
                #    day=tr.startTime().day
                #    mon=tr.startTime().month
                #    #self.logtime("   day %d/%d"%(mon,day))
                cnt=cnt+1
                #begin=tr.startTime().unixTime()
                #ending=tr.endTime().unixTime()
                #cachekey="%s-%d-%d"%("T",begin,ending)
                #if cachekey in self.mycache:
                #   T=self.mycache[cachekey]
                #   #self.logtime("Got T for %s from cache"%tr)
                #else:
                T=self.getGrids(dbid,"T","SFC",tr,noDataError=0)
                if T is None:
                    continue
                #
                #   Windchill only valid when T<51 so see if all
                #   T are above this level
                #
                minT=minimum.reduce(minimum.reduce(T))
                if minT>=50.95:
                    WindChill=T
                    clip(WindChill, -120.0, 120.0, WindChill)
                    self.createGrid(dbid,"WindChill","SCALAR",WindChill,tr)
                    continue
                Wind=self.getGrids(dbid,"Wind","SFC",tr,noDataError=0)
                if Wind is None:
                    continue
                #
                mag = Wind[0] * 1.152
                WindChill = where(less_equal(mag, 3.0), T, 35.74 + (0.6215 * T) - 
                           (35.75 * (mag ** 0.16)) + (0.4275 * T * (mag ** 0.16)))
                #
                # clip values where WindChill > T
                #
                WindChill = where(greater(WindChill, T), T, WindChill)
                #
                # substitute the temperature if WindChill >= 51 degrees
                #
                WindChill = where(greater_equal(T, 50.95), T, WindChill)
                #
                #  Clip to WindChill range
                #
                clip(WindChill, -120.0, 120.0, WindChill)
                #
                #
                #
                self.createGrid(dbid,"WindChill","SCALAR",WindChill,tr)
        return
    #=========================================================================
    #
    #  makeAT - make all the Apparent T grids
    #
    def makeAT(self,dbid,fcsttr):
        allTInfo=self.getGridInfo(dbid,"T","SFC",fcsttr)
        if len(allTInfo)>0:
            cnt=0
            for info in allTInfo:
                tr=info.gridTime()
                cnt=cnt+1
                Wind=self.getGrids(dbid,"Wind","SFC",tr,noDataError=0)
                T=self.getGrids(dbid,"T","SFC",tr,noDataError=0)
                #Td=self.getGrids(dbid,"Td","SFC",tr,noDataError=0)
                RH=self.getGrids(dbid,"RH","SFC",tr,noDataError=0)
                if T is None and RH is None:
                    continue
                mag = Wind[0] * 1.15
                WindChillValue = np.copy(T)
                mask = mag > 3
                WindChillValue[mask] = (35.74 + (0.6215 * T) - (35.75 * (mag ** 0.16)) + (0.4275 * T * (mag ** 0.16)))[mask]
                mask = WindChillValue > T
                WindChillValue[mask] = T[mask]
        
                # heat index
                A = -42.379
                B = 2.04901523 * T
                C = 10.14333127 * RH
                D = -0.22475541 * T * RH
                E = -0.00683783 * T ** 2
                F = -0.05481717 * RH ** 2
                G = 0.00122874 * T ** 2 * RH
                H = 0.00085282 * T * RH ** 2
                I = -0.00000199 * T ** 2 * RH ** 2
        
                HeatIndexValue = A + B + C + D + E + F + G + H + I
        
                # make the adjustments for low humidity
                rhLessThan13 = RH < 13.0
                T80to112 = ((T >= 80) & (T <= 112))
                downMask = (rhLessThan13 & T80to112)
        
                # make array that is T where conditions are true and 100, otherwise
                adjustT = np.copy(T)
                adjustT[downMask] = 80.0
                HeatIndexValue[downMask] = (HeatIndexValue - (((13.0 - RH) / 4.0) * np.sqrt((17.0 - abs(adjustT - 95.0)) / 17)))[downMask]
        
                # make the adjustments for high humidity
                rhGreater85 = RH > 85.0
                T80to87 = ((T >= 80.0) & (T <= 87.0))
                mask = (rhGreater85 & T80to87)
                HeatIndexValue[mask] = (HeatIndexValue + (((RH - 85.0) / 10.0) * ((87 - T) / 5.0)))[mask]
        
                ApparentT = np.copy(T)
                mask = T < 51
                ApparentT[mask] = WindChillValue[mask]
                mask = T > 79
                ApparentT[mask] = HeatIndexValue[mask]
        
                # ensure ApparentT does not exceed the bounds of the grid. Important in service backup
                ApparentT[ApparentT < -120] = -120
                ApparentT[ApparentT > 140] = 140

                self.createGrid(dbid,"ApparentT","SCALAR",ApparentT,tr)
        return
    
    def maxminAT(self,dbid,fcsttr):   
        mintrList = self.getGridInfo(dbid,"MinApT","SFC",fcsttr)
        maxtrList = self.getGridInfo(dbid,"MaxApT","SFC",fcsttr)
        for info in mintrList:
            tr=info.gridTime()
            minAT=self.getGrids(dbid,"ApparentT","SFC",tr,"Min",noDataError=0)
            self.createGrid(dbid,"MinApT","SCALAR",minAT,tr)
        for info in maxtrList:
            tr=info.gridTime()
            maxAT=self.getGrids(dbid,"ApparentT","SFC",tr,"Max",noDataError=0)
            self.createGrid(dbid,"MaxApT","SCALAR",maxAT,tr)
            
    def solarFlux(self, lats, startYear, startMonth, startDay):
        fEot, fR0r, fRdecl, fDeclsc1, fDeclsc2, nJulianDate = self.eqnOfTime(startYear, startMonth, startDay, lats)
        
        fSF = (fDeclsc1+fDeclsc2)*fR0r

        # In the case of a negative declination, solar flux is null
        fCoeff = self.empty()
        fMask = fSF >= 0
        fCoeff[fMask] = -1.56e-12*fSF[fMask]**4 + 5.972e-9*fSF[fMask]**3 - 8.364e-6*fSF[fMask]**2  + 5.183e-3*fSF[fMask] - 0.435
        
        fSFT = self.empty()
        fSFT = fSF * fCoeff 
        fSFTMask = fSFT < 0
        fSFT[fSFTMask] = 0

        return fSFT,fEot,fR0r,fRdecl, fDeclsc1,fDeclsc2, nJulianDate
    
    def eqnOfTime(self, startYear, startMonth, startDay, lats):
        #  date 
        nJulianDate, leap = self.Julian(startYear, startMonth, startDay)
        
        # Check if it is a leap year
        if(leap):
            fDivide = 366.0
        else:
            fDivide = 365.0
        # Correction for "equation of time"
        fA = nJulianDate/fDivide*2*pi
        
        # Calculates solar constant as a function of the Julian day
        fR0r = (1.0/(1.0-9.464e-4*np.sin(fA)-0.01671*np.cos(fA) - 1.489e-4*np.cos(2.0*fA)-2.917e-5*np.sin(3.0*fA) - 3.438e-4*np.cos(4.0*fA))**2) * 0.1367e4
        # Declination stuff
        fRdecl = 0.412*np.cos((nJulianDate+10.0)*2.0*pi/fDivide-pi)
        fDeclsc1 = np.sin(np.deg2rad(lats))*np.sin(fRdecl)
        fDeclsc2 = np.cos(np.deg2rad(lats))*np.cos(fRdecl)

        # in minutes
        fEotMin = 0.002733 -7.343*np.sin(fA)+ .5519*np.cos(fA) - 9.47*np.sin(2.0*fA) - 3.02*np.cos(2.0*fA) - 0.3289*np.sin(3.*fA) -0.07581*np.cos(3.0*fA) - 0.1935*np.sin(4.0*fA) - 0.1245*np.cos(4.0*fA)
        # Express in fraction of hour
        fEot = fEotMin/60.0
        # Express in radians
        fEot = fEot*15*pi/180.0

        return fEotMin, fR0r, fRdecl, fDeclsc1, fDeclsc2, nJulianDate
    
    def Julian(self, startYear, startMonth, startDay):
    
        if calendar.isleap(startYear): # Bissextil year, 366 days
            lMonth = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
            leap = 1
        else: # Normal year, 365 days
            lMonth = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
            leap = 0

        nJulian = lMonth[startMonth-1] + startDay
        return nJulian, leap

    def sunrSuns(self,startYear,startMonth,startDay,lons,lats):
        # Compute day of mean solar time
        daysSince2000Jan0 = (367*(startYear)-((7*((startYear)+(((startMonth)+9)/12)))/4)+((275*(startMonth))/9)+(startDay)-730530)
        d = daysSince2000Jan0 + 0.5 - (lons/360.0)        
        # Compute sidereal time of this moment
        revInGMT = ((180.0 + 356.0470 + 282.9404) + (0.9856002585 + 4.70935E-5) * d) + 180.0 + lons
        sidTimeGMT = revInGMT - 360.0 * np.floor(revInGMT * (1.0 / 360.0)) 
        # Compute sun's RA and Declination at this moment
        # Compute mean elements
        revIn2 = 356.0470 + 0.9856002585 * d
        M = revIn2 - 360.0 * np.floor(revIn2 * (1.0/360.0))
        w = 282.9404 + 4.70935E-5 * d
        e = 0.016709 - 1.151E-9 * d
        
        # Compute true longitude and radius vector 
        E = M + e * np.rad2deg(1) * np.sin(np.deg2rad(M)) * (1.0 + e * np.cos(np.deg2rad(M)))
        x = np.cos(np.deg2rad(E)) - e
        y = np.sqrt(1.0 - np.square(e)) * np.sin(np.deg2rad(E))
        # Solar distance
        r = np.sqrt(np.square(x) + np.square(y))
        # True Anomaly 
        v = np.arctan2(y,x) * np.rad2deg(1)
        # True solar longitude 
        solLons = v + w
        # Make it 0 - 360 degrees.
        solLonMask = solLons >= 360.0
        solLons[solLonMask] = solLons[solLonMask] - 360.0
        
        # Compute ecliptic rectangular coordinates (z=0) 
        x2 = r * np.cos(np.deg2rad(solLons))
        y2 = r * np.sin(np.deg2rad(solLons))
    
        # Compute obliquity of ecliptic (inclination of Earth's axis) 
        obl_ecl = 23.4393 - 3.563E-7 * d
        
        # Convert to equatorial rectangular coordinates - x is unchanged 
        z = y2 * np.sin(np.deg2rad(obl_ecl))
        y2 = y2 * np.cos(np.deg2rad(obl_ecl))

        # Convert to spherical coordinates 
        RA = np.arctan2(y2,x2) * np.rad2deg(1)
        dec = np.arctan2(z, np.sqrt(np.square(x2) + np.square(y2))) * np.rad2deg(1)
        
        # Compute time when Sun is at south - in hours UT 
        tsouth = 12.0-((sidTimeGMT-RA)-360.0*np.floor((sidTimeGMT-RA)*(1/360.0)+0.5))/15.0
        # Compute the Sun's apparent radius, degrees 
        sradius = 0.2666 / r
        # Do correction to upper limb, if necessary
        upper_limb = 1 
        altit = -35.0/60.0
        if upper_limb:
            altit = altit - sradius
    
        # Compute the diurnal arc that the Sun traverses to reach the specified altitude altit: 
        cost = (np.sin(np.deg2rad(altit))-np.sin(np.deg2rad(lats))*np.sin(np.deg2rad(dec)))/(np.cos(np.deg2rad(lats))*np.cos(np.deg2rad(dec)))
        # Sun always below altit
        costMaskA = cost >= 1.0
        # Sun Always above altit
        costMaskB = cost <= -1.0
        costMaskC = ~costMaskA & ~costMaskB
        t = self.empty()
        t[costMaskB] = 12.0
        t[costMaskC] = (np.arccos(cost[costMaskC])*np.rad2deg(1))/15.0
        # Sunrise and sunset time in hours utc
        sunriseArr = self.empty()
        sunsetArr = self.empty()
        sunriseArr = tsouth-t
        sunsetArr = tsouth+t
        return sunriseArr,sunsetArr
    
    ###################################################################################################################################################################################################
    
    #######################################
    ######### Main WBGT Methods ########### 
    #######################################
    
    def localSolarFluxAdj(self,sunsetArrGMT,sunriseArrGMT,solarnoon_GMT,startHour,daylenArr,maxflux):
        # Turn the start hour into an gridded array for comparison
        startHourArr = self.newGrid(startHour)
        # boolean to check if nighttime
        is_night = self.empty(dtype=np.bool)
        # Copy the sunsetArr for adjusting
        coSunset = copy.deepcopy(sunsetArrGMT)
        # Subtract 24 hours from the sunset array to get small positive values past 24 or 00z
        coSunset -= 24.
        # Create a mask where sunset is past 00z
        past00zMask = coSunset > 0.
        # Create a start hour mask where sunset is past 00z and the start hour is prior to sunset
        startHourMask = past00zMask & (startHourArr < coSunset)
        # Add 24 hours to the start hour in this particular case (i.e. 00z start hour is 24)
        startHourArr[startHourMask] += 24.
        # If it's past 00z compare the single start hour to the small past 00z sunset hour and the sunrise hour
        is_night[past00zMask] = (float(startHour) > coSunset[past00zMask]) & (float(startHour) < sunriseArrGMT[past00zMask])
        # If it's not past 00z compare the single start hour to the <24 sunset value and the sunrise hour
        is_night[~past00zMask] = (float(startHour) > sunsetArrGMT[~past00zMask]) | (float(startHour) < sunriseArrGMT[~past00zMask])
        # if not nighttime, find distance in hours from solar noon
        stdDist = self.empty()
        distFromMax = self.empty()
        # Used the starthourArr here because if sunset is past 00z then we need to use the start hour higher value (i.e. instead of using the 00z start hour if sunset is past 00z we need to use 24)
        distFromMax[~is_night] = np.absolute(solarnoon_GMT[~is_night] - startHourArr[~is_night])
        # standardize distance from peak to current hour as percent of distance between solar noon and sunrise/sunset
        stdDist[~is_night] = 1.0 - (distFromMax[~is_night] / (daylenArr[~is_night]/2.0))        
        # find location of stdDist on normal distribution as % of peak value
        distDict = {
            'A':{'crit':stdDist>=0.9,'ep1':1.0,'ep2':0.9,'pcnt1':1.00,'pcnt2':0.97},
            'B':{'crit':(stdDist>=0.8)&(stdDist<0.9),'ep1':0.9,'ep2':0.8,'pcnt1':0.97,'pcnt2':0.86},
            'C':{'crit':(stdDist>=0.7)&(stdDist<0.8),'ep1':0.8,'ep2':0.7,'pcnt1':0.86,'pcnt2':0.66},
            'D':{'crit':(stdDist>=0.6)&(stdDist<0.7),'ep1':0.7,'ep2':0.6,'pcnt1':0.66,'pcnt2':0.48},
            'E':{'crit':(stdDist>=0.5)&(stdDist<0.6),'ep1':0.6,'ep2':0.5,'pcnt1':0.48,'pcnt2':0.31},
            'F':{'crit':(stdDist>=0.4)&(stdDist<0.5),'ep1':0.5,'ep2':0.4,'pcnt1':0.31,'pcnt2':0.17},
            'G':{'crit':(stdDist>=0.3)&(stdDist<0.4),'ep1':0.4,'ep2':0.3,'pcnt1':0.17,'pcnt2':0.10},
            'H':{'crit':(stdDist>=0.2)&(stdDist<0.3),'ep1':0.3,'ep2':0.2,'pcnt1':0.10,'pcnt2':0.05},
            'I':{'crit':(stdDist>=0.1)&(stdDist<0.2),'ep1':0.2,'ep2':0.1,'pcnt1':0.05,'pcnt2':0.03},
            'J':{'crit':(stdDist>=0.0)&(stdDist<0.1),'ep1':0.1,'ep2':0.0,'pcnt1':0.03,'pcnt2':0.00}
        }
        sortedDistKeys = sorted(list(distDict.keys()))
        ep1 = self.empty()
        ep2 = self.empty()
        pcnt1 = self.empty()
        pcnt2 = self.empty()
        for distKey in sortedDistKeys:
            distMask = distDict[distKey]['crit']
            ep1[distMask] = distDict[distKey]['ep1']
            ep2[distMask] = distDict[distKey]['ep2']
            pcnt1[distMask] = distDict[distKey]['pcnt1']
            pcnt2[distMask] = distDict[distKey]['pcnt2']
        # determine approximate distance between endpoints (ep1 and ep2) with the following linear interpolation: y2= (x2-x1)(y3-y1)/(x3-x1) + y1
        pctOfMax = self.newGrid(1)
        pctOfMax[~is_night] = ((stdDist[~is_night] - ep1[~is_night]) * (pcnt2[~is_night] - pcnt1[~is_night]) / (ep2[~is_night] - ep1[~is_night])) + pcnt1[~is_night]
        if self.__pointDiagnostics:
            self.__testFile.write("sunsetArrGMT, sunriseArrGMT, daylenArr, distFromMax, stdDist, is_night\n")
            self.__testFile.write("{}, {}, {}, {}, {}, {}\n".format(sunsetArrGMT[self.__testCoords[1]][self.__testCoords[0]],sunriseArrGMT[self.__testCoords[1]][self.__testCoords[0]],daylenArr[self.__testCoords[1]][self.__testCoords[0]],distFromMax[self.__testCoords[1]][self.__testCoords[0]],stdDist[self.__testCoords[1]][self.__testCoords[0]],is_night[self.__testCoords[1]][self.__testCoords[0]]))
            self.__testFile.write("pcnt1,pcnt2,pctOfMax\n")
            self.__testFile.write("{}, {}, {}\n".format(pcnt1[self.__testCoords[1]][self.__testCoords[0]],pcnt2[self.__testCoords[1]][self.__testCoords[0]],pctOfMax[self.__testCoords[1]][self.__testCoords[0]]))
        # if nighttime, assign no max flux. otherwise multiply by pctOfMax to obtain final flux adjustment            
        maxflux[is_night] = 0
        maxflux[~is_night] = maxflux[~is_night] * pctOfMax[~is_night]
        return maxflux
    
    def getZenAng(self, fRdecl, fEot, lats, lons, startHour, thisSPressMB, thisTinC):
        # solar declination from the sunrise sunset stuff
        decl = fRdecl*(180./np.pi)
        # solar equation of time in minutes from sunrise sunset
        eqt = fEot
        # True solar time including a 'fix' using longitude and then converting to minutes
        toff = eqt+(4*lons)
        tst = startHour*60+toff
        # Hour Angle
        hang = (tst/4)-180
                
        # Cosine of the solar zenith angle
        cz =  np.cos(np.deg2rad(1) * decl) * np.cos(np.deg2rad(1) * lats) * np.cos(np.deg2rad(1) * hang) + np.sin(np.deg2rad(1) * decl) * np.sin(np.deg2rad(1) * lats)
        # handle round off errors
        czMask = np.absolute(cz) > 1.0
        if czMask.any():
            czMaskPos = cz[czMask] >= 0
            czMaskNeg = cz[czMask] < 0
            cz[czMaskPos] = 1.0
            cz[czMaskNeg] = -1.0
        zenetr = np.arccos(cz) * np.rad2deg(1)
        # Limit the degrees below the horizon to 9
        zenMask = zenetr > 99.0
        zenetr[zenMask] = 99.0
        elevetr = 90.0 - zenetr
        
        # If the sun is near zenith the algorithm bombs; refraction near 0
        elMaskA = elevetr > 85.0
        # Otherwise we have refraction
        elMaskB = ~elMaskA & (elevetr >= 5.0)
        elMaskC = ~elMaskA & ~elMaskB & (elevetr >= -0.575)
        tanelev = np.tan(np.deg2rad(1)*elevetr)
        refcor=self.empty()
        refcor[elMaskB] = (58.1 / tanelev[elMaskB]) - (0.07 / tanelev[elMaskB] ** 3) + (0.000086 / tanelev[elMaskB] ** 5)
        refcor[elMaskC] = 1735.0 + ( elevetr[elMaskC] * (-518.2 + elevetr[elMaskC] * ( 103.4 + elevetr[elMaskC] * ( -12.79 + elevetr[elMaskC] * 0.711 ) ) ) )
        refcor[~elMaskA & ~elMaskB & ~ elMaskC] = -20.774 / tanelev[~elMaskA & ~elMaskB & ~ elMaskC]
        prestemp = ( thisSPressMB * 283.0 ) / ( 1013.0 * ( 273.0 + thisTinC ) )
        refcor *= prestemp / 3600.0
        # Refracted solar elevation angle
        elevref = elevetr + refcor
        # Limit the degrees below the horizon to 9
        elMaskD = elevref < -9.0
        elevref[elMaskD] = -9.0
        # Refracted solar zenith angle
        zenref = 90.0 - elevref
        # To prevent incredibly small values for cosine result and thus a crazy value for max flux, constrain the zenith angle so that it never gets too close to 90 degrees
        zenMaskA = (zenref > 89.5) & (zenref <= 90.0)
        zenMaskB = (zenref > 90.0) & (zenref < 90.5)
        zenref[zenMaskA] = 89.5
        zenref[zenMaskB] = 90.5
        # Cosine of refracted solar zenith angle
        cos_zenith = np.cos(np.deg2rad(1)*zenref)
        if self.__pointDiagnostics:
            self.__testFile.write("zenith\n")
            self.__testFile.write("{}\n".format(zenref[self.__testCoords[1]][self.__testCoords[0]]))
            self.__testFile.write("cos_zenith\n")
            self.__testFile.write("{}\n".format(cos_zenith[self.__testCoords[1]][self.__testCoords[0]]))
        return zenref, cos_zenith
        
    def wbgtMainCalc(self,maxflux,thisSky,nJulianDate,startHour,lons,lats,thisTdinC,thisTinC,thisSPressMB,this2mWindSpd_mh,thisWindSpd_ms,thisT,eachGridTime,threatKeys,albedoGrid,fEot,fR0r,fRdecl,fDeclsc1,fDeclsc2):
        if self.__pointDiagnostics:
            self.__testFile.write("thisSky, thisTdinC, thisTinC, thisSPressMB, thisT, albedoGrid\n")
            self.__testFile.write("{}, {}, {}, {}, {}, {}\n".format(thisSky[self.__testCoords[1]][self.__testCoords[0]],thisTdinC[self.__testCoords[1]][self.__testCoords[0]],thisTinC[self.__testCoords[1]][self.__testCoords[0]],thisSPressMB[self.__testCoords[1]][self.__testCoords[0]],thisT[self.__testCoords[1]][self.__testCoords[0]],albedoGrid[self.__testCoords[1]][self.__testCoords[0]]))
        
        # Calculate the zenith angle
        zenref, cos_zenith = self.getZenAng(fRdecl,fEot,lats, lons, startHour, thisSPressMB, thisTinC)
        
        
        # estimate actual solar flux from max solar flux and sky cover - Kasten and Czeplak, 1980
        adjFlux = maxflux * (1. - (0.75*(np.power(thisSky/100.0, 3.4))))
        diffuse = thisSky/100.0 # diffuse rad
        direct = 1.0 - diffuse # direct rad
        if self.__pointDiagnostics:
            self.__testFile.write("adjFlux\n")
            self.__testFile.write("{}\n".format(adjFlux[self.__testCoords[1]][self.__testCoords[0]]))
        # Modify the original computation of direct vs diffuse proportions cap direct rad @ 0.75 (T. Blaine)
        directMask = direct > 0.75
        direct[directMask] = 0.75
        diffuse[directMask] = 0.25
        
        # Stefan-Boltzmann Constant
        sb = 5.67*(10**-8)
        
        # Blaine Thomas suggestion, vary the convective heat transfer 
        # coefficient between night and day, using the zenith angle to differentiate
        # Using a value of 0 at night and and 0.228 during the daytime seems to give the
        # closest approximation to the black globe temp    
        h = self.newGrid(0.228)
        zenithMaskC = zenref > 87.0
        h[zenithMaskC] = 0.0
        a0 = (17.67*(thisTdinC-thisTinC))/(thisTdinC+243.5)
        a = np.exp(a0)
        b0 = (17.502*(thisTinC))/(thisTinC + 240.97)
        b = np.exp(b0)
        ea = a * (1.0007 + 0.00000346 * thisSPressMB) * (6.112 * b)
        ea2 = 0.575 * np.power(ea,0.14285714)
        
        # Uses hires albedo data to improve upon equation from Dimiceli and Piltz
        b1 = direct / (4 * cos_zenith * sb) + (1 + albedoGrid / sb) * diffuse
        b2 = ea2 * np.power(thisTinC, 4)
        b = adjFlux * b1 + b2
        c = (this2mWindSpd_mh**(0.58)) * h / (5.3865 * (10**-8))
        
        # black globe temperature
        hMask = h == 0
        tg = self.empty()
        tg[hMask] = (b2[hMask]**0.25) * 1.8 + 32
        tg[~hMask] = ((b[~hMask] + c[~hMask] * thisTinC[~hMask] + 7680000) / (c[~hMask] + 256000)) * 1.8 + 32
                  
        # wet bulb temperature constants
        c1 = 0.0091
        c2 = 6106.40
        
        fp = 0.0006355 * thisSPressMB
        es = 6.11 * np.power(10, ((((thisTinC) * 7.5)) / ((thisTinC + 237.3))))
        ed = 6.11 * np.power(10, ((((thisTdinC) * 7.5)) / ((thisTdinC + 237.3))))
        s1 = (es-ed)
        s2 = (thisTinC - thisTdinC)
        
        s2Mask = s2 == 0
        Twc = self.empty()
        Twc[s2Mask] = thisTinC[s2Mask]
        Twc[~s2Mask] = ((thisTinC[~s2Mask] * fp[~s2Mask]) + (thisTdinC[~s2Mask] * (s1[~s2Mask]) / (s2[~s2Mask]) )) / (fp[~s2Mask] + (s1[~s2Mask]) / (s2[~s2Mask]))
        
        Tw = self.empty()    
        Tw = (1.8 * Twc) + 32
        
        # There's been some controversy wrt this method of computing the nwb temp Lilijgren uses true nwb while Dimiceli uses psychometric wb.
        # After testing the Lilijgen method, it was found to be too slow for the large size grids in NDFD (at least for now)
        for i in range(5):
            Twk = Twc + 273.15
            ew = 6.11 * np.power(10,((Twc * 7.5)/(Twc + 237.3)))
            de1 = fp * (thisTinC-Twc)
            de = de1-(ew-ed)
            der = ((ew*(c1-(c2/(Twk**2)))-fp))
            Twk = Twk-(de/der)
            Twc = Twk - 273.15
            Tw = (1.8 * Twc) + 32
        
        # Wet bulb depression in C
        Twbd = thisTinC - Twc
        
        # Use Twbd in a new way to calculate the natural wet bulb, per Blaine Thomas
        nwbAlt = self.empty()
        nwbAlt = Twc + (0.001651 * adjFlux) - (0.09555 * thisWindSpd_ms) + (0.13235 * Twbd) + 0.20249
        nwbAltMask = nwbAlt < Twc
        nwbAlt[nwbAltMask] = Twc[nwbAltMask]
        # Convert to F
        nwbAltinF = (1.8 * nwbAlt) + 32
        
        if self.__pointDiagnostics:
            self.__testFile.write("tg, nwbAltinF, thisT\n")
            self.__testFile.write("{}, {}, {}\n".format(tg[self.__testCoords[1]][self.__testCoords[0]],nwbAltinF[self.__testCoords[1]][self.__testCoords[0]],thisT[self.__testCoords[1]][self.__testCoords[0]]))
            self.__testFile.write("************************************************************\n")
        # WBGT final calculation
        wbgt = self.empty()
        wbgt =  (0.2 * tg) + (0.7 * nwbAltinF) + (0.1 * thisT)
        wbgtNaN = np.isnan(wbgt)
        wbgt[wbgtNaN] = -120.0
        np.clip(wbgt, -120.0, 140.0, out=wbgt) # because we're dumb
        self.createGrid("Fcst","WBGT","SCALAR",wbgt,eachGridTime)
        self.eaRegionParsing(wbgt,eachGridTime,threatKeys)
        return
        
    def eaRegionParsing(self,wbgt,eachGridTime,threatKeys):
        idxGrid = self.empty()
        idxes = ["low","elevated","moderate","high"]
        for eachTuple in self.validEAs:
            if eachTuple[1] is not None:
                setattr(self,eachTuple[0],self.encodeEditArea(eachTuple[1]))
                for num, eachIdx in enumerate(idxes):
                    idxGrid[(wbgt>self.regionsDict[eachTuple[0]][eachIdx])&(getattr(self,eachTuple[0]))] = num+1
            else:
                setattr(self,eachTuple[0],self.empty(dtype=np.bool))
        self.createGrid("Fcst","WBGTRisk","DISCRETE",(idxGrid,threatKeys),eachGridTime,discreteKeys=threatKeys,discreteOverlap=0,discreteAuxDataLength=2,defaultColorTable="Hazards")
        return
    
    def _findLongestModelDB(self,model,field,level,gfeTR):
        alldbs = self.availableDatabases()
        modeldbs = [x for x in alldbs if model in x.__str__()]
        lenCheck = 0
        for modeldb in modeldbs:
            ginfo = self.getGridInfo(modeldb,field,level,gfeTR)
            thisLen = len(ginfo) 
            if thisLen > lenCheck:
                lenCheck = thisLen
                theDB = modeldb
                dbTime = modeldb.modelTime()
                dbHour = int(dbTime.stringFmt("%H"))
                dbGInfo = ginfo
        return theDB, dbTime, dbHour, dbGInfo
    
    def _stitchDBs(self,db1,db1Info,db2,db2Info,field,level):
        db1End = db1Info[len(db1Info)-1].gridTime()
        db2Times = [x.gridTime() for x in db2Info if x.gridTime().__gt__(db1End)]
        db1Times = [x.gridTime() for x in db1Info]
        allTimes = db1Times + db2Times
        db1Grids = self.getGrids(db1,field,level,db1Times,noDataError=0)
        db2Grids = self.getGrids(db2,field,level,db2Times,noDataError=0)
        stitchedGrids = db2Grids.copy()
        stitchedGrids.update(db1Grids)
        return stitchedGrids,allTimes
    
    def _stretchedDB(self,db1,db1Info,field,level):
        extraData = {}
        db1Times = [x.gridTime() for x in db1Info]
        db1Grids = self.getGrids(db1,field,level,db1Times,noDataError=0)
        sortedTrs = sorted(db1Grids)
        for idx, validTr in enumerate(sortedTrs):
            if idx < len(sortedTrs)-1:
                thisEnd = validTr.endTime()
                nextStart = sortedTrs[idx+1].startTime()
                numGrids = (nextStart - thisEnd) / 3600
                extraCount = 0
                while extraCount < numGrids:
                    newStart = thisEnd
                    newEnd = thisEnd.__add__(3600)
                    newTR = TimeRange.TimeRange(newStart,newEnd)
                    extraData[newTR] = db1Grids[validTr]
                    thisEnd = newEnd
                    extraCount += 1
        hourlyGrids = db1Grids.copy()
        hourlyGrids.update(extraData)
        allTimes = sorted(list(hourlyGrids.keys()))
        fullTr = TimeRange.TimeRange(allTimes[0].startTime(),allTimes[len(allTimes)-1].endTime())
        return hourlyGrids,allTimes,fullTr
    
    def _createMaxGrids(self,fullTr,gfeTRList, threatKeys):
        threatDict = {
            "Low":0,
            "Elevated":1,
            "Moderate":2,
            "High":3,
            "Extreme":4,
            }
        self.createFromScratchCmd(["MaxWBGTRisk"],fullTr)
        self.fragmentCmd(["MaxWBGTRisk"],fullTr)
        indivRiskGrids = self.getGrids("Fcst", "WBGTRisk", "SFC", gfeTRList, noDataError=0)
        maxGridInfo = self.getGridInfo("Fcst","MaxWBGTRisk","SFC",fullTr)
        for eachMaxGrid in maxGridInfo:
            anyGrids = False
            compositeGrid = self.empty()
            eachGridTr = eachMaxGrid.gridTime()
            for gridKeyTR in list(indivRiskGrids.keys()):
                gridDict = {}
                if eachGridTr.overlaps(gridKeyTR):
                    if indivRiskGrids[gridKeyTR] is not None:
                        anyGrids = True
                        for fakeIdx,eachKey in reversed(list(enumerate(indivRiskGrids[gridKeyTR][1]))):
                            idxMask = indivRiskGrids[gridKeyTR][0] == fakeIdx
                            indivRiskGrids[gridKeyTR][0][idxMask] = threatDict[eachKey]
                        maxMask = indivRiskGrids[gridKeyTR][0] > compositeGrid
                        compositeGrid[maxMask] = indivRiskGrids[gridKeyTR][0][maxMask]
            if anyGrids:
                self.createGrid("Fcst","MaxWBGTRisk","DISCRETE",(compositeGrid,threatKeys),eachGridTr,discreteKeys=threatKeys,discreteOverlap=0,discreteAuxDataLength=2,defaultColorTable="Hazards")
            else:
                self.deleteGrid("Fcst", "MaxWBGTRisk", "SFC", eachGridTr)
        return
    
    def solarFlux(self, lats, startYear, startMonth, startDay):
        fEot, fR0r, fRdecl, fDeclsc1, fDeclsc2, nJulianDate = self.eqnOfTime(startYear, startMonth, startDay, lats)
        
        fSF = (fDeclsc1+fDeclsc2)*fR0r

        # In the case of a negative declination, solar flux is null
        fCoeff = self.empty()
        fMask = fSF >= 0
        fCoeff[fMask] = -1.56e-12*fSF[fMask]**4 + 5.972e-9*fSF[fMask]**3 - 8.364e-6*fSF[fMask]**2  + 5.183e-3*fSF[fMask] - 0.435
        
        fSFT = self.empty()
        fSFT = fSF * fCoeff 
        fSFTMask = fSFT < 0
        fSFT[fSFTMask] = 0

        return fSFT,fEot,fR0r,fRdecl, fDeclsc1,fDeclsc2, nJulianDate
    
    def eqnOfTime(self, startYear, startMonth, startDay, lats):
        #  date 
        nJulianDate, leap = self.Julian(startYear, startMonth, startDay)
        
        # Check if it is a leap year
        if(leap):
            fDivide = 366.0
        else:
            fDivide = 365.0
        # Correction for "equation of time"
        fA = nJulianDate/fDivide*2*pi
        
        # Calculates solar constant as a function of the Julian day
        fR0r = (1.0/(1.0-9.464e-4*np.sin(fA)-0.01671*np.cos(fA) - 1.489e-4*np.cos(2.0*fA)-2.917e-5*np.sin(3.0*fA) - 3.438e-4*np.cos(4.0*fA))**2) * 0.1367e4
        # Declination stuff
        fRdecl = 0.412*np.cos((nJulianDate+10.0)*2.0*pi/fDivide-pi)
        fDeclsc1 = np.sin(np.deg2rad(lats))*np.sin(fRdecl)
        fDeclsc2 = np.cos(np.deg2rad(lats))*np.cos(fRdecl)

        # in minutes
        fEotMin = 0.002733 -7.343*np.sin(fA)+ .5519*np.cos(fA) - 9.47*np.sin(2.0*fA) - 3.02*np.cos(2.0*fA) - 0.3289*np.sin(3.*fA) -0.07581*np.cos(3.0*fA) - 0.1935*np.sin(4.0*fA) - 0.1245*np.cos(4.0*fA)
        # Express in fraction of hour
        fEot = fEotMin/60.0
        # Express in radians
        fEot = fEot*15*pi/180.0

        return fEotMin, fR0r, fRdecl, fDeclsc1, fDeclsc2, nJulianDate
    
    def Julian(self, startYear, startMonth, startDay):
    
        if calendar.isleap(startYear): # Bissextil year, 366 days
            lMonth = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
            leap = 1
        else: # Normal year, 365 days
            lMonth = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
            leap = 0

        nJulian = lMonth[startMonth-1] + startDay
        return nJulian, leap

    def sunrSuns(self,startYear,startMonth,startDay,lons,lats):
        # Compute day of mean solar time
        daysSince2000Jan0 = (367*(startYear)-((7*((startYear)+(((startMonth)+9)/12)))/4)+((275*(startMonth))/9)+(startDay)-730530)
        d = daysSince2000Jan0 + 0.5 - (lons/360.0)        
        # Compute sidereal time of this moment
        revInGMT = ((180.0 + 356.0470 + 282.9404) + (0.9856002585 + 4.70935E-5) * d) + 180.0 + lons
        sidTimeGMT = revInGMT - 360.0 * np.floor(revInGMT * (1.0 / 360.0)) 
        # Compute sun's RA and Declination at this moment
        # Compute mean elements
        revIn2 = 356.0470 + 0.9856002585 * d
        M = revIn2 - 360.0 * np.floor(revIn2 * (1.0/360.0))
        w = 282.9404 + 4.70935E-5 * d
        e = 0.016709 - 1.151E-9 * d
        
        # Compute true longitude and radius vector 
        E = M + e * np.rad2deg(1) * np.sin(np.deg2rad(M)) * (1.0 + e * np.cos(np.deg2rad(M)))
        x = np.cos(np.deg2rad(E)) - e
        y = np.sqrt(1.0 - np.square(e)) * np.sin(np.deg2rad(E))
        # Solar distance
        r = np.sqrt(np.square(x) + np.square(y))
        # True Anomaly 
        v = np.arctan2(y,x) * np.rad2deg(1)
        # True solar longitude 
        solLons = v + w
        # Make it 0 - 360 degrees.
        solLonMask = solLons >= 360.0
        solLons[solLonMask] = solLons[solLonMask] - 360.0
        
        # Compute ecliptic rectangular coordinates (z=0) 
        x2 = r * np.cos(np.deg2rad(solLons))
        y2 = r * np.sin(np.deg2rad(solLons))
    
        # Compute obliquity of ecliptic (inclination of Earth's axis) 
        obl_ecl = 23.4393 - 3.563E-7 * d
        
        # Convert to equatorial rectangular coordinates - x is unchanged 
        z = y2 * np.sin(np.deg2rad(obl_ecl))
        y2 = y2 * np.cos(np.deg2rad(obl_ecl))

        # Convert to spherical coordinates 
        RA = np.arctan2(y2,x2) * np.rad2deg(1)
        dec = np.arctan2(z, np.sqrt(np.square(x2) + np.square(y2))) * np.rad2deg(1)
        
        # Compute time when Sun is at south - in hours UT 
        tsouth = 12.0-((sidTimeGMT-RA)-360.0*np.floor((sidTimeGMT-RA)*(1/360.0)+0.5))/15.0
        # Compute the Sun's apparent radius, degrees 
        sradius = 0.2666 / r
        # Do correction to upper limb, if necessary
        upper_limb = 1 
        altit = -35.0/60.0
        if upper_limb:
            altit = altit - sradius
    
        # Compute the diurnal arc that the Sun traverses to reach the specified altitude altit: 
        cost = (np.sin(np.deg2rad(altit))-np.sin(np.deg2rad(lats))*np.sin(np.deg2rad(dec)))/(np.cos(np.deg2rad(lats))*np.cos(np.deg2rad(dec)))
        # Sun always below altit
        costMaskA = cost >= 1.0
        # Sun Always above altit
        costMaskB = cost <= -1.0
        costMaskC = ~costMaskA & ~costMaskB
        t = self.empty()
        t[costMaskB] = 12.0
        t[costMaskC] = (np.arccos(cost[costMaskC])*np.rad2deg(1))/15.0
        # Sunrise and sunset time in hours utc
        sunriseArr = self.empty()
        sunsetArr = self.empty()
        sunriseArr = tsouth-t
        sunsetArr = tsouth+t
        return sunriseArr,sunsetArr
                
    def execute(self, timeRange, editArea, varDict):

        ###################   Configuration section   ######################
        #
        # Note: These config items are here so the tool can run without
        # a configuration file. If you wish ot override one of these
        # values, you need to install the config file and define your
        # overrides in that file.  Changes made to this section will only
        # take effect if there is no DiurnalConfig file present
        #
        ####################################################################
        self.configDict = {}

        self.configDict["sourceList"] = ["NAM12", "GFS",
                        "Obs", "MSAS", "Official", "F-table (Climo)"]
                
        # Sources that are based on observations must be listed here
        # since they are processed differently than forecast models 
        self.configDict["obsModelList"] = ["Obs", "MSAS"]

        # These models will be selected by default in the GUI
        self.configDict["primaryModel"] = "NBM"
        self.configDict["secondaryModel"] = "GFS"
        
        # if no other model guidance is found, use lastResort as guidance
        self.configDict["lastResort"] = "GFS"

        # This flag controls whether the tool derives the model min/max from
        # the hourly grids (True) or the min/max grids fetched from the model
        # If set to False, the tool will use the model min/max grids in the
        # database where/when it is lower/higher than the hourly min/max.
        self.configDict["useHourlyMinMax"] = True

        # Previous versions and model substitutes will be color coded
        # based on this list.  If the first sub was made the first color
        # is used, second sub, second color, etc.
        self.configDict["highlightColors"] = ["yellow", # second choice / previous version
                                              "orange", # third choice
                                              "red",    # and so on
                                              "red",
                                              "red",
                                             ]

        self.configDict["DEBUG"] = False  # set this to True to see the intermediate grids

        self.configDict["weList"] = ["T"]  # selected elements by default

        self.configDict["startTime"] = "Today"
        self.configDict["endTime"] = "Day 7"
        
        self.configDict["displayDaysToExclude"] = False

        # the number of days that will appear in the menu to exclude
        self.configDict["daysToExclude"] = 14

        self.configDict["adjustDST"] = 0
        
        # A True value will check if there is an active edit area defined
        # and if so, inform the user and ask if they wish to continue
        self.configDict["checkForEditArea"] = False

        # The number of grid points below which a dialog will appear
        # informing the user that an edit area is selected.  If the
        # number of edit area grid points is above this number, no
        # dialog will appear.
        self.configDict["editAreaPoints"] = 50
        ###################   End Configuration section   ##################

        # Replace any config items with those found in the DiurnalConfig module
        self.replaceConfigValues()

        allowedWEs = ["T", "RH"]  # list of allowable elements. This must never change.

        # make a mask for the specified editArea
        mask = self.encodeEditArea(editArea)
        
        # Run from the GFE, so make the variableList to prompt for info
        if varDict is None:
            # Set fetch time period name and set up GUI
            trNames = self.getTimePeriods()
    
            variableList = []
            trNames.insert(0, "Selected TimeRange")
            variableList.append(("Beginning", self.configDict["startTime"], "radio", trNames))
            variableList.append(("Ending", self.configDict["endTime"], "radio", trNames))
            variableList.append(("Element", self.configDict["weList"], "check", allowedWEs))
    
            variableList.append(("Primary Model for Diurnal:", self.configDict["primaryModel"],
                                 "radio", self.configDict["sourceList"]))
            variableList.append(("Secondary Model for Diurnal:", self.configDict["secondaryModel"],
                                 "radio", self.configDict["sourceList"]))

            # See if we need to make an excludeList
            if self.configDict["displayDaysToExclude"]:
                # build a list of timeRanges for possible exclusion
                oneDay = 24 * 3600
                today00Z = int(time.time() / (oneDay)) * oneDay
                GUIList = []
                self._timeDict = {}
                excludeEnd = today00Z - self.configDict["daysToExclude"] * oneDay
                for t in range(today00Z, excludeEnd, -oneDay):
                    tr = self.makeTimeRange(t, t + oneDay)
                    dateStr = self.makeDateStr(t)
                    GUIList.append(dateStr)
                    self._timeDict[dateStr] = tr

                variableList.append(("Exclude Obs:", [], "check", GUIList))

            # Display the GUI and check for cancel
            varDict = {}
            processVarList = ProcessVariableList.ProcessVariableList(
                "Diurnal", variableList, varDict) 
            status = processVarList.status()
            if status.upper() != "OK": 
                self.cancel()

            # check for a small number of points
            numPoints = np.sum(np.sum(mask))
            if numPoints == 0:
                mask = self.newGrid(True, np.bool)
            elif (
                self.configDict["checkForEditArea"]
                and numPoints < self.configDict["editAreaPoints"]
            ):
                if not self.userWantsToContinue():
                    return
                
        # Create a copy of this mask for use with the Diurnal Smoothing
        fullMaskArea = mask

            # end if varDict is None

        self._excludeTRList = []
        if self.configDict["displayDaysToExclude"]:
            excludedDays = varDict["Exclude Obs:"]
            self._excludeTRList = []
            for exDay in excludedDays:
                self._excludeTRList.append(self._timeDict[exDay])

        ### end if tool was run from GFE  #####

        # Check for info in varDict
        infoList = ["Beginning", "Ending", "Element", "Primary Model for Diurnal:",
                    "Secondary Model for Diurnal:",
                    ]

        for info in infoList:
            if info not in varDict.keys():
                self.statusBarMsg(info + " not found in varDict. Make sure this is " + \
                                  "defined before calling Diurnal.", "U")
                self.cancel()


        # Determine model and timeRange over which the tool will execute
        startName = varDict["Beginning"]
        endName = varDict["Ending"]
        
        primary = varDict["Primary Model for Diurnal:"]
        secondary = varDict["Secondary Model for Diurnal:"]
        diurnalModelList = [primary, secondary, self.configDict["lastResort"]]
        
        climoOnlyAreas = varDict.get("Climo-Only Areas", [])

        self._smoothAreas = varDict.get("Custom Smooth Areas", [])

        # Make sure the parms we're editing are loaded in the GFE
        calcWEs = []  # make a copy
        for weName in varDict["Element"]:
            calcWEs.append(weName)

        # remove any WEs that are loaded
        parmList = self.loadedParms()
        for parmName, level, dbID in parmList:
            if dbID == self.mutableID() and parmName in calcWEs:
                calcWEs.remove(parmName)

        # load any unloaded WEs
        for weName in calcWEs:
            self.loadParm("Fcst", weName, "SFC")
            
        weList = varDict["Element"]  # get the list of weather elements
        for weName in weList:
            executeTR = self.makeExecuteTR(startName, endName, timeRange)
            #print "executeTR is now ",executeTR

            # This tool access and processes the same set of grids multiple
            # times. So, they are represented as global variables, rather
            # than passing them around with each step in the process.
            self._modelMaxGridDict = {}  # model max grids
            self._modelMinGridDict = {}  # model min grids
            self._fcstMinGridDict = {}   # Fcst min grids
            self._fcstMaxGridDict = {}   # Fcst max grids
            self._modelHrGridDict = {}  # hourly model grids
            self._fTrendDict = {}   # hourly fTrends that determine which min
                                    # or max grid will be used
            self._highlightDict = {}

            # Calculating from the f-table is implemented in another tool
            if primary == "F-table (Climo)":
                if weName == "T":
                    self.callProcedure("DiurnalTFromClimo", varDict=varDict, timeRange=timeRange,
                                       editArea=editArea)
                    continue
                else:
                    self.statusBarMsg("Climo option is only valid for Temperature.", "S")
                    return
            
            # Fetch model grids and interpolate to one hour
            # Expand the timeRange 12 hours so we can get good mins and maxes
            start = executeTR.startTime().unixTime() - (12 * 3600)
            end = executeTR.endTime().unixTime() + (12 * 3600)
            hourlyTR = self.makeTimeRange(start, end)
            #print "hourlyTR is ",hourlyTR
            self._modelHrGridDict = self.makeHourlyGrids(diurnalModelList, weName, hourlyTR)

            if self.configDict["DEBUG"]:
                trList = list(self._modelHrGridDict.keys())
                trList.sort(key=lambda a: a.startTime().unixTime())
                ##print len(trList), "hourly grids generated"
                for tr in trList:
                    self.createGrid("Fcst", weName + "HourlyGrid", "SCALAR",
                                    self._modelHrGridDict[tr], tr)


            # Fetch the timeranges for the Fcst min and max
            minTRList, maxTRList = self.getMinMaxTRLists(weName, executeTR)

            # Fetch the model min and max grids
            self._modelMinGridDict, self._modelMaxGridDict = self.getMinMaxGrids(
                diurnalModelList[0], weName, minTRList, maxTRList)

            if self.configDict["DEBUG"]:
                trList = list(self._modelMinGridDict.keys())
                trList.sort(key=lambda a: a.startTime().unixTime())
                for tr in trList:
                    self.createGrid("Fcst", weName + "RawModelMin", "SCALAR",
                                    self._modelMinGridDict[tr], tr,
                                    minAllowedValue=-100.0, maxAllowedValue=150.0)
                trList = list(self._modelMaxGridDict.keys())
                trList.sort(key=lambda a: a.startTime().unixTime())
                for tr in trList:
                    self.createGrid("Fcst",  weName + "RawModelMax", "SCALAR",
                                    self._modelMaxGridDict[tr], tr,
                                    minAllowedValue=-100.0, maxAllowedValue=150.0)

                
            # Fetch the Fcst min and max grids
            self._fcstMinGridDict, self._fcstMaxGridDict = self.getMinMaxGrids(
                "Fcst", weName, minTRList, maxTRList)

            # QC the set of min and max grids we found
            if not self.checkFcstMinMaxGrids(weName, executeTR,
                                             self._fcstMinGridDict,
                                             self._fcstMaxGridDict):
                self.cancel()
            # Calculate the trend for each hour -1 decreasing, 1 increasing
            self.calcFTrendGrids()

            # Finally calculate hourly grids based on the Fcst min/max and model trends
            self.calcDiurnalGrids(weName, executeTR, mask,fullMaskArea)

            # check to see if we should replace grids in the latter time frame
            # with climatology.
            if secondary == "F-table (Climo)" and weName == "T":
                self.overwriteWithClimo(primary, executeTR, editArea)
            elif secondary == "F-table (Climo)" and weName == "RH":
                self.statusBarMsg("The Climo option only works for T grids.", "S")

            # Highlight and report time periods for missing grids
            highlightDict = self._highlightDict
            self.highlightTRs(weName, highlightDict, executeTR)
            self.reportMissingGrids(weName, highlightDict, executeTR)

        #if "Td" in includeParms:
        self.loadIfNot("Td")
        self.makeTd("Fcst",executeTR)
        #:
        #self.loadIfNot("HeatIndex")
        #self.makeHI("Fcst",executeTR)
        #
        #self.loadIfNot("WindChill")
        #self.makeWC("Fcst",executeTR)
        #
        self.loadIfNot("ApparentT")
        self.makeAT("Fcst",executeTR)
        
        self.loadIfNot("MaxApT")
        self.loadIfNot("MinApT")
        self.maxminAT("Fcst",executeTR)
        
        ############# WBGT
        self.__pointDiagnostics = 0
        threatKeys = self.getDiscreteKeys("WBGTRisk")
        lats,lons = self.getLatLonGrids()
        wbgtGrid = self.empty()
        
        if self.__pointDiagnostics:
            # Get coords of a test point
            self.__testCoords = self.getGridCell(26.0,-81.0)
            # Open diagnostic file for writing
            self.__testFile = open("/tmp/WBGT.txt", "w")

        # Delete existing grids
        gfeTR = TimeRange.allTimes().toJavaObj()        
        self.deleteCmd(["WBGT","WBGTRisk","MaxWBGTRisk"],gfeTR)
        self.saveElements(["WBGT","WBGTRisk","MaxWBGTRisk"])
        
#         print("*********************************************************************************")
        # Import albedo and roughness grids
        try:
            albedoGrid = self.getGrids("WBGTdb","Albedo","SFC",TimeRange.allTimes())
            roughnessGrid = self.getGrids("WBGTdb","Roughness","SFC",TimeRange.allTimes())
        except:
            print("Albedo or Roughness Grid not found. Exiting")
            return
        # Get the rap database that goes out at least 36 hours
#         rapDB, rapTime, rapHour, rapGInfo = self._findLongestModelDB("D2D_RAP13_","p","SFC",gfeTR)
#         print(("RAP Database Used for first 36 hours of PMSL Data: %s"%(rapDB)))
        # Get the fullest GFS database (i.e. doesn't have a lot of skipped data)
        gfsDB, gfsTime, gfsHour, gfsGInfo = self._findLongestModelDB("D2D_GFS_","p","SFC",gfeTR)
#         print(("GFS Database Used for remainder of PMSL Data: %s"%(gfsDB)))
        # Make a dictionary of grids that are hourly through the first 3 periods with the RAP and 3/6 hourly further in time using the GFS
#         stitchedDB,gfeTRList = self._stitchDBs(rapDB,rapGInfo,gfsDB,gfsGInfo,"p","SFC")
        # Make a dictionary of grids that are hourly through the whole period by holding pmsl constant until the next grid.
        stitchedDB,gfeTRList,fullTr = self._stretchedDB(gfsDB, gfsGInfo, "p", "SFC")
        
        # Get terrain grid for station pressure calculation
        terrain = self.getTopo()
        terrain_meters = terrain * 0.3048

        # Get forecast grids
        tGrids = self.getGrids("Fcst", "T", "SFC", gfeTRList,mode="TimeWtAverage",noDataError=0)
        tdGrids = self.getGrids("Fcst", "Td", "SFC", gfeTRList,mode="TimeWtAverage",noDataError=0)
        windGrids = self.getGrids("Fcst", "Wind", "SFC", gfeTRList,mode="TimeWtAverage",noDataError=0)
        skyGrids = self.getGrids("Fcst", "Sky", "SFC", gfeTRList,mode="TimeWtAverage",noDataError=0)
        
        # Set up WBGT edit areas and thresholds
        self.regionsDict = {
            "Region1":{"editArea":'ea1',
                       "low":72.3,
                       "elevated":76.1,
                       "moderate":80.1,
                       "high":84.0},
            "Region2":{"editArea":'ea2',
                       "low":75.9,
                       "elevated":78.7,
                       "moderate":83.7,
                       "high":87.6},
            "Region3":{"editArea":'ea3',
                       "low":78.3,
                       "elevated":82.0,
                       "moderate":86.0,
                       "high":90.0},
            }
        thisSite = self.getSiteID()
        if thisSite in ["SJU","HFO","GUM"]:
            self.validEAs = [("Region3","Entire_Domain")]
        else:
            self.eaList = self.editAreaList("WBGT Regions CONUS")
            print(self.eaList)
            self.validEAs = [(region,self.regionsDict[region]["editArea"]) if self.regionsDict[region]["editArea"] in self.eaList else (region, None) for region in list(self.regionsDict.keys())]
        print(("Valid Edit Areas Are: %s"%(self.validEAs)))
        print("*********************************************************************************")
        
        # Loop through each grid and do the calculations
        for eachGridTime in gfeTRList:
            # Be sure none of these are None
            thisP = stitchedDB[eachGridTime]
            thisT = tGrids[eachGridTime]
            thisTd = tdGrids[eachGridTime]
            thisWind = windGrids[eachGridTime]
            thisSky = skyGrids[eachGridTime]
            testArr = [x for x in (thisT,thisTd,thisWind,thisSky,thisP) if x is None]
            # If all grids exist
            if not testArr:
                # Convert pressure grid to station pressure in MSLP
#                 thisPinHg = thisP * 0.0002953
#                 thisSPressMB = (thisPinHg * (((288 - (0.0065 * terrain_meters))/288)**5.2561)) * 33.8639
                thisSPressMB = thisP / 100
                # Get grid info for time reference
                gridStartUTC = eachGridTime.startTime()
                startYear = gridStartUTC.year
                startMonth = gridStartUTC.month
                startDay = gridStartUTC.day
                startHour = gridStartUTC.hour
                startMinute = gridStartUTC.minute
                # Get sunrise and sunset information
                daylenArr = self.empty()
                maxflux = self.empty()
                solarnoon_GMT = self.empty()
                sunriseArrGMT, sunsetArrGMT = self.sunrSuns(startYear, startMonth, startDay, lons, lats)
                daylenArr = sunsetArrGMT - sunriseArrGMT
                # Get various solar fields
                maxflux,fEot,fR0r,fRdecl,fDeclsc1,fDeclsc2, nJulianDate = self.solarFlux(lats, startYear, startMonth, startDay)
                solarnoon_GMT = (sunriseArrGMT + (0.5*daylenArr))
                if self.__pointDiagnostics:
                    self.__testFile.write("************************************************************\n")
                    self.__testFile.write("Grid Time\n")
                    self.__testFile.write("{}\n".format(eachGridTime))
                    self.__testFile.write("startMonth, startDay, startHour, maxflux, solarnoon_GMT\n")
                    self.__testFile.write("{}, {}, {}, {}, {}\n".format(startMonth, startDay, startHour, maxflux[self.__testCoords[1]][self.__testCoords[0]], solarnoon_GMT[self.__testCoords[1]][self.__testCoords[0]]))
                #######################
                ### WBGT PARAMETERS ###
                #######################
                
                # temps convert to degC
                thisTinC = (thisT-32)*5/9
                thisTdinC = (thisTd-32)*5/9
                
                # wind convert from 10m to 2m, and to m/s
                # using logarithmic wind profile from Allen et. al 1998,
                # "Crop evapotranspiration - Guidelines for computing crop water requirements"
                # Eq. 47
                # convert to m/s
                thisWindSpd_ms = thisWind[0] * 0.5144
                #calculate approx wind at 2m height
                this2mWindSpd_ms = (thisWindSpd_ms)*( np.log(2/roughnessGrid) / np.log(10/roughnessGrid) )
                #adjust for calm winds
                calmMask = this2mWindSpd_ms < 1.788
                this2mWindSpd_ms[calmMask] = 1.788
                #convert to m/hr for WBGT calc
                this2mWindSpd_mh = this2mWindSpd_ms*3600.0
                if self.__pointDiagnostics:
                    self.__testFile.write("this2mWindSpd_ms, this2mWindSpd_mh\n")
                    self.__testFile.write("{}, {}\n".format(this2mWindSpd_ms[self.__testCoords[1]][self.__testCoords[0]],this2mWindSpd_mh[self.__testCoords[1]][self.__testCoords[0]])) 
                #adjust maxflux from above for day/night
                maxflux = self.localSolarFluxAdj(sunsetArrGMT, sunriseArrGMT, solarnoon_GMT, startHour, daylenArr, maxflux)
                if self.__pointDiagnostics:
                    self.__testFile.write("maxflux\n")
                    self.__testFile.write("{}\n".format(maxflux[self.__testCoords[1]][self.__testCoords[0]]))
                # Calculate WBGT
                self.wbgtMainCalc(maxflux, thisSky, nJulianDate, startHour, lons, lats, thisTdinC, thisTinC, thisSPressMB, this2mWindSpd_mh, thisWindSpd_ms, thisT, eachGridTime, threatKeys, albedoGrid,fEot,fR0r,fRdecl,fDeclsc1,fDeclsc2)
        self._createMaxGrids(fullTr,gfeTRList,threatKeys)
        #
        #  Save the auto-generated parms - this was added because ForecastBuilder
        #    runs diurnal and does it's own save - but doesn't know that our
        #    version makes some other fields - so doesn't save them.
        #
        self.saveElements(["T","Td","RH","HeatIndex","WindChill","ApparentT","MaxApT","MinApT","WBGT","WBGTRisk","MaxWBGTRisk"])

        self.statusBarMsg("Diurnal is FINISHED.", "R")
        return
        
