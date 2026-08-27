# ----------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# obsdump - part of MatchObsAll routines - version 0.214
#
#    Searches through METAR and mesonet netcdf files, and creates ASCII
#    comma-delimeted files of Temp, Dewpoint, Wind Direction, Wind
#    speed, Wind Gust.  If there are many obs in a single hour, the one
#    nearest the "top of the hour" is saved.
#
# Author: Tim Barker - SOO BOI
#   2021/04/21 - version 0.207  MBC   Added support for WaveHeight.
#   2020/12/18 - version 0.210. JRW.  Updated to Python3 using updatePython
#   2018/06/04 - version 0.206. JRW.  Added DUMPAREA to remove obs outside a defined
#                                     edit area.  If no DUMPAREA is selected it will
#                                     default to the entire domain.
#   2016/11/04 - version 0.205. JRW.  Fixed Low Vis problems
#   2016/08/10 - version 0.204. JRW.  Fixed WindGust/pkwnd computation
#   2016/08/10 - version 0.203. JRW.  Removed need to import ObsdumpConfig.  Will
#                                     now import MatchObsAllConfigCR
#   2016/07/01 - version 0.202. JRW.  Fixed bugs related to cloud heights.
#                                     Fixed formatting from 3 spaces to 4
#   2016/03/10 - version 0.117. JRW.  Modified to run using the DAF instead of uengine.
#   2016/01/16 - version 0.116. JRW.  Fixed Sky cover obs to use the lowest coverage available.
#                                     Fixed CloudBasePrimary to set NonCeiling values when there
#                                     is no OVC, BKN, or VV in the metar. It will set to lowest
#                                     FEW, or SCT layer.
#   2015/12/22 - version 0.115. JRW.  Modified Visibility to output in 100's of miles.
#                                     This fixes a problem when Vis was less then 1 mile.
#   2015/07/15 - version 0.114. JRW.  Fixed Ceiling problem with missing heights.
#   2015/05/05 - version 0.112. JRW.  Modified to use new Uengine Structure.
#                                     Fixed time issue with Maritime Obs.
#   2014/12/05 - version 0.111. JRW.  Added Aviation elements.
#   2012/01/16 - version 0.110. Change mesonet access because data is now
#                in HDF5 files rather than in postgres database.
#   2011/07/20 - version 0.108. AWIPS-2 port.
#                  Changed how configuration is read in.
#   2008/06/10 - version 0.105. Changed into its own procedure rather
#                than a tool - so that there is no dependence on WG1
#                parm and no problem with locked grids if the procedure
#                unexpectedly fails.
#   2008/05/27 - version 0.103 fix issue with moving ships not having
#                right lat/lon.
#   2008/04/11 - version 0.100.  Added type checks for getGridCell
#                because OB8.3 is picky about type checking.
#   2006/10/10 - version 0.99. Better handle character type of netCDF
#                data.  Fix problems with checking times outside
#                window.
#   2005/02/14 - version 0.98. Update version number for rest of
#                MatchObsAll system
#   2003/05/27 - version 0.95. Fix typo in getRoughLimits.  Fix area
#                of obs dumped for non-square domains.
#   2003/04/08 - version 0.93. Non-polygon edit areas will no longer
#                cause obsdump to fail. Also will flush stdout when
#                writing log messages
#   2003/02/26 - version 0.92. Version update only for MatchObsAll
#   2003/02/07 - add maritime data. Better log messages when files
#                not found.  Better log messages for start/stop and
#                writing datafiles.
#   2003/02/06 - remove any commas from site ID or Name, since that
#                will mess up other programs that split lines based
#                on commas and expect values in certain columns.
#   2003/01/20 - First version - used to support MatchObsAll tool
#--------------------------------------------------------------------
#
#  Normally don't want Obsdump to show up in menus - it gets run
#  via cron every hour.  But you could set MenuItems to "Populate"
#  (or something similar) for testing if you really need it.
#
MenuItems = None
import numpy as np

from collections import defaultdict
import datetime
import os
import sys
import calendar
import re
import time

import SmartScript
from ufpy.dataaccess import DataAccessLayer as DAL
from shapely.geometry import MultiPoint
import AbsTime
import TimeRange
from java.util import Date
import dateutil.tz
GMT_ZONE = dateutil.tz.gettz('GMT')

#-----------------------------------------------------------------
#
# Get site configuration
#
import MatchObsAllConfigCR as MOAC
#SITEID=self.getSiteID()
#
#  Set defauls if not set in the site's configuration file
#
if "DEBUG" not in MOAC.Config:
    MOAC.Config["DEBUG"]=0
if "PREVHOURS" not in MOAC.Config:
    MOAC.Config["PREVHOURS"]=8
if "OBSDIR" not in MOAC.Config:
    MOAC.Config["OBSDIR"]="/localapps/runtime/MatchObsAll/data/XXX"
if "sortareas" not in MOAC.Config:
    MOAC.Config["sortareas"]=[]
if "DUMPAREA" not in MOAC.Config:
    MOAC.Config["DUMPAREA"]=None
#-----------------------------------------------------------------

#------------------------------------------------------------------
# logtime - write a string with a date/time stamp
#
def logtime(message, priority=0):
    if priority > MOAC.Config["DEBUG"]:
        return
    gt = datetime.datetime.now(GMT_ZONE)
    ts = gt.strftime("%Y/%m/%d %H:%M:%S")
    print("%s| %s" % (ts, message))
    sys.stdout.flush()
    return

class StationData(object):
    def __init__(self, sid, name, lat, lon, elev):
        self._sid = sid
        self._name = name
        self._lat = lat
        self._lon = lon
        if elev > 20000:
            self._elev = 20000
        elif elev < -1000:
            logtime("Woops!  elev = {}".format(elev))
            self._elev = -1000
        else:
            self._elev = int(elev)

    def get_sid(self):
        return self._sid

    def get_name(self):
        return self._name

    def get_lat(self):
        return self._lat

    def get_lon(self):
        return self._lon

    def get_elev(self):
        return self._elev

class Observation(object):
    def __init__(self, stnid, obtime):
        self._stnid = stnid
        self._obtime = obtime
        self._offset = float("inf")
        self._temp = None
        self._dpt = None
        self._wdir = None
        self._wspd = None
        self._pkwd = None
        self._gust = None
        self._vsby = None
        self._wvehght = None
        self._sky_cover = []
        self._sky_layer_base = []

    def get_stnid(self):
        return self._stnid
    def set_stnid(self, t):
        try:
            self._stnid = t
        except:
            self._stnid = None

    def get_obtime(self):
        return self._obtime
    def set_obtime(self, t):
        try:
            self._obtime = t
        except:
            self._obtime = None

    def get_offset(self):
        return self._offset
    def set_offset(self, t):
        try:
            self._offset = int(t)
        except:
            self._offset = float("inf")

    def get_temp(self):
        return self._format_return_val(self._temp)
    def set_temp(self, t, unit='C'):
        if t is not None and t > -9000:
            try:
                self._temp = self._update_t_td('T', t, unit)
            except Exception as e:
                logtime(e.args)
                self._temp = None

    def get_dpt(self):
        return self._format_return_val(self._dpt)
    def set_dpt(self, t, unit='C'):
        if t is not None and t > -9000:
            try:
                self._dpt = self._update_t_td('Td', t, unit)
            except Exception as e:
                logtime(e.args)
                self._dpt = None

    def get_wdir(self):
        return self._format_return_val(self._wdir)
    def set_wdir(self, t):
        if t is not None and t > -9000:
            try:
                if t >= 0 and t <= 360:
                    self._wdir = int(t)
                else:
                    raise ValueError("Wind direction out of range: {}".format(t))
            except Exception as e:
                logtime(e.args)
                self._wdir = None

    def get_wspd(self):
        return self._format_return_val(self._wspd)
    def set_wspd(self, t, unit='kt'):
        #logtime("Input wind speed: {}".format(t))
        if t is not None and t > -9000:
            try:
                if unit == 'm/s':
                    self._wspd = int(round(self._MStoKT(t)))
                elif unit == 'kt':
                    self._wspd = int(t)
                else:
                    raise ValueError("Unknown unit: {}".format(t))
                if self._wspd < 0 or self._wspd > 150:
                    raise ValueError("Wind speed out of range: {}".format(self._wspd))
            except Exception as e:
                logtime(e.args)
                self._wspd = None

    def get_pkwd(self):
        return self._format_return_val(self._pkwd)
    def set_pkwd(self, t, unit='kt'):
        #logtime("Input peak wind speed: {}".format(t))
        if t is not None and t > -9000:
            try:
                if unit == 'm/s':
                    self._pkwd = int(round(self._MStoKT(t)))
                elif unit == 'kt':
                    self._pkwd = int(t)
                else:
                    raise ValueError("Unknown unit: {}".format(t))
                if self._pkwd < 0 or self._pkwd > 150:
                    raise ValueError("Peak Wind out of range: {}".format(self._pkwd))
            except Exception as e:
                logtime(e.args)
                self._pkwd = None

    def get_gust(self):
        if self.get_pkwd() != self._format_return_val(None):
           if self._format_return_val(self._gust) != self._format_return_val(None):
              if self._pkwd > self._gust:
                 return self.get_pkwd()

        if self._format_return_val(self._gust) != self._format_return_val(None):
           if self._gust > 0:
              return self._format_return_val(self._gust)
        return self._format_return_val(None)
    def set_gust(self, t, unit='kt'):
        #logtime("Input wind gust: {}".format(t))
        if t is not None and t > -9000:
            try:
                if unit == 'm/s':
                    self._gust = int(round(self._MStoKT(t)))
                elif unit == 'kt':
                    self._gust = int(t)
                else:
                    raise ValueError("Unknown unit: {}".format(t))
                if self._gust < 0 or self._gust > 150:
                    raise ValueError("Wind gust out of range: {}".format(self._gust))
            except Exception as e:
                logtime(e.args)
                self._gust = None

    def get_vsby(self):
        return self._format_return_val(self._vsby, "vis")
    def set_vsby(self, t):
        if t is not None and t > -9000:
            try:
                if int(t) < 0 or int(t) > 10:
                    raise ValueError("Visibility out of range: {}".format(t))
                self._vsby = t
            except Exception as e:
                logtime(e.args)
                self._vsby = None


    def get_wvehght(self):
        return self._format_return_val(self._wvehght)
    def set_wvehght(self, t):
        if t is not None and t > -9000:
           self._wvehght = int(self.MtoFT(t))
           #print t
        #return
    def reset_sky(self):
        self._sky_cover = []
        self._sky_layer_base = []

    def get_sky_cover(self):
        return self._sky_cover
    def add_sky_cover(self, s):
        if s and s in ["VV", "OVC", "BKN", "SCT", "FEW", "CLR"]:
            self._sky_cover.append(s)
        elif s:
            print("sky cover: %s is not setup"%s)
        #if s and s > -9000:
        #    self._sky_cover.append(s)

    def get_sky_layer_base(self):
        return self._sky_layer_base
    def add_sky_layer_base(self, s):
        if s and s > -9000:
            self._sky_layer_base.append(int(s))

    def get_skyamt(self):
        coverList = self.get_sky_cover()
        if coverList:
            if "VV" in coverList:
                return self._format_return_val(100)
            if "OVC" in coverList:
                return self._format_return_val(95)
            if "BKN" in coverList:
                return self._format_return_val(75)
            if "SCT" in coverList:
                return self._format_return_val(50)
            if "FEW" in coverList:
                return self._format_return_val(25)
            if "CLR" in coverList:
                return self._format_return_val(5)
            return self._format_return_val(0)
        else:
            return self._format_return_val(None)

    def get_cig(self):
        coverList = self.get_sky_cover()
        hghts = self.get_sky_layer_base()
        try:
           z = sorted(zip(hghts, coverList))
        except:
           return self._format_return_val(None)
        try:
           hghts_sorted, cover_sorted = list(zip(*z))
        except:
           return self._format_return_val(None)
        index = -1
        if cover_sorted:
#            if self.get_stnid() == "KHFJ":
#               logtime("COVER LIST = {}".format(repr(cover_sorted)))
#               logtime("HEIGHTS = {}".format(repr(hghts_sorted)))
            if "BKN" in cover_sorted or "OVC" in cover_sorted or "VV" in cover_sorted:
                try:
                    return self._format_return_val(hghts_sorted[cover_sorted.index("BKN")] / 100)
                except:
                    pass
                try:
                    return self._format_return_val(hghts_sorted[cover_sorted.index("OVC")] / 100)
                except:
                    pass
                try:
                    return self._format_return_val(hghts_sorted[cover_sorted.index("VV")] / 100)
                except:
                    pass
            elif "SCT" in cover_sorted or "FEW" in cover_sorted:
                try:
                    return self._format_return_val(hghts_sorted[cover_sorted.index("SCT")] / 100)
                except:
                    pass
                try:
                    return self._format_return_val(hghts_sorted[cover_sorted.index("FEW")] / 100)
                except:
                    pass
            else:
                return self._format_return_val(250)
        else:
            return self._format_return_val(None)

    def _format_return_val(self, val, parm=None):
        if parm == "vis":
            retstr = "    "
            if (val is not None):
                val *= 100
                rval = int(round(val))
                retstr = "{:4d}".format(rval)
        else:
            retstr = "   "
            if (val is not None):
                rval = int(round(val))
                retstr = "{:3d}".format(rval)
        return retstr

    def _update_t_td(self, parm, t, unit):
        #logtime("Input value for {} = {}".format(parm, t))
        if unit == 'C':
            rtn = int(round(self._CtoF(t)))
        elif unit == 'K':
            rtn = int(round(self._KtoF(t)))
        elif unit == 'F':
            rtn = int(round(t))
        else:
            raise ValueError("Unknown unit: {}".format(unit))
        if rtn < -80 or rtn > 130:
            raise ValueError("{} out of range: {}".format(parm, rtn))

        return rtn

    def _CtoF(self, t):
        return (t * (9.0 / 5.0)) + 32.0

    def _KtoF(self, t):
        return (t - 273.15) * (9.0 / 5.0) + 32.0

    def _KtoC(self, t):
        return t - 273.15

    def _MStoKT(self, sms):
        return sms * 1.943

    def MtoFT(self,sms):
        return sms * 3.281


#
#  The Main Procedure
#
class Procedure (SmartScript.SmartScript):
    def __init__(self, dbss):
        SmartScript.SmartScript.__init__(self, dbss)
    def execute(self):
        "Generate files with obs data nearest top of hour"
        t0 = time.time()
        #self.read_config()
        logtime("Starting obsdump")
        #
        #  setup sorting edit areas
        #
        self.sortmasks = []
        self.sortnames = []
        for eaname in MOAC.Config["sortareas"]:
            logtime("Processing sort area: {}".format(eaname))
            mask = self.encodeEditArea(eaname)
            if (mask is not None):
                self.sortmasks.append(mask)
                self.sortnames.append(eaname)
            else:
                msg = "Edit area %s is not a valid polygon in IFP Server" % eaname
                logtime(msg)

        if MOAC.Config["DUMPAREA"] is None:
            self.dumpMask=self.newGrid(1)
        else:
            mask=self.encodeEditArea(MOAC.Config["DUMPAREA"])
            if (mask is not None):
                self.dumpMask = mask.copy()
            else:
                self.dumpMask=self.newGrid(1)
        #
        #  get rough lat/lon box of grid
        #
        ny, nx = self.getGridShape()
        #e = self._empty
        #self.xmax = e.shape[1] - 1
        #self.ymax = e.shape[0] - 1
        (self.LATMIN, self.LATMAX, self.LONMIN, self.LONMAX) = self.getRoughLimits(nx-1, ny-1)
        #projData = self.getGridLoc().getProjection()
        #geom = projData.getGeometry()
        #self.LATMIN = projData.getLatLonLL().y
        #self.LONMIN = projData.getLatLonLL().x
        #self.LATMAX = projData.getLatLonUR().y
        #self.LONMAX = projData.getLatLonUR().x
        #
        # Get station naem dictionary
        #
        self.stationNameDict = {}
        with open(stnfile) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                data = line.split("|")
                if len(data) < 3:
                    continue
                sid = data[0].strip()
                name = data[2].replace(",", " ")
                if sid not in self.stationNameDict:
                    self.stationNameDict[sid] = name
        #self.stationNameDict = self.get_station_names()
        self.Obs = {}
        self.StnData = {}
        logtime("Latitude Limits : %8.3f -- %8.3f" % (self.LATMIN, self.LATMAX), 5)
        logtime("Longitude Limits: %8.3f -- %8.3f" % (self.LONMIN, self.LONMAX), 5)
        self.DTNOW = datetime.datetime.now(GMT_ZONE)
        self.DTNOW = self.DTNOW.replace(minute=0, second=0, microsecond=0)
        self.DELTA_1HOUR = datetime.timedelta(hours=1)
        self.DELTA_30MIN = datetime.timedelta(minutes=30)
        self.DT_LIST = []
        for offset in range(MOAC.Config["PREVHOURS"], -1, -1):
            #logtime("OFFSET: {}".format(offset))
            dt = self.DTNOW - (self.DELTA_1HOUR * offset)
            #logtime("DT = {}".format(repr(dt)))
            self.Obs[dt] = defaultdict(dict)
            self.DT_LIST.append(dt)
        self.get_metars()
        logtime("Execution time after metars= {}s".format(time.time() - t0))
        self.get_mesonet()
        logtime("Execution time after mesonet= {}s".format(time.time() - t0))
        self.get_maritime()
        logtime("Execution time after maritime= {}s".format(time.time() - t0))
        self.write_output()
        logtime("Execution time = {}s".format(time.time() - t0))

    def get_station_names(self):
        rtn = {}
        req = DAL.newDataRequest("common_obs_spatial")
        req.setParameters("stationid", "name")
        envelope = MultiPoint(((self.LONMAX, self.LATMIN), (self.LONMIN, self.LATMAX)))
        req.setEnvelope(envelope)
        geoms = DAL.getGeometryData(req)
        for geom in geoms:
            sid = geom.getString("stationid")
            if sid not in rtn:
                rtn[sid] = geom.getString("name")
        return rtn

    def get_metars(self):
        req = DAL.newDataRequest("obs")
        req.setParameters("elevation", "tempFromTenths", "temperature",
                          "dpFromTenths", "dewpoint", "windSpeed",
                          "windDir", "windGust", "pkwndSpeed", "visibility",
                          "skyCover", "skyLayerBase")
        envelope = MultiPoint(((self.LONMIN, self.LATMIN), (self.LONMAX, self.LATMAX)))
        req.setEnvelope(envelope)
        for dt in self.DT_LIST:
            logtime("Working on METARs for datetime = {}".format(repr(dt)))
            beginDt = dt - self.DELTA_30MIN
            endDt = dt + self.DELTA_30MIN
            start = AbsTime.AbsTime(self.get_epoch(beginDt))
            end = AbsTime.AbsTime(self.get_epoch(endDt))
            tr = TimeRange.TimeRange(start, end)
            observations = DAL.getGeometryData(req, tr)
            #logtime("Working on TR: {}".format(repr(tr)))
            if not observations:
                logtime("MAJOR PROBLEM.  WAS NOT ABLE TO FETCH ANY METAR OBS!")
                continue
            logtime("{} METARs retrieved".format(len(observations)), 10)
            for ob in observations:
                sid = ob.getLocationName()
                if sid not in self.StnData:
                    if "elevation" in ob.getParameters():
                        if sid in self.stationNameDict:
                            stnName = self.stationNameDict[sid]
                        else:
                            stnName = sid
                        self.StnData[sid] = StationData(sid, stnName,
                                                        ob.getGeometry().y,
                                                        ob.getGeometry().x,
                                                        ob.getNumber("elevation")*3.2808)
                dtThisOb = AbsTime.AbsTime(ob.getDataTime().getRefTime())
                if self.Obs[dt][sid]:
                    # A record already exists for this time/station
                    logtime("An observation already exists for SID={} and DT={}".format(sid, repr(dt)), 100)
                    thisOffset = abs(int((dt - dtThisOb).total_seconds()))
                    logtime("thisOffset = {}, previous offset = {}".format(thisOffset, self.Obs[dt][sid].get_offset()), 100)
                    if thisOffset <= self.Obs[dt][sid].get_offset():
                        # This ob is the same or closer to the top of the hour
                        if thisOffset < self.Obs[dt][sid].get_offset():
                            # This is a new ob, closer to the top of the hour
                            # Reset the sky cover
                            self.Obs[dt][sid].reset_sky()
                        self.Obs[dt][sid].set_obtime(dtThisOb)
                        self.Obs[dt][sid].set_offset(abs(dt - dtThisOb).total_seconds())
                        self.set_ob_data("metar", ob, dt, sid)
                else:
                    # This is a new Observation
                    self.Obs[dt][sid] = Observation(sid, dtThisOb)
                    self.Obs[dt][sid].set_offset(abs(dt - dtThisOb).total_seconds())
                    self.set_ob_data("metar", ob, dt, sid)

    def get_mesonet(self):
        req = DAL.newDataRequest("ldadmesonet")
        req.setParameters("elevation", "temperature", "dewpoint", "windDir",
                          "windSpeed", "windGust")
        envelope = MultiPoint(((self.LONMIN, self.LATMIN), (self.LONMAX, self.LATMAX)))
        req.setEnvelope(envelope)
        for dt in self.DT_LIST:
            beginDt = dt - self.DELTA_30MIN
            endDt = dt + self.DELTA_30MIN
            start = AbsTime.AbsTime(self.get_epoch(beginDt))
            end = AbsTime.AbsTime(self.get_epoch(endDt))
            tr = TimeRange.TimeRange(start, end)
            observations = DAL.getGeometryData(req, tr)
            if not observations:
                logtime("MAJOR PROBLEM.  WAS NOT ABLE TO FETCH ANY MESONET OBS!")
                continue
            #dtStr = dt.strftime('%Y%m%d%H')
            logtime("{} MESONETs retrieved".format(len(observations)))
            for ob in observations:
                sid = ob.getLocationName()
                if sid not in self.StnData:
                    if sid in self.stationNameDict:
                        stnName = self.stationNameDict[sid]
                    else:
                        stnName = sid
                    self.StnData[sid] = StationData(sid, stnName,
                                                    ob.getGeometry().y,
                                                    ob.getGeometry().x,
                                                    ob.getNumber("elevation")*3.2808)
                dtThisOb = AbsTime.AbsTime(ob.getDataTime().getRefTime())
                if self.Obs[dt][sid]:
                    # A record already exists for this time/station
                    thisOffset = abs((dt - dtThisOb).total_seconds())
                    if thisOffset <= self.Obs[dt][sid].get_offset():
                        # This ob is closer to the top of the hour
                        self.Obs[dt][sid].set_obtime(dtThisOb)
                        self.Obs[dt][sid].set_offset(abs(dt - dtThisOb).total_seconds())
                        self.set_ob_data("mesonet", ob, dt, sid)
                else:
                    # This is a new Observation
                    self.Obs[dt][sid] = Observation(sid, dtThisOb)
                    self.Obs[dt][sid].set_offset(abs(dt - dtThisOb).total_seconds())
                    self.set_ob_data("mesonet", ob, dt, sid)

    def get_maritime(self):
        req = DAL.newDataRequest("sfcobs")
        req.setParameters("elevation", "temperature", "dewpoint", "windDir",
                          "windSpeed", "windGust","waveHeight")
        envelope = MultiPoint(((self.LONMIN, self.LATMIN), (self.LONMAX, self.LATMAX)))
        req.setEnvelope(envelope)
        for dt in self.DT_LIST:
            beginDt = dt - self.DELTA_30MIN
            endDt = dt + self.DELTA_30MIN
            start = AbsTime.AbsTime(self.get_epoch(beginDt))
            end = AbsTime.AbsTime(self.get_epoch(endDt))
            tr = TimeRange.TimeRange(start, end)
            observations = DAL.getGeometryData(req, tr)
            if not observations:
                logtime("MAJOR PROBLEM.  WAS NOT ABLE TO FETCH ANY MARITIME OBS!")
                continue
            #dtStr = dt.strftime('%Y%m%d%H')
            logtime("{} sfcobs (maritime) retrieved".format(len(observations)))
            for ob in observations:
                sid = ob.getLocationName()
                if sid not in self.StnData:
                    if sid in self.stationNameDict:
                        stnName = self.stationNameDict[sid]
                    else:
                        stnName = sid
                    self.StnData[sid] = StationData(sid, stnName,
                                                    ob.getGeometry().y,
                                                    ob.getGeometry().x,
                                                    ob.getNumber("elevation")*3.2808)
                dtThisOb = AbsTime.AbsTime(ob.getDataTime().getRefTime())
                if self.Obs[dt][sid]:
                    # A record already exists for this time/station
                    thisOffset = abs((dt - dtThisOb).total_seconds())
                    if thisOffset <= self.Obs[dt][sid].get_offset():
                        # This ob is closer to the top of the hour
                        self.Obs[dt][sid].set_obtime(dtThisOb)
                        self.Obs[dt][sid].set_offset(abs(dt - dtThisOb).total_seconds())
                        self.set_ob_data("maritime", ob, dt, sid)
                else:
                    # This is a new Observation
                    self.Obs[dt][sid] = Observation(sid, dtThisOb)
                    self.Obs[dt][sid].set_offset(abs(dt - dtThisOb).total_seconds())
                    self.set_ob_data("maritime", ob, dt, sid)

    #-------------------------------------------------------------------
    # set_ob_data - Set the data parameters in an Observation record
    #-------------------------------------------------------------------
    def set_ob_data(self, obType, ob, dt, sid):
        if 'skyCover' in ob.getParameters() or 'skyLayerBase' in ob.getParameters():
            # Handle sky cover
            self.Obs[dt][sid].add_sky_cover(ob.getString('skyCover'))
            if ob.getString('skyCover') == "CLR":
                self.Obs[dt][sid].add_sky_layer_base(25000)
            else:
                self.Obs[dt][sid].add_sky_layer_base(ob.getNumber('skyLayerBase'))
        else:
            # Handle everything else
            if obType == "mesonet" or obType == "maritime":
                tempUnit = 'K'
                windUnit = "m/s"
            else:
                tempUnit = 'C'
                windUnit = "kt"
            if 'tempFromTenths' in ob.getParameters() and ob.getNumber('tempFromTenths') > -900:
                # tempFromTenths not missing
                self.Obs[dt][sid].set_temp(ob.getNumber('tempFromTenths'), unit=tempUnit)
            else:
                # tempFromTenths is missing, use temperature
                self.Obs[dt][sid].set_temp(ob.getNumber('temperature'), unit=tempUnit)
            if 'dpFromTenths' in ob.getParameters() and ob.getNumber('dpFromTenths') > -900:
                # dpFromTenths not missing
                self.Obs[dt][sid].set_dpt(ob.getNumber('dpFromTenths'), unit=tempUnit)
            else:
                # dpFromTenths is missing, use dewpoint
                self.Obs[dt][sid].set_dpt(ob.getNumber('dewpoint'), unit=tempUnit)
            self.Obs[dt][sid].set_wvehght(ob.getNumber('waveHeight'))
            self.Obs[dt][sid].set_wdir(ob.getNumber('windDir'))
            self.Obs[dt][sid].set_wspd(ob.getNumber('windSpeed'), unit=windUnit)
            self.Obs[dt][sid].set_gust(ob.getNumber('windGust'), unit=windUnit)
            if 'pkwndSpeed' in ob.getParameters():
                self.Obs[dt][sid].set_pkwd(ob.getNumber('pkwndSpeed'), unit=windUnit)
            if 'visibility' in ob.getParameters():
                self.Obs[dt][sid].set_vsby(ob.getNumber('visibility'))

    #-------------------------------------------------------------------
    # sort_obs - Sort the observations by the sort areas
    #
    def sort_obs(self, dt):
        order = []
        outlines = []
        for id in self.Obs[dt]:

            (x, y) = self.getGridCell(self.StnData[id].get_lat(),
                                      self.StnData[id].get_lon())
            if x is None or y is None:
                # Ob is outside the GFE grid
                continue
            #
            if self.dumpMask[y, x]<0.5:
                continue
            #
            sortval = 9999
            #logtime("{} location is {}, {}".format(id, x, y))
            for i in range(len(self.sortmasks)):
                #logtime("MASK TYPE: {}".format(type(self.sortmasks[i][y, x])))
                #logtime("MASK SHAPE: {}".format(self.sortmasks[i][y, x].shape))
                if self.sortmasks[i][y, x] > 0:
                    sortval = i
                    break
            sorttup = (sortval, 20000 - self.StnData[id].get_elev(), id)
            order.append(sorttup)
        order.sort()

        return order

    #-------------------------------------------------------------------
    # write_output - Write the observations out to the text file for
    #                MOA to read.
    #
    def write_output(self):
        for dt in self.DT_LIST:
            outlines = []
            # Added sitename to the end of obsdump list
            siteName=''
            sortedStations = self.sort_obs(dt)
            if sortedStations is not None:
               lastarea = None
               for key in sortedStations:
                   (areanum, elev, idkey) = key

                   if areanum != lastarea:
                       outlines.append("#\n")
                       if (areanum == 9999):
                           outlines.append("#   Misc stations\n")
                           siteName='MSC'
                       else:
                           outlines.append("#  %s\n" % self.sortnames[areanum])
                           siteName=self.sortnames[areanum]
                       outlines.append("#\n")
                   try:
                       ############################################################
                       # Append siteName at the end of each line to output CWA
                       # This will allow sorting by CWA obsqc_py
                       # Added , {} and siteName at the end.
                       ############################################################
                       outlines.append("{:>5},{:40},{:7.4f},{:9.4f},{:5d},{},{},{},{},{},{},{},{},{},{}\n".format(idkey,
                                                                                                        self.StnData[idkey].get_name(),
                                                                                                        self.StnData[idkey].get_lat(),
                                                                                                        self.StnData[idkey].get_lon(),
                                                                                                        self.StnData[idkey].get_elev(),
                                                                                                        self.Obs[dt][idkey].get_temp(),
                                                                                                        self.Obs[dt][idkey].get_dpt(),
                                                                                                        self.Obs[dt][idkey].get_wdir(),
                                                                                                        self.Obs[dt][idkey].get_wspd(),
                                                                                                        self.Obs[dt][idkey].get_gust(),
                                                                                                        self.Obs[dt][idkey].get_skyamt(),
                                                                                                        self.Obs[dt][idkey].get_cig(),
                                                                                                        self.Obs[dt][idkey].get_vsby(),
                                                                                                        self.Obs[dt][idkey].get_wvehght(),
                                                                                                        siteName))
                   except:
                       pass
                       #logtime("vsby is %s"%self.Obs[dt][idkey].get_vsby())
                   lastarea = areanum
               #outfile = "/tmp/test_obsdump_{}.dat".format(dt.strftime('%Y%m%d%H'))
               COMMANDS=["umask 002",
                         "mkdir -p %s"%(MOAC.Config["OBSDIR"] )
                        ]
               for command in COMMANDS:
                   os.system(command)
               outfile = os.path.join(MOAC.Config["OBSDIR"], "{}00.dat".format(dt.strftime('%Y%m%d%H')))
               with open(outfile, 'w') as f:
                   f.writelines(outlines)

    #-------------------------------------------------------------------
    # getRoughLimits - based on the corners of the grid - get a rough
    #                  idea of the lat/lon limits for points we want
    #                  to keep
    #
    def getRoughLimits(self, xmax, ymax):
       (lat1, lon1) = self.getLatLon(0, 0)
       (lat2, lon2) = self.getLatLon(0, ymax)
       (lat3, lon3) = self.getLatLon(xmax, 0)
       (lat4, lon4) = self.getLatLon(xmax, ymax)
       LATMIN = min((lat1, lat2, lat3, lat4))
       LATMAX = max((lat1, lat2, lat3, lat4))
       LONMIN = min((lon1, lon2, lon3, lon4))
       LONMAX = max((lon1, lon2, lon3, lon4))
       return (float(LATMIN), float(LATMAX), float(LONMIN), float(LONMAX))

    def get_epoch(self, dt):
        dt0 = datetime.datetime(1970, 1, 1, tzinfo=GMT_ZONE)
        return (dt - dt0).total_seconds()
