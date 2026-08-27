HideTool = 0
ToolType = "numeric"
WeatherElementEdited = "Wx"
ScreenList = ["Wx"]

import numpy as np
import re
import OKDialog
import sys
import tkinter as tk

import TimeRange
import copy

import AppDialog
import time


################################################################################
# WxTool:
#   Author: Mark Loeffelbein WFO Missoula, MT
#       Altered previous WxTool to keep variables and to limit t-storm wording.
#       Unable to find the author of the original tool.
# Add variables to KeepList to keep them in the Wx grid.
# Make sure to add a slash in between the wx identifiers.
#   example:
#       KeepList = "F/ZR/K/ZR"
#   to add R
#       KeepList = "F/ZR/K/ZR/R"
# This will only work for wx type and not intensity or coverage. Adding R to
# the list will keep R, -R, +R, but adding -R will remove R from the wx grid
# since this script doesn't understand intensities.
#
#  Modified by Ron Miller (WFO Spokane) for the following:
#     Removed Precip Type from GUI since it's determined by SnowLevel
#     Added Precip Intensity from QPF/SnowAmt option
#     Set Rain/Snow mix to above SnowLevel (per Western Region SnowLevel team)
#     Added Dialog box for errors (e.g. Specifying Thunder and Sprinkles)
#
#  Modified by Ron Miller in Aug 2011
#     Put all of the widgets in a self-constructed GUI
#     Added sliders with a skew-t picture to depict the mixed rain/snow layer
#
#  Requires OKDialog (utility) installed on the server.
#
#  Requires at least one snowlevel grid within WX grid time range.
#
#
#  2014/08/05 - Modified so that a specified T coverage was selected - it would
#               apply to all grids, not just the first (previously on subsequent
#               grids it would switch back to 'match pop').  Tim Barker
#  2014/08/06 = Modified to add ability to match T to PoP, but only up to a
#               certain limit of coverage/prob. Barker/BOI
#  2015/08/20 - Removed all where statements.  Barker/BOI
#  2015/09/02 - Added colorization to  thunder/freezing/sleet tabs when they
#               are non-None, so that even if the panel is minimized the
#               forecaster gets a hint that the tool will (potentially) add
#               it (red) or keep it (yellow). Barker/BOI
#  2022/08/01 - Updated for py3 and HCI.  KA/BOI
#  2024/05/20 - Added bidirectional PotThunder support & modernizations.  KA/BOI
#  2026/04/20 - Added more PotThunder fixes. Updated latest python to shorten code. KA/BOI
#
################################################################################
#
# configuration section

# the adjustable rain/snow mix layer must have limits on the thickness.

maxThicknessRainSnowLayer = 2000

# the default top and bottom of the mixed layer, with respect to the snow
# level.  Positive is above snow level, negative is below.

defaultTopRainSnow = +400
defaultBotRainSnow = -400

# default Keep Wx list

KeepList = "WP,L,F,ZF,IF,IC,H,BS,BN,K,BD,FR,ZY,VA"

#  Thresholds for rain and snow intensity.  These are based on 6-hr QPF and
#  SnowAmt grids.  So, for example, if the grid point has 0.25" of rain in
#  6 hours, it will be assigned "moderate" rainfall.  Or a grid that has
#  5" of SnowAmt in 6 hours will be assigned "heavy" snowfall.
moderateRthresh = 0.25
heavyRthresh = 0.50
moderateSthresh = 2.5
heavySthresh = 5.0

#  Elevation to specify the difference between valleys and mountains for precip
#  intensity.  For grid points above this elevation, the above intensity thresholds
#  are doubled (i.e. need greater amounts of rain or snow in the mountains to
#  be labeled as "moderate" or "heavy".
valleyMtnThreshold = 15000

# Thunder coverage to PotThunder value mapping
THUNDER_TO_POTTHUNDER = {
    "SChc": 15,
    "Chc": 25,
    "Lkly": 55,
    "Def": 75,
    "Iso": 15,
    "Sct": 25,
    "Num": 55,
    "Wide": 75
}

#
#  TOOLNAME is the title at the top of the dialog box
#
TOOLNAME="WxFromPoP"
#
#
#  End Configuration Section
################################################################################
#
#  Global regular expression for window geometry used in multiple classes
#
REGEO=re.compile(r"(\d+)x(\d+)([-+])(\d+)([-+])(\d+)")
#
################################################################################

# Set up Class
import SmartScript

class Tool (SmartScript.SmartScript):
    def __init__(self, dbss):
        self._dbss=dbss
        SmartScript.SmartScript.__init__(self, dbss)
        #
        self.GuiSettings={}
        self.GuiDefaults()

    def preProcessTool(self, ToolTimeRange, varDict):
        # Compile a regular expression to look for weather type and intensity
        # combinations that can be changed to moderate (m) or heavy (+)
        # intensities safely.  This will be used later in the tool.
        self._weatherInten = re.compile(':([RW?|SW?|IP|ZR]):-:')
        self.SelectTR=ToolTimeRange

        self.dialog=WxTool(title=TOOLNAME, callbackMethod=self.startWxTool)
        self.dialog.wm_attributes("-topmost", 1)
        self.dialog.wait_visibility()
        self.readStatus()
        self.setupGui()
        self.dialog.mainloop()

        #  When the GUI closes, call the SmartScript cancel routine within
        #  this pre-processTool routine - so that if the user clicked on
        #  the cancel button without ever modifying anything - then no
        #  grid will be marked as edited
        self.cancel()
        return

    #  Cancel within preProcessGrid so that execute does not get run
    #  and mark the grid as edited.  Cannot put this in preProcessTool
    #  because it would then post a cancel dialog box.
    def preProcessGrid(self):
        self.cancel()

    #  Dummy execute routine.  All actions done in startWxTool
    def execute(self):
        "Make Rain/Snow match PoP and SnowLevel"
        return

##########################################################################################################

    def startWxTool(self, buttonType):

        if buttonType=="Cancel":
            self.saveStatus()
            return

        self.mutableName=self.mutableID().modelName()
        selectTR=self.SelectTR

        # get the various settings from the GUI
        self.coverage=self.dialog.coverage.get()
        self.convective=self.dialog.convective.get()
        self.ptype=self.dialog.ptype.get()
        self.intensity=self.dialog.intensity.get()
        self.keepWx=self.dialog.keepWx.get()
        self.keepWxList=self.dialog.keepWxList.get()

        self.thunderAnswer=self.dialog.thunder.get()
        if self.convective=="Stratiform": # if stratiform - no thunder will be allowed
            self.thunderAnswer="None"

        if self.thunderAnswer=="Keep":
            self.keepWxList += ",T"

        self.thunderType=self.dialog.thunderType.get()
        self.forceThunderLevel=self.dialog.tforceLev.get()
        self.matchThunderLevel=self.dialog.tmatchLev.get()

        self.botRainSnow=self.dialog.botRainSnow.get()
        self.topRainSnow=self.dialog.topRainSnow.get()

        self.freezing=self.dialog.freezing.get()
        if self.freezing=="Keep":
            self.keepWxList += ",ZR"
        self.fmmatchLev=self.dialog.fmmatchLev.get()

        self.sleet=self.dialog.sleet.get()
        if self.sleet=="Keep":
            self.keepWxList += ",IP"
        self.immatchLev=self.dialog.immatchLev.get()

        self.popThresh=self.dialog.popThresh.get()

        if self.intensity == "from QPF/SnowAmt":
            self.intensityFromGrids = "Yes"
        else:
            self.intensityFromGrids = "No"

        self.Topo=self.getTopo()

        # get the selectionTimeRange, This may correspond
        # to a single grid, or a range of grids. Then get
        # the list of grids within that time period
        javaselectTR=self._dbss.getParmOp().getSelectionTimeRange()
        selectTR=TimeRange.TimeRange(javaselectTR)
        infoList = self.getGridInfo(self.mutableName, "Wx", "SFC", selectTR)

        # Now loop through each Wx grid that's in this time range.  Call the
        # runWxTool routine for each Wx grid selected
        for info in infoList:
             self.GridTimeRange = info.gridTime()
             self.runWxTool()
        self.saveStatus()
        return


    def runWxTool(self):
        "Assign coverage/probability of Wx to match completed PoP grids.  More flexibility built into the tool to account for precip variables."

        # get the PoP, SnowLevel, and Wx grids for the selected time range
        fPoP = self.getGrids(self.mutableName, "PoP", "SFC", self.GridTimeRange, mode="Max", cache=0, noDataError=0)
        if fPoP is None:
            self.statusBarMsg("No PoP available", "U")
            return

        #  Round PoP to nearest 10th - which is what it would be if it were saved.
        PoP10= np.add( np.multiply(fPoP, 10.0), 0.5)
        fPoP= np.divide(PoP10.astype(int), 10.0)
        #  Round PoP to nearest integer
        fPoP= np.add(fPoP, 0.5)
        self.PoP=fPoP.astype(int)

        self.SnowLevel = self.getGrids(self.mutableName, "SnowLevel", "SFC", self.GridTimeRange, noDataError=0)
        if self.SnowLevel is None:
            self.statusBarMsg("No SnowLevel available", "U")
            return

        self.Wx = self.getGrids(self.mutableName, "Wx", "SFC", self.GridTimeRange, noDataError=0)
        if self.Wx is None:
            self.statusBarMsg("No Wx available", "U")
            return

        editAreaMask=self.getEditAreaMask()
        self.maxT=self.getGrids(self.mutableName, "T", "SFC", self.GridTimeRange, mode="Max", noDataError=0)
        self.minT=self.getGrids(self.mutableName, "T", "SFC", self.GridTimeRange, mode="Min", noDataError=0)

        self.potThunder=self.getGrids(self.mutableName, "PotThunder", "SFC", self.GridTimeRange, noDataError=0)
        if self.potThunder is None:
            self.potThunder = self.newGrid(0.0)

        # Initialize tracking grid for thunder assignments
        self.thunderAssignments = self.newGrid(0)

        # Create spatial masks based on PotThunder thresholds for 'From PotThunder' mode
        # 0: None, 1: SChc/Iso, 2: Chc/Sct, 3: Lkly/Num, 4: Def/Wide
        potMasks = [
            self.potThunder < 14.5,
            (self.potThunder >= 14.5) & (self.potThunder < 24.5),
            (self.potThunder >= 24.5) & (self.potThunder < 54.5),
            (self.potThunder >= 54.5) & (self.potThunder < 74.5),
            self.potThunder >= 74.5
        ]

        # if the user selected keepWx=Y, then call routine to remove all of the
        # Wx that's not in the keepWxList.  This will be stored in oldWxValues/Keys
        if self.keepWx =="Y":
            (oldWxValues, oldWxStrings) = self.keepWxInGrid(self.Wx, editAreaMask)
        else:
            # otherwise, create an empty grid with "<NoWx>"
            oldWxStrings = ["<NoCov>:<NoWx>:<NoInten>:<NoVis>:",]
            oldWxValues = self.newGrid(0)

        # append levels of thunderstorm chances to ltgcov
        self.limitT = 0
        self.ltgcat = 1
        self.thunder="N"

        wxStrings=[]
        wxStrings.append("<NoCov>:<NoWx>:<NoInten>:<NoVis>:")
        wxValues=self.newGrid(0)

        # Assign weather types
        if self.convective == "Stratiform":
            wxval = "R"
            rain =  "R"
            snow =  "S"

        elif self.convective == "Showers":
            wxval = "RW"
            rain =  "RW"
            snow =  "SW"

        if self.ptype == "Rain/Snow" and self.convective == "Stratiform":
            wxval = "R"
            rain =  "R"
            snow =  "S"

        elif self.ptype == "Rain/Snow" and self.convective == "Showers":
            wxval = "RW"
            rain =  "RW"
            snow =  "SW"

        if self.ptype == "Snow" and self.convective == "Stratiform":
            wxval = "S"
            snow =  "S"
            rain =  "R"

        elif self.ptype == "Snow" and self.convective == "Showers":
            wxval = "SW"
            snow =  "SW"
            rain =  "RW"

        # For convective ZR, there is no ZRW
        if self.ptype == "Freezing Rain" and self.convective == "Stratiform":
            wxval = "ZR"
        elif self.ptype == "Freezing Rain" and self.convective == "Showers":
            wxval = "ZR"

        # For convective IP, there is no IPW
        if self.ptype == "Sleet" and self.convective == "Stratiform":
            wxval = "IP"
        elif self.ptype == "Sleet" and self.convective == "Showers":
            wxval = "IP"

        # Set the max and min elevations for rain/snow mix
        # if rsmix = 0, then rsmin = rsmax = snlvl no rain/snow mix gets assigned
        rs_bot = self.SnowLevel + self.botRainSnow
        rs_top = self.SnowLevel + self.topRainSnow

        # set masks for rain, mix and snow
        rainMask    = self.Topo <= rs_bot
        snowMask    = self.Topo >= rs_top
        mixMask     = (self.Topo >= rs_bot) & (self.Topo < rs_top)
        rainButCold = ((rainMask | mixMask)) & (self.maxT < 30.0)
        snowMask |= rainButCold
        rainMask &= np.logical_not(rainButCold)
        mixMask &= np.logical_not(rainButCold)
        snowButWarm = ((snowMask | mixMask)) & (self.minT > 45.0)
        rainMask |= snowButWarm
        snowMask &= np.logical_not(snowButWarm)
        mixMask &= np.logical_not(snowButWarm)
        rainMask &= (editAreaMask > 0.5)
        snowMask &= (editAreaMask > 0.5)
        mixMask &= (editAreaMask > 0.5)

        # Assign intensity
        if self.intensity == "Light":
            inten = "-"
        elif self.intensity == "Moderate":
            inten = "m"
        elif self.intensity == "Heavy":
            inten = "+"
        elif self.intensity == "Sprinkles/Flurries":
            inten = "--"
        else:
            inten = "-"

        #set the visibility
        vis = ":<NoVis>:"

        # If thunder is selected for rain or snow event and is convective type
        if self.thunder == "Y" and self.convective == "Showers" and inten != "--":
            snow = "SW,T"
            rain = "RW,T"

        # If thunder is selected in a stratiform event, turn thunder off
        if self.thunder == "Y" and self.convective == "Stratiform":
            dialog = OKDialog.OKDialog(None, "Urgent", "  Can't have Thunder with stratiform precip.  ")
            self.thunder = "N"

        # If thunder is selected in a sprinkle/flurry event, turn thunder off
        if self.thunder == "Y" and inten == "--":
            dialog = OKDialog.OKDialog(None, "Urgent", "  Can't have Thunder with sprinkles.  ")
            self.thunder = "N"

        #  Masks for PoP categories 0-4
        pmask=[]
        pmask.append(self.PoP < self.popThresh)
        pmask.append( (self.PoP >= self.popThresh) & (self.PoP <= 24))
        pmask.append( (self.PoP > 24) & (self.PoP <= 54))
        pmask.append( (self.PoP > 54) & (self.PoP <= 74))
        pmask.append(self.PoP > 74)

        covs=["", "Iso", "Sct", "Num", "Wide"]
        probs=["", "SChc", "Chc", "Lkly", "Def"]

        ltgtype=covs[self.ltgcat]
        if self.coverage=="Probability":
            ltgtype=probs[self.ltgcat]

        descriptors=probs
        if ((self.coverage=="Coverage")and(self.convective=="Showers")):
            descriptors=covs

        #  Loop over the 5 coverage classes (none, slight, chance, likely, categorical)
        for cat in range(5):

            #  Determine the thunder (lightning) category.
            lcat=0
            if self.thunderAnswer=="Force":
                lcat=self.forceThunderLevel
            elif self.thunderAnswer=="Match":
                lcat=min(cat, self.matchThunderLevel)

            ltgtype=covs[lcat]
            if self.thunderType=="Probability":
                ltgtype=probs[lcat]

            fcat=0
            if self.freezing=="Force":
                fcat=cat
            elif self.freezing=="Mix":
                fcat=min(cat, self.fmmatchLev)
            fztype=probs[fcat]

            icat=0
            if self.sleet=="Force":
                icat=cat
            elif self.sleet=="Mix":
                icat=min(cat, self.immatchLev)
            iptype=probs[icat]

            #  For none category - could still have thunder
            if cat==0:
                if self.thunderAnswer == "FromPot":
                    for pot_cat in range(1, 5):
                        addthunder = pmask[cat] & (editAreaMask > 0.5) & potMasks[pot_cat]
                        if addthunder.any():
                            pot_ltgtype = covs[pot_cat] if self.thunderType == "Coverage" else probs[pot_cat]
                            wxstr = f"{pot_ltgtype}:T:<NoInten>:<NoVis>:"
                            idx = self.getIndex(wxstr, wxStrings)
                            wxValues[addthunder] = idx
                            if pot_ltgtype in THUNDER_TO_POTTHUNDER:
                                self.thunderAssignments[addthunder] = THUNDER_TO_POTTHUNDER[pot_ltgtype]
                elif (self.thunderAnswer=="Force") and (lcat>0):
                   addthunder = pmask[cat] & (editAreaMask > 0.5)
                   wxstr = f"{ltgtype}:T:<NoInten>:<NoVis>:"
                   idx = self.getIndex(wxstr, wxStrings)
                   wxValues[addthunder] = idx
                   if ltgtype in THUNDER_TO_POTTHUNDER:
                       self.thunderAssignments[addthunder] = THUNDER_TO_POTTHUNDER[ltgtype]
                continue

            #  Rain only areas
            prain= pmask[cat] & rainMask
            if prain.any():
                wxtype = "RW" if self.convective == "Showers" else "R"
                base_wxstr = f"{descriptors[cat]}:{wxtype}:{inten}:<NoVis>:"

                if self.freezing in ("Force", "Mix"):
                    fzstring = f"{fztype}:ZR:{inten}:<NoVis>:"
                    base_wxstr = fzstring if self.freezing == "Force" else f"{base_wxstr}^{fzstring}"
                if self.sleet in ("Force", "Mix"):
                    ipstring = f"{iptype}:IP:{inten}:<NoVis>:"
                    base_wxstr = ipstring if self.sleet == "Force" else f"{base_wxstr}^{ipstring}"

                if self.thunderAnswer == "FromPot":
                    for pot_cat in range(1, 5):
                        pot_prain = prain & potMasks[pot_cat]
                        if pot_prain.any():
                            pot_ltgtype = covs[pot_cat] if self.thunderType == "Coverage" else probs[pot_cat]
                            wxstr = f"{base_wxstr}^{pot_ltgtype}:T:<NoInten>:<NoVis>:"
                            idx = self.getIndex(wxstr, wxStrings)
                            wxValues[pot_prain] = idx
                            if pot_ltgtype in THUNDER_TO_POTTHUNDER:
                                self.thunderAssignments[pot_prain] = THUNDER_TO_POTTHUNDER[pot_ltgtype]

                    no_pot_prain = prain & potMasks[0]
                    if no_pot_prain.any():
                        idx = self.getIndex(base_wxstr, wxStrings)
                        wxValues[no_pot_prain] = idx
                else:
                    wxstr = base_wxstr
                    if self.thunderAnswer in ("Force", "Match"):
                        wxstr += f"^{ltgtype}:T:<NoInten>:<NoVis>:"
                        if ltgtype in THUNDER_TO_POTTHUNDER:
                            self.thunderAssignments[prain] = THUNDER_TO_POTTHUNDER[ltgtype]

                    idx = self.getIndex(wxstr, wxStrings)
                    wxValues[prain] = idx

            #  Snow only areas
            psnow= pmask[cat] & snowMask
            if psnow.any():
                wxtype = "SW" if self.convective == "Showers" else "S"
                base_wxstr = f"{descriptors[cat]}:{wxtype}:{inten}:<NoVis>:"

                if self.freezing in ("Force", "Mix"):
                    fzstring = f"{fztype}:ZR:{inten}:<NoVis>:"
                    base_wxstr = fzstring if self.freezing == "Force" else f"{base_wxstr}^{fzstring}"
                if self.sleet in ("Force", "Mix"):
                    ipstring = f"{iptype}:IP:{inten}:<NoVis>:"
                    base_wxstr = ipstring if self.sleet == "Force" else f"{base_wxstr}^{ipstring}"

                if self.thunderAnswer == "FromPot":
                    for pot_cat in range(1, 5):
                        pot_psnow = psnow & potMasks[pot_cat]
                        if pot_psnow.any():
                            pot_ltgtype = covs[pot_cat] if self.thunderType == "Coverage" else probs[pot_cat]
                            wxstr = f"{base_wxstr}^{pot_ltgtype}:T:<NoInten>:<NoVis>:"
                            idx = self.getIndex(wxstr, wxStrings)
                            wxValues[pot_psnow] = idx
                            if pot_ltgtype in THUNDER_TO_POTTHUNDER:
                                self.thunderAssignments[pot_psnow] = THUNDER_TO_POTTHUNDER[pot_ltgtype]

                    no_pot_psnow = psnow & potMasks[0]
                    if no_pot_psnow.any():
                        idx = self.getIndex(base_wxstr, wxStrings)
                        wxValues[no_pot_psnow] = idx
                else:
                    wxstr = base_wxstr
                    if self.thunderAnswer in ("Force", "Match"):
                        wxstr += f"^{ltgtype}:T:<NoInten>:<NoVis>:"
                        if ltgtype in THUNDER_TO_POTTHUNDER:
                            self.thunderAssignments[psnow] = THUNDER_TO_POTTHUNDER[ltgtype]

                    idx = self.getIndex(wxstr, wxStrings)
                    wxValues[psnow] = idx

            #  Rain/Snow mix areas
            pmix= pmask[cat] & mixMask
            if pmix.any():
                if self.convective == "Showers":
                    base_wxstr = f"{descriptors[cat]}:RW:{inten}:<NoVis>:^{descriptors[cat]}:SW:{inten}:<NoVis>:"
                else:
                    base_wxstr = f"{descriptors[cat]}:R:{inten}:<NoVis>:^{descriptors[cat]}:S:{inten}:<NoVis>:"

                if self.freezing in ("Force", "Mix"):
                    fzstring = f"{fztype}:ZR:{inten}:<NoVis>:"
                    base_wxstr = fzstring if self.freezing == "Force" else f"{base_wxstr}^{fzstring}"
                if self.sleet in ("Force", "Mix"):
                    ipstring = f"{iptype}:IP:{inten}:<NoVis>:"
                    base_wxstr = ipstring if self.sleet == "Force" else f"{base_wxstr}^{ipstring}"

                if self.thunderAnswer == "FromPot":
                    for pot_cat in range(1, 5):
                        pot_pmix = pmix & potMasks[pot_cat]
                        if pot_pmix.any():
                            pot_ltgtype = covs[pot_cat] if self.thunderType == "Coverage" else probs[pot_cat]
                            wxstr = f"{base_wxstr}^{pot_ltgtype}:T:<NoInten>:<NoVis>:"
                            idx = self.getIndex(wxstr, wxStrings)
                            wxValues[pot_pmix] = idx
                            if pot_ltgtype in THUNDER_TO_POTTHUNDER:
                                self.thunderAssignments[pot_pmix] = THUNDER_TO_POTTHUNDER[pot_ltgtype]

                    no_pot_pmix = pmix & potMasks[0]
                    if no_pot_pmix.any():
                        idx = self.getIndex(base_wxstr, wxStrings)
                        wxValues[no_pot_pmix] = idx
                else:
                    wxstr = base_wxstr
                    if self.thunderAnswer in ("Force", "Match"):
                        wxstr += f"^{ltgtype}:T:<NoInten>:<NoVis>:"
                        if ltgtype in THUNDER_TO_POTTHUNDER:
                            self.thunderAssignments[pmix] = THUNDER_TO_POTTHUNDER[ltgtype]

                    idx = self.getIndex(wxstr, wxStrings)
                    wxValues[pmix] = idx

        ########################################################################
        #  Determine Wx intensities based upon QPF and SnowAmt if requested by user
        wxValues, wxStrings = self.WxIntensity(snowMask, wxValues, wxStrings)

        #Now merge the oldWx and newWx grids
        (newWxValues, newWxStrings) = self.mergeWx(oldWxValues, oldWxStrings, wxValues, wxStrings)

        (newWxValues, newWxStrings) = self.compactWx(newWxValues, newWxStrings)

        # Save the new grid
        self.Wx = [newWxValues, newWxStrings]
        self.createGrid(self.mutableName, "Wx", "WEATHER", self.Wx, self.GridTimeRange)

        # Update PotThunder grid based on thunder assignments
        self.updatePotThunder(editAreaMask)

        return

    def updatePotThunder(self, editAreaMask):
        """Update PotThunder grid to ensure it supports the Wx thunder coverage."""
        if self.potThunder is None:
            self.potThunder = self.newGrid(0.0)

        newPotThunder = self.potThunder.copy()

        # Mask where thunder was assigned in the Wx grid inside the edit area
        thunderMask = (self.thunderAssignments > 0) & (editAreaMask > 0.5)

        if np.any(thunderMask):
            # Take the maximum of the current PotThunder or the baseline required for the new Wx
            current_vals = newPotThunder[thunderMask]
            required_vals = self.thunderAssignments[thunderMask]
            newPotThunder[thunderMask] = np.maximum(current_vals, required_vals)

        # Clear PotThunder from regions that had it removed by Wx updates
        noThunderMask = (self.thunderAssignments == 0) & (editAreaMask > 0.5)
        newPotThunder[noThunderMask & (newPotThunder >= 14.5)] = 14

        # Save the updated PotThunder grid
        self.createGrid(self.mutableName, "PotThunder", "SCALAR", newPotThunder, self.GridTimeRange)


# Merges the oldWx and newWx grids
    def mergeWx(self, oldWxValues, oldWxKeys, newWxValues, newWxKeys):
        # create the final Wx grid
        finalWxKeys = []
        finalWxValues = self.newGrid(0)

        emptyWxList=["<NoCov>:<NoWx>:<NoInten>:<NoVis>:",]
        # loop through each index in the newWxKeys
        for newIndex in range(len(newWxKeys)):
            # loop through each index in the oldWxKeys
            for oldIndex in range(len(oldWxKeys)):
                # see if there are any gridpoints that the newWx in the new grid and the
                # oldWx in the old grid.  If so, we'll need to combine them.
                valuesBoth= (oldWxValues == oldIndex) & (newWxValues == newIndex)
                needBoth=np.zeros(valuesBoth.shape)
                needBoth[valuesBoth]=1
                existBoth=valuesBoth.any()
                numneeded=np.sum(needBoth)

                if existBoth:
                    oldString=oldWxKeys[oldIndex]
                    newString=newWxKeys[newIndex]
                    if newString==oldString:
                       newUglyStr = oldString
                    else:
                       # handle the case if BOTH are NoWx
                       if ((oldString in emptyWxList) and (newString in emptyWxList)):
                          newUglyStr = emptyWxList[0]
                       else:
                          # handle the cases if EITHER are NoWx
                          if oldString in emptyWxList:
                             newUglyStr = newString
                          elif newString in emptyWxList:
                             newUglyStr = oldString
                          else:
                             newUglyStr = newString+"^"+oldString

                    index = self.getIndex(newUglyStr, finalWxKeys)
                    finalWxValues[valuesBoth]=index
                else:
                    pass
        return finalWxValues, finalWxKeys

    #==========================================================================
    #  Get rid of any wx keys that are not used (but always keep <NoWx>)
    #
    def compactWx(self, WxValues, WxStrings):
        NewWxStrings = []
        NewWxValues = self.newGrid(0)

        numremoved=0
        emptyWx="<NoCov>:<NoWx>:<NoInten>:<NoVis>:"
        # loop through each index in the newWxKeys
        for idx in range(len(WxStrings)):
            wxString=WxStrings[idx]
            gotidx=WxValues == idx
            pts=np.sum(gotidx)
            if ((pts>0)or(wxString==emptyWx)):
                newIdx=self.getIndex(wxString, NewWxStrings)
                NewWxValues[gotidx]=newIdx
            else:
                numremoved+=1
        return NewWxValues, NewWxStrings

    #=====================================================================
    #  getEditAreaMask - return mask of 0/1 where current editArea exists
    #
    def getEditAreaMask(self):
        editArea=self.getActiveEditArea()
        editAreaObj=editArea.getGrid()
        if not editAreaObj.isAnyBitsSet():
            editArea.invert()
        editAreaMask=editArea.getGrid().getNDArray()
        return editAreaMask

    def getEmptyValue(self, keys):
        uglyString = "<NoCov>:<NoWx>:<NoInten>:<NoVis>:"
        return self.getIndex(uglyString, keys)

    def getByteValue(self, cov, wxval, inten, keys):
        if not cov:
            uglyString = "<NoCov>:<NoWx>:<NoInten>:<NoVis>:"
        else:
            uglyString = cov + ":" + wxval + ":" + inten + ":<NoVis>:"
            if self.thunder == "Y":
                uglyString += "^" + cov + ":T:<NoInten>:<NoVis>:"
        return self.getIndex(uglyString, keys)


    def keepWxInGrid(self, WxGrid, editAreaMask):
        wxList = self.keepWxList.split(",")

        wxValues, keys = WxGrid
        WxOut = []
        emptyWx="<NoCov>:<NoWx>:<NoInten>:<NoVis>:"

        numkeys=len(keys)
        trans=[]
        for i in range(numkeys):
            trans.append(i) # default is no change
            keyInArea= (wxValues == i) & (editAreaMask > 0.5)
            usedInArea=np.sum(keyInArea)
            if usedInArea<1:
                continue
            changeneeded=0
            codes=keys[i].split('^')
            newcodes=[]
            for j in range(len(codes)):
                (cover, type, inten, vis, attr)=codes[j].split(":")
                found = 0
                for k in wxList:
                    if (type==k):
                        found = 1
                        newcodes.append(codes[j])
                        break
                if found==0:
                    changeneeded=1
            if (changeneeded==1):
                if (not newcodes):
                    trans[i]=self.getIndex(emptyWx, keys)
                else:
                   newkey="^".join(newcodes)
                   trans[i]=self.getIndex(newkey, keys)
            else:
                trans[i]=i

        newwxValues=wxValues
        for i in range(numkeys):
           keyInArea= (wxValues == i) & (editAreaMask > 0.5)
           newwxValues[keyInArea]=trans[i]
        return (newwxValues, keys)

    ############################################################################
    #  modify Wx intensity based upon QPF and/or SnowAmt

    def WxIntensity(self, snowMask, wxValues, wxKeys):
        PoP = self.getGrids(self.mutableName, 'PoP', 'SFC', self.GridTimeRange,
                            noDataError=0)
        QPF = self.getGrids(self.mutableName, 'QPF', 'SFC', self.GridTimeRange, mode="Sum",
                            noDataError=0)

        if QPF is None:
            return (wxValues, wxKeys)

        try:
            QPF[QPF<0.01]=0.00
        except Exception:
            QPF=self.newGrid(0)

        SnowAmt = self.getGrids(self.mutableName, 'SnowAmt', 'SFC', self.GridTimeRange, mode="Sum",
                            noDataError=0)
        if SnowAmt is None:
            return (wxValues, wxKeys)

        valleyMask = self.Topo <= valleyMtnThreshold
        mtnMask = self.Topo > valleyMtnThreshold

        startQPF = curGridStart = self.GridTimeRange.startTime()
        duration = float(self.GridTimeRange.duration()) / 3600.0

        if duration > 6.0:
            QPF *= (6.0/duration)
            SnowAmt *= (6.0/duration)

        modPcpnDict = {}
        hvyPcpnDict = {}

        modMask = self.newGrid(0)
        hvyMask = self.newGrid(0)

        Amount=QPF.copy()
        Amount[snowMask]=SnowAmt[snowMask]

        modRthresh = self.newGrid(moderateRthresh*2)
        modRthresh[valleyMask]=moderateRthresh
        hvyRthresh = self.newGrid(heavyRthresh*2)
        hvyRthresh[valleyMask]=heavyRthresh
        modSthresh = self.newGrid(moderateSthresh*2)
        modSthresh[valleyMask]=moderateSthresh
        hvySthresh = self.newGrid(heavySthresh*2)
        hvySthresh[valleyMask]=heavySthresh

        modThresh = modRthresh.copy()
        modThresh[snowMask]=modSthresh[snowMask]
        hvyThresh = hvyRthresh.copy()
        hvyThresh[snowMask]=hvySthresh[snowMask]

        modMask = ((self.PoP >= self.popThresh) & (Amount >= modThresh) & (Amount < hvyThresh))
        hvyMask = ((self.PoP >= self.popThresh) & (Amount >= hvyThresh))

        for wx in wxKeys:
            modPcpn = wx
            hvyPcpn = wx

            if ':-:' in wx:
                intenMatch = self._weatherInten.findall(wx)

                for group in intenMatch:
                    modPcpn = modPcpn.replace(':%s:-:' % (group),
                                              ':%s:m:' % (group))
                    hvyPcpn = hvyPcpn.replace(':%s:-:' % (group),
                                              ':%s:+:' % (group))

                modPcpn = wx.replace(':-:', ':m:')
                hvyPcpn = wx.replace(':-:', ':+:')

                modPcpnDict[wx] = modPcpn
                hvyPcpnDict[wx] = hvyPcpn

            wxMask = wxValues == self.getIndex(wx, wxKeys)

            if wx in modPcpnDict:
                modWxMask = (wxMask & modMask)
                numModPoints = np.sum(np.sum(modWxMask, 1))
                if numModPoints > 0:
                    wxValues[modWxMask]=self.getIndex(modPcpnDict[wx], wxKeys)

            if wx in hvyPcpnDict:
                hvyWxMask = wxMask & hvyMask
                numHvyPoints = np.sum(np.sum(hvyWxMask, 1))
                if numHvyPoints > 0:
                    wxValues[hvyWxMask]=self.getIndex(hvyPcpnDict[wx], wxKeys)

        return wxValues, wxKeys

    def ComplexWxString(self, description):
        coverages, wxTypes, intens, vis, attrs = \
                   complexConvertWxDescription(description)

        if "TRW" in wxTypes:
            wxTypes.remove("TRW")
            wxTypes.append("T")
            wxTypes.append("RW")

        if coverages[0] == "*":
            coverages[0] = "<NoCov>"
        if wxTypes[0] == "*":
            wxTypes[0] = "<NoWx>"
        if intens[0] == "*":
            intens[0] = "<NoInten>"
        if vis[0] == "*":
            vis[0] = "<NoVis>"

        if attrs[0] == '*':
            attrStr = ""
        else:
            first = 1
            attrStr = ""
            for attr in attrs:
                if not first:
                    attrStr += ","
                first = 0
                attrStr += attr

        Wx = ""
        first = 1
        for wxType in wxTypes:
            if not first:
                Wx += "^"
            first = 0
            intensity = intens[0]
            if wxType == "T" and not intensity == "+":
                intensity = "<NoInten>"
            if self.limitT == 1 and wxType == "T":
                if coverages[0] not in self.ltgcov:
                    Wx = Wx + self.ltgcov[-1]+":"+wxType+":"+intensity+":"+vis[0]+":"\
                         +attrStr
                else:
                    Wx = Wx + coverages[0]+":"+wxType+":"+intensity+":"+vis[0]+":"\
                         +attrStr
            else:
                Wx = Wx + coverages[0]+":"+wxType+":"+intensity+":"+vis[0]+":"\
                 +attrStr

        return Wx

    def GuiDefaults(self):
        self.GuiSettings["coverage"]    ="Probability"  # Coverage or Probability
        self.GuiSettings["convective"]  ="Showers"      # Showers or Stratiform
        self.GuiSettings["mixState"]    =0              # 0=closed 1=open
        self.GuiSettings["botRainSnow"] =-400           # int 0 - -1000
        self.GuiSettings["topRainSnow"] =400            # int 0 - 1000
        self.GuiSettings["thunderState"]=0              # 0=closed 1=open
        self.GuiSettings["thunderType"] ="Probability"  # Coverage of Probability
        self.GuiSettings["thunder"]     ="Keep"         # None, Keep, Force, Match, FromPot
        self.GuiSettings["tforceLev"]   =1              # int 1-4 1=iso 4=wide
        self.GuiSettings["tmatchLev"]   =4              # int 1-4 1=iso 4=wide
        self.GuiSettings["fzState"]     =0              # 0=closed 1=open
        self.GuiSettings["freezing"]    ="None"         # None, Keep, Force, Mix
        self.GuiSettings["fmmatchLev"]  =3              # int 1-4 1=schc 4-def
        self.GuiSettings["ipState"]     =0              # 0=closed 1=open
        self.GuiSettings["sleet"]       ="None"         # None, Keep, Force, Mix
        self.GuiSettings["immatchLev"]  =3              # int 1-4 1=schc 4-def
        return

    def setupGui(self):
        if self.GuiSettings["coverage"]=="Coverage":
            self.dialog.covbut.invoke()
        else:
            self.dialog.prbbut.invoke()

        if self.GuiSettings["mixState"]==0:
            self.dialog.mixState.set(1)
            self.dialog.changeMixState()
        else:
            self.dialog.mixState.set(0)
            self.dialog.changeMixState()

        self.dialog.botRainSnow.set(self.GuiSettings["botRainSnow"])
        self.dialog.bottomMove(None)

        self.dialog.topRainSnow.set(self.GuiSettings["topRainSnow"])
        self.dialog.topMove(None)

        if self.GuiSettings["thunderState"]==0:
            self.dialog.thunderState.set(1)
            self.dialog.changeThunderState()
        else:
            self.dialog.thunderState.set(0)
            self.dialog.changeThunderState()

        if self.GuiSettings["thunderType"]=="Coverage":
            self.dialog.tcb.invoke()
        else:
            self.dialog.tpb.invoke()

        thunderValue=self.GuiSettings["thunder"]
        self.dialog.thunder.set(thunderValue)
        if thunderValue=="None":
            self.dialog.th1.invoke()
        elif thunderValue=="Keep":
            self.dialog.th2.invoke()
        elif thunderValue=="Force":
            self.dialog.th3.invoke()
        elif thunderValue=="Match":
            self.dialog.th4.invoke()
        elif thunderValue=="FromPot":
            self.dialog.th5.invoke()

        tf=self.GuiSettings["tforceLev"]
        self.dialog.tforceLev.set(tf)
        if tf==1:
            self.dialog.tf1.select()
        elif tf==2:
            self.dialog.tf2.select()
        elif tf==3:
            self.dialog.tf3.select()
        elif tf==4:
            self.dialog.tf4.select()

        tm=self.GuiSettings["tmatchLev"]
        self.dialog.tmatchLev.set(tm)
        if tm==1:
            self.dialog.tm1.select()
        elif tm==2:
            self.dialog.tm2.select()
        elif tm==3:
            self.dialog.tm3.select()
        elif tm==4:
            self.dialog.tm4.select()

        if self.GuiSettings["convective"]=="Showers":
            self.dialog.shobut.invoke()
        else:
            self.dialog.strbut.invoke()

        if self.GuiSettings["fzState"]==0:
            self.dialog.fzState.set(1)
            self.dialog.changeFzState()
        else:
            self.dialog.fzState.set(0)
            self.dialog.changeFzState()

        self.dialog.freezing.set(self.GuiSettings["freezing"])
        self.dialog.changeFreezing()

        fm=self.GuiSettings["fmmatchLev"]
        self.dialog.fmmatchLev.set(fm)
        if fm==1:
            self.dialog.fmm1.select()
        elif fm==2:
            self.dialog.fmm2.select()
        elif fm==3:
            self.dialog.fmm3.select()
        elif fm==4:
            self.dialog.fmm4.select()
        self.dialog.changeFreezing()

        if self.GuiSettings["ipState"]==0:
            self.dialog.ipState.set(1)
            self.dialog.changeIpState()
        else:
            self.dialog.ipState.set(0)
            self.dialog.changeIpState()

        self.dialog.sleet.set(self.GuiSettings["sleet"])
        self.dialog.changeSleet()

        im=self.GuiSettings["immatchLev"]
        self.dialog.immatchLev.set(im)
        if im==1:
            self.dialog.imm1.select()
        elif im==2:
            self.dialog.imm2.select()
        elif im==3:
            self.dialog.imm3.select()
        elif im==4:
            self.dialog.imm4.select()
        self.dialog.changeSleet()
        return

    def saveStatus(self):
        self.GuiSettings["coverage"]    =self.dialog.coverage.get()
        self.GuiSettings["convective"]  =self.dialog.convective.get()
        self.GuiSettings["mixState"]    =self.dialog.mixState.get()
        self.GuiSettings["botRainSnow"] =self.dialog.botRainSnow.get()
        self.GuiSettings["topRainSnow"] =self.dialog.topRainSnow.get()
        self.GuiSettings["thunderState"]=self.dialog.thunderState.get()
        self.GuiSettings["thunderType"] =self.dialog.thunderType.get()
        self.GuiSettings["thunder"]     =self.dialog.thunder.get()
        self.GuiSettings["tforceLev"]   =self.dialog.tforceLev.get()
        self.GuiSettings["tmatchLev"]   =self.dialog.tmatchLev.get()
        self.GuiSettings["fzState"]     =self.dialog.fzState.get()
        self.GuiSettings["freezing"]    =self.dialog.freezing.get()
        self.GuiSettings["fmmatchLev"]  =self.dialog.fmmatchLev.get()
        self.GuiSettings["ipState"]     =self.dialog.ipState.get()
        self.GuiSettings["sleet"]       =self.dialog.sleet.get()
        self.GuiSettings["immatchLev"]  =self.dialog.immatchLev.get()

        objectDirName="%sSettings"%TOOLNAME
        self.saveObject("GuiSettings", self.GuiSettings, objectDirName)
        return

    def readStatus(self):
        objectDirName="%sSettings"%TOOLNAME
        lf=self.getUserFile("GuiSettings", objectDirName)
        if lf.exists():
            tempstate=self.getObject("GuiSettings", objectDirName)
            for key in tempstate:
                if key in self.GuiSettings:
                    self.GuiSettings[key]=tempstate[key]
        return

    def getUserFile(self, name, category):
        from com.raytheon.uf.common.localization import PathManagerFactory, LocalizationContext
        LocalizationType=LocalizationContext.LocalizationType
        LocalizationLevel=LocalizationContext.LocalizationLevel
        pathMgr = PathManagerFactory.getPathManager()
        path = 'gfe/userPython/' + category + '/' + name
        lc = pathMgr.getContext(LocalizationType.valueOf('CAVE_STATIC'), LocalizationLevel.valueOf('USER'))
        lf = pathMgr.getLocalizationFile(lc, path)
        return lf

class WxTool(AppDialog.AppDialog):
    def __init__(self, title="nonModal Dialog", callbackMethod=None, **kwargs):
        self.__callbackMethod = callbackMethod
        AppDialog.AppDialog.__init__(self, **kwargs)
        self.title(title)
        return

    def buttonbox(self):
        buttonFrame = tk.Frame(self)
        tk.Button(buttonFrame, text="Run",
            command=self.__runCB, width=10, state=tk.NORMAL).pack(\
            side=tk.LEFT, pady=5, padx=10)
        tk.Button(buttonFrame, text="Run/Dismiss",
            command=self.__okCB, width=12, state=tk.NORMAL).pack(\
            side=tk.LEFT, pady=5, padx=10)
        tk.Button(buttonFrame, text="Cancel", width=10,
            command=self.cancelCB).pack(\
            side=tk.LEFT, pady=5, padx=10)
        buttonFrame.pack(side=tk.BOTTOM)

    def __runCB(self):
        self.validate()
        self.__callbackMethod("Run")

    def __okCB(self):
        self.validate()
        self.__callbackMethod("OK")
        self.ok()

    def cancelCB(self):
        self.__callbackMethod("Cancel")
        self.destroy()

    def cancel(self, event=None):
        self.__callbackMethod("Cancel")
        self.destroy()

    def body(self, master):
        cfont = 'Arial 14'
        dfont = 'Arial 6'
        bgcol = "#B4B4B4"

        mainFrame=tk.Frame(master)

        #################################
        # Coverage or Probability

        coverageFrame=tk.Frame(master, relief=tk.GROOVE, borderwidth=3, background=bgcol)
        self.coverage=tk.StringVar()
        g=tk.Label(coverageFrame, text="Coverage or Probability?",
                        font=cfont, bg=bgcol, width=25, height=3)
        g.pack(side=tk.LEFT)
        self.covbut=tk.Radiobutton(coverageFrame, text="Coverage", variable=self.coverage, value="Coverage", bg=bgcol,
                                        highlightthickness=0, height=1, command=self.changeCoverage)
        self.covbut.pack(side=tk.LEFT)
        self.prbbut=tk.Radiobutton(coverageFrame, text="Probability", variable=self.coverage, value="Probability", bg=bgcol,
                                        highlightthickness=0, height=1, command=self.changeCoverage)
        self.prbbut.pack(side=tk.LEFT)
        coverageFrame.pack(side=tk.TOP, fill=tk.BOTH, expand=0)

        #################################
        # Stratiform or Convective

        convectiveFrame=tk.Frame(master, relief=tk.GROOVE, borderwidth=3, background=bgcol)
        self.convective=tk.StringVar()
        g=tk.Label(convectiveFrame, text="Showers or Stratiform?",
                        font=cfont, bg=bgcol, width=25, height=3)
        g.pack(side=tk.LEFT)
        self.shobut=tk.Radiobutton(convectiveFrame, text="Showers", variable=self.convective, value="Showers", bg=bgcol,
                                        highlightthickness=0, height=1, command=self.changeConvective)
        self.shobut.pack(side=tk.LEFT)
        self.strbut=tk.Radiobutton(convectiveFrame, text="Stratiform", variable=self.convective, value="Stratiform", bg=bgcol,
                                        highlightthickness=0, height=1, command=self.changeConvective)
        self.strbut.pack(side=tk.LEFT)
        convectiveFrame.pack(side=tk.TOP, fill=tk.BOTH, expand=0)

        self.ptype=tk.StringVar()
        self.ptype.set("Rain/Snow")
        self.popThresh=tk.IntVar()
        self.popThresh.set(15)
        self.intensity=tk.StringVar()
        self.intensity.set("Light")
        self.keepWx=tk.StringVar()
        self.keepWx.set("Y")
        self.keepWxList=tk.StringVar()
        self.keepWxList.set(KeepList)

        sliderFrame=tk.Frame(master, relief=tk.GROOVE, borderwidth=3, background=bgcol)
        mixOpen=tk.Frame(sliderFrame, background=bgcol)
        self.mixState=tk.IntVar()
        self.mixState.set(0)
        self.mixButton=tk.Button(mixOpen, text=">", padx=0, pady=0, command=self.changeMixState)
        self.mixButton.pack(side=tk.LEFT, anchor=tk.W)
        g=tk.Label(mixOpen, text="Rain-Snow-Mix around Snowlevel", font=cfont, bg=bgcol)
        g.pack(side=tk.LEFT, anchor=tk.W)
        mixOpen.pack(side=tk.TOP, anchor=tk.W)

        self.mixChoice=tk.Frame(sliderFrame, background=bgcol)

        botSliderFrame=tk.Frame(self.mixChoice, background=bgcol)
        b=tk.Label(botSliderFrame, text="Bottom",
                        font=cfont, bg=bgcol, height=1)
        b.pack(side=tk.TOP)
        self.botRainSnow=tk.IntVar()
        botSliderObject=tk.Scale(botSliderFrame, from_=maxThicknessRainSnowLayer, to=maxThicknessRainSnowLayer*-1.0,
                                      orient=tk.VERTICAL, showvalue=1,
                                      variable=self.botRainSnow, bg=bgcol, resolution=100, command=self.bottomMove,
                          )
        self.botRainSnow.set(defaultBotRainSnow)
        botSliderObject.pack(side=tk.LEFT)
        botSliderFrame.pack(side=tk.RIGHT, fill=tk.BOTH)

        soundingFrame=tk.Frame(self.mixChoice, background=bgcol)
        self.soundingObject=tk.Canvas(soundingFrame, width=100, height=100, background="black",
                          )
        self.soundingObject.pack(side=tk.BOTTOM)
        soundingFrame.pack(side=tk.RIGHT, fill=tk.BOTH)

        topSliderFrame=tk.Frame(self.mixChoice, background=bgcol)
        t=tk.Label(topSliderFrame, text="Top",
                        font=cfont, bg=bgcol, height=1)
        t.pack(side=tk.TOP)
        self.topRainSnow=tk.IntVar()
        topSliderObject=tk.Scale(topSliderFrame, from_=maxThicknessRainSnowLayer, to=maxThicknessRainSnowLayer*-1.0,
                                      orient=tk.VERTICAL, showvalue=1,
                                      variable=self.topRainSnow, bg=bgcol, resolution=100, command=self.topMove,
                          )
        self.topRainSnow.set(defaultTopRainSnow)
        topSliderObject.pack(side=tk.BOTTOM)
        topSliderFrame.pack(side=tk.RIGHT, fill=tk.BOTH)

        labelSliderFrame=tk.Frame(self.mixChoice, background=bgcol)
        l3=tk.Label(labelSliderFrame, text="", font=dfont, bg=bgcol, width=5, height=2)
        l3.pack(side=tk.BOTTOM)
        rainlabel=tk.Label(labelSliderFrame, text="Rain",
                        font=cfont, bg="green", width=5, height=1)
        rainlabel.pack(side=tk.BOTTOM)
        mixlabel=tk.Label(labelSliderFrame, text="Mix",
                        font=cfont, bg="#FF00FF", width=5, height=1)
        mixlabel.pack(side=tk.BOTTOM)
        snowlabel=tk.Label(labelSliderFrame, text="Snow",
                        font=cfont, bg="cyan", width=5, height=1)
        snowlabel.pack(side=tk.BOTTOM)
        labelSliderFrame.pack(side=tk.RIGHT, anchor=tk.E, ipadx=10, fill=tk.BOTH)

        padFrame=tk.Frame(self.mixChoice, background=bgcol)
        padFrame.pack(side=tk.LEFT, fill=tk.BOTH, expand=tk.YES)

        self.mixChoice.pack(side=tk.TOP, anchor=tk.W, fill=tk.BOTH, expand=tk.YES)
        if self.mixState.get()==0:
            self.mixChoice.pack_forget()

        sliderFrame.pack(side=tk.TOP, fill=tk.BOTH, expand=0)

        #################################
        # Thunder
        thunderFrame=tk.Frame(master, relief=tk.GROOVE, borderwidth=3, background=bgcol)

        self.thunderOpen=tk.Frame(thunderFrame, background=bgcol)
        self.thunderState=tk.IntVar()
        self.thunderState.set(0)
        self.thunderButton=tk.Button(self.thunderOpen, text=">", padx=0, pady=0, command=self.changeThunderState)
        self.thunderButton.pack(side=tk.LEFT, anchor=tk.W)
        self.thunderLabel=tk.Label(self.thunderOpen, text="Thunder", anchor=tk.W, font=cfont, bg=bgcol)
        self.thunderLabel.pack(side=tk.LEFT, anchor=tk.W, fill=tk.X, expand=1)
        self.thunderOpen.pack(side=tk.TOP, anchor=tk.W, fill=tk.X, expand=1)

        self.thunderChoice=tk.Frame(thunderFrame, background=bgcol)
        self.thunder=tk.StringVar()

        self.tSpacer=tk.Frame(self.thunderChoice, background=bgcol, width=50)
        self.tSpacer.pack_propagate(0)
        self.tSpacer.pack(side=tk.LEFT, fill=tk.Y, expand=1)

        self.thunderType=tk.StringVar()
        self.trow1=tk.Frame(self.thunderChoice, background=bgcol)
        self.tcb=tk.Radiobutton(self.trow1, variable=self.thunderType, text="Coverage", value="Coverage", bg=bgcol,
                                     highlightthickness=0, height=1, command=self.changeTType)
        self.tcb.pack(side=tk.LEFT, anchor=tk.W)
        self.tpb=tk.Radiobutton(self.trow1, variable=self.thunderType, text="Probability", value="Probability", bg=bgcol,
                                     highlightthickness=0, height=1, command=self.changeTType)
        self.tpb.pack(side=tk.RIGHT, anchor=tk.E)
        self.trow1.pack(side=tk.TOP)

        self.trow2=tk.Frame(self.thunderChoice, background=bgcol)
        self.th1=tk.Radiobutton(self.trow2, variable=self.thunder, text="None", value="None", bg=bgcol, highlightthickness=0, height=1, command=self.changeThunder)
        self.th1.pack(side=tk.LEFT, anchor=tk.W)
        self.trow2.pack(side=tk.TOP, fill=tk.X, expand=tk.YES)

        self.trow3=tk.Frame(self.thunderChoice, background=bgcol)
        self.th2=tk.Radiobutton(self.trow3, variable=self.thunder, text="Keep Existing", value="Keep", bg=bgcol, highlightthickness=0, height=1, command=self.changeThunder)
        self.th2.pack(side=tk.LEFT, anchor=tk.W)
        self.trow3.pack(side=tk.TOP, fill=tk.X, expand=tk.YES)

        self.trow4=tk.Frame(self.thunderChoice, background=bgcol)
        self.th3=tk.Radiobutton(self.trow4, variable=self.thunder, text="Force:", value="Force", bg=bgcol, highlightthickness=0, height=1, command=self.changeThunder)
        self.th3.pack(side=tk.LEFT, anchor=tk.W)
        self.tforceLev=tk.IntVar()
        self.tf4=tk.Radiobutton(self.trow4, variable=self.tforceLev, text="Def", value=4, bg=bgcol, highlightthickness=0, height=1)
        self.tf4.pack(side=tk.RIGHT, anchor=tk.E)
        self.tf3=tk.Radiobutton(self.trow4, variable=self.tforceLev, text="Lkly", value=3, bg=bgcol, highlightthickness=0, height=1)
        self.tf3.pack(side=tk.RIGHT, anchor=tk.E)
        self.tf2=tk.Radiobutton(self.trow4, variable=self.tforceLev, text="Chc", value=2, bg=bgcol, highlightthickness=0, height=1)
        self.tf2.pack(side=tk.RIGHT, anchor=tk.E)
        self.tf1=tk.Radiobutton(self.trow4, variable=self.tforceLev, text="SChc", value=1, bg=bgcol, highlightthickness=0, height=1)
        self.tf1.pack(side=tk.RIGHT, anchor=tk.E)
        self.tf1.select()
        self.trow4.pack(side=tk.TOP, fill=tk.X, expand=tk.YES)

        self.trow5=tk.Frame(self.thunderChoice, background=bgcol)
        self.th4=tk.Radiobutton(self.trow5, variable=self.thunder, text="Match PoP:  Max of:", value="Match", bg=bgcol, highlightthickness=0, height=1, command=self.changeThunder)
        self.th4.pack(side=tk.LEFT, anchor=tk.W)
        self.tmatchLev=tk.IntVar()
        self.tm4=tk.Radiobutton(self.trow5, variable=self.tmatchLev, text="Def", value=4, bg=bgcol, highlightthickness=0, height=1)
        self.tm4.pack(side=tk.RIGHT, anchor=tk.E)
        self.tm3=tk.Radiobutton(self.trow5, variable=self.tmatchLev, text="Lkly", value=3, bg=bgcol, highlightthickness=0, height=1)
        self.tm3.pack(side=tk.RIGHT, anchor=tk.E)
        self.tm2=tk.Radiobutton(self.trow5, variable=self.tmatchLev, text="Chc", value=2, bg=bgcol, highlightthickness=0, height=1)
        self.tm2.pack(side=tk.RIGHT, anchor=tk.E)
        self.tm1=tk.Radiobutton(self.trow5, variable=self.tmatchLev, text="SChc", value=1, bg=bgcol, highlightthickness=0, height=1)
        self.tm1.pack(side=tk.RIGHT, anchor=tk.E)
        self.tm4.select()
        self.trow5.pack(side=tk.TOP, fill=tk.X, expand=tk.YES)

        self.trow6=tk.Frame(self.thunderChoice, background=bgcol)
        self.th5=tk.Radiobutton(self.trow6, variable=self.thunder, text="From PotThunder Grid", value="FromPot", bg=bgcol, highlightthickness=0, height=1, command=self.changeThunder)
        self.th5.pack(side=tk.LEFT, anchor=tk.W)
        self.trow6.pack(side=tk.TOP, fill=tk.X, expand=tk.YES)

        self.th2.invoke()

        self.thunderChoice.pack(side=tk.TOP, anchor=tk.W)
        if self.thunderState.get()==0:
            self.thunderChoice.pack_forget()

        thunderFrame.pack(side=tk.TOP, fill=tk.BOTH, expand=0)

        ###############################################################
        #  Freezing Rain
        fzFrame=tk.Frame(master, relief=tk.GROOVE, borderwidth=3, background=bgcol)
        fzOpen=tk.Frame(fzFrame, background=bgcol)
        self.fzState=tk.IntVar()
        self.fzState.set(0)
        self.fzButton=tk.Button(fzOpen, text=">", padx=0, pady=0, command=self.changeFzState)
        self.fzButton.pack(side=tk.LEFT, anchor=tk.W)
        self.fzLabel=tk.Label(fzOpen, text="Freezing Rain", anchor=tk.W, font=cfont, bg=bgcol)
        self.fzLabel.pack(side=tk.LEFT, anchor=tk.W, fill=tk.X, expand=1)
        fzOpen.pack(side=tk.TOP, anchor=tk.W, fill=tk.X, expand=1)

        self.freezing=tk.StringVar()
        self.fzChoice=tk.Frame(fzFrame, background=bgcol)

        self.fSpacer=tk.Frame(self.fzChoice, background=bgcol, width=50)
        self.fSpacer.pack_propagate(0)
        self.fSpacer.pack(side=tk.LEFT, anchor=tk.W, fill=tk.Y, expand=1)

        self.frow1=tk.Frame(self.fzChoice, background=bgcol)
        fz1=tk.Radiobutton(self.frow1, variable=self.freezing, text="None", value="None", bg=bgcol, highlightthickness=0, height=1, command=self.changeFreezing)
        fz1.pack(side=tk.LEFT, anchor=tk.W)
        self.frow1.pack(side=tk.TOP, anchor=tk.W, fill=tk.X, expand=1)

        self.frow2=tk.Frame(self.fzChoice, background=bgcol)
        fz2=tk.Radiobutton(self.frow2, variable=self.freezing, text="Keep Existing", value="Keep", bg=bgcol, highlightthickness=0, height=1, command=self.changeFreezing)
        fz2.pack(side=tk.LEFT, anchor=tk.W)
        self.frow2.pack(side=tk.TOP, anchor=tk.W, fill=tk.X, expand=1)

        self.frow3=tk.Frame(self.fzChoice, background=bgcol)
        fz3=tk.Radiobutton(self.frow3, variable=self.freezing, text="Force Only Freezing Rain", value="Force", bg=bgcol, highlightthickness=0, height=1, command=self.changeFreezing)
        fz3.pack(side=tk.LEFT, anchor=tk.W)
        self.frow3.pack(side=tk.TOP, anchor=tk.W, fill=tk.X, expand=1)

        self.frow4=tk.Frame(self.fzChoice, background=bgcol)
        fz4=tk.Radiobutton(self.frow4, variable=self.freezing, text="Mix with Rain/Snow:", value="Mix", bg=bgcol, highlightthickness=0, height=1, command=self.changeFreezing)
        fz4.pack(side=tk.LEFT, anchor=tk.W)
        self.frow4.pack(side=tk.TOP, anchor=tk.W, fill=tk.X, expand=1)

        self.fmSpacer=tk.Frame(self.fzChoice, background=bgcol, width=30)
        self.fmSpacer.pack_propagate(0)
        self.fmSpacer.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

        self.frow5=tk.Frame(self.fzChoice, background=bgcol)
        self.fmmatchLev=tk.IntVar()
        self.fmm4=tk.Radiobutton(self.frow5, variable=self.fmmatchLev, text="Def", value=4, bg=bgcol, highlightthickness=0, height=1)
        self.fmm4.pack(side=tk.RIGHT, anchor=tk.E)
        self.fmm3=tk.Radiobutton(self.frow5, variable=self.fmmatchLev, text="Lkly", value=3, bg=bgcol, highlightthickness=0, height=1)
        self.fmm3.pack(side=tk.RIGHT, anchor=tk.E)
        self.fmm2=tk.Radiobutton(self.frow5, variable=self.fmmatchLev, text="Chc", value=2, bg=bgcol, highlightthickness=0, height=1)
        self.fmm2.pack(side=tk.RIGHT, anchor=tk.E)
        self.fmm1=tk.Radiobutton(self.frow5, variable=self.fmmatchLev, text="SChc", value=1, bg=bgcol, highlightthickness=0, height=1)
        self.fmm1.pack(side=tk.RIGHT, anchor=tk.E)
        self.fm2=tk.Label(self.frow5, text="             Max of:", bg=bgcol)
        self.fm2.pack(side=tk.LEFT, anchor=tk.W, fill=tk.X, expand=1)
        self.fmm4.select()
        self.frow5.pack(side=tk.TOP, anchor=tk.E, fill=tk.X, expand=1)

        self.fzChoice.pack(side=tk.TOP, anchor=tk.N, fill=tk.X, expand=1)
        if self.fzState.get()==0:
            self.fzChoice.pack_forget()
        fzFrame.pack(side=tk.TOP, anchor=tk.N, fill=tk.X, expand=1)

        ###############################################################
        #  Sleet
        ipFrame=tk.Frame(master, relief=tk.GROOVE, borderwidth=3, background=bgcol)
        ipOpen=tk.Frame(ipFrame, background=bgcol)
        self.ipState=tk.IntVar()
        self.ipState.set(0)
        self.ipButton=tk.Button(ipOpen, text=">", padx=0, pady=0, command=self.changeIpState)
        self.ipButton.pack(side=tk.LEFT, anchor=tk.W)
        self.ipLabel=tk.Label(ipOpen, text="Sleet", anchor=tk.W, font=cfont, bg=bgcol)
        self.ipLabel.pack(side=tk.LEFT, anchor=tk.W, fill=tk.X, expand=1)
        ipOpen.pack(side=tk.TOP, anchor=tk.W, fill=tk.X, expand=1)

        self.sleet=tk.StringVar()
        self.ipChoice=tk.Frame(ipFrame, background=bgcol)

        self.iSpacer=tk.Frame(self.ipChoice, background=bgcol, width=50)
        self.iSpacer.pack_propagate(0)
        self.iSpacer.pack(side=tk.LEFT, anchor=tk.W, fill=tk.Y, expand=1)

        self.irow1=tk.Frame(self.ipChoice, background=bgcol)
        ip1=tk.Radiobutton(self.irow1, variable=self.sleet, text="None", value="None", bg=bgcol, highlightthickness=0, height=1, command=self.changeSleet)
        ip1.pack(side=tk.LEFT, anchor=tk.W)
        self.irow1.pack(side=tk.TOP, anchor=tk.W, fill=tk.X, expand=1)

        self.irow2=tk.Frame(self.ipChoice, background=bgcol)
        ip2=tk.Radiobutton(self.irow2, variable=self.sleet, text="Keep Existing", value="Keep", bg=bgcol, highlightthickness=0, height=1, command=self.changeSleet)
        ip2.pack(side=tk.LEFT, anchor=tk.W)
        self.irow2.pack(side=tk.TOP, anchor=tk.W, fill=tk.X, expand=1)

        self.irow3=tk.Frame(self.ipChoice, background=bgcol)
        ip3=tk.Radiobutton(self.irow3, variable=self.sleet, text="Force Only Sleet", value="Force", bg=bgcol, highlightthickness=0, height=1, command=self.changeSleet)
        ip3.pack(side=tk.LEFT, anchor=tk.W)
        self.irow3.pack(side=tk.TOP, anchor=tk.W, fill=tk.X, expand=1)

        self.irow4=tk.Frame(self.ipChoice, background=bgcol)
        ip4=tk.Radiobutton(self.irow4, variable=self.sleet, text="Mix with Rain/Snow:", value="Mix", bg=bgcol, highlightthickness=0, height=1, command=self.changeSleet)
        ip4.pack(side=tk.LEFT, anchor=tk.W)
        self.irow4.pack(side=tk.TOP, anchor=tk.W, fill=tk.X, expand=1)

        self.imSpacer=tk.Frame(self.ipChoice, background=bgcol, width=30)
        self.imSpacer.pack_propagate(0)
        self.imSpacer.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

        self.irow5=tk.Frame(self.ipChoice, background=bgcol)
        self.immatchLev=tk.IntVar()
        self.imm4=tk.Radiobutton(self.irow5, variable=self.immatchLev, text="Def", value=4, bg=bgcol, highlightthickness=0, height=1)
        self.imm4.pack(side=tk.RIGHT, anchor=tk.E)
        self.imm3=tk.Radiobutton(self.irow5, variable=self.immatchLev, text="Lkly", value=3, bg=bgcol, highlightthickness=0, height=1)
        self.imm3.pack(side=tk.RIGHT, anchor=tk.E)
        self.imm2=tk.Radiobutton(self.irow5, variable=self.immatchLev, text="Chc", value=2, bg=bgcol, highlightthickness=0, height=1)
        self.imm2.pack(side=tk.RIGHT, anchor=tk.E)
        self.imm1=tk.Radiobutton(self.irow5, variable=self.immatchLev, text="SChc", value=1, bg=bgcol, highlightthickness=0, height=1)
        self.imm1.pack(side=tk.RIGHT, anchor=tk.E)
        self.im2=tk.Label(self.irow5, text="             Max of:", bg=bgcol)
        self.im2.pack(side=tk.LEFT, anchor=tk.W, fill=tk.X, expand=1)
        self.imm4.select()
        self.irow5.pack(side=tk.TOP, anchor=tk.E, fill=tk.X, expand=1)

        self.ipChoice.pack(side=tk.TOP, anchor=tk.N, fill=tk.X, expand=1)
        if self.ipState.get()==0:
            self.ipChoice.pack_forget()
        ipFrame.pack(side=tk.TOP, anchor=tk.N, fill=tk.X, expand=1)

        self.prbbut.invoke()
        fz1.invoke()
        mainFrame.pack(side=tk.TOP, fill=tk.BOTH, expand=tk.YES)
        return mainFrame

    def topMove(self, event):
        top=self.topRainSnow.get()
        bot=self.botRainSnow.get()
        if top<0:
            self.topRainSnow.set(0)
        elif top<bot:
            self.botRainSnow.set(top)
        elif (top-bot) > maxThicknessRainSnowLayer:
            self.botRainSnow.set(top-maxThicknessRainSnowLayer)
        self.drawSounding()
        return

    def bottomMove(self, event):
        top=self.topRainSnow.get()
        bot=self.botRainSnow.get()
        if bot>0:
            self.botRainSnow.set(0)
        elif top<bot:
            self.topRainSnow.set(bot)
        elif (top-bot) > maxThicknessRainSnowLayer:
            self.topRainSnow.set(bot+maxThicknessRainSnowLayer)
        self.drawSounding()
        return

    def changeCoverage(self):
        cov=self.coverage.get()
        if cov=="Coverage":
            self.tcb.invoke()
        else:
            self.tpb.invoke()
        return

    def changeConvective(self):
        convective=self.convective.get()
        if convective=="Stratiform":
            self.coverage.set("Probability")
            self.covbut.configure(state=tk.DISABLED)
            self.tcb.configure(state=tk.DISABLED)
            self.tpb.configure(state=tk.DISABLED)
            self.th1.configure(state=tk.DISABLED)
            self.th2.configure(state=tk.DISABLED)
            self.th3.configure(state=tk.DISABLED)
            self.th4.configure(state=tk.DISABLED)
            self.th5.configure(state=tk.DISABLED)
            self.tf1.configure(state=tk.DISABLED)
            self.tf2.configure(state=tk.DISABLED)
            self.tf3.configure(state=tk.DISABLED)
            self.tf4.configure(state=tk.DISABLED)
            self.tm1.configure(state=tk.DISABLED)
            self.tm2.configure(state=tk.DISABLED)
            self.tm3.configure(state=tk.DISABLED)
            self.tm4.configure(state=tk.DISABLED)
        else:
            self.covbut.configure(state=tk.NORMAL)
            self.tcb.configure(state=tk.NORMAL)
            self.tpb.configure(state=tk.NORMAL)
            self.th1.configure(state=tk.NORMAL)
            self.th2.configure(state=tk.NORMAL)
            self.th3.configure(state=tk.NORMAL)
            self.th4.configure(state=tk.NORMAL)
            self.th5.configure(state=tk.NORMAL)

            thunderValue=self.thunder.get()
            if thunderValue=="Force":
                self.tf1.configure(state=tk.NORMAL)
                self.tf2.configure(state=tk.NORMAL)
                self.tf3.configure(state=tk.NORMAL)
                self.tf4.configure(state=tk.NORMAL)
                self.tm1.configure(state=tk.DISABLED)
                self.tm2.configure(state=tk.DISABLED)
                self.tm3.configure(state=tk.DISABLED)
                self.tm4.configure(state=tk.DISABLED)
            elif thunderValue=="Match":
                self.tm1.configure(state=tk.NORMAL)
                self.tm2.configure(state=tk.NORMAL)
                self.tm3.configure(state=tk.NORMAL)
                self.tm4.configure(state=tk.NORMAL)
                self.tf1.configure(state=tk.DISABLED)
                self.tf2.configure(state=tk.DISABLED)
                self.tf3.configure(state=tk.DISABLED)
                self.tf4.configure(state=tk.DISABLED)
            else:
                self.tf1.configure(state=tk.DISABLED)
                self.tf2.configure(state=tk.DISABLED)
                self.tf3.configure(state=tk.DISABLED)
                self.tf4.configure(state=tk.DISABLED)
                self.tm1.configure(state=tk.DISABLED)
                self.tm2.configure(state=tk.DISABLED)
                self.tm3.configure(state=tk.DISABLED)
                self.tm4.configure(state=tk.DISABLED)
        return

    def changeMixState(self):
        if self.mixState.get()==1:
            self.mixChoice.pack_forget()
            self.mixButton.configure(text=">")
            self.mixState.set(0)
        else:
            self.mixChoice.pack(side=tk.TOP, anchor=tk.W, fill=tk.BOTH, expand=tk.YES)
            self.mixButton.configure(text="V")
            self.mixState.set(1)
        self.updateSize()
        return

    def changeThunderState(self):
        if self.thunderState.get()==1:
            self.thunderChoice.pack_forget()
            self.thunderButton.configure(text=">")
            self.thunderState.set(0)
        else:
            self.thunderChoice.pack(side=tk.TOP, anchor=tk.W)
            self.thunderButton.configure(text="V")
            self.thunderState.set(1)
        self.updateSize()
        return

    def changeFzState(self):
        if self.fzState.get()==1:
            self.fzChoice.pack_forget()
            self.fzButton.configure(text=">")
            self.fzState.set(0)
        else:
            self.fzChoice.pack(side=tk.TOP, anchor=tk.W)
            self.fzButton.configure(text="V")
            self.fzState.set(1)
        self.updateSize()
        return

    def changeIpState(self):
        if self.ipState.get()==1:
            self.ipChoice.pack_forget()
            self.ipButton.configure(text=">")
            self.ipState.set(0)
        else:
            self.ipChoice.pack(side=tk.TOP, anchor=tk.W)
            self.ipButton.configure(text="V")
            self.ipState.set(1)
        self.updateSize()
        return

    def updateSize(self):
        geom=self.geometry()
        mobj=REGEO.match(geom)
        if (mobj):
            ofx1=int(mobj.group(4))
            ofy1=int(mobj.group(6))
        else:
            ofx1=0
            ofy1=0

        self.update_idletasks()
        wh=self.winfo_reqheight()
        ww=self.winfo_reqwidth()
        newgeom="%dx%d+%d+%d"%(ww, wh, ofx1, ofy1)
        self.geometry(newgeom)
        return

    def changeThunder(self):
        bgcol = "#B4B4B4"
        state=self.thunder.get()
        if state=="Match":
            self.tm1.configure(state=tk.NORMAL)
            self.tm2.configure(state=tk.NORMAL)
            self.tm3.configure(state=tk.NORMAL)
            self.tm4.configure(state=tk.NORMAL)
            self.thunderOpen.configure(background="#D24040")
            self.thunderLabel.configure(background="#D24040")
        else:
            self.tm1.configure(state=tk.DISABLED)
            self.tm2.configure(state=tk.DISABLED)
            self.tm3.configure(state=tk.DISABLED)
            self.tm4.configure(state=tk.DISABLED)

        if state=="Force":
            self.tf1.configure(state=tk.NORMAL)
            self.tf2.configure(state=tk.NORMAL)
            self.tf3.configure(state=tk.NORMAL)
            self.tf4.configure(state=tk.NORMAL)
            self.thunderOpen.configure(background="#D24040")
            self.thunderLabel.configure(background="#D24040")
        else:
            self.tf1.configure(state=tk.DISABLED)
            self.tf2.configure(state=tk.DISABLED)
            self.tf3.configure(state=tk.DISABLED)
            self.tf4.configure(state=tk.DISABLED)

        if state=="FromPot":
            self.thunderOpen.configure(background="#D24040")
            self.thunderLabel.configure(background="#D24040")

        if state=="Keep":
            self.thunderOpen.configure(background="#E1E100")
            self.thunderLabel.configure(background="#E1E100")

        if state=="None":
            self.thunderOpen.configure(background=bgcol)
            self.thunderLabel.configure(background=bgcol)

        self.update_idletasks()
        return

    def changeTType(self):
        state=self.thunderType.get()
        if state=="Probability":
            self.tm1.configure(text="SChc")
            self.tm2.configure(text="Chc")
            self.tm3.configure(text="Lkly")
            self.tm4.configure(text="Def")
            self.tf1.configure(text="SChc")
            self.tf2.configure(text="Chc")
            self.tf3.configure(text="Lkly")
            self.tf4.configure(text="Def")
        else:
            self.tm1.configure(text="Iso")
            self.tm2.configure(text="Sct")
            self.tm3.configure(text="Num")
            self.tm4.configure(text="Wide")
            self.tf1.configure(text="Iso")
            self.tf2.configure(text="Sct")
            self.tf3.configure(text="Num")
            self.tf4.configure(text="Wide")

    def changeFreezing(self):
        bgcol = "#B4B4B4"
        state=self.freezing.get()
        if state=="Mix":
            self.fm2.configure(state=tk.NORMAL)
            self.fmm1.configure(state=tk.NORMAL)
            self.fmm2.configure(state=tk.NORMAL)
            self.fmm3.configure(state=tk.NORMAL)
            self.fmm4.configure(state=tk.NORMAL)
            self.fzLabel.configure(background="#D24040")
        else:
            self.fm2.configure(state=tk.DISABLED)
            self.fmm1.configure(state=tk.DISABLED)
            self.fmm2.configure(state=tk.DISABLED)
            self.fmm3.configure(state=tk.DISABLED)
            self.fmm4.configure(state=tk.DISABLED)
            if state=="Force":
                self.fzLabel.configure(background="#D24040")
            elif state=="Keep":
                self.fzLabel.configure(background="#E1E100")
            else:
                self.fzLabel.configure(background=bgcol)

        return

    def changeSleet(self):
        bgcol = "#B4B4B4"
        state=self.sleet.get()
        if state=="Mix":
            self.im2.configure(state=tk.NORMAL)
            self.imm1.configure(state=tk.NORMAL)
            self.imm2.configure(state=tk.NORMAL)
            self.imm3.configure(state=tk.NORMAL)
            self.imm4.configure(state=tk.NORMAL)
            self.ipLabel.configure(background="#D24040")
        else:
            self.im2.configure(state=tk.DISABLED)
            self.imm1.configure(state=tk.DISABLED)
            self.imm2.configure(state=tk.DISABLED)
            self.imm3.configure(state=tk.DISABLED)
            self.imm4.configure(state=tk.DISABLED)
            if state=="Force":
                self.ipLabel.configure(background="#D24040")
            elif state=="Keep":
                self.ipLabel.configure(background="#E1E100")
            else:
                self.ipLabel.configure(background=bgcol)

        return

    def drawSounding(self):
        top=self.topRainSnow.get()
        bot=self.botRainSnow.get()

        topLevel = (-1.0*top/(maxThicknessRainSnowLayer*2.0))*100 + 50
        botLevel = (-1.0*bot/(maxThicknessRainSnowLayer*2.0))*100 + 50

        self.soundingObject.delete(tk.ALL)
        self.soundingObject.create_rectangle(0, 0, 100, topLevel, fill="cyan")
        self.soundingObject.create_rectangle(0, topLevel, 100, botLevel, fill="#FF00FF")
        self.soundingObject.create_rectangle(0, botLevel, 100, 100, fill="green")
        self.soundingObject.create_line(40, 0, 50, 50, fill="red", width=2.0)
        self.soundingObject.create_line(50, 50, 70, 100, fill="red", width=2.0)
        self.soundingObject.create_line(35, 0, 45, 50, fill="darkgreen", width=2.0)
        self.soundingObject.create_line(45, 50, 50, 100, fill="darkgreen", width=2.0)
        self.soundingObject.create_line(0, 50, 100, 50, fill="black", width=4.0)
        self.soundingObject.create_text(20, 45, fill="black", text="SnowLevel", font="arial")
