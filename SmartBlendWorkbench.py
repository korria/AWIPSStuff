
##
# ----------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# SmartBlendWorkbench - version 1.0
#
# Operational model blending workbench for GFE
#
# Major capabilities:
#   * Multi-model weighted blending for scalar and vector grids
#   * Preview diagnostics before apply
#   * Blend vs NBM / Blend vs Fcst awareness
#   * "Preserve my edits" mask based on Fcst departures from NBM
#   * Apply only high-confidence area
#   * Apply only where blend differs from Fcst by threshold
#   * BOIVerify recent-skill suggested weights using URMA by default
#   * Presets + saved custom presets
#   * Persisted last-used settings per weather element
#
# Notes:
#   * This is designed as an office-ready first version. It is intentionally
#     conservative and transparent rather than overly automatic.
#   * The preview summarizes what the blend would do over the currently selected
#     mutable grids in the selected time range.
#   * The "Suggest Skill" button can take some time since it searches recent
#     BOIVerify cases.
# ----------------------------------------------------------------------------
##

ToolType = "numeric"
WeatherElementEdited = "variableElement"
ScreenList = ["SCALAR", "VECTOR"]

# ----------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------

PERSIST_CATEGORY = "SmartBlendWorkbench"

# model|versions|group|color|parm list
MODEL_CONFIG = [
    "Fcst|1|BASE|#000000|ALL",
    "Official|1|BASE|#555555|ALL",
    "NBM|1|NBM|#3030ff|ALL",
    "NBMEXP|1|NBM|#5050ff|ALL",
    "NBM10Pct|1|NBM|#8080ff|MaxT,MinT,T,Td,Sky,Wind,WindGust",
    "NBM25Pct|1|NBM|#6a6aff|MaxT,MinT,T,Td,Sky,Wind,WindGust",
    "NBM50Pct|1|NBM|#5252ff|MaxT,MinT,T,Td,Sky,Wind,WindGust",
    "NBM75Pct|1|NBM|#3a3aff|MaxT,MinT,T,Td,Sky,Wind,WindGust",
    "NBM90Pct|1|NBM|#2222ff|MaxT,MinT,T,Td,Sky,Wind,WindGust",
    "SuperBlend|1|CONS|#7f007f|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Sky,QPF,SnowLevel,SnowRatio,Wind,WindGust",
    "CONSAll|1|CONS|#666666|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,QPF",
    "BCCONSAll|1|CONS|#999999|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH",
    "CONSShort|1|CONS|#808080|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Sky,QPF,SnowLevel,SnowRatio,Wind,WindGust",
    "BCCONSShort|1|CONS|#b0b0b0|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Wind,WindGust",
    "WPCGuide|2|WPC|#3f7fff|ALL",
    "GFS|2|GFS|#7f0000|ALL",
    "GFSBC|2|GFS|#ff0000|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Wind",
    "ADJMAV|2|GFS|#bfbf00|ALL",
    "ADJMAVBC|2|GFS|#ffff00|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Wind",
    "ADJMEX|2|GFS|#7f3f00|ALL",
    "ADJMEXBC|2|GFS|#ff7f00|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Wind",
    "ECMWF|2|EC|#00007f|ALL",
    "ECMWFBC|2|EC|#0000ff|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Wind",
    "ECAIFS|1|EC|#11117f|ALL",
    "ECAIFSBC|1|EC|#1111ff|ALL",
    "ADJECE|2|EC|#3f007f|ALL",
    "ADJECEBC|2|EC|#7f007f|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Wind",
    "ADJECS|1|EC|#bf00bf|PoP",
    "ADJECM|2|EC|#7f7fff|MaxT,MinT,PoP",
    "ADJECMBC|2|EC|#9999ff|MaxT,MinT,PoP",
    "CMCnh|2|CAN|#7fbfff|ALL",
    "CMCnhBC|2|CAN|#3f7fff|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Wind",
    "CMCreg|2|CAN|#4fa0ff|ALL",
    "CMCregBC|2|CAN|#2f90ff|ALL",
    "HRDPS|2|CAN|#0f7fff|ALL",
    "HRDPSBC|2|CAN|#077fff|ALL",
    "NAM12|2|NAM|#007f00|ALL",
    "NAM12BC|2|NAM|#00ff00|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Wind",
    "ADJMET|2|NAM|#007f7f|ALL",
    "ADJMETBC|2|NAM|#00bfbf|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Wind",
    "SREF|1|NAM|#003f3f|ALL",
    "SREFBC|1|NAM|#007f7f|MaxT,MinT,MaxRH,MinRH,TdMrn,TdAft,T,Td,RH,PoP,Wind",
    "NAMNest|2|NAM|#77773f|ALL",
    "NAMNestBC|2|NAM|#66663f|ALL",
    "RAP13|2|RAP|#aa5500|MaxT,MinT,QPF,Wind",
    "HRRR|3|HRRR|#ff5500|ALL",
    "HRRRBC|3|HRRR|#ff8800|ALL",
    "HREF|1|HRRR|#ff00aa|ALL",
    "HREFBC|1|HRRR|#ff55cc|ALL",
    "MPASRT|2|MPA|#2f2f2f|ALL",
    "MPASRTBC|2|MPA|#3f3f3f|ALL",
    "LRRFS|1|MPA|#77ff3f|ALL",
    "LRRFSBC|1|MPA|#66ff3f|ALL",
    "UWWRF1|2|WRF|#5f5f5f|ALL",
    "UWWRF1BC|2|WRF|#4f4f4f|ALL",
    "HIRESWarw|2|HIRES|#1f7f3f|MaxT,MinT,MixHgt,MixHgtMax,MixHgtMin,TransWind,Wind,Haines,Sky,QPF,SnowLevel",
    "HIRESWarwBC|1|HIRES|#3fff7f|MaxT,MinT",
    "HIRESWFV3|2|HIRES|#2f9f5f|MaxT,MinT,MixHgt,MixHgtMax,MixHgtMin,TransWind,Wind,Haines,Sky,QPF,SnowLevel",
    "HIRESWFV3BC|1|HIRES|#4fff8f|MaxT,MinT",
]

PARM_ALTERNATES = {
    "QPF": ["QPFRaw", "QPF50Prcntl"],
    "SnowLevel": ["FzLevel"],
    "MixHgt": ["MixHgtMax", "MixHgtMin"],
    "PoP": ["PoP6", "PPI06", "ProbCDBZ20", "ProbCDBZ30", "ProbCDBZ40", "ProbCDBZ50"],
}

DEFAULT_PRESETS = {
    "DEFAULT": [("Fcst", 2), ("NBM", 5)],
    "MaxT": [("Fcst", 2), ("NBM", 5), ("ECMWF", 3), ("GFS", 2), ("ADJMAVBC", 2), ("ADJECEBC", 2)],
    "MinT": [("Fcst", 2), ("NBM", 5), ("ECMWF", 3), ("GFS", 2), ("ADJMAVBC", 2), ("ADJECEBC", 2)],
    "T": [("Fcst", 2), ("NBM", 4), ("HRRR", 3), ("RAP13", 2), ("ECMWF", 2)],
    "Td": [("Fcst", 2), ("NBM", 4), ("HRRR", 3), ("RAP13", 2), ("ECMWF", 2)],
    "RH": [("Fcst", 2), ("NBM", 4), ("HRRR", 3), ("RAP13", 2), ("ECMWF", 2)],
    "MaxRH": [("Fcst", 2), ("NBM", 4), ("HRRR", 3), ("RAP13", 2), ("ECMWF", 2)],
    "MinRH": [("Fcst", 2), ("NBM", 4), ("HRRR", 3), ("RAP13", 2), ("ECMWF", 2)],
    "PoP": [("Fcst", 2), ("NBM", 5), ("SuperBlend", 3), ("HREF", 2), ("ECMWF", 2)],
    "QPF": [("Fcst", 2), ("NBM", 4), ("HREF", 3), ("HRRR", 3), ("ECMWF", 2), ("WPCGuide", 2)],
    "Sky": [("Fcst", 2), ("NBM", 5), ("HRRR", 3), ("RAP13", 2), ("ECMWF", 2)],
    "Wind": [("Fcst", 2), ("NBM", 4), ("HRRR", 4), ("RAP13", 3), ("NAMNest", 2), ("HREF", 2)],
    "WindGust": [("Fcst", 2), ("NBM", 4), ("HRRR", 4), ("RAP13", 3), ("HREF", 2)],
    "SnowLevel": [("Fcst", 2), ("NBM", 4), ("ECMWF", 3), ("GFS", 2), ("HREF", 2)],
}

EDGE_STYLES = ["Flat", "Edge", "Taper"]
VECTOR_MODES = ["Both", "Magnitude Only", "Direction Only"]
SCALAR_MODES = ["TimeWtAverage", "Max"]

DEFAULT_EDGE_STYLE = "Taper"
DEFAULT_EDGE_WIDTH = 5
EDGE_WIDTH_MAX = 30
DEFAULT_VECTOR_MODE = "Both"
DEFAULT_SCALAR_MODE = "TimeWtAverage"
USE_NEGATIVE_WEIGHTS = 0

DEFAULT_NORMALIZE = 1
DEFAULT_PRESERVE_EDITS = 1
DEFAULT_HIGH_CONF_ONLY = 0
DEFAULT_CHANGE_ONLY = 0
DEFAULT_CHANGE_THRESHOLD = 2.0
DEFAULT_PRESERVE_THRESHOLD = 3.0
DEFAULT_PRE_FLAG = 0
DEFAULT_FIRST_FLAG = 1
DEFAULT_LAST_FLAG = 1
DEFAULT_POST_FLAG = 0

CHECK_SMALL_EDIT_AREA = True
MIN_EDIT_AREA_POINTS = 20

SKILL_LOOKBACK_DAYS = 30
SKILL_MAX_CASES_PER_MODEL = 10
SKILL_OBS_SOURCE = "URMA"

OBS_ELEMENT_MAP = {
    "MaxT": "T",
    "MinT": "T",
    "MaxRH": "RH",
    "MinRH": "RH",
    "TdMrn": "Td",
    "TdAft": "Td",
    "QPF": "QPE06",
}

ELEMENT_THRESHOLDS = {
    "DEFAULT": {"depart": 2.0, "spread_good": 2.0, "spread_warn": 5.0, "change": 3.0},
    "MaxT": {"depart": 2.0, "spread_good": 2.0, "spread_warn": 5.0, "change": 4.0},
    "MinT": {"depart": 2.0, "spread_good": 2.0, "spread_warn": 5.0, "change": 4.0},
    "T": {"depart": 2.0, "spread_good": 2.0, "spread_warn": 4.0, "change": 3.0},
    "Td": {"depart": 2.0, "spread_good": 2.0, "spread_warn": 4.0, "change": 3.0},
    "RH": {"depart": 5.0, "spread_good": 5.0, "spread_warn": 12.0, "change": 8.0},
    "MaxRH": {"depart": 5.0, "spread_good": 5.0, "spread_warn": 12.0, "change": 8.0},
    "MinRH": {"depart": 5.0, "spread_good": 5.0, "spread_warn": 12.0, "change": 8.0},
    "PoP": {"depart": 10.0, "spread_good": 10.0, "spread_warn": 25.0, "change": 15.0},
    "QPF": {"depart": 0.05, "spread_good": 0.05, "spread_warn": 0.20, "change": 0.10},
    "Sky": {"depart": 10.0, "spread_good": 10.0, "spread_warn": 25.0, "change": 15.0},
    "Wind": {"depart": 3.0, "spread_good": 3.0, "spread_warn": 8.0, "change": 5.0},
    "WindGust": {"depart": 4.0, "spread_good": 4.0, "spread_warn": 10.0, "change": 6.0},
    "SnowLevel": {"depart": 200.0, "spread_good": 200.0, "spread_warn": 500.0, "change": 300.0},
}

# ----------------------------------------------------------------------------
# IMPORTS
# ----------------------------------------------------------------------------

import numpy as np

import AppDialog
import SmartScript
import TimeRange
import AbsTime
import BOIVerifyUtility

try:
    import tkinter as tk
except ImportError:
    import Tkinter as tk

try:
    import tkinter.messagebox as msg
except ImportError:
    import tkMessageBox as msg

try:
    import tkinter.simpledialog as simpledialog
except ImportError:
    import tkSimpleDialog as simpledialog

# ----------------------------------------------------------------------------
# DIALOG
# ----------------------------------------------------------------------------

class ToolDialog(AppDialog.AppDialog):
    def __init__(self, title="SmartBlend Workbench", callbackMethod=None,
                 sources=None, presetNames=None, initialPreset="Default", **kwargs):
        self.__callbackMethod = callbackMethod
        self.sources = sources if sources is not None else []
        self.presetNames = presetNames if presetNames is not None else ["Default"]

        self.weightVars = []
        self.percentVars = []
        self.statusVars = []
        self.rowFrames = []
        self.nameLabels = []
        self.scaleWidgets = []
        self.groupVars = {}

        self.summaryVars = {}

        # _buildVariables MUST run before AppDialog.__init__ because
        # the parent __init__ calls body() -> _buildTopControls() which
        # reads self.presetVar etc.  We pass no master to the Tk vars
        # because self is not yet a valid Tk widget.
        self._buildVariables(initialPreset)

        AppDialog.AppDialog.__init__(self, **kwargs)
        self.title(title)
        self.wm_attributes("-topmost", 1)
        self.protocol("WM_DELETE_WINDOW", self.closeCB)

    def _buildVariables(self, initialPreset):
        self.presetVar = tk.StringVar()
        self.presetVar.set(initialPreset)

        self.edgeStyleVar = tk.StringVar()
        self.edgeStyleVar.set(DEFAULT_EDGE_STYLE)

        self.edgeWidthVar = tk.IntVar()
        self.edgeWidthVar.set(DEFAULT_EDGE_WIDTH)

        self.vectorModeVar = tk.StringVar()
        self.vectorModeVar.set(DEFAULT_VECTOR_MODE)

        self.scalarModeVar = tk.StringVar()
        self.scalarModeVar.set(DEFAULT_SCALAR_MODE)

        self.normalizeVar = tk.IntVar()
        self.normalizeVar.set(DEFAULT_NORMALIZE)

        self.preserveEditsVar = tk.IntVar()
        self.preserveEditsVar.set(DEFAULT_PRESERVE_EDITS)

        self.highConfOnlyVar = tk.IntVar()
        self.highConfOnlyVar.set(DEFAULT_HIGH_CONF_ONLY)

        self.changeOnlyVar = tk.IntVar()
        self.changeOnlyVar.set(DEFAULT_CHANGE_ONLY)

        self.changeThreshVar = tk.DoubleVar()
        self.changeThreshVar.set(DEFAULT_CHANGE_THRESHOLD)

        self.preserveThreshVar = tk.DoubleVar()
        self.preserveThreshVar.set(DEFAULT_PRESERVE_THRESHOLD)

        self.preFlagVar = tk.IntVar()
        self.preFlagVar.set(DEFAULT_PRE_FLAG)

        self.firstFlagVar = tk.IntVar()
        self.firstFlagVar.set(DEFAULT_FIRST_FLAG)

        self.lastFlagVar = tk.IntVar()
        self.lastFlagVar.set(DEFAULT_LAST_FLAG)

        self.postFlagVar = tk.IntVar()
        self.postFlagVar.set(DEFAULT_POST_FLAG)

        for key in ["availability", "comparison", "confidence", "verdict"]:
            self.summaryVars[key] = tk.StringVar()
            self.summaryVars[key].set("")

    def buttonbox(self):
        buttonFrame = tk.Frame(self)

        tk.Button(buttonFrame, text="Preview", width=10, command=self.previewCB).pack(side=tk.LEFT, padx=4, pady=5)
        tk.Button(buttonFrame, text="Apply", width=10, command=self.applyCB).pack(side=tk.LEFT, padx=4, pady=5)
        tk.Button(buttonFrame, text="Apply/Close", width=12, command=self.applyCloseCB).pack(side=tk.LEFT, padx=4, pady=5)
        tk.Button(buttonFrame, text="Suggest Skill", width=12, command=self.skillCB).pack(side=tk.LEFT, padx=4, pady=5)
        tk.Button(buttonFrame, text="Load Preset", width=11, command=self.presetCB).pack(side=tk.LEFT, padx=4, pady=5)
        tk.Button(buttonFrame, text="Save Preset", width=11, command=self.savePresetCB).pack(side=tk.LEFT, padx=4, pady=5)
        tk.Button(buttonFrame, text="Zero All", width=10, command=self.zeroAllCB).pack(side=tk.LEFT, padx=4, pady=5)
        tk.Button(buttonFrame, text="Close", width=10, command=self.closeCB).pack(side=tk.LEFT, padx=4, pady=5)

        buttonFrame.pack(side=tk.BOTTOM, fill=tk.X)

    def body(self, master):
        outer = tk.Frame(master)

        self._buildTopControls(outer)
        self._buildGroupControls(outer)
        self._buildSourcePane(outer)
        self._buildSummaryPane(outer)

        outer.pack(side=tk.TOP, fill=tk.BOTH, expand=tk.TRUE)
        return outer

    def _buildTopControls(self, master):
        top = tk.Frame(master, relief=tk.GROOVE, borderwidth=2)

        row1 = tk.Frame(top)
        tk.Label(row1, text="Preset:").pack(side=tk.LEFT, padx=4)
        self.presetMenu = tk.OptionMenu(row1, self.presetVar, *self.presetNames)
        self.presetMenu.pack(side=tk.LEFT, padx=4)

        tk.Label(row1, text="Edge:").pack(side=tk.LEFT, padx=(12, 4))
        tk.OptionMenu(row1, self.edgeStyleVar, *EDGE_STYLES).pack(side=tk.LEFT, padx=4)

        tk.Label(row1, text="Edge Width:").pack(side=tk.LEFT, padx=(12, 4))
        tk.Scale(row1, from_=1, to=EDGE_WIDTH_MAX, orient=tk.HORIZONTAL,
                 variable=self.edgeWidthVar, length=120).pack(side=tk.LEFT, padx=4)

        tk.Label(row1, text="Scalar Calc:").pack(side=tk.LEFT, padx=(12, 4))
        tk.OptionMenu(row1, self.scalarModeVar, *SCALAR_MODES).pack(side=tk.LEFT, padx=4)

        tk.Label(row1, text="Vector:").pack(side=tk.LEFT, padx=(12, 4))
        tk.OptionMenu(row1, self.vectorModeVar, *VECTOR_MODES).pack(side=tk.LEFT, padx=4)
        row1.pack(side=tk.TOP, fill=tk.X, pady=2)

        row2 = tk.Frame(top)
        tk.Checkbutton(row2, text="Normalize weights", variable=self.normalizeVar,
                       command=self.setPercents).pack(side=tk.LEFT, padx=6)
        tk.Checkbutton(row2, text="Preserve prior Fcst edits vs NBM",
                       variable=self.preserveEditsVar).pack(side=tk.LEFT, padx=6)
        tk.Checkbutton(row2, text="Apply only high-confidence area",
                       variable=self.highConfOnlyVar).pack(side=tk.LEFT, padx=6)
        tk.Checkbutton(row2, text="Apply only where |Blend-Fcst| exceeds threshold",
                       variable=self.changeOnlyVar).pack(side=tk.LEFT, padx=6)

        tk.Label(row2, text="Change Thresh:").pack(side=tk.LEFT, padx=(12, 4))
        tk.Entry(row2, textvariable=self.changeThreshVar, width=6).pack(side=tk.LEFT, padx=2)

        tk.Label(row2, text="Preserve Thresh:").pack(side=tk.LEFT, padx=(12, 4))
        tk.Entry(row2, textvariable=self.preserveThreshVar, width=6).pack(side=tk.LEFT, padx=2)
        row2.pack(side=tk.TOP, fill=tk.X, pady=2)

        row3 = tk.Frame(top)
        tk.Checkbutton(row3, text="Use previous hour", variable=self.preFlagVar,
                       command=self.prePush).pack(side=tk.LEFT, padx=6)
        tk.Checkbutton(row3, text="Use first hour", variable=self.firstFlagVar).pack(side=tk.LEFT, padx=6)
        tk.Checkbutton(row3, text="Use last hour", variable=self.lastFlagVar).pack(side=tk.LEFT, padx=6)
        tk.Checkbutton(row3, text="Use next hour", variable=self.postFlagVar,
                       command=self.postPush).pack(side=tk.LEFT, padx=6)
        row3.pack(side=tk.TOP, fill=tk.X, pady=2)

        top.pack(side=tk.TOP, fill=tk.X, pady=4)

    def _buildGroupControls(self, master):
        groups = []
        for src in self.sources:
            if src["group"] not in groups:
                groups.append(src["group"])

        fr = tk.Frame(master, relief=tk.GROOVE, borderwidth=2)
        tk.Label(fr, text="Groups:", font=("Helvetica", 9, "bold")).pack(side=tk.LEFT, padx=4)

        for grp in groups:
            var = tk.IntVar(self)
            var.set(1)
            self.groupVars[grp] = var
            tk.Checkbutton(fr, text=grp, variable=var,
                           command=self.updateGroupVisibility).pack(side=tk.LEFT, padx=3)

        fr.pack(side=tk.TOP, fill=tk.X, pady=3)

    def _buildSourcePane(self, master):
        wrapper = tk.Frame(master, relief=tk.GROOVE, borderwidth=2)
        hdr = tk.Frame(wrapper)
        tk.Label(hdr, text="Model / Run", width=36, anchor=tk.W,
                 font=("Helvetica", 10, "bold")).pack(side=tk.LEFT, padx=2)
        tk.Label(hdr, text="Status", width=16, anchor=tk.W,
                 font=("Helvetica", 10, "bold")).pack(side=tk.LEFT, padx=2)
        tk.Label(hdr, text="Wt", width=5, anchor=tk.E,
                 font=("Helvetica", 10, "bold")).pack(side=tk.LEFT, padx=2)
        tk.Label(hdr, text="Pct", width=6, anchor=tk.E,
                 font=("Helvetica", 10, "bold")).pack(side=tk.LEFT, padx=2)
        hdr.pack(side=tk.TOP, fill=tk.X)

        canvas = tk.Canvas(wrapper, height=420, width=980, bd=0, highlightthickness=0)
        sb = tk.Scrollbar(wrapper, orient=tk.VERTICAL, command=canvas.yview)
        self.listFrame = tk.Frame(canvas)

        self.listFrame.bind(
            "<Configure>",
            lambda event, canv=canvas: canv.configure(scrollregion=canv.bbox("all"))
        )

        canvas.create_window((0, 0), window=self.listFrame, anchor="nw")
        canvas.configure(yscrollcommand=sb.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=tk.TRUE)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        smallerFont = ("Helvetica", 9)
        origin = -10 if USE_NEGATIVE_WEIGHTS else 0

        for src in self.sources:
            fr = tk.Frame(self.listFrame, bg="#d9d9d9")
            self.rowFrames.append(fr)

            wv = tk.IntVar(self)
            wv.set(0)
            self.weightVars.append(wv)

            pv = tk.StringVar(self)
            pv.set("  0%")
            self.percentVars.append(pv)

            sv = tk.StringVar(self)
            sv.set(src.get("status", "unchecked"))
            self.statusVars.append(sv)

            fgc = src.get("color", "#000000")
            if src.get("oldrun", 0) == 1:
                fgc = "#606060"

            nameLabel = tk.Label(fr, text=src["display"], width=36, anchor=tk.W,
                                 fg=fgc, bg="#d9d9d9", font=smallerFont)
            self.nameLabels.append(nameLabel)
            nameLabel.pack(side=tk.LEFT, padx=2)

            tk.Label(fr, textvariable=sv, width=16, anchor=tk.W, bg="#d9d9d9",
                     font=smallerFont).pack(side=tk.LEFT, padx=2)

            scale = tk.Scale(fr, orient=tk.HORIZONTAL, from_=origin, to=10,
                             resolution=1, bd=0, highlightthickness=0,
                             variable=wv, length=110, sliderlength=16, width=10,
                             command=self.setPercents, font=smallerFont)
            self.scaleWidgets.append(scale)
            scale.pack(side=tk.LEFT, padx=2)

            tk.Label(fr, textvariable=wv, width=5, anchor=tk.E, bg="#d9d9d9",
                     font=smallerFont).pack(side=tk.LEFT, padx=2)

            tk.Label(fr, textvariable=pv, width=6, anchor=tk.E, bg="#d9d9d9",
                     font=smallerFont).pack(side=tk.LEFT, padx=2)

            fr.pack(side=tk.TOP, fill=tk.X, padx=1, pady=1)

        wrapper.pack(side=tk.TOP, fill=tk.BOTH, expand=tk.TRUE, pady=4)

    def _buildSummaryPane(self, master):
        box = tk.Frame(master, relief=tk.GROOVE, borderwidth=2)

        tk.Label(box, text="Availability / Sources", anchor=tk.W,
                 font=("Helvetica", 10, "bold")).pack(side=tk.TOP, anchor=tk.W, padx=6, pady=(4, 0))
        tk.Label(box, textvariable=self.summaryVars["availability"],
                 justify=tk.LEFT, anchor=tk.W).pack(side=tk.TOP, anchor=tk.W, padx=6)

        tk.Label(box, text="Blend vs NBM / Fcst", anchor=tk.W,
                 font=("Helvetica", 10, "bold")).pack(side=tk.TOP, anchor=tk.W, padx=6, pady=(6, 0))
        tk.Label(box, textvariable=self.summaryVars["comparison"],
                 justify=tk.LEFT, anchor=tk.W).pack(side=tk.TOP, anchor=tk.W, padx=6)

        tk.Label(box, text="Confidence / Disagreement", anchor=tk.W,
                 font=("Helvetica", 10, "bold")).pack(side=tk.TOP, anchor=tk.W, padx=6, pady=(6, 0))
        tk.Label(box, textvariable=self.summaryVars["confidence"],
                 justify=tk.LEFT, anchor=tk.W).pack(side=tk.TOP, anchor=tk.W, padx=6)

        tk.Label(box, text="Verdict", anchor=tk.W,
                 font=("Helvetica", 10, "bold")).pack(side=tk.TOP, anchor=tk.W, padx=6, pady=(6, 0))
        tk.Label(box, textvariable=self.summaryVars["verdict"],
                 justify=tk.LEFT, anchor=tk.W, fg="#7f0000").pack(side=tk.TOP, anchor=tk.W, padx=6, pady=(0, 4))

        box.pack(side=tk.TOP, fill=tk.X, pady=4)

    def updatePresetChoices(self, presetNames):
        self.presetNames = presetNames[:]
        menu = self.presetMenu["menu"]
        menu.delete(0, "end")
        for name in self.presetNames:
            menu.add_command(label=name, command=lambda val=name: self.presetVar.set(val))

    def prePush(self):
        if self.preFlagVar.get() == 1:
            self.firstFlagVar.set(1)

    def postPush(self):
        if self.postFlagVar.get() == 1:
            self.lastFlagVar.set(1)

    def updateGroupVisibility(self):
        for i, src in enumerate(self.sources):
            vis = self.groupVars[src["group"]].get()
            if vis == 1:
                self.rowFrames[i].pack(side=tk.TOP, fill=tk.X, padx=1, pady=1)
            else:
                self.rowFrames[i].pack_forget()

    def setPercents(self, *args):
        total = 0.0
        for i, wv in enumerate(self.weightVars):
            val = float(wv.get())
            if self.normalizeVar.get() == 1:
                total += abs(val)
            else:
                total += max(0.0, val)

            bg = "#ffffd9" if val != 0 else "#d9d9d9"
            self.rowFrames[i].config(bg=bg)
            self.nameLabels[i].config(bg=bg)
            self.scaleWidgets[i].config(bg=bg)

        if total <= 0.0:
            for pv in self.percentVars:
                pv.set("  0%")
        else:
            for i, pv in enumerate(self.percentVars):
                w = float(self.weightVars[i].get())
                if self.normalizeVar.get() == 1:
                    frac = 100.0 * abs(w) / total
                else:
                    frac = 100.0 * max(w, 0.0) / total
                pv.set("%3d%%" % int(round(frac)))

    def setWeights(self, weights):
        for i, val in enumerate(weights):
            if i < len(self.weightVars):
                self.weightVars[i].set(int(val))
        self.setPercents()

    def getState(self):
        return {
            "weights": [int(v.get()) for v in self.weightVars],
            "edgeStyle": self.edgeStyleVar.get(),
            "edgeWidth": int(self.edgeWidthVar.get()),
            "vectorMode": self.vectorModeVar.get(),
            "scalarMode": self.scalarModeVar.get(),
            "normalize": int(self.normalizeVar.get()),
            "preserveEdits": int(self.preserveEditsVar.get()),
            "highConfOnly": int(self.highConfOnlyVar.get()),
            "changeOnly": int(self.changeOnlyVar.get()),
            "changeThreshold": float(self.changeThreshVar.get()),
            "preserveThreshold": float(self.preserveThreshVar.get()),
            "preFlag": int(self.preFlagVar.get()),
            "firstFlag": int(self.firstFlagVar.get()),
            "lastFlag": int(self.lastFlagVar.get()),
            "postFlag": int(self.postFlagVar.get()),
            "presetName": self.presetVar.get(),
        }

    def loadState(self, state):
        if state is None:
            return
        if "weights" in state:
            self.setWeights(state["weights"])
        if "edgeStyle" in state:
            self.edgeStyleVar.set(state["edgeStyle"])
        if "edgeWidth" in state:
            self.edgeWidthVar.set(state["edgeWidth"])
        if "vectorMode" in state:
            self.vectorModeVar.set(state["vectorMode"])
        if "scalarMode" in state:
            self.scalarModeVar.set(state["scalarMode"])
        if "normalize" in state:
            self.normalizeVar.set(state["normalize"])
        if "preserveEdits" in state:
            self.preserveEditsVar.set(state["preserveEdits"])
        if "highConfOnly" in state:
            self.highConfOnlyVar.set(state["highConfOnly"])
        if "changeOnly" in state:
            self.changeOnlyVar.set(state["changeOnly"])
        if "changeThreshold" in state:
            self.changeThreshVar.set(state["changeThreshold"])
        if "preserveThreshold" in state:
            self.preserveThreshVar.set(state["preserveThreshold"])
        if "preFlag" in state:
            self.preFlagVar.set(state["preFlag"])
        if "firstFlag" in state:
            self.firstFlagVar.set(state["firstFlag"])
        if "lastFlag" in state:
            self.lastFlagVar.set(state["lastFlag"])
        if "postFlag" in state:
            self.postFlagVar.set(state["postFlag"])
        self.setPercents()

    def updateSummary(self, preview):
        self.summaryVars["availability"].set(preview.get("availabilityText", ""))
        self.summaryVars["comparison"].set(preview.get("comparisonText", ""))
        self.summaryVars["confidence"].set(preview.get("confidenceText", ""))
        self.summaryVars["verdict"].set(preview.get("verdictText", ""))

        srcStatuses = preview.get("sourceStatuses", {})
        for i, src in enumerate(self.sources):
            if src["key"] in srcStatuses:
                self.statusVars[i].set(srcStatuses[src["key"]])

    def previewCB(self):
        self.__callbackMethod("Preview")

    def applyCB(self):
        self.__callbackMethod("Apply")

    def applyCloseCB(self):
        self.__callbackMethod("ApplyClose")

    def skillCB(self):
        self.__callbackMethod("SuggestSkill")

    def presetCB(self):
        self.__callbackMethod("LoadPreset")

    def savePresetCB(self):
        self.__callbackMethod("SavePreset")

    def zeroAllCB(self):
        for wv in self.weightVars:
            wv.set(0)
        self.setPercents()

    def askPresetName(self):
        return simpledialog.askstring("Save Preset", "Preset name:")

    def closeCB(self):
        self.__callbackMethod("Close")
        self.cancel()

# ----------------------------------------------------------------------------
# TOOL
# ----------------------------------------------------------------------------

class Tool(SmartScript.SmartScript):
    def __init__(self, dbss):
        SmartScript.SmartScript.__init__(self, dbss)
        self._dbss = dbss
        self.setToolType("numeric")

    # ----------------------------------------------------------------------
    # lifecycle
    # ----------------------------------------------------------------------

    def preProcessTool(self):
        self.listAllParms = self.availableParms()
        self.VU = BOIVerifyUtility.BOIVerifyUtility(self._dbss, None)

        self.customPresets = {}
        self.lastStateByParm = {}
        self._readPersistentState()

        self.currentSources = []
        self.currentParmName = None
        self.currentParmLevel = None
        self.currentWxType = None
        self.currentSelectedTR = None
        self.currentGridInfos = []
        self.previewCache = None

    def preProcessGrid(self, WEname=None, ToolTimeRange=None):
        selectedParms = self.selectedParms()
        if not selectedParms:
            self.statusBarMsg("SmartBlendWorkbench: no selected parm", "R")
            self.cancel()
            return

        mutID = self.mutableID().modelIdentifier()
        parmName = None
        parmLevel = None

        for pName, pLevel, pDb in selectedParms:
            if pDb.modelIdentifier() == mutID:
                parmName = pName
                parmLevel = pLevel
                break

        if parmName is None:
            parmName, parmLevel, pDb = selectedParms[0]

        self.currentParmName = parmName
        self.currentParmLevel = parmLevel

        parmObj = self.getParm(self.mutableID(), parmName, parmLevel)
        if parmObj is None:
            self.statusBarMsg("SmartBlendWorkbench: could not find selected parm", "R")
            self.cancel()
            return

        parmInfo = parmObj.getGridInfo()
        self.currentWxType = str(parmInfo.getGridType())

        javaTR = self._dbss.getParmOp().getSelectionTimeRange()
        self.currentSelectedTR = TimeRange.TimeRange(javaTR)

        self.currentGridInfos = self.getGridInfo(self.mutableID().modelIdentifier(),
                                                 parmName, parmLevel, self.currentSelectedTR)
        if not self.currentGridInfos:
            self.statusBarMsg("SmartBlendWorkbench: no mutable grids in selected time range", "R")
            self.cancel()
            return

        self.currentSources = self._buildSources(parmName, parmLevel)

        presetNames = self._presetNamesForParm(parmName)
        initialPreset = "Default"

        self.dlg = ToolDialog(
            title="SmartBlend Workbench - %s" % parmName,
            callbackMethod=self._dialogAction,
            sources=self.currentSources,
            presetNames=presetNames,
            initialPreset=initialPreset
        )

        state = self.lastStateByParm.get(parmName)
        if state is not None:
            # remap model-keyed weights to current positional order
            modelWeights = state.get("modelWeights")
            if modelWeights is not None:
                newWeights = [0] * len(self.currentSources)
                for i, src in enumerate(self.currentSources):
                    if src["model"] in modelWeights:
                        newWeights[i] = modelWeights[src["model"]]
                state["weights"] = newWeights
            self.dlg.loadState(state)
        else:
            self._applyPreset("Default")

        self.dlg.setPercents()
        self.dlg.mainloop()
        self.cancel()

    def execute(self, variableElement):
        return variableElement

    # ----------------------------------------------------------------------
    # dialog actions
    # ----------------------------------------------------------------------

    def _dialogAction(self, action):
        if action == "Close":
            self._saveDialogState()
            self._writePersistentState()
            return

        if action == "LoadPreset":
            self._applyPreset(self.dlg.presetVar.get())
            return

        if action == "SavePreset":
            self._savePresetFromDialog()
            return

        if action == "SuggestSkill":
            self._suggestSkillWeights()
            return

        if action == "Preview":
            preview = self._computePreview(self.dlg.getState())
            self.previewCache = preview
            self.dlg.updateSummary(preview)
            return

        if action == "Apply":
            self._applyBlend(self.dlg.getState())
            self._saveDialogState()
            self._writePersistentState()
            return

        if action == "ApplyClose":
            self._applyBlend(self.dlg.getState())
            self._saveDialogState()
            self._writePersistentState()
            self.dlg.closeCB()
            return

    def _saveDialogState(self):
        if self.currentParmName is not None and self.dlg is not None:
            state = self.dlg.getState()
            # persist a model-keyed weight map so weights survive source-list changes
            modelWeights = {}
            for i, src in enumerate(self.currentSources):
                if i < len(state["weights"]):
                    modelWeights[src["model"]] = state["weights"][i]
            state["modelWeights"] = modelWeights
            self.lastStateByParm[self.currentParmName] = state

    # ----------------------------------------------------------------------
    # source setup
    # ----------------------------------------------------------------------

    def _buildSources(self, parmName, parmLevel):
        sources = []
        seen = {}

        for entry in MODEL_CONFIG:
            model, versions, group, color, parmlist = entry.split("|")
            model = model.strip()
            versions = abs(int(versions.strip()))
            group = group.strip()
            color = color.strip()
            parmlist = parmlist.strip()

            if not self._acceptParmList(parmName, parmlist):
                continue

            for run in range(0, -versions, -1):
                if model == "Fcst":
                    db = self.mutableID()
                elif model == "Official":
                    db = self.findDatabase("Official")
                else:
                    db = self.findDatabase(model, run)

                if db is None:
                    continue

                modid = db.modelIdentifier()
                if modid is None or modid == "":
                    continue

                sourceParm = self._findAvailableParmForModel(db, parmName, parmLevel)
                if sourceParm is None:
                    continue

                modtime = db.modelTime()
                year = modtime.year
                oldrun = 1 if run < 0 else 0

                if model == "Fcst":
                    display = "Forecast:"
                elif model == "Official":
                    display = "Official:"
                elif year == 1970:
                    display = "%s:" % model
                else:
                    display = "%s %02d/%02d %02dZ:" % (model, modtime.month, modtime.day, modtime.hour)

                key = "%s|%s|%s|%s|%d" % (model, modid, sourceParm, parmLevel, run)
                if key in seen:
                    continue
                seen[key] = 1

                src = {
                    "key": key,
                    "model": model,
                    "group": group,
                    "dbId": modid,
                    "parmName": sourceParm,
                    "parmLevel": parmLevel,
                    "display": display,
                    "oldrun": oldrun,
                    "color": color,
                    "status": "ready",
                }
                sources.append(src)

        return sources

    def _findAvailableParmForModel(self, db, parmName, parmLevel):
        tries = [parmName]
        if parmName in PARM_ALTERNATES:
            tries.extend(PARM_ALTERNATES[parmName])

        modid = db.modelIdentifier()

        # exact name + level first
        for tryName in tries:
            for pName, pLevel, pDb in self.listAllParms:
                if pName == tryName and pLevel == parmLevel and pDb.modelIdentifier() == modid:
                    return tryName

        # relaxed level if exact level missing
        for tryName in tries:
            for pName, pLevel, pDb in self.listAllParms:
                if pName == tryName and pDb.modelIdentifier() == modid:
                    return tryName

        return None

    def _acceptParmList(self, parmName, parmlist):
        parts = parmlist.split(",")
        if parts[0] == "ALL":
            return True
        invert = 0
        if parts[0].startswith("^"):
            invert = 1
            parts[0] = parts[0][1:]
        found = parmName in parts
        return (not found) if invert else found

    # ----------------------------------------------------------------------
    # presets / persistence
    # ----------------------------------------------------------------------

    def _presetNamesForParm(self, parmName):
        names = ["Default"]
        if parmName in DEFAULT_PRESETS:
            names.append(parmName)
        customNames = []
        for name, rec in self.customPresets.items():
            pe = rec.get("element", "ALL")
            if pe in ["ALL", parmName]:
                customNames.append(name)
        customNames.sort()
        names.extend(customNames)
        names.append("Skill30d")
        return names

    def _applyPreset(self, presetName):
        weights = [0] * len(self.currentSources)
        state = None

        if presetName == "Skill30d":
            self._suggestSkillWeights()
            return

        if presetName in self.customPresets:
            rec = self.customPresets[presetName]
            state = rec.get("state")
            if state is not None:
                self.dlg.loadState(state)
                return

        key = presetName if presetName in DEFAULT_PRESETS else "DEFAULT"
        recipe = DEFAULT_PRESETS.get(key, DEFAULT_PRESETS["DEFAULT"])
        modelToWeight = {}
        for model, wt in recipe:
            modelToWeight[model] = wt

        for i, src in enumerate(self.currentSources):
            if src["model"] in modelToWeight:
                weights[i] = modelToWeight[src["model"]]

        self.dlg.setWeights(weights)

    def _savePresetFromDialog(self):
        pname = self.dlg.askPresetName()
        if pname is None:
            return
        pname = pname.strip()
        if pname == "":
            return

        self.customPresets[pname] = {
            "element": self.currentParmName,
            "state": self.dlg.getState(),
        }
        self.dlg.updatePresetChoices(self._presetNamesForParm(self.currentParmName))
        self.dlg.presetVar.set(pname)
        self._writePersistentState()
        self.statusBarMsg("SmartBlendWorkbench: preset saved", "S")

    def _writePersistentState(self):
        try:
            self.saveObject("customPresets", self.customPresets, PERSIST_CATEGORY)
            self.saveObject("lastStateByParm", self.lastStateByParm, PERSIST_CATEGORY)
        except Exception as e:
            self.statusBarMsg("SmartBlendWorkbench: could not save state: %s" % str(e), "A")

    def _readPersistentState(self):
        try:
            lf = self._getUserFile("customPresets", PERSIST_CATEGORY)
            if lf.exists():
                self.customPresets = self.getObject("customPresets", PERSIST_CATEGORY)
        except Exception:
            self.customPresets = {}
            self.statusBarMsg("SmartBlendWorkbench: could not load saved presets", "A")

        try:
            lf = self._getUserFile("lastStateByParm", PERSIST_CATEGORY)
            if lf.exists():
                self.lastStateByParm = self.getObject("lastStateByParm", PERSIST_CATEGORY)
        except Exception:
            self.lastStateByParm = {}
            self.statusBarMsg("SmartBlendWorkbench: could not load last-used settings", "A")

    def _getUserFile(self, name, category):
        from com.raytheon.uf.common.localization import PathManagerFactory, LocalizationContext
        pathMgr = PathManagerFactory.getPathManager()
        lc = pathMgr.getContext(
            LocalizationContext.LocalizationType.valueOf('CAVE_STATIC'),
            LocalizationContext.LocalizationLevel.valueOf('USER')
        )
        return pathMgr.getLocalizationFile(lc, 'gfe/userPython/' + category + '/' + name)

    # ----------------------------------------------------------------------
    # preview computation
    # ----------------------------------------------------------------------

    def _computePreview(self, settings):
        parmName = self.currentParmName
        parmLevel = self.currentParmLevel
        wxType = self.currentWxType
        thresh = self._thresholds(parmName)
        mask = self._getActiveMask()

        agg = {
            "requestedWeightedCount": 0.0,
            "usedCount": 0.0,
            "missingCount": 0.0,
            "droppedWeightFraction": 0.0,
            "dominantSourceFraction": 0.0,
            "dominantSourceText": "none",
            "spreadMean": 0.0,
            "spreadP90": 0.0,
            "highConfFrac": 0.0,
            "lowConfFrac": 0.0,
            "diffNBMMean": None,
            "diffNBMP90": None,
            "diffNBMFrac": None,
            "diffFcstMean": None,
            "diffFcstP90": None,
            "diffFcstFrac": None,
            "supportFrac": None,
            "revertFrac": None,
        }

        n = 0
        sourceStatusMap = {}
        domCount = {}

        for src in self.currentSources:
            sourceStatusMap[src["key"]] = []

        for gridInfo in self.currentGridInfos:
            gridTR = gridInfo.gridTime()

            fcstGrid = self.getGrids(self.mutableID().modelIdentifier(), parmName, parmLevel,
                                     gridTR, noDataError=0, cache=0)
            if fcstGrid is None:
                continue

            nbmGrid = None
            nbmSrc = self._latestModelSource("NBM")
            if nbmSrc is not None:
                nbmRes, nbmMeta = self._buildBlendForTR(
                    parmName, parmLevel, gridTR,
                    weights=[1], customSources=[nbmSrc],
                    settings=settings, allowNormalize=0
                )
                if nbmRes is not None:
                    nbmGrid = nbmRes["result"]

            blendRes, meta = self._buildBlendForTR(
                parmName, parmLevel, gridTR,
                weights=settings["weights"], customSources=None,
                settings=settings, allowNormalize=1
            )
            if blendRes is None:
                continue

            n += 1

            for key, status in meta["sourceStatuses"].items():
                if key in sourceStatusMap:
                    sourceStatusMap[key].append(status)

            blendGrid = blendRes["result"]
            spreadGrid = blendRes["spread"]

            magBlend = self._gridMagnitude(blendGrid, wxType)
            magFcst = self._gridMagnitude(fcstGrid, wxType)
            magNBM = self._gridMagnitude(nbmGrid, wxType) if nbmGrid is not None else None

            agg["requestedWeightedCount"] += meta["requestedWeightedCount"]
            agg["usedCount"] += meta["usedCount"]
            agg["missingCount"] += meta["missingCount"]
            agg["droppedWeightFraction"] += meta["droppedWeightFraction"]
            agg["dominantSourceFraction"] += meta["dominantSourceFraction"]
            agg["spreadMean"] += self._maskedMean(spreadGrid, mask)
            agg["spreadP90"] += self._maskedPercentile(spreadGrid, mask, 90)
            agg["highConfFrac"] += self._maskedFrac(spreadGrid <= thresh["spread_good"], mask)
            agg["lowConfFrac"] += self._maskedFrac(spreadGrid > thresh["spread_warn"], mask)

            dname = meta["dominantSourceText"]
            domCount[dname] = domCount.get(dname, 0) + 1

            diffFcst = np.abs(magBlend - magFcst)
            if agg["diffFcstMean"] is None:
                agg["diffFcstMean"] = 0.0
                agg["diffFcstP90"] = 0.0
                agg["diffFcstFrac"] = 0.0
            agg["diffFcstMean"] += self._maskedMean(diffFcst, mask)
            agg["diffFcstP90"] += self._maskedPercentile(diffFcst, mask, 90)
            agg["diffFcstFrac"] += self._maskedFrac(diffFcst > thresh["change"], mask)

            if magNBM is not None:
                diffNBM = np.abs(magBlend - magNBM)
                if agg["diffNBMMean"] is None:
                    agg["diffNBMMean"] = 0.0
                    agg["diffNBMP90"] = 0.0
                    agg["diffNBMFrac"] = 0.0
                    agg["supportFrac"] = 0.0
                    agg["revertFrac"] = 0.0

                agg["diffNBMMean"] += self._maskedMean(diffNBM, mask)
                agg["diffNBMP90"] += self._maskedPercentile(diffNBM, mask, 90)
                agg["diffNBMFrac"] += self._maskedFrac(diffNBM > thresh["depart"], mask)

                fcstDepartMask = np.abs(magFcst - magNBM) > thresh["depart"]
                supportMask = fcstDepartMask & (np.abs(magBlend - magFcst) <= np.abs(magBlend - magNBM))
                revertMask = fcstDepartMask & (np.abs(magBlend - magNBM) < np.abs(magBlend - magFcst))
                agg["supportFrac"] += self._maskedFrac(supportMask, mask)
                agg["revertFrac"] += self._maskedFrac(revertMask, mask)

        if n == 0:
            return {
                "availabilityText": "No preview available.\nEither all weights were zero or no source grids were found.",
                "comparisonText": "",
                "confidenceText": "",
                "verdictText": "No blend computed.",
                "sourceStatuses": {},
            }

        for key in ["requestedWeightedCount", "usedCount", "missingCount",
                    "droppedWeightFraction", "dominantSourceFraction",
                    "spreadMean", "spreadP90", "highConfFrac", "lowConfFrac"]:
            agg[key] /= float(n)

        for key in ["diffNBMMean", "diffNBMP90", "diffNBMFrac",
                    "diffFcstMean", "diffFcstP90", "diffFcstFrac",
                    "supportFrac", "revertFrac"]:
            if agg[key] is not None:
                agg[key] /= float(n)

        # aggregate statuses
        finalStatuses = {}
        for src in self.currentSources:
            key = src["key"]
            vals = sourceStatusMap[key]
            if len(vals) == 0:
                finalStatuses[key] = "weight=0"
            elif all([v == "used" for v in vals]):
                finalStatuses[key] = "used"
            elif all([v == "missing" for v in vals]):
                finalStatuses[key] = "missing"
            elif all([v == "weight=0" for v in vals]):
                finalStatuses[key] = "weight=0"
            else:
                finalStatuses[key] = "partial"

        dominantText = "none"
        if len(domCount) > 0:
            dominantText = sorted(domCount.items(), key=lambda x: x[1], reverse=True)[0][0]

        availabilityText = (
            "Requested weighted sources: %d\n"
            "Contributing sources: %d\n"
            "Dropped weighted sources: %d\n"
            "Dropped requested weight: %d%%\n"
            "Dominant effective source: %s"
        ) % (
            int(round(agg["requestedWeightedCount"])),
            int(round(agg["usedCount"])),
            int(round(agg["missingCount"])),
            int(round(100.0 * agg["droppedWeightFraction"])),
            dominantText
        )

        comparisonLines = []
        if agg["diffNBMMean"] is not None:
            comparisonLines.append(
                "Blend vs NBM: mean |diff|=%s, p90=%s, area > thresh=%d%%" % (
                    self._fmt(agg["diffNBMMean"]),
                    self._fmt(agg["diffNBMP90"]),
                    int(round(100.0 * agg["diffNBMFrac"]))
                )
            )

        if agg["diffFcstMean"] is not None:
            comparisonLines.append(
                "Blend vs Fcst: mean |diff|=%s, p90=%s, area > thresh=%d%%" % (
                    self._fmt(agg["diffFcstMean"]),
                    self._fmt(agg["diffFcstP90"]),
                    int(round(100.0 * agg["diffFcstFrac"]))
                )
            )

        if agg["supportFrac"] is not None:
            comparisonLines.append(
                "Where Fcst already departs from NBM: support=%d%%, trends back to NBM=%d%%" % (
                    int(round(100.0 * agg["supportFrac"])),
                    int(round(100.0 * agg["revertFrac"]))
                )
            )

        comparisonText = "\n".join(comparisonLines)

        confidenceText = (
            "Mean spread=%s\n"
            "P90 spread=%s\n"
            "High-confidence area=%d%%\n"
            "Low-confidence area=%d%%"
        ) % (
            self._fmt(agg["spreadMean"]),
            self._fmt(agg["spreadP90"]),
            int(round(100.0 * agg["highConfFrac"])),
            int(round(100.0 * agg["lowConfFrac"]))
        )

        verdict, reasons = self._gradePreview(agg, thresh)
        verdictText = "%s\n%s" % (verdict, "\n".join(["- " + r for r in reasons]))

        return {
            "availabilityText": availabilityText,
            "comparisonText": comparisonText,
            "confidenceText": confidenceText,
            "verdictText": verdictText,
            "sourceStatuses": finalStatuses,
        }

    def _gradePreview(self, agg, thresh):
        score = 0
        reasons = []

        dropped = agg["droppedWeightFraction"]
        if dropped > 0.30:
            score += 30
            reasons.append("Large fraction of requested weight is unavailable")
        elif dropped > 0.15:
            score += 15
            reasons.append("Some weighted sources are unavailable")

        dom = agg["dominantSourceFraction"]
        if dom > 0.75:
            score += 20
            reasons.append("Effective blend is dominated by one source")
        elif dom > 0.55:
            score += 10
            reasons.append("One source carries a large share of the blend")

        lowConfFrac = agg["lowConfFrac"]
        if lowConfFrac > 0.30:
            score += 25
            reasons.append("Large low-confidence / high-spread area")
        elif lowConfFrac > 0.15:
            score += 12
            reasons.append("Moderate low-confidence area present")

        if agg["diffFcstFrac"] is not None:
            if agg["diffFcstFrac"] > 0.40:
                score += 20
                reasons.append("Blend makes a broad large change from Fcst")
            elif agg["diffFcstFrac"] > 0.20:
                score += 10
                reasons.append("Blend changes a meaningful portion of the area")

        if agg["revertFrac"] is not None and agg["revertFrac"] > 0.25:
            score += 15
            reasons.append("Blend trends back toward NBM over existing manual forecast departures")

        if score >= 50:
            verdict = "Blend Quality: WARNING"
        elif score >= 25:
            verdict = "Blend Quality: CAUTION"
        else:
            verdict = "Blend Quality: GOOD"

        if not reasons:
            reasons.append("No major warning signals detected")

        return verdict, reasons

    # ----------------------------------------------------------------------
    # apply
    # ----------------------------------------------------------------------

    def _applyBlend(self, settings):
        parmName = self.currentParmName
        parmLevel = self.currentParmLevel
        wxType = self.currentWxType
        thresh = self._thresholds(parmName)

        weights = settings["weights"]
        if max([abs(w) for w in weights] + [0]) == 0:
            self.statusBarMsg("SmartBlendWorkbench: all weights are zero", "R")
            return

        mask = self._getActiveMask()
        numPoints = int(np.sum(mask))
        if CHECK_SMALL_EDIT_AREA and numPoints < MIN_EDIT_AREA_POINTS:
            if not msg.askyesno("Continue?", "Small edit area detected. Continue?"):
                return

        numApplied = 0

        # these don't change per grid -- look up once
        parmObj = self.getParm(self.mutableID(), parmName, parmLevel)
        parmInfo = parmObj.getGridInfo()
        minLimit = parmInfo.getMinValue()
        maxLimit = parmInfo.getMaxValue()
        nbmSrc = self._latestModelSource("NBM")

        for gridInfo in self.currentGridInfos:
            gridTR = gridInfo.gridTime()

            oldGrid = self.getGrids(self.mutableID().modelIdentifier(), parmName, parmLevel,
                                    gridTR, noDataError=0, cache=0)
            if oldGrid is None:
                continue

            blendRes, meta = self._buildBlendForTR(
                parmName, parmLevel, gridTR,
                weights=weights, customSources=None,
                settings=settings, allowNormalize=1
            )
            if blendRes is None:
                continue

            newGrid = blendRes["result"]
            spreadGrid = blendRes["spread"]
            extraMask = self.newGrid(1.0)
            extraMask[:] = 1.0

            magNew = self._gridMagnitude(newGrid, wxType)
            magOld = self._gridMagnitude(oldGrid, wxType)

            nbmGrid = None
            if nbmSrc is not None:
                nbmRes, nbmMeta = self._buildBlendForTR(
                    parmName, parmLevel, gridTR,
                    weights=[1], customSources=[nbmSrc],
                    settings=settings, allowNormalize=0
                )
                if nbmRes is not None:
                    nbmGrid = nbmRes["result"]

            if settings["preserveEdits"] == 1 and nbmGrid is not None:
                magNBM = self._gridMagnitude(nbmGrid, wxType)
                # "overwrite allowed" mask: 1.0 where Fcst is close to NBM (ok to blend),
                # 0.0 where forecaster already departed (protect those edits)
                overwriteMask = np.where(np.abs(magOld - magNBM) <= settings["preserveThreshold"], 1.0, 0.0)
                extraMask *= overwriteMask

            if settings["highConfOnly"] == 1:
                highMask = np.where(spreadGrid <= thresh["spread_good"], 1.0, 0.0)
                extraMask *= highMask

            if settings["changeOnly"] == 1:
                changeMask = np.where(np.abs(magNew - magOld) >= settings["changeThreshold"], 1.0, 0.0)
                extraMask *= changeMask

            finalGrid = self._applyInEditArea(newGrid, oldGrid, extraMask,
                                              settings["edgeStyle"], settings["edgeWidth"], wxType)

            if wxType == "SCALAR":
                np.clip(finalGrid, minLimit, maxLimit, out=finalGrid)
                self.createGrid(self.mutableID().modelIdentifier(), parmName, wxType, finalGrid, gridTR)
            else:
                mag, direc = finalGrid
                np.clip(mag, minLimit, maxLimit, out=mag)
                self.createGrid(self.mutableID().modelIdentifier(), parmName, wxType, (mag, direc), gridTR)

            numApplied += 1

        if numApplied > 0:
            self.statusBarMsg("SmartBlendWorkbench: blend applied", "S")
        else:
            self.statusBarMsg("SmartBlendWorkbench: nothing applied", "A")

    # ----------------------------------------------------------------------
    # blend engine
    # ----------------------------------------------------------------------

    def _buildBlendForTR(self, parmName, parmLevel, targetTR, weights,
                         customSources=None, settings=None, allowNormalize=1):
        if settings is None:
            settings = self.dlg.getState()

        wxType = self.currentWxType
        sources = customSources if customSources is not None else self.currentSources
        readTR = self._adjustTimeRange(targetTR, settings)

        requestedWeightedCount = 0
        requestedAbsWeight = 0.0
        missingWeight = 0.0
        usedCount = 0
        sourceStatuses = {}

        if wxType == "SCALAR":
            numer = None
            denom = 0.0
            diagWeights = []
            diagGrids = []
            diagNames = []

            for i, src in enumerate(sources):
                wt = weights[i]
                if wt == 0:
                    sourceStatuses[src["key"]] = "weight=0"
                    continue

                requestedWeightedCount += 1
                requestedAbsWeight += abs(float(wt))

                scalarMode = self._modeForParm(parmName, settings["scalarMode"])
                grid = self.getGrids(src["dbId"], src["parmName"], src["parmLevel"],
                                     readTR, mode=scalarMode, noDataError=0, cache=0)
                if grid is None:
                    sourceStatuses[src["key"]] = "missing"
                    missingWeight += abs(float(wt))
                    continue

                if numer is None:
                    numer = self.empty()

                numer += (grid * float(wt))
                denom += float(wt)

                diagWeights.append(abs(float(wt)))
                diagGrids.append(grid.copy())
                diagNames.append(src["display"])
                usedCount += 1
                sourceStatuses[src["key"]] = "used"

            if numer is None or abs(denom) < 1.0e-6:
                return None, self._emptyMeta(sourceStatuses, requestedWeightedCount, usedCount, missingWeight, requestedAbsWeight)

            result = numer / denom if (allowNormalize == 1 and settings["normalize"] == 1) else numer
            spread = self._weightedSpread(diagGrids, diagWeights, result)

            dominantName = "none"
            dominantFrac = 0.0
            if len(diagWeights) > 0:
                sw = float(sum(diagWeights))
                maxw = max(diagWeights)
                idx = diagWeights.index(maxw)
                dominantName = diagNames[idx]
                dominantFrac = maxw / max(sw, 1.0e-6)

            return {"result": result, "spread": spread}, {
                "sourceStatuses": sourceStatuses,
                "requestedWeightedCount": requestedWeightedCount,
                "usedCount": usedCount,
                "missingCount": requestedWeightedCount - usedCount,
                "droppedWeightFraction": (missingWeight / requestedAbsWeight) if requestedAbsWeight > 0 else 0.0,
                "dominantSourceText": dominantName,
                "dominantSourceFraction": dominantFrac,
            }

        elif wxType == "VECTOR":
            vecMode = settings["vectorMode"]

            uNumer = None
            vNumer = None
            magNumer = None
            denom = 0.0

            diagWeights = []
            diagMags = []
            diagNames = []

            for i, src in enumerate(sources):
                wt = weights[i]
                if wt == 0:
                    sourceStatuses[src["key"]] = "weight=0"
                    continue

                requestedWeightedCount += 1
                requestedAbsWeight += abs(float(wt))

                grid = self.getGrids(src["dbId"], src["parmName"], src["parmLevel"],
                                     readTR, noDataError=0, cache=0)
                if grid is None:
                    sourceStatuses[src["key"]] = "missing"
                    missingWeight += abs(float(wt))
                    continue

                # guard: some models may return scalar for a vector parm
                if not isinstance(grid, (tuple, list)) or len(grid) != 2:
                    sourceStatuses[src["key"]] = "missing"
                    missingWeight += abs(float(wt))
                    continue

                mag, direc = grid
                if vecMode == "Both":
                    u, v = self.MagDirToUV(mag, direc)
                    if uNumer is None:
                        uNumer = self.empty()
                        vNumer = self.empty()
                    uNumer += (u * float(wt))
                    vNumer += (v * float(wt))
                    diagMags.append(mag.copy())
                elif vecMode == "Magnitude Only":
                    if magNumer is None:
                        magNumer = self.empty()
                    magNumer += (mag * float(wt))
                    diagMags.append(mag.copy())
                else:
                    u, v = self.MagDirToUV(10.0, direc)
                    if uNumer is None:
                        uNumer = self.empty()
                        vNumer = self.empty()
                    uNumer += (u * float(wt))
                    vNumer += (v * float(wt))
                    diagMags.append(mag.copy())

                denom += float(wt)
                diagWeights.append(abs(float(wt)))
                diagNames.append(src["display"])
                usedCount += 1
                sourceStatuses[src["key"]] = "used"

            if usedCount == 0 or abs(denom) < 1.0e-6:
                return None, self._emptyMeta(sourceStatuses, requestedWeightedCount, usedCount, missingWeight, requestedAbsWeight)

            if vecMode == "Magnitude Only":
                magResult = magNumer / denom if (allowNormalize == 1 and settings["normalize"] == 1) else magNumer
                fcstGrid = self.getGrids(self.mutableID().modelIdentifier(), parmName, parmLevel,
                                         targetTR, noDataError=0, cache=0)
                if fcstGrid is not None:
                    direcResult = fcstGrid[1]
                else:
                    direcResult = self.empty()
                    direcResult[:] = 0.0
                result = (magResult, direcResult)
                spread = self._weightedSpread(diagMags, diagWeights, magResult)
            elif vecMode == "Direction Only":
                doNorm = (allowNormalize == 1 and settings["normalize"] == 1)
                uFinal = uNumer / denom if doNorm else uNumer
                vFinal = vNumer / denom if doNorm else vNumer
                magBlend, dirBlend = self.UVToMagDir(uFinal, vFinal)
                fcstGrid = self.getGrids(self.mutableID().modelIdentifier(), parmName, parmLevel,
                                         targetTR, noDataError=0, cache=0)
                if fcstGrid is not None:
                    magUse = fcstGrid[0]
                else:
                    magUse = magBlend
                result = (magUse, dirBlend)
                spread = self._weightedSpread(diagMags, diagWeights, magUse)
            else:
                doNorm = (allowNormalize == 1 and settings["normalize"] == 1)
                uFinal = uNumer / denom if doNorm else uNumer
                vFinal = vNumer / denom if doNorm else vNumer
                result = self.UVToMagDir(uFinal, vFinal)
                spread = self._weightedSpread(diagMags, diagWeights, result[0])

            dominantName = "none"
            dominantFrac = 0.0
            if len(diagWeights) > 0:
                sw = float(sum(diagWeights))
                maxw = max(diagWeights)
                idx = diagWeights.index(maxw)
                dominantName = diagNames[idx]
                dominantFrac = maxw / max(sw, 1.0e-6)

            return {"result": result, "spread": spread}, {
                "sourceStatuses": sourceStatuses,
                "requestedWeightedCount": requestedWeightedCount,
                "usedCount": usedCount,
                "missingCount": requestedWeightedCount - usedCount,
                "droppedWeightFraction": (missingWeight / requestedAbsWeight) if requestedAbsWeight > 0 else 0.0,
                "dominantSourceText": dominantName,
                "dominantSourceFraction": dominantFrac,
            }

        return None, {}

    def _emptyMeta(self, sourceStatuses, requestedWeightedCount, usedCount, missingWeight, requestedAbsWeight):
        return {
            "sourceStatuses": sourceStatuses,
            "requestedWeightedCount": requestedWeightedCount,
            "usedCount": usedCount,
            "missingCount": requestedWeightedCount - usedCount,
            "droppedWeightFraction": (missingWeight / requestedAbsWeight) if requestedAbsWeight > 0 else 1.0,
            "dominantSourceText": "none",
            "dominantSourceFraction": 0.0,
        }

    def _weightedSpread(self, grids, weights, blend):
        out = self.empty()
        out[:] = 0.0
        if len(grids) <= 1:
            return out

        denom = max(float(sum(weights)), 1.0e-6)
        for grid, wt in zip(grids, weights):
            out += (wt * ((grid - blend) ** 2.0))
        out = np.sqrt(out / denom)
        return out

    def _applyInEditArea(self, newGrid, oldGrid, extraMask, edgeStyle, edgeWidth, wxType):
        editArea = self.getActiveEditArea()
        editAreaMask = editArea.getGrid()
        if not editAreaMask.isAnyBitsSet():
            editArea.invert()

        if edgeStyle == "Flat":
            edgegrid = self.encodeEditArea(editArea)
        elif edgeStyle == "Taper":
            edgegrid = self.taperGrid(editArea, edgeWidth)
        else:
            # "Edge" -- minimal taper at boundary
            edgegrid = self.taperGrid(editArea, 0)

        factor = edgegrid * extraMask

        if wxType == "SCALAR":
            return oldGrid + ((newGrid - oldGrid) * factor)

        newMag, newDir = newGrid
        oldMag, oldDir = oldGrid
        newU, newV = self.MagDirToUV(newMag, newDir)
        oldU, oldV = self.MagDirToUV(oldMag, oldDir)
        uFinal = oldU + ((newU - oldU) * factor)
        vFinal = oldV + ((newV - oldV) * factor)
        return self.UVToMagDir(uFinal, vFinal)

    # ----------------------------------------------------------------------
    # BOIVerify skill weighting
    # ----------------------------------------------------------------------

    def _suggestSkillWeights(self):
        parmName = self.currentParmName
        parmLevel = self.currentParmLevel
        mask = self._getActiveMask()

        if not self.VU.checkFile(OBS_ELEMENT_MAP.get(parmName, parmName), SKILL_OBS_SOURCE):
            self.statusBarMsg("SmartBlendWorkbench: no %s obs file for %s" % (SKILL_OBS_SOURCE, parmName), "A")
            return

        results = {}
        modelsDone = {}

        for src in self.currentSources:
            model = src["model"]
            if model in modelsDone:
                continue
            modelsDone[model] = 1

            if model in ["Fcst", "Official"]:
                results[model] = None
                continue

            score = self._recentModelSkill(model, parmName, parmLevel, mask)
            results[model] = score

        valid = [v for v in results.values() if v is not None and v > 0]
        if not valid:
            self.statusBarMsg("SmartBlendWorkbench: no BOIVerify skill data found", "A")
            return

        inv = {}
        for model, mae in results.items():
            if mae is not None and mae > 0:
                inv[model] = 1.0 / max(mae, 1.0e-6)

        vmax = max(inv.values())
        weights = [0] * len(self.currentSources)
        latestSeen = {}

        for i, src in enumerate(self.currentSources):
            model = src["model"]
            if model not in inv:
                continue
            if model in latestSeen:
                weights[i] = 0
                continue
            latestSeen[model] = 1
            weights[i] = max(1, min(10, int(round(10.0 * inv[model] / vmax))))

        for i, src in enumerate(self.currentSources):
            if src["model"] == "Fcst":
                weights[i] = max(weights[i], 2)
            elif src["model"] == "NBM":
                weights[i] = max(weights[i], 4)

        self.dlg.setWeights(weights)
        self.statusBarMsg("SmartBlendWorkbench: loaded BOIVerify skill-suggested weights", "S")

    def _recentModelSkill(self, model, parmName, parmLevel, mask):
        obsElement = OBS_ELEMENT_MAP.get(parmName, parmName)
        archiveModel = "Official" if model == "Fcst" else model
        modeVal = self._archiveModeForParm(parmName)

        count = 0
        errsum = 0.0

        startUnix = self.currentSelectedTR.startTime().unixTime()
        endUnix = self.currentSelectedTR.endTime().unixTime()
        duration = endUnix - startUnix

        if not self.VU.checkFile(obsElement, SKILL_OBS_SOURCE):
            return None

        parmToUse = parmName
        if not self.VU.checkFile(parmToUse, archiveModel):
            if parmName in PARM_ALTERNATES and len(PARM_ALTERNATES[parmName]) > 0:
                altparm = PARM_ALTERNATES[parmName][0]
                if self.VU.checkFile(altparm, archiveModel):
                    parmToUse = altparm
                else:
                    return None
            else:
                return None

        for daysBack in range(1, SKILL_LOOKBACK_DAYS + 1):
            histStart = startUnix - (daysBack * 86400)
            histEnd = histStart + duration

            try:
                obsGrid = self.VU.getVerGrids(SKILL_OBS_SOURCE, 0, obsElement,
                                              histStart, histEnd, mode=modeVal)
            except Exception:
                obsGrid = None

            if obsGrid is None:
                continue
            if not isinstance(obsGrid, np.ndarray) and len(obsGrid) == 2:
                obsGrid = obsGrid[0]

            cases = self._getModelCases(parmToUse, archiveModel, histStart, histEnd)
            if not cases:
                continue

            casekeys = sorted(list(cases.keys()))
            casekeys.reverse()

            base = None
            recs = None
            for ck in casekeys:
                bstr, sstr, estr = ck.split(",")
                base = int(bstr)
                recs = cases[ck]
                break

            if base is None or recs is None:
                continue

            try:
                modGrid = self.VU.getVerGrids(archiveModel, base, parmToUse,
                                              histStart, histEnd, mode=modeVal, recList=recs)
            except Exception:
                modGrid = None

            if modGrid is None:
                continue
            if not isinstance(modGrid, np.ndarray) and len(modGrid) == 2:
                modGrid = modGrid[0]

            diff = np.abs(modGrid - obsGrid)
            errsum += self._maskedMean(diff, mask)
            count += 1

            if count >= SKILL_MAX_CASES_PER_MODEL:
                break

        if count == 0:
            return None
        return errsum / float(count)

    def _getModelCases(self, parm, model, perStart, perEnd):
        cases = {}
        if not self.VU.checkFile(parm, model):
            return cases

        try:
            endsin = (self.VU.fncEtime[:] > perStart) & (self.VU.fncEtime[:] <= perEnd)
            startsin = (self.VU.fncStime[:] < perEnd) & (self.VU.fncStime[:] >= perStart)
            crosses = (self.VU.fncStime[:] < perStart) & (self.VU.fncEtime[:] > perEnd)
            recmatch = endsin | startsin | crosses

            if recmatch.any():
                recnumberList = list(np.compress(recmatch, self.VU.fncRecs[:]))
                baselist = []
                for rec in recnumberList:
                    recnum = int(rec)
                    base = self.VU.fncBtime[recnum]
                    if base not in baselist:
                        baselist.append(base)

                for base in baselist:
                    reclist = [int(rec) for rec in recnumberList if self.VU.fncBtime[int(rec)] == base]
                    cases["%d,%d,%d" % (base, perStart, perEnd)] = reclist
        except (AttributeError, IndexError, TypeError):
            pass

        return cases

    # ----------------------------------------------------------------------
    # time / thresholds / modes
    # ----------------------------------------------------------------------

    def _adjustTimeRange(self, targetTR, settings):
        start = targetTR.startTime().unixTime()
        end = targetTR.endTime().unixTime()

        if settings["preFlag"] == 1:
            start -= 3600
        if settings["firstFlag"] == 0:
            start += 3600
        if settings["lastFlag"] == 0:
            end -= 3600
        if settings["postFlag"] == 1:
            end += 3600

        # prevent zero or negative duration from collapsing the range
        if end <= start:
            start = targetTR.startTime().unixTime()
            end = targetTR.endTime().unixTime()

        return TimeRange.TimeRange(AbsTime.AbsTime(start), AbsTime.AbsTime(end))

    def _modeForParm(self, parmName, scalarMode):
        if parmName in ["QPF", "SnowAmt"]:
            return "Sum"
        if parmName in ["MaxT", "MaxRH"]:
            return "Max"
        if parmName in ["MinT", "MinRH"]:
            return "Min"
        return scalarMode

    def _archiveModeForParm(self, parmName):
        if parmName in ["QPF", "SnowAmt"]:
            return "Sum"
        if parmName in ["MaxT", "MaxRH"]:
            return "Max"
        if parmName in ["MinT", "MinRH"]:
            return "Min"
        return "TimeWtAverage"

    def _thresholds(self, parmName):
        if parmName in ELEMENT_THRESHOLDS:
            return ELEMENT_THRESHOLDS[parmName]
        return ELEMENT_THRESHOLDS["DEFAULT"]

    # ----------------------------------------------------------------------
    # references / utilities
    # ----------------------------------------------------------------------

    def _latestModelSource(self, modelName):
        for src in self.currentSources:
            if src["model"] == modelName:
                return src
        return None

    def _getActiveMask(self):
        editArea = self.getActiveEditArea()
        mask = self.encodeEditArea(editArea)
        if np.sum(mask) == 0:
            mask = self.newGrid(1.0)
            mask[:] = 1.0
        return mask

    def _gridMagnitude(self, grid, wxType):
        if grid is None:
            return None
        if wxType == "SCALAR":
            return grid
        return grid[0]

    def _maskedValues(self, grid, mask):
        flat = np.ravel(grid)
        mflat = np.ravel(mask)
        return np.compress(mflat > 0, flat)

    def _maskedMean(self, grid, mask):
        vals = self._maskedValues(grid, mask)
        if vals.size == 0:
            return 0.0
        return float(np.add.reduce(vals)) / float(vals.size)

    def _maskedPercentile(self, grid, mask, pct):
        vals = self._maskedValues(grid, mask)
        if vals.size == 0:
            return 0.0
        svals = np.sort(vals)
        idx = int(round((pct / 100.0) * (len(svals) - 1)))
        idx = max(0, min(idx, len(svals) - 1))
        return float(svals[idx])

    def _maskedFrac(self, boolGrid, mask):
        vals = self._maskedValues(boolGrid.astype(np.float32), mask)
        if vals.size == 0:
            return 0.0
        return float(np.add.reduce(vals > 0.5)) / float(vals.size)

    def _fmt(self, val):
        aval = abs(float(val))
        if aval >= 100:
            return "%.0f" % val
        elif aval >= 10:
            return "%.1f" % val
        elif aval >= 1:
            return "%.2f" % val
        else:
            return "%.3f" % val
