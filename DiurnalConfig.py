# ----------------------------------------------------------------------------
# SVN: $Revision: 32073 $  $Date: 2022-05-01 01:34:25 +0000 (Sun, 01 May 2022) $
# $URL: https://vlab.noaa.gov/svn/nwsscp/Gfe/Procedures/NwsDiurnal/tags/latest_stable/DiurnalConfig-example.Procedure $
#
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# DiurnalConfig
#
# Author: Tom LeFebvre
# Maintainer: Jonathan Lamb (WFO CHS)
#
# 2026/09/02 - version 3.0 KA. Keep the hourly model curve anchored
#          to Official MinT/MaxT; defer cycle-time defaults to the GFE clock.
#         2.6.1 - 27 Apr 2022
#
# ----------------------------------------------------------------------------

# The MenuItems list defines the GFE menu item(s) under which the
# Procedure is to appear.
# Possible items are: Populate, Edit, Consistency, Verify, Hazards

# Giving MenuItems a dummy menu name causes it to not appear in any menu
MenuItems = ["Nowhere"]

import SmartScript

## Local configuration definitions.  These values will override the
## values defined in the Diurnal procedure.

configDict = {}

configDict["sourceList"] = ["LLRFS", "ECAM","UWWRF1","HRDPS", "NAM12","NAMNest", "NamDNG", "HIRESWarw", 'HIRESWFV3', "GFS1hr", "SREFBC", "GFS", "ECMWF",
                            "NBM","NBMEXP", "SuperBlend", "CONSAll","CONSRaw", "Obs", "Official", "F-table (Climo)"]
        
# Sources that are based on observations must be listed here
# since they are processed differently than forecast models 
configDict["obsModelList"] = ["Obs", "RTMA"]

# These models will be selected by default in the GUI
configDict["primaryModel"] = "NBM"
configDict["secondaryModel"] = "SuperBlend"


# if no other model guidance is found, use lastResort as guidance
configDict["lastResort"] = "GFS"

# This flag controls whether the tool derives the model min/max from
# the hourly grids (True) or the min/max grids fetched from the model
# If set to False, the tool will use the model min/max grids in the
# database where/when it is lower/higher than the hourly min/max.
configDict["useHourlyMinMax"] = True

# Previous versions and model substitutes will be color coded
# based on this list.  If the first sub was made the first color
# is used, second sub, second color, etc.
configDict["highlightColors"] = ["yellow", # second choice / previous version
                                 "orange", # third choice
                                 "red",    # and so on
                                 "red",
                                 "red",
                                ]

configDict["DEBUG"] = False  # set this to True to see the intermediate grids

# The list of selected elements when the GUI appears
configDict["weList"] = ["T", "RH"]  

configDict["startTime"] = "Current Cycle"
configDict["endTime"] = "M15 D14"

# Set this to 1 if you want Diurnal to adjust for Daylight Saving Time
# Set to 0 for no adjustment.
configDict["adjustDST"] = 0
configDict["timeZone"] = "America/Boise"

# Set this to True if you want the option to exclude particular
# days of Obs-like data.
configDict["displayDaysToExclude"] = False

# The number of days that will appear in the menu to exclude
configDict["daysToExclude"] = 14

# A True value will check if there is an active edit area defined
# and if so, inform the user and ask if they wish to continue
configDict["checkForEditArea"] = True

# The number of grid points below which a dialog will appear
# informing the user that an edit area is selected.  If the
# number of edit area grid points is above this number, no
# dialog will appear.
configDict["editAreaPoints"] = 50

# F-table for the climatology method
configDict["fTable"] = \
    [[ 68, 81, 86, 89, 87, 89, 89, 86, 80, 69, 62, 61], # 0
     [ 53, 62, 65, 72, 72, 75, 75, 69, 57, 51, 50, 50], # 1
     [ 44, 50, 53, 54, 53, 55, 55, 50, 45, 41, 42, 40], # 2
     [ 37, 42, 44, 44, 42, 43, 42, 40, 37, 35, 35, 34], # 3
     [ 31, 36, 36, 36, 34, 35, 34, 33, 30, 30, 30, 29], # 4
     [ 26, 31, 30, 30, 27, 28, 28, 27, 24, 25, 25, 25], # 5
     [ 21, 25, 24, 25, 21, 22, 22, 22, 20, 20, 20, 21], # 6
     [ 17, 20, 18, 20, 16, 16, 17, 19, 16, 17, 16, 18], # 7
     [ 13, 16, 15, 15, 12, 11, 13, 14, 12, 14, 13, 15], # 8
     [ 10, 12, 11, 10,  7,  8,  8, 10,  9, 11, 10, 11], # 9
     [  7,  9,  7,  6,  3,  3,  4,  6,  5,  7,  7,  7], # 10
     [  5,  5,  3,  3,  0,  0,  1,  2,  2,  4,  5,  5], # 11
     [  2,  2,  1,  0,  0,  1,  0,  0,  0,  1,  1,  1], # 12
     [  0,  0,  0,  6, 12, 15, 13, 10,  3,  0,  0,  0], # 13
     [  1,  2, 12, 24, 28, 31, 29, 28, 21, 16,  7,  1], # 14
     [ 15, 20, 30, 41, 44, 47, 46, 46, 41, 37, 26, 19], # 15
     [ 36, 41, 48, 56, 59, 60, 60, 62, 58, 56, 47, 40], # 16
     [ 56, 59, 63, 69, 71, 72, 74, 76, 72, 71, 65, 60], # 17
     [ 73, 75, 75, 81, 80, 81, 84, 86, 83, 83, 80, 77], # 18
     [ 86, 87, 87, 89, 89, 89, 92, 94, 91, 92, 91, 89], # 19
     [ 96, 96, 94, 96, 95, 95, 97, 97, 97, 97, 97, 97], # 20
     [100,100, 98,100, 98,100,100,100,100,100,100,100], # 21
     [ 98,100,100,100,100,100,100,100, 98, 98, 97, 97], # 22
     [ 89, 96, 96, 97, 96, 96, 96, 95, 94, 91, 85, 84]] # 23
##     J   F   M   A   M   J   J   A   S   O   N   D
                 


## The following empty code is here to fool the ifpServer into
## thinking it's a procedure.  This is so it will appear right
## next to the primary procedure.

class Procedure (SmartScript.SmartScript):
    def __init__(self, dbss):
        SmartScript.SmartScript.__init__(self, dbss)

    def execute(self):

        return
