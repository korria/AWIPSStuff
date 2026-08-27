#
#  Configuration for MatchObsAllQC
#
#=============================================================================
#
#  Directory where 'dumped' hourly datafile are kept
#
DATADIR="/localapps/runtime/MatchObsAll/data"
#
#  Base directory of this script - so that we can find the data directory
#  for the location cache
#
BASEDIR="/localapps/runtime/MatchObsAllQC"
#
#
# Path of the geom.dat file
#
GEOMPATH="/localapps/runtime/MatchObsAll/data/geom.dat"
#
#
# CWA ORDER (Order surrounding CWAs)
#
CWAORDER = {"BOI":1,"MFR":2,"PDT":3,"OTX":4,"MSO":5,"PIH":6,"SLC":7,"LKN":8,"REV":9,"TFX":"A"}
#
#
# ZONE ORDER (Same convention as CWA ORDER)
#
ZONEORDER = {"IDZ011":0,"IDZ013":1,"IDZ028":2,"IDZ033":3,"IDZ012":4,"ORZ064":5,"IDZ014":6,"IDZ016":7,"IDZ029":8,"IDZ015":9,"IDZ030":"A","ORZ062":"B","ORZ063":"C","ORZ061":"D"}
#
#
# MatchObsAllDir
#
MOADIR = "/localapps/runtime/MatchObsAll"
