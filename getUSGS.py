#!/awips2/python/bin/python
#
#  getUSGS
#
#  Quick program to get XML formatted 15-minute river data
#  from USGS website.  Decode it, turn it into SHEF, and put
#  it into a product to send over to AWIPS for decoding and
#  storing.
#
#  HISTORY:
#   2024-08-06 - KA - Modernized code, and rewrote to ingest json from USGS
#   2022-07-12 - KA - Updated sites and for HCI/py3.
#   2017-12-11 - Tim Barker - SOO BOI - Original Version
#                             Basically a copy of the getIdahoPower
#                             program
#
#==============================================================================
#
#
#
import requests
import json
from datetime import datetime, timezone
import urllib.request, urllib.parse, urllib.error,time,re,os,sys,calendar
import xml.etree.ElementTree as ET
from io import StringIO
from xml.etree.ElementTree import tostring
#
#
#
#==============================================================================
#
#  logdate - write a message to stdout with a time/date stamp and flush the
#            output so it gets written right away
#
def logdate(message):
    times = datetime.utcnow()
    formatted_time = times.strftime("%Y/%m/%d %H:%M:%S")
    print(f"{formatted_time} : {message}")
    sys.stdout.flush()

#==============================================================================
def getLastOb():
   #
   #  Put in missing times for all stations - in case we cant read the file
   #
   lastob={}
   for stnid in list(STNIDS.keys()):
      lastob[stnid]=0
   #
   #  read info about last observation time at each site
   #
   try:
      inobj=open(LASTFILE,"r")
   except:
      logdate("Could not read last observation times from %s"%LASTFILE)
      return lastob
   inlines=inobj.readlines()
   inobj.close()
   for fullline in inlines:
      line=fullline.strip()
      (stnid,obtime)=line.split(":",2)
      if stnid in list(STNIDS.keys()):
         lastob[stnid]=int(obtime)
   return lastob

# Function to convert local time to UTC and extract components
def convert_to_utc_components(local_time_str):
    # Parse the local time string
    local_time = datetime.strptime(local_time_str, "%Y-%m-%dT%H:%M:%S.%f%z")

    # Convert to UTC
    utc_time = local_time.astimezone(timezone.utc)

    # Extract year, month, day, hour, minute
    year = utc_time.year
    month = utc_time.month
    day = utc_time.day
    hour = utc_time.hour
    minute = utc_time.minute

    # Get seconds since epoch
    epoch_seconds = int(utc_time.timestamp())

    return year, month, day, hour, minute, epoch_seconds
#==============================================================================
#
#  M A I N
#
logdate("Starting getUSGS")
starttime=time.time()
#
#  Associate the Handbook5 ID with the USGS Station Number
#
STNIDS={
         "HAGI1":(13108150,),
         "BWKI1":(13135500,),
         "HALI1":(13139510,),
         "JBGN2":(13162225,),
         "BSCI1":(13206305,),
         "TFAI1":(13092747,),
         "CJSI1":(13171500,),
         "CCXI1":(13159800,),
         "CDWI1":(13211205,),
         "CAMI1":(13258500,),
         "SNYI1":(13213100,),
         "BUHI1":(13094000,),
         "SRPI1":(13093383,),
         "MLDI1":(13152500,),
         "STFI1":(13090500,),
         "SCBI1":(13140800,),
         "WEII1":(13269000,),
         "PARI1":(13213000,),
         "SFRI1":(13106500,),
         "SKHI1":(13154500,),
         "MFPI1":(13237920,),
         "OWYO3":(13183000,),
         "SFMO3":(13215480,),
         "DWSI1":(13089500,),
        }
LASTFILE="/ldad/localapps/getUSGS/lastobs.dat"
OUTFILE="/ldad/localapps/getUSGS/PRODBOIRR5BOI"
FULLFILE="/data/Incoming/PRODBOIRR5BOI"

#
#  Get last obs times
#
lastob=getLastOb()
#
#  Setup blank output
#
outlines=[]
outdata=0
#
#  Loop over configured stations
#
for stnid in list(STNIDS.keys()):
   logdate("Getting data for station %s"%stnid)
   stnnum=STNIDS[stnid][0]
   lastdata=0 # last time of new data for this station
   #
   #  period=P1D for 1 day
   #  period PT6H for 6 hours
   #
   url=f"https://waterservices.usgs.gov/nwis/iv/?format=json&period=PT6H&parameterCd=00065&siteStatus=all&sites={stnnum}"

   # Fetch the data from the URL
   response = requests.get(url)
   data = response.json()

   # Parse the date and values
   parsed_data = []

   # Navigating to the time series data in the JSON structure
   time_series = data['value']['timeSeries'][0]['values'][0]['value']

   # Extracting the date and value
   for entry in time_series:
      date = entry['dateTime']

      year, month, day, hour, minute, epoch_seconds = convert_to_utc_components(date)
      if epoch_seconds>lastob[stnid]:
         if(epoch_seconds>lastdata):
            lastdata=epoch_seconds
         if entry['value'] is not None:
            value = float(entry['value'])  # Convert value to float
            formatted_str = f".A {stnid} {year % 100:02d}{month:02d}{day:02d} Z DH{hour:02d}{minute:02d} /HG {value:8.2f}"
            outlines.append(formatted_str)
            outdata=1

   #
   #  If we got data from this station later than the previous data for
   #  this station - then take the lastest of the new data - and set it
   #  as the last data time for this station
   #
   if lastdata>0:
      lastob[stnid]=lastdata
#
#  If nothing needs to be output - then log it and leave everything as it was
#
if outdata==0:
   logdate("No new data")
   logdate("Finished getUSGS")
   sys.exit(0)
#
#  Put data in OUTFILE
#
(gyea,gmon,gday,ghou,gmin,gsec,gwda,gyda,gdst)=time.gmtime()
logdate("Data sent to AWIPS:")
outobj=open(OUTFILE,"w")
outobj.write("SRUS55 KBOI DDHHMM\n")
outobj.write("RR5BOI\n")
for outline in outlines:
   outobj.write("%s\n"%outline)
   logdate("  %s"%outline)
outobj.close()
os.chmod(OUTFILE,0o666)
time.sleep(1)
#
#  Copy file to FULLFILE
#
os.system("/bin/cp %s %s"%(OUTFILE,FULLFILE))
logdate(" ")
#
#  Store new info about last observation at each site
#
outobj=open(LASTFILE,"w")
for stnid in list(STNIDS.keys()):
   outobj.write("%s:%d\n"%(stnid,lastob[stnid]))
outobj.close()
#
#  Finished
#
logdate("Finished getUSGS")
sys.exit(0)
