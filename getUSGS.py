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
#   2026-11-23 - KA - Moved to the USGS Water Data API (waterservices.usgs.gov
#                     is being retired). One request for all stations instead
#                     of one per station, with a timeout and retries; a failed
#                     run exits 1 without touching lastobs, so the next run
#                     fills the gap (one connection reset used to crash it).
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
from datetime import datetime, timezone
import os,sys,time
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

#==============================================================================
#
#  getObs - one request for the last 6 hours of gage height (00065) at every
#           station. Returns a list of {monitoring_location_id, time, value}.
#           If USGS can't be reached it exits 1 and leaves lastobs alone, so
#           the next run picks up the gap.
#
def getObs():
   params={
      "f":"json",
      "monitoring_location_id":",".join("USGS-%d"%v[0] for v in STNIDS.values()),
      "parameter_code":"00065",
      "time":"PT6H",
      "skipGeometry":"true",
      "properties":"monitoring_location_id,time,value",
      "limit":10000,
   }
   for attempt in range(1,TRIES+1):
      try:
         response=requests.get(API_URL,params=params,timeout=60)
         response.raise_for_status()   # a 429 here means USGS wants an API key
         data=response.json()
         if any(link.get("rel")=="next" for link in data.get("links",[])):
            logdate("WARNING: USGS response has another page - some data not read")
         return [feature["properties"] for feature in data["features"]]
      except (requests.RequestException,ValueError,KeyError) as e:
         logdate("Attempt %d of %d failed: %s"%(attempt,TRIES,e))
         if attempt<TRIES:
            time.sleep(30*attempt)
   logdate("Could not get data from USGS - will try again next run")
   sys.exit(1)
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
API_URL="https://api.waterdata.usgs.gov/ogcapi/v1/collections/continuous/items"
TRIES=3

#
#  Get last obs times
#
lastob=getLastOb()
#
#  Get the obs, keeping only times newer than the last ob at each station
#
logdate("Getting data for %d stations"%len(STNIDS))
SITES={"USGS-%d"%v[0]:stnid for stnid,v in STNIDS.items()}  # USGS id -> Handbook5 id
newobs={stnid:{} for stnid in STNIDS}  # stnid -> {epoch seconds: (UTC time, value)}
seen=set()
for ob in getObs():
   stnid=SITES.get(ob["monitoring_location_id"])
   if stnid is None:
      continue
   seen.add(stnid)
   if ob["value"] is None:
      continue
   obtime=datetime.fromisoformat(ob["time"]).astimezone(timezone.utc)
   epoch=int(obtime.timestamp())
   if epoch>lastob[stnid]:
      newobs[stnid].setdefault(epoch,(obtime,float(ob["value"])))
missing=[stnid for stnid in STNIDS if stnid not in seen]
if missing:
   logdate("No data from USGS for: %s"%" ".join(missing))
#
#  SHEF lines in station order, oldest first; the newest becomes the last ob
#
outlines=[]
for stnid,obs in newobs.items():
   for epoch in sorted(obs):
      (obtime,value)=obs[epoch]
      outlines.append(f".A {stnid} {obtime:%y%m%d} Z DH{obtime:%H%M} /HG {value:8.2f}")
   if obs:
      lastob[stnid]=max(obs)
#
#  If nothing needs to be output - then log it and leave everything as it was
#
if not outlines:
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
