AWIPS simulated satellite display for RRFS (GOES-West ABI bands 13, 8, 9, 10).
Source fields: SBTA1813, SBTA188, SBTA189, SBTA1810 at level NTAT, decoded from the
NCEP table (center 7, 4.2.3.192 codes 82, 77, 78, 79). getRRFS.py already sends them.

First check they are in AWIPS: Product Browser > Grid > RRFS model > look for SBTA1813.
If they are missing, the grib table needs a site entry; this package will not help.

Copy derivedParameters/definitions/*.xml to SITE.
Paste the *_ADD.xml snippets inside the existing SITE styleRuleset and FieldsMenus files.
In the Volume Browser pick the RRFS model, the Sfc plane, and load as Image.
Fields: SimIR13 (clean IR), SimWV08 (upper WV), SimWV09 (mid WV), SimWV10 (low WV).
IR uses the real GOES IR colormap and range; WV uses the gray WV colormap with the
legacy GOES stretch (145 to -110 C), so both look like real GOES-West imagery.
If RRFS fields are ever trimmed in getRRFS.py, keep SBTA1813, SBTA188, SBTA189, SBTA1810.
