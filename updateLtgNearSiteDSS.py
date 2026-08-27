#!/awips2/python/bin/python
#
# Program to copy DSS Events to LtgNearSite
#
# KA - 2026/05/26

import re
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime


# ------------------------------------------------------------
# Hardcoded file locations
# ------------------------------------------------------------

textFilePath = Path("/localapps/runtime/LtgNearSite/etc/events.txt")
xmlFilePath = Path("/localapps/runtime/LtgNearSite/etc/LtgNearSite.xml")

makeBackup = True
pruneMissingDynamic = False


# Static sites that should always stay in the XML
staticSiteNames = {
    "BOI",
    "EUL",
    "BNO",
    "BKE",
    "MYL",
    "ONO",
    "JER",
    "TWF",
}


fieldRegex = re.compile(r"^([A-Za-z]+):\s*(.*)$")


def normalizeName(name: str) -> str:
    return " ".join(name.strip().split())


def parseEventTime(value: str) -> str:
    value = value.strip()

    for dateFormat in ("%m/%d/%Y %H:%M", "%Y/%m/%d %H:%M"):
        try:
            eventDateTime = datetime.strptime(value, dateFormat)
            return eventDateTime.strftime("%Y/%m/%d %H:%M")
        except ValueError:
            pass

    raise ValueError(f"Could not parse date/time: {value}")


def parseTextEvents(textFilePath: Path) -> list[dict]:
    events = []
    currentEvent = {}

    with textFilePath.open("r", encoding="utf-8") as textFile:
        for rawLine in textFile:
            line = rawLine.strip()
            match = fieldRegex.match(line)

            if not match:
                continue

            key, value = match.groups()
            key = key.strip()
            value = value.strip()

            if key == "eventName":
                if currentEvent:
                    events.append(currentEvent)

                currentEvent = {"eventName": value}
            else:
                if currentEvent:
                    currentEvent[key] = value

    if currentEvent:
        events.append(currentEvent)

    parsedEvents = []

    for event in events:
        requiredFields = [
            "eventName",
            "eventBegin",
            "eventEnd",
            "eventLat",
            "eventLon",
        ]

        missingFields = [
            fieldName for fieldName in requiredFields if not event.get(fieldName)
        ]

        if missingFields:
            print(
                f"Skipping incomplete event: "
                f"{event.get('eventName', 'UNKNOWN')} missing {missingFields}"
            )
            continue

        parsedEvents.append(
            {
                "name": normalizeName(event["eventName"]),
                "distance": "10",
                "latitude": event["eventLat"].strip(),
                "longitude": event["eventLon"].strip(),
                "start": parseEventTime(event["eventBegin"]),
                "end": parseEventTime(event["eventEnd"]),
            }
        )

    return parsedEvents


def getChildText(siteElement: ET.Element, tagName: str) -> str:
    childElement = siteElement.find(tagName)

    if childElement is None or childElement.text is None:
        return ""

    return childElement.text.strip()


def isStaticSite(siteElement: ET.Element) -> bool:
    siteName = normalizeName(getChildText(siteElement, "name"))

    if siteName in staticSiteNames:
        return True

    hasStart = siteElement.find("start") is not None
    hasEnd = siteElement.find("end") is not None

    # Any site without start/end is treated as static and preserved.
    if not hasStart and not hasEnd:
        return True

    return False


def setChild(siteElement: ET.Element, tagName: str, value: str) -> None:
    childElement = siteElement.find(tagName)

    if childElement is None:
        childElement = ET.SubElement(siteElement, tagName)

    childElement.text = value


def indentXml(element: ET.Element, level: int = 0) -> None:
    indentText = "\n" + level * "    "

    if len(element):
        if not element.text or not element.text.strip():
            element.text = indentText + "    "

        for childElement in element:
            indentXml(childElement, level + 1)

        if not childElement.tail or not childElement.tail.strip():
            childElement.tail = indentText
    else:
        if level and (not element.tail or not element.tail.strip()):
            element.tail = indentText


def updateXml(events: list[dict]) -> None:
    xmlTree = ET.parse(xmlFilePath)
    rootElement = xmlTree.getroot()

    if rootElement.tag != "data":
        raise ValueError(f"Expected root tag <data>, found <{rootElement.tag}>")

    existingSitesByName = {}

    for siteElement in rootElement.findall("site"):
        siteName = normalizeName(getChildText(siteElement, "name"))

        if siteName:
            existingSitesByName[siteName] = siteElement

    eventNamesFromText = {event["name"] for event in events}

    addedCount = 0
    updatedCount = 0
    skippedStaticCount = 0
    removedDynamicCount = 0

    for event in events:
        eventName = event["name"]
        existingSite = existingSitesByName.get(eventName)

        if existingSite is not None:
            if isStaticSite(existingSite):
                skippedStaticCount += 1
                continue

            setChild(existingSite, "name", event["name"])
            setChild(existingSite, "distance", event["distance"])
            setChild(existingSite, "latitude", event["latitude"])
            setChild(existingSite, "longitude", event["longitude"])
            setChild(existingSite, "start", event["start"])
            setChild(existingSite, "end", event["end"])

            updatedCount += 1

        else:
            siteElement = ET.SubElement(rootElement, "site")

            ET.SubElement(siteElement, "name").text = event["name"]
            ET.SubElement(siteElement, "distance").text = event["distance"]
            ET.SubElement(siteElement, "latitude").text = event["latitude"]
            ET.SubElement(siteElement, "longitude").text = event["longitude"]
            ET.SubElement(siteElement, "start").text = event["start"]
            ET.SubElement(siteElement, "end").text = event["end"]

            addedCount += 1

    if pruneMissingDynamic:
        for siteElement in list(rootElement.findall("site")):
            siteName = normalizeName(getChildText(siteElement, "name"))

            if isStaticSite(siteElement):
                continue

            if siteName not in eventNamesFromText:
                rootElement.remove(siteElement)
                removedDynamicCount += 1

    indentXml(rootElement)
    xmlTree.write(xmlFilePath, encoding="utf-8", xml_declaration=False)

    print(f"Updated XML: {xmlFilePath}")
    print(f"Added: {addedCount}")
    print(f"Updated: {updatedCount}")
    print(f"Skipped static: {skippedStaticCount}")
    print(f"Removed dynamic: {removedDynamicCount}")


def main() -> None:
    if not textFilePath.exists():
        raise FileNotFoundError(f"Text file not found: {textFilePath}")

    if not xmlFilePath.exists():
        raise FileNotFoundError(f"XML file not found: {xmlFilePath}")

    if makeBackup:
        backupFilePath = xmlFilePath.with_suffix(xmlFilePath.suffix + ".bak")
        shutil.copy2(xmlFilePath, backupFilePath)
        print(f"Backup created: {backupFilePath}")

    events = parseTextEvents(textFilePath)
    updateXml(events)


if __name__ == "__main__":
    main()