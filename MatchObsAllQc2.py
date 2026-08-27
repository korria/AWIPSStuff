#!/awips2/python/bin/python
"""
MatchObsAllQC - Weather Observation Quality Control

A GUI program to mark observations as bad in MatchObsAll data obs dump files, based on
Tim Barker's MatchObsAllQC program.

Author: Korri Anderson

Version history, newest first.  (Releases before 2.4 were not dated.)

2026/08/27 - version 2.5. KA. Concurrent-write safety.  Every save rewrites the
             whole file from memory, so a save from another workstation between
             this session's last read and its next write was silently discarded.
             The mtimes were already recorded on every read and write but never
             consulted before overwriting; writeAllFile and writeManFile now
             check them and ask before clobbering, reloading instead if the
             operator declines.  Both files are also written through a temp file
             and swapped in with os.replace, so a crash or a full disk can no
             longer leave a truncated file for other workstations to read.
2026/08/27 - version 2.4. KA. Data integrity: "*" placeholders no longer count
             as overrides, and the .manual readers drop entries the files no
             longer contain, so another operator's deletions stop being written
             back.  CWA filter: "Select None" hides every station, and a rebuild
             keeps the forecaster's selection.  Display: readShapes fills in
             stations that appear after the lat/lon limits settle, geom.dat rows
             are parsed individually, table spans are cleared, "Issues only"
             respects the displayed hour, and sort keys carry the station ID
             beside them so hyphenated IDs survive.  Lifecycle: RunnerThread no
             longer shadows QThread.finished, has a timeout, and is waited on at
             close; the selectors rebuild without firing signals mid-load.
             Overrides must be numeric and non-blank.  Added Stations > Clear
             Ignore-All for Non-Reporting Stations.  Per-station debug output is
             behind MOAQC_DEBUG.
  (undated) - version 2.3. KA. Added support for Adwaita-dark theme.
  (undated) - version 2.2. KA. PySide6 - AWIPS 23.4.2.
  (undated) - version 2.1. KA. Uses PySide2 instead of Tkinter, also allows for
             better sorting features.

"""

import os
import re
import sys
import copy
import time
import shutil
import calendar
import subprocess
import traceback
from pathlib import Path

try:
    from PySide6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QPushButton, QLabel, QScrollArea, QSplitter, QMenu,
        QTableWidget, QTableWidgetItem, QHeaderView, QTabWidget,
        QMessageBox, QComboBox, QCheckBox, QGroupBox, QTextEdit,
        QDialog, QLineEdit, QRadioButton, QButtonGroup, QSpinBox, QGridLayout
    )
    from PySide6.QtCore import Qt, QTimer, Signal, Slot, QSize, QThread
    from PySide6.QtGui import QAction, QFont, QCursor, QColor, QIcon, QPalette
except ImportError as e:
    print(f"Error importing PySide6: {e}")
    sys.exit(1)

# Import from configuration - with robust fallbacks
try:
    from MatchObsAllQC_config import *
    print(f"Loaded configuration. DATADIR={DATADIR}")
except ImportError as e:
    print(f"Warning: Could not import MatchObsAllQc_config: {e}")
    sys.exit(1)

#  Per-station chatter on the hot paths (every redisplay, every QC write) is
#  synchronous stdout that blocks the GUI.  Set MOAQC_DEBUG=1 to turn it back on.
DEBUG = os.environ.get("MOAQC_DEBUG", "") not in ("", "0")


def dbg(*args, **kwargs):
    """print() that only speaks when MOAQC_DEBUG is set."""
    if DEBUG:
        print(*args, **kwargs)


# Make sure the data directory exists
if not os.path.exists(DATADIR):
    print(f"Warning: Data directory {DATADIR} does not exist")
    try:
        os.makedirs(DATADIR)
        print(f"Created data directory: {DATADIR}")
    except Exception as e:
        print(f"Could not create data directory: {e}")

class OverrideDialog(QDialog):
    """Dialog for entering override values"""

    def __init__(self, parent=None, currentValue=""):
        super().__init__(parent)
        self.setWindowTitle("Override Value")
        self.setWindowFlags(self.windowFlags() | Qt.WindowStaysOnTopHint)
        self.resize(250, 120)

        # Setup layout
        layout = QVBoxLayout()
        self.setLayout(layout)

        # Current value display
        currentLayout = QHBoxLayout()
        currentLabel = QLabel("Current:")
        self.currentDisplay = QLabel(currentValue)
        self.currentDisplay.setFrameStyle(QLabel.Panel | QLabel.Sunken)
        currentLayout.addWidget(currentLabel)
        currentLayout.addWidget(self.currentDisplay)
        layout.addLayout(currentLayout)

        # Override value input
        overrideLayout = QHBoxLayout()
        overrideLabel = QLabel("Override:")
        self.overrideInput = QLineEdit()
        self.overrideInput.setMaxLength(6)
        overrideLayout.addWidget(overrideLabel)
        overrideLayout.addWidget(self.overrideInput)
        layout.addLayout(overrideLayout)

        # Buttons
        buttonLayout = QHBoxLayout()
        self.okButton = QPushButton("OK")
        self.cancelButton = QPushButton("Cancel")
        buttonLayout.addWidget(self.okButton)
        buttonLayout.addWidget(self.cancelButton)
        layout.addLayout(buttonLayout)

        # Connect signals
        self.okButton.clicked.connect(self.accept)
        self.cancelButton.clicked.connect(self.reject)

        #  OK stays disabled until something is typed.  An empty string is the
        #  same sentinel the tool uses for "ignore", so OK-on-blank used to turn
        #  an override into an ignore and drop the ob entirely.
        self.okButton.setEnabled(False)
        self.overrideInput.textChanged.connect(
            lambda text: self.okButton.setEnabled(bool(text.strip()))
        )

    def accept(self):
        """Refuse anything MatchObsAll could not read back as a number."""
        text = self.overrideInput.text().strip()
        if not text:
            return

        try:
            float(text)
        except ValueError:
            QMessageBox.warning(
                self, "Invalid Override",
                f"'{text}' is not a number.\n\n"
                f"Override values are written straight into the .manual file "
                f"for MatchObsAll to parse, so they have to be numeric."
            )
            return

        super().accept()

    def getValue(self):
        """Return the entered override value"""
        return self.overrideInput.text().strip()

    @staticmethod
    def getOverride(parent=None, currentValue=""):
        """Static method to create the dialog and return entered value"""
        dialog = OverrideDialog(parent, currentValue)
        result = dialog.exec()

        if result == QDialog.Accepted:
            return dialog.getValue()
        else:
            return None

class StationSearchComboBox(QComboBox):
    """Custom combobox with search functionality for stations"""

    # Add a custom signal for Enter key press
    returnPressed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)

        # Enable search/filter functionality
        self.setEditable(True)
        self.setInsertPolicy(QComboBox.NoInsert)
        self.setMaxVisibleItems(20)  # Show more items in the dropdown

        # Create a line edit with custom behavior to maintain focus
        self.customLineEdit = QLineEdit(self)
        self.setLineEdit(self.customLineEdit)

        # Completely disable automatic completion
        self.setCompleter(None)

        # Store the full list of stations
        self.allStations = []

        # Use a timer to delay filtering to improve typing experience
        self.filterTimer = QTimer()
        self.filterTimer.setSingleShot(True)
        self.filterTimer.setInterval(1000)  # 1000ms delay
        self.filterTimer.timeout.connect(self.executeFilter)

        # Store the current text for delayed filtering
        self.currentFilterText = ""

        # Connect the text changes to the delayed filter
        self.lineEdit().textEdited.connect(self.delayedFilter)

        # Flag to prevent recursive signal handling
        self.ignoreSelectionChanges = False

    def filterItems(self, text):
        """Filter the dropdown items based on the entered text"""
        # Apply an extreme approach to prevent interruptions
        self.ignoreSelectionChanges = True

        # Remember current text and cursor position
        current_text = text
        cursor_pos = self.lineEdit().cursorPosition()

        # First, clear the existing items without updating the view
        self.blockSignals(True)
        self.clear()

        # Check if we should show all items or filter
        if not text:
            # Show all stations
            for station in self.allStations:
                self.addItem(station)
        else:
            # Filter stations that contain the search text
            filteredStations = [station for station in self.allStations if text.upper() in station.upper()]
            for station in filteredStations:
                self.addItem(station)

        # Important: force-restore our text content and cursor position before unblocking signals
        # This prevents the widget from changing the text when items are added
        self.lineEdit().setText(current_text)
        self.lineEdit().setCursorPosition(cursor_pos)

        # Now it's safe to unblock signals
        self.blockSignals(False)

        # Show the dropdown if we have results and text isn't empty
        if self.count() > 0 and text:
            # Showing the popup without stealing focus
            self.showPopup()

        # Reset our flag
        self.ignoreSelectionChanges = False

    def setAllStations(self, stations):
        """Set the full list of stations.

        Rebuilt with signals blocked: clear() plus the first addItem() used to
        emit currentIndexChanged mid-reload and clobber the selected station.
        """
        previous = self.currentText()
        self.allStations = stations

        self.blockSignals(True)
        try:
            self.clear()
            for station in stations:
                self.addItem(station)

            if previous:
                index = self.findText(previous)
                if index >= 0:
                    self.setCurrentIndex(index)
        finally:
            self.blockSignals(False)

    def delayedFilter(self, text):
        """Store the text and restart the filter timer"""
        self.currentFilterText = text
        self.filterTimer.start()

    def executeFilter(self):
        """Execute the filtering after the timer expires"""
        self.filterItems(self.currentFilterText)

    def keyPressEvent(self, event):
        """Handle key press events to improve typing experience"""
        if event.key() == Qt.Key_Return or event.key() == Qt.Key_Enter:
            # Handle Enter/Return key to select a station
            current_text = self.lineEdit().text().strip()

            # Case 1: Exact match to a station ID
            if current_text in self.allStations:
                # Find index if it's in the current dropdown
                for i in range(self.count()):
                    if self.itemText(i) == current_text:
                        self.setCurrentIndex(i)
                        # Emit both signals
                        self.activated.emit(i)
                        self.returnPressed.emit()
                        self.hidePopup()
                        return

                # If we get here, it's not in current dropdown but is a valid station
                # Emit return pressed with the text as is
                self.returnPressed.emit()
                self.hidePopup()
                return

            # Case 2: Match to a currently filtered item
            if self.count() > 0:
                # Select the first item
                self.setCurrentIndex(0)
                self.activated.emit(0)
                self.returnPressed.emit()
                self.hidePopup()
                return

        # Handle arrow keys differently to avoid selection issues during typing
        elif event.key() in (Qt.Key_Up, Qt.Key_Down):
            # Only handle these if the popup is visible
            if self.view().isVisible():
                super().keyPressEvent(event)
                return

        # For all other keys, use default behavior
        super().keyPressEvent(event)

    def focusOutEvent(self, event):
        """Override focus out to check if we should select a station"""
        # Before losing focus, check if current text matches a station
        current_text = self.lineEdit().text().strip()
        if current_text in self.allStations:
            # Find the index and trigger selection
            for i in range(self.count()):
                if self.itemText(i) == current_text:
                    self.setCurrentIndex(i)
                    self.activated.emit(i)
                    break

        super().focusOutEvent(event)

    def showPopup(self):
        """Override showPopup to show popup without affecting typing"""
        # Save current text and position
        current_text = self.lineEdit().text()
        cursor_pos = self.lineEdit().cursorPosition()

        # Show popup
        super().showPopup()

        # Restore text and cursor - prevents popup from changing text
        self.lineEdit().setText(current_text)
        self.lineEdit().setCursorPosition(cursor_pos)

class DataTable(QTableWidget):
    """Custom table widget for displaying observation data"""

    cellRightClicked = Signal(QTableWidgetItem, str)
    cellLeftClicked = Signal(int, int)  # Row, column signals for left clicks

    # Status colors
    COLORS = {
        "plain": QColor(221, 221, 221),       # light gray
        "ignore": QColor(255, 150, 255),      # Light magenta
        "ignoreall": QColor(255, 0, 255),     # Darker magenta
        "override": QColor(255, 0, 0),        # Light red
        "overrideall": QColor(200, 0, 0),     # Darker red
        "header": QColor(200, 200, 200),      # Gray
        "label": QColor(131, 56, 236),        # Blue-gray
        "nodata": QColor(255, 255, 0),        # Yellow - for stations with no data
        "unchanged": QColor(255, 165, 0),     # Orange - for stations with unchanged values
        "sanity": QColor(0, 200, 200),        # Cyan - physically impossible value
        "spike": QColor(255, 200, 0)          # Amber - temporally suspicious value
    }

    def __init__(self, parent=None):
        super().__init__(parent)

        # Set up table appearance
        self.setShowGrid(True)
        self.setAlternatingRowColors(True)

        # Set table to use fixed width columns rather than stretching
        self.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.horizontalHeader().setStretchLastSection(False)

        # Don't show vertical header (row numbers)
        self.verticalHeader().setVisible(False)

        palette = QPalette()
        palette.setColor(QPalette.Base,            QColor("#ccc"))
        palette.setColor(QPalette.AlternateBase,   QColor("#ddd"))
        palette.setColor(QPalette.Text,            QColor("#000"))
        palette.setColor(QPalette.Shadow,          QColor("#222"))
        palette.setColor(QPalette.Highlight,       QColor("#09f"))
        palette.setColor(QPalette.HighlightedText, QColor("#fff"))
 
        self.setPalette(palette)

        # Setup context menu
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.showContextMenu)

        # Allow shift/ctrl click range selection of cells
        self.setSelectionMode(QTableWidget.ExtendedSelection)
        self.setSelectionBehavior(QTableWidget.SelectItems)

        # Connect the cell clicked signal
        self.cellClicked.connect(self.handleCellClick)

        # Set cursor to hand when hovering over rows
        self.setCursor(Qt.PointingHandCursor)

    def showContextMenu(self, point):
        """Show context menu on right-click"""
        item = self.itemAt(point)
        if item is None:
            return

        # Extract field type and status from item data
        fieldType = item.data(Qt.UserRole)
        if not fieldType or not isinstance(fieldType, dict):
            return

        # Emit signal with clicked item and field type
        self.cellRightClicked.emit(item, fieldType.get("field", ""))

    def handleCellClick(self, row, column):
        """Handle regular cell click (left mouse button)"""
        # Emit signal with row and column for the main window to handle
        self.cellLeftClicked.emit(row, column)

class RunnerThread(QThread):
    """Thread for running background processes"""

    #  NOT named "finished": QThread already emits a zero-argument finished()
    #  when run() returns, and shadowing it bound moaFinished to both, so the
    #  slot fired twice or was called with no arguments at all.
    runFinished = Signal(bool, str)  # Success flag and message

    #  Ceiling on the child process so a hung "ssh pv2" cannot pin the thread
    #  forever with the Run button disabled.
    TIMEOUT_SECONDS = 30 * 60

    def __init__(self, command):
        super().__init__()
        self.command = command

    def run(self):
        """Run the command in a separate thread"""
        process = None
        try:
            process = subprocess.Popen(
                self.command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )

            try:
                stdout, stderr = process.communicate(timeout=self.TIMEOUT_SECONDS)
            except subprocess.TimeoutExpired:
                process.kill()
                process.communicate()
                self.runFinished.emit(
                    False,
                    f"Timed out after {self.TIMEOUT_SECONDS // 60} minutes"
                )
                return

            if process.returncode == 0:
                self.runFinished.emit(True, "Success")
            else:
                self.runFinished.emit(False, stderr)

        except Exception as e:
            if process is not None and process.poll() is None:
                process.kill()
            self.runFinished.emit(False, str(e))

class MatchObsAllQC(QMainWindow):
    """Main application window for MatchObsAllQC"""

    def __init__(self):
        super().__init__()

        # Window setup
        self.setWindowTitle("MatchObsAllQC 2.5")
        self.resize(1200, 1200)

        # Initialize data structures
        self.raw = {}          # Raw observation data
        self.man = {}          # Manual corrections
        self.all = {}          # All-station overrides
        self.static = {}       # Static station info
        self.geom = {}         # Geometry data for stations

        # Configuration
        self.currentDatetime = None
        self.currentStation = None
        self.displayMode = "Time"  # 'Time' or 'Station'
        self.dataFiles = []
        self.manFiles = []
        self.allFile = "all.manual"
        self.dates = []
        self.allStations = []
        self.fileTimes = {}
        self.manTimes = {}
        self.allTimes = 0

        # Fixed name length
        self.nameLengthValue = 30

        # CWA filtering
        self.availableCwas = []
        self.selectedCwas = []  # All CWAs selected by default

        # Undo stack: list of (description, restore_callable)
        self.undoStack = []
        self.MAX_UNDO = 50

        # Caches for auto-detected issue flags, recomputed on each redisplay.
        # sanityFlags: {(stationId, datetime): set_of_field_id_strings}
        # spikeFlags:  {(stationId, datetime, field_id): True}
        self.sanityFlags = {}
        self.spikeFlags = {}
        # {datetime: set_of_station_ids} - rebuilt alongside the flags above
        self.flaggedStationsByTime = {}

        # Physical sanity check limits (Fahrenheit / knots / degrees).
        # A value outside [min, max] gets flagged. Keys are field ids.
        self.SANITY_LIMITS = {
            "tmp": {"min": -80,  "max": 130},
            "dew": {"min": -100, "max": 90},
            "spd": {"min": 0,    "max": 200},
            "gst": {"min": 0,    "max": 250},
            "dir": {"min": 0,    "max": 360},
        }

        # Spike thresholds: max plausible change between consecutive obs.
        # A value differing by MORE than this from BOTH temporal neighbors
        # is flagged as a spike (so genuine fronts/wind shifts aren't flagged
        # because their new value persists into the next ob). DIR is omitted
        # because of the 0/360 wraparound.
        self.SPIKE_THRESHOLDS = {
            "tmp": 25,   # °F
            "dew": 25,   # °F
            "spd": 50,   # kt
            "gst": 60,   # kt
        }

        # Current operation context
        self.tagId = None
        self.tagTime = None
        self.tagVal = None

        # Thread for background processing
        self.runnerThread = None

        # Setup for field handling
        self.setupFields()
        self.setupUi()

        # Set up location limits - to be used when reading files
        self.minLat = 100.0
        self.maxLat = -100.0
        self.minLon = 190.0
        self.maxLon = -190.0
        self.changedLatLonLimits = False

        # Try to load data
        try:
            self.statusBar.showMessage("Loading geography information...")
            self.getGeography()

            self.statusBar.showMessage("Loading observation data...")
            self.readAllDataFiles()

            # Display initial data if we have dates
            if self.dates:
                self.displayTime(self.dates[0])
        except Exception as e:
            print(f"Error loading data: {e}")
            print(traceback.format_exc())
            QMessageBox.warning(self, "Data Loading Error",
                                f"Error loading data:\n{str(e)}.")

        self.statusBar.showMessage("Ready")

    def setupFields(self):
        """Setup field definitions for the application"""
        # Define static fields (station metadata)
        self.staticFields = [
            {"id": "id", "label": "ID", "shortLabel": "ID"},
            {"id": "name", "label": "Name", "shortLabel": "Name"},
            {"id": "lat", "label": "Latitude", "shortLabel": "Lat"},
            {"id": "lon", "label": "Longitude", "shortLabel": "Lon"},
            {"id": "elev", "label": "Elevation", "shortLabel": "Elev"}
        ]

        # Define variable fields (observation data)
        self.variableFields = [
            {"id": "tmp", "label": "Temperature", "shortLabel": "TMP"},
            {"id": "dew", "label": "Dewpoint", "shortLabel": "DEW"},
            {"id": "dir", "label": "Wind Direction", "shortLabel": "DIR"},
            {"id": "spd", "label": "Wind Speed", "shortLabel": "SPD"},
            {"id": "gst", "label": "Wind Gust", "shortLabel": "GST"}
        ]

        # Set counts for convenience
        self.numStaticFields = len(self.staticFields)
        #  staticFields[0] is the ID, which is kept as the dictionary key rather
        #  than in the value list, so one fewer attribute is actually stored.
        #  Reading numStaticFields of them collected vals[5] - the first
        #  observation - into self.static as a phantom fifth "static" field.
        self.numStaticAttrs = self.numStaticFields - 1
        self.numVariableFields = len(self.variableFields)
        self.dataFieldStart = 5  # Index where variable data starts in data files

    def setupUi(self):
        """Set up the user interface"""
        # Create central widget
        centralWidget = QWidget()
        self.setCentralWidget(centralWidget)

        # Main layout
        mainLayout = QVBoxLayout(centralWidget)

        # Create toolbar area
        toolbarLayout = QHBoxLayout()

        # Mode selector
        self.modeSelector = QComboBox()
        self.modeSelector.addItem("Time View")
        self.modeSelector.addItem("Station View")
        self.modeSelector.currentIndexChanged.connect(self.changeDisplayMode)
        toolbarLayout.addWidget(QLabel("View Mode:"))
        toolbarLayout.addWidget(self.modeSelector)

        # Time selector (for Time View)
        self.timeLayout = QHBoxLayout()
        self.prevTimeBtn = QPushButton("←")
        self.prevTimeBtn.setToolTip("Previous Time")
        self.prevTimeBtn.clicked.connect(self.prevTime)

        self.timeSelector = QComboBox()
        self.timeSelector.currentIndexChanged.connect(self.timeSelected)

        self.nextTimeBtn = QPushButton("→")
        self.nextTimeBtn.setToolTip("Next Time")
        self.nextTimeBtn.clicked.connect(self.nextTime)

        self.timeLayout.addWidget(self.prevTimeBtn)
        self.timeLayout.addWidget(self.timeSelector)
        self.timeLayout.addWidget(self.nextTimeBtn)
        toolbarLayout.addLayout(self.timeLayout)

        # Station selector (for Station View)
        self.stationLayout = QHBoxLayout()
        self.prevStationBtn = QPushButton("←")
        self.prevStationBtn.setToolTip("Previous Station")
        self.prevStationBtn.clicked.connect(self.prevStation)

        self.stationSelector = StationSearchComboBox()
        self.stationSelector.currentIndexChanged.connect(self.stationSelected)
        self.stationSelector.returnPressed.connect(self.handleStationEnterPressed)

        self.nextStationBtn = QPushButton("→")
        self.nextStationBtn.setToolTip("Next Station")
        self.nextStationBtn.clicked.connect(self.nextStation)

        self.stationInfoLabel = QLabel()

        self.stationLayout.addWidget(self.prevStationBtn)
        self.stationLayout.addWidget(self.stationSelector)
        self.stationLayout.addWidget(self.nextStationBtn)
        self.stationLayout.addWidget(self.stationInfoLabel)
        toolbarLayout.addLayout(self.stationLayout)

        # Hide station controls initially
        self.hideStationControls()

        # Spacer
        toolbarLayout.addStretch(1)

        # Run MOA button
        self.runMoaBtn = QPushButton("Run MatchObsAll")
        self.runMoaBtn.clicked.connect(self.runMoa)
        toolbarLayout.addWidget(self.runMoaBtn)

        mainLayout.addLayout(toolbarLayout)

        # Main options area
        optionsArea = QVBoxLayout()

        # Display options layout
        displayLayout = QHBoxLayout()

        # Display checkboxes
        displayGroup = QGroupBox("Display Options")
        displayInnerLayout = QHBoxLayout(displayGroup)

        self.showNameCb = QCheckBox("Name")
        self.showNameCb.setChecked(True)
        self.showNameCb.stateChanged.connect(self.redisplay)

        self.showLatCb = QCheckBox("Latitude")
        self.showLatCb.stateChanged.connect(self.redisplay)

        self.showLonCb = QCheckBox("Longitude")
        self.showLonCb.stateChanged.connect(self.redisplay)

        self.showElevCb = QCheckBox("Elevation")
        self.showElevCb.setChecked(True)
        self.showElevCb.stateChanged.connect(self.redisplay)

        self.showCwaCb = QCheckBox("CWA")
        self.showCwaCb.stateChanged.connect(self.redisplay)

        # Field checkboxes
        self.fieldCheckboxes = {}
        for field in self.variableFields:
            cb = QCheckBox(field["shortLabel"])
            cb.setChecked(True)
            cb.stateChanged.connect(self.redisplay)
            self.fieldCheckboxes[field["id"]] = cb
            displayInnerLayout.addWidget(cb)

        displayInnerLayout.addWidget(self.showNameCb)
        displayInnerLayout.addWidget(self.showLatCb)
        displayInnerLayout.addWidget(self.showLonCb)
        displayInnerLayout.addWidget(self.showElevCb)
        displayInnerLayout.addWidget(self.showCwaCb)

        # Issues-only filter (Time view only): hide stations that have no
        # overrides, no manual ignores, no unchanged-value flags, and no
        # auto-detected sanity or spike flags.
        self.issuesOnlyCb = QCheckBox("Issues only")
        self.issuesOnlyCb.setToolTip(
            "Time view: show only stations with overrides, ignores, "
            "unchanged values, sanity-check failures, or spikes."
        )
        self.issuesOnlyCb.stateChanged.connect(self.redisplay)
        displayInnerLayout.addWidget(self.issuesOnlyCb)

        displayInnerLayout.addStretch(1)

        displayLayout.addWidget(displayGroup)

        # Sort options
        sortGroup = QGroupBox("Sort By")
        sortLayout = QHBoxLayout(sortGroup)

        self.sortOptions = QButtonGroup()

        sortId = QRadioButton("ID")
        self.sortOptions.addButton(sortId, 6)  # Use 6 as the ID value
        sortLayout.addWidget(sortId)

        sortElev = QRadioButton("Elevation")
        sortElev.setChecked(True)
        self.sortOptions.addButton(sortElev, 5)
        sortLayout.addWidget(sortElev)

        for i, field in enumerate(self.variableFields):
            rb = QRadioButton(field["shortLabel"])
            self.sortOptions.addButton(rb, i)
            sortLayout.addWidget(rb)

        self.sortOptions.buttonClicked.connect(self.redisplay)
        displayLayout.addWidget(sortGroup)

        # Create a group box for sort direction
        sortDirGroup = QGroupBox("Sort Direction")
        sortDirLayout = QHBoxLayout(sortDirGroup)

        self.sortDirOptions = QButtonGroup()

        sortAscending = QRadioButton("Ascending")
        sortDescending = QRadioButton("Descending")
        sortDescending.setChecked(True)  # Default to descending (higher values first)

        self.sortDirOptions.addButton(sortAscending, 0)
        self.sortDirOptions.addButton(sortDescending, 1)

        sortDirLayout.addWidget(sortAscending)
        sortDirLayout.addWidget(sortDescending)

        self.sortDirOptions.buttonClicked.connect(self.redisplay)
        displayLayout.addWidget(sortDirGroup)

        # Font size control
        fontGroup = QGroupBox("Font Size")
        fontLayout = QHBoxLayout(fontGroup)

        self.fontSize = QSpinBox()
        self.fontSize.setRange(10, 24)
        self.fontSize.setValue(18)
        self.fontSize.setSingleStep(1)
        self.fontSize.valueChanged.connect(self.changeFont)
        fontLayout.addWidget(self.fontSize)

        displayLayout.addWidget(fontGroup)

        # Add display options to main options area
        optionsArea.addLayout(displayLayout)

        # CWA Filtering - UPDATED SECTION
        cwaLayout = QHBoxLayout()

        # Create the CWA filter group with exactly the same style as Display Options
        self.cwaFilterGroup = QGroupBox("CWA Filter")

        # Use the same style that Qt applies to the Display Options GroupBox by default
        # This ensures visual consistency between the two sections

        # Create a horizontal layout for the CWA filter group - match Display Options margins
        cwaInnerLayout = QHBoxLayout(self.cwaFilterGroup)

        # Get the default margins that are used in the Display Options group
        # This ensures consistent spacing across both groups
        margins = cwaInnerLayout.contentsMargins()

        # Button container layout for Select All/None buttons
        cwaButtonsLayout = QHBoxLayout()
        cwaButtonsLayout.setContentsMargins(0, 0, 0, 0)

        # Create buttons with the same effective height as checkbox controls in Display Options
        self.selectAllCwasBtn = QPushButton("All")
        self.selectAllCwasBtn.setMaximumWidth(60)
        self.selectAllCwasBtn.clicked.connect(self.selectAllCwas)

        self.selectNoneCwasBtn = QPushButton("None")
        self.selectNoneCwasBtn.setMaximumWidth(60)
        self.selectNoneCwasBtn.clicked.connect(self.selectNoneCwas)

        cwaButtonsLayout.addWidget(self.selectAllCwasBtn)
        cwaButtonsLayout.addWidget(self.selectNoneCwasBtn)
        cwaInnerLayout.addLayout(cwaButtonsLayout)

        # Create a container for checkboxes with better spacing for visibility
        self.cwaCheckboxContainer = QWidget()
        self.cwaCheckboxLayout = QHBoxLayout(self.cwaCheckboxContainer)
        self.cwaCheckboxLayout.setContentsMargins(0, 0, 0, 0)
        self.cwaCheckboxLayout.setSpacing(15)  # Increased spacing between checkboxes
        # Align checkboxes to match height of those in Display Options
        self.cwaCheckboxLayout.setAlignment(Qt.AlignVCenter)

        # The checkboxes will be added dynamically in updateCwaFilter()
        self.cwaCheckboxes = {}

        # Add checkbox container directly to layout without scroll area
        cwaInnerLayout.addWidget(self.cwaCheckboxContainer)

        # Add the CWA filter group to the CWA layout
        cwaLayout.addWidget(self.cwaFilterGroup)

        # Add CWA filter to main options area
        optionsArea.addLayout(cwaLayout)

        # Ensure that the Display Options and CWA Filter have the same height policy
        self.cwaFilterGroup.setSizePolicy(displayGroup.sizePolicy())

        # Add the options area to the main layout
        mainLayout.addLayout(optionsArea)

        # Create a container widget for the table with fixed width
        tableContainer = QWidget()
        tableLayout = QHBoxLayout(tableContainer)
        tableLayout.setContentsMargins(0, 0, 0, 0)

        # Create data table
        self.dataTable = DataTable()
        self.dataTable.cellRightClicked.connect(self.handleCellRightClick)
        self.dataTable.cellLeftClicked.connect(self.handleCellLeftClick)

        # Add table to layout with spacers on either side to center it
        tableLayout.addStretch(1)
        tableLayout.addWidget(self.dataTable)
        tableLayout.addStretch(1)

        # Add the table container to the main layout
        mainLayout.addWidget(tableContainer)

        # Create status bar
        self.statusBar = self.statusBar()

        # Create menus
        self.createMenus()

        # Apply initial font
        self.changeFont()

    def updateCwaFilter(self):
        """Update the CWA filter controls with available CWAs"""
        # Get all available CWAs from geometry data
        availableCwas = set()
        for stationId, geomData in self.geom.items():
            if len(geomData) > 5:  # Make sure we have CWA data
                cwa = geomData[5]
                if cwa and cwa.strip():
                    availableCwas.add(cwa)

        # Convert to sorted list
        self.availableCwas = sorted(list(availableCwas))

        # If we have new CWAs, rebuild the filter
        if set(self.availableCwas) != set(self.cwaCheckboxes.keys()):
            print(f"Updating CWA filter with {len(self.availableCwas)} CWAs")

            #  Remember what was checked.  A rebuild used to tick every box and
            #  overwrite selectedCwas, so the first station arriving from a new
            #  CWA threw away the forecaster's filter mid-QC.
            previousCwas = set(self.cwaCheckboxes.keys())
            previouslyChecked = set(
                cwa for cwa, box in self.cwaCheckboxes.items() if box.isChecked()
            )

            #  Empty the layout properly.  deleteLater() on widgets alone left the
            #  addStretch() spacer behind and a new one was appended each rebuild.
            while self.cwaCheckboxLayout.count():
                item = self.cwaCheckboxLayout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()

            self.cwaCheckboxes = {}

            # Create new checkboxes in a horizontal layout
            for cwa in self.availableCwas:
                checkbox = QCheckBox(cwa)
                #  Keep the previous state.  A CWA that is new since the last
                #  build (including every CWA on the first build) starts checked.
                checkbox.setChecked(cwa not in previousCwas
                                    or cwa in previouslyChecked)
                checkbox.stateChanged.connect(self.onCwaFilterChanged)

                # Make the checkbox more compact
                checkbox.setStyleSheet("QCheckBox { padding: 0px; margin: 0px; }")

                # Add to the horizontal layout
                self.cwaCheckboxLayout.addWidget(checkbox)
                self.cwaCheckboxes[cwa] = checkbox

            # Add a stretch at the end to prevent spreading checkboxes too far apart
            self.cwaCheckboxLayout.addStretch(1)

            # Update selected CWAs from the boxes we just restored
            self.updateSelectedCwas()

    def onCwaFilterChanged(self):
        """Handle changes to CWA filter checkboxes"""
        self.updateSelectedCwas()
        self.redisplay()

    def updateSelectedCwas(self):
        """Update the list of selected CWAs based on checkbox states"""
        self.selectedCwas = []
        for cwa, checkbox in self.cwaCheckboxes.items():
            if checkbox.isChecked():
                self.selectedCwas.append(cwa)

        # Print the selected CWAs for debugging
        print(f"Selected CWAs: {self.selectedCwas}")

    def selectAllCwas(self):
        """Select all CWAs in the filter"""
        for cwa, checkbox in self.cwaCheckboxes.items():
            checkbox.setChecked(True)
        self.updateSelectedCwas()
        self.redisplay()

    def selectNoneCwas(self):
        """Deselect all CWAs in the filter"""
        for cwa, checkbox in self.cwaCheckboxes.items():
            checkbox.setChecked(False)
        self.updateSelectedCwas()
        self.redisplay()

    def createMenus(self):
        """Create application menus"""
        # File menu
        fileMenu = self.menuBar().addMenu("&File")
        exitAction = QAction("E&xit", self)
        exitAction.triggered.connect(self.close)
        fileMenu.addAction(exitAction)

        # Edit menu
        editMenu = self.menuBar().addMenu("&Edit")
        self.undoAction = QAction("&Undo Last QC Action", self)
        self.undoAction.setShortcut("Ctrl+Z")
        self.undoAction.setToolTip(
            "Revert the most recent ignore/override/cancel action."
        )
        self.undoAction.triggered.connect(self.undoLast)
        self.undoAction.setEnabled(False)
        editMenu.addAction(self.undoAction)

        # View menu
        viewMenu = self.menuBar().addMenu("&View")
        refreshAction = QAction("&Refresh Data", self)
        refreshAction.triggered.connect(self.refreshData)
        viewMenu.addAction(refreshAction)

        # Stations menu
        stationsMenu = self.menuBar().addMenu("&Stations")
        ignoreNonReportingAction = QAction(
            "Ignore Non-Reporting Stations…", self
        )
        ignoreNonReportingAction.setToolTip(
            "Find stations with zero raw observations across the loaded "
            "period and mark every variable as ignore-all."
        )
        ignoreNonReportingAction.triggered.connect(
            self.ignoreNonReportingStations
        )
        stationsMenu.addAction(ignoreNonReportingAction)

        clearIgnoreAllAction = QAction(
            "Clear Ignore-All for Non-Reporting Stations...", self
        )
        clearIgnoreAllAction.setToolTip(
            "Remove the all.manual ignore-all entry for stations with zero raw "
            "observations across the loaded period.  Per-time QC is untouched."
        )
        clearIgnoreAllAction.triggered.connect(
            self.clearIgnoreAllForNonReporting
        )
        stationsMenu.addAction(clearIgnoreAllAction)

        # Run menu
        runMenu = self.menuBar().addMenu("&Run")
        runMoaAction = QAction("&Run MatchObsAll", self)
        runMoaAction.triggered.connect(self.runMoa)
        runMenu.addAction(runMoaAction)

    def changeFont(self):
        """Change the font size in the data table"""
        font = QFont("Carlito", self.fontSize.value())
        self.dataTable.setFont(font)
        self.redisplay()

    def changeDisplayMode(self, index):
        """Change between Time and Station view modes"""
        if index == 0:  # Time View
            self.displayMode = "Time"
            self.showTimeControls()
            self.hideStationControls()
            if self.dates:
                self.displayTime(self.currentDatetime or self.dates[0])
        else:  # Station View
            self.displayMode = "Station"
            self.hideTimeControls()
            self.showStationControls()
            if self.allStations:
                self.displayStation(self.currentStation or self.allStations[0])

    def showTimeControls(self):
        """Show time selection controls"""
        for i in range(self.timeLayout.count()):
            item = self.timeLayout.itemAt(i)
            if item and item.widget():
                item.widget().show()

    def hideTimeControls(self):
        """Hide time selection controls"""
        for i in range(self.timeLayout.count()):
            item = self.timeLayout.itemAt(i)
            if item and item.widget():
                item.widget().hide()

    def showStationControls(self):
        """Show station selection controls"""
        for i in range(self.stationLayout.count()):
            item = self.stationLayout.itemAt(i)
            if item and item.widget():
                item.widget().show()

    def hideStationControls(self):
        """Hide station selection controls"""
        for i in range(self.stationLayout.count()):
            item = self.stationLayout.itemAt(i)
            if item and item.widget():
                item.widget().hide()

    def timeSelected(self, index):
        """Handle time selection from dropdown"""
        if index >= 0 and index < self.timeSelector.count():
            date = self.timeSelector.itemData(index)
            self.displayTime(date)

    def stationSelected(self, index):
        """Handle station selection from dropdown"""
        if index >= 0 and index < self.stationSelector.count():
            station = self.stationSelector.itemText(index)
            self.displayStation(station)

    def nextTime(self):
        """Move to next (more recent) time"""
        if not self.dates:
            return

        if self.currentDatetime in self.dates:
            index = self.dates.index(self.currentDatetime)
            if index > 0:
                self.displayTime(self.dates[index - 1])
        else:
            self.displayTime(self.dates[0])

    def prevTime(self):
        """Move to previous (less recent) time"""
        if not self.dates:
            return

        if self.currentDatetime in self.dates:
            index = self.dates.index(self.currentDatetime)
            if index < len(self.dates) - 1:
                self.displayTime(self.dates[index + 1])
        else:
            self.displayTime(self.dates[0])

    def nextStation(self):
        """Move to next station alphabetically"""
        if not self.allStations:
            return

        sortedStations = sorted(self.allStations)
        if self.currentStation in sortedStations:
            index = sortedStations.index(self.currentStation)
            if index < len(sortedStations) - 1:
                nextStation = sortedStations[index + 1]
                self.displayStation(nextStation)

                # Update the station selector to match
                selectorIndex = self.stationSelector.findText(nextStation)
                if selectorIndex >= 0:
                    self.stationSelector.blockSignals(True)
                    self.stationSelector.setCurrentIndex(selectorIndex)
                    self.stationSelector.blockSignals(False)

                # Update the text in the line edit directly
                self.stationSelector.lineEdit().setText(nextStation)

                # Maintain focus on the selector
                self.stationSelector.lineEdit().setFocus()
        else:
            firstStation = sortedStations[0]
            self.displayStation(firstStation)

            # Update selector and maintain focus
            selectorIndex = self.stationSelector.findText(firstStation)
            if selectorIndex >= 0:
                self.stationSelector.blockSignals(True)
                self.stationSelector.setCurrentIndex(selectorIndex)
                self.stationSelector.blockSignals(False)

            self.stationSelector.lineEdit().setText(firstStation)
            self.stationSelector.lineEdit().setFocus()

    def prevStation(self):
        """Move to previous station alphabetically"""
        if not self.allStations:
            return

        sortedStations = sorted(self.allStations)
        if self.currentStation in sortedStations:
            index = sortedStations.index(self.currentStation)
            if index > 0:
                prevStation = sortedStations[index - 1]
                self.displayStation(prevStation)

                # Update the station selector to match
                selectorIndex = self.stationSelector.findText(prevStation)
                if selectorIndex >= 0:
                    self.stationSelector.blockSignals(True)
                    self.stationSelector.setCurrentIndex(selectorIndex)
                    self.stationSelector.blockSignals(False)

                # Update the text in the line edit directly
                self.stationSelector.lineEdit().setText(prevStation)

                # Maintain focus on the selector
                self.stationSelector.lineEdit().setFocus()
        else:
            # If no current station, start with the first one
            if sortedStations:
                firstStation = sortedStations[0]
                self.displayStation(firstStation)

                # Update selector and maintain focus
                selectorIndex = self.stationSelector.findText(firstStation)
                if selectorIndex >= 0:
                    self.stationSelector.blockSignals(True)
                    self.stationSelector.setCurrentIndex(selectorIndex)
                    self.stationSelector.blockSignals(False)

                self.stationSelector.lineEdit().setText(firstStation)
                self.stationSelector.lineEdit().setFocus()

    def handleStationEnterPressed(self):
        """Handle when Enter is pressed in the station selector"""
        # Get the current text
        entered_text = self.stationSelector.lineEdit().text().strip()

        # Check if it's a valid station
        if entered_text in self.allStations:
            print(f"Enter pressed with valid station: {entered_text}")
            self.displayStation(entered_text)
        else:
            # If dropdown has items, select the first one
            if self.stationSelector.count() > 0:
                station = self.stationSelector.itemText(0)
                print(f"Enter pressed, selecting first match: {station}")
                self.displayStation(station)

                # Update the text field to show what was selected
                self.stationSelector.lineEdit().setText(station)

        # Return focus to the station selector
        self.stationSelector.lineEdit().setFocus()

    def redisplay(self):
        """Redisplay current data with updated settings"""
        if self.displayMode == "Time":
            self.displayTime(self.currentDatetime)
        else:
            self.displayStation(self.currentStation)

    def refreshData(self):
        """Reload data from files"""
        try:
            self.statusBar.showMessage("Refreshing data...")
            self.readAllDataFiles()
            self.redisplay()
            self.statusBar.showMessage("Data refreshed", 3000)
        except Exception as e:
            print(f"Error refreshing data: {e}")
            print(traceback.format_exc())
            QMessageBox.warning(self, "Data Refresh Error",
                               f"Error refreshing data:\n{str(e)}")
            self.statusBar.showMessage("Error refreshing data", 3000)

    def handleCellLeftClick(self, row, column):
        """Handle left-click on a cell to switch views"""
        try:
            # Skip if clicked on a header row or empty row
            if row < 0 or column < 0:
                return

            # If Shift/Ctrl is held, the user is extending a selection — don't
            # hijack the click to switch views. Let Qt handle the selection.
            mods = QApplication.keyboardModifiers()
            if mods & (Qt.ShiftModifier | Qt.ControlModifier):
                return

            # Different behavior based on current display mode
            if self.displayMode == "Time":
                # In time view, clicking a row switches to station view for that station
                stationIdItem = self.dataTable.item(row, 0)
                if stationIdItem and not stationIdItem.text().startswith(" "):  # Avoid header rows
                    stationId = stationIdItem.text()

                    # Check if this is actually a station ID (not a header)
                    if stationId in self.allStations:
                        # Store the current time before switching views

                        # Switch to station view
                        self.modeSelector.setCurrentIndex(1)  # Trigger change to station view
                        self.displayStation(stationId)
            else:
                # In station view, clicking a row switches to time view for that time
                timeItem = self.dataTable.item(row, 0)
                if timeItem:
                    datetime = timeItem.data(Qt.UserRole)
                    if datetime and datetime in self.dates:
                        # Switch to time view
                        self.modeSelector.setCurrentIndex(0)  # Trigger change to time view
                        self.displayTime(datetime)
        except Exception as e:
            print(f"Error handling cell click: {e}")

    def handleCellRightClick(self, item, fieldType):
        """Handle right-click on a cell to show context menu"""
        if not item or not fieldType:
            return

        # Extract row and column
        row = item.row()
        col = item.column()

        # Get the id and datetime based on display mode
        if self.displayMode == "Time":
            # In time view, we need the station ID from the first column
            stationIdItem = self.dataTable.item(row, 0)
            if not stationIdItem:
                return

            self.tagId = stationIdItem.text()
            self.tagTime = self.currentDatetime
        else:
            # In station view, we need the datetime from the first column
            timeItem = self.dataTable.item(row, 0)
            if not timeItem:
                return

            self.tagTime = timeItem.data(Qt.UserRole)
            self.tagId = self.currentStation

        # Set the field type
        self.tagVal = fieldType

        # Determine the status
        status = item.data(Qt.UserRole).get("status", "plain")

        # Create the context menu
        contextMenu = QMenu(self)

        # --- Multi-cell selection actions (when 2+ variable cells selected) ---
        # These appear at the top so they don't get lost below per-cell actions.
        selectedCells = self._selectedVariableCells()
        if len(selectedCells) > 1:
            ignoreSelectedAction = QAction(
                f"Ignore {len(selectedCells)} selected cell(s)", self
            )
            ignoreSelectedAction.triggered.connect(self.ignoreSelectedCells)
            contextMenu.addAction(ignoreSelectedAction)

            # Count how many of the selected cells actually have a per-time entry
            clearable = sum(
                1 for _i, sid, dt, fid in selectedCells
                if dt in self.getList(fid, sid)
            )
            if clearable > 0:
                cancelSelectedAction = QAction(
                    f"Cancel ignore/override on {clearable} of "
                    f"{len(selectedCells)} selected cell(s)", self
                )
                cancelSelectedAction.triggered.connect(
                    self.cancelIgnoreSelectedCells
                )
                contextMenu.addAction(cancelSelectedAction)

            contextMenu.addSeparator()

        if status == "plain" or status == "unchanged" or status == "sanity" or status == "spike":
            # For plain data, offer ignore and override options
            ignoreAction = QAction("Ignore this value", self)
            ignoreAction.triggered.connect(self.ignore)
            contextMenu.addAction(ignoreAction)

            ignoreAllAction = QAction("Ignore all values for this station (this variable)", self)
            ignoreAllAction.triggered.connect(self.ignoreAll)
            contextMenu.addAction(ignoreAllAction)

            contextMenu.addSeparator()

            overrideAction = QAction("Override this value", self)
            overrideAction.triggered.connect(self.override)
            contextMenu.addAction(overrideAction)

        elif status == "ignore":
            # For ignored data
            cancelAction = QAction("Cancel ignore", self)
            cancelAction.triggered.connect(self.cancelIgnore)
            contextMenu.addAction(cancelAction)

        elif status == "ignoreall":
            # For ignoreall data
            cancelAction = QAction("Cancel ignore all (this variable)", self)
            cancelAction.triggered.connect(self.cancelIgnoreAll)
            contextMenu.addAction(cancelAction)

        elif status == "override":
            # For overridden data
            cancelAction = QAction("Cancel override", self)
            cancelAction.triggered.connect(self.cancelOverride)
            contextMenu.addAction(cancelAction)

        elif status == "overrideall":
            # For overrideall data
            cancelAction = QAction("Cancel override all", self)
            cancelAction.triggered.connect(self.cancelOverrideAll)
            contextMenu.addAction(cancelAction)

        # --- Station-wide actions (always available, regardless of cell status) ---
        contextMenu.addSeparator()

        ignoreAllObsAction = QAction("Ignore ALL obs for this station (every variable)", self)
        ignoreAllObsAction.triggered.connect(lambda: self.ignoreAllObsForStation())
        contextMenu.addAction(ignoreAllObsAction)

        # Only offer the cancel if the station already has any ignore-all set
        if self.hasAnyIgnoreAll(self.tagId):
            cancelAllObsAction = QAction("Cancel ignore-all for this station (every variable)", self)
            cancelAllObsAction.triggered.connect(lambda: self.cancelIgnoreAllObsForStation())
            contextMenu.addAction(cancelAllObsAction)

        # Show the context menu
        contextMenu.exec(QCursor.pos())

    def ignore(self):
        """Handle the 'Ignore this value' action"""
        dbg(f"Ignoring value for field {self.tagVal}, station {self.tagId}, "
            f"time {self.tagTime}")

        #  Validate the field id before touching anything.  This used to build
        #  an index that nothing downstream read; only the guard mattered.
        if self.tagVal not in [field["id"] for field in self.variableFields]:
            print(f"Error: Could not find field {self.tagVal}")
            return

        # Snapshot for undo BEFORE mutating
        sid, fid, dt = self.tagId, self.tagVal, self.tagTime
        snap = self._snapshotManCell(sid, fid, dt)

        def restore():
            self._restoreManCell(snap, sid, fid, dt)
            self.writeManFile(dt)

        self._pushUndo(f"Ignore {fid} at {sid}", restore)

        # Get existing field list or create a new one
        fieldList = self.getList(self.tagVal, self.tagId)

        # Set the value to empty string (ignore)
        fieldList[self.tagTime] = ""

        # Update the data structure
        self.setList(fieldList, self.tagVal, self.tagId)

        # Write to file and redisplay
        self.writeManFile(self.tagTime)
        self.redisplay()

    def ignoreAll(self):
        """Handle the 'Ignore all values for this station' action"""
        # Snapshot for undo BEFORE mutating
        prev = self.getAll(self.tagVal, self.tagId)
        sid, fid = self.tagId, self.tagVal

        def restore():
            self.setAll(prev, fid, sid)
            self.writeAllFile()

        self._pushUndo(f"Ignore all {fid} for {sid}", restore)

        self.setAll("", self.tagVal, self.tagId)
        self.writeAllFile()
        self.redisplay()

    def override(self):
        """Handle the 'Override this value' action"""
        # Get the current value
        fieldList = self.getList(self.tagVal, self.tagId)
        currentValue = ""
        if self.tagTime in fieldList:
            currentValue = fieldList[self.tagTime]

        # Ask for the override value
        ask = OverrideDialog.getOverride(self, currentValue)
        if ask is not None:
            sid, fid, dt = self.tagId, self.tagVal, self.tagTime

            # Snapshot BEFORE mutating
            snap = self._snapshotManCell(sid, fid, dt)

            def restore():
                self._restoreManCell(snap, sid, fid, dt)
                self.writeManFile(dt)

            self._pushUndo(f"Override {fid} at {sid} → {ask}", restore)

            override = ask
            fieldList = self.getList(self.tagVal, self.tagId)
            fieldList[self.tagTime] = override
            self.setList(fieldList, self.tagVal, self.tagId)
            self.writeManFile(self.tagTime)
        self.redisplay()

    def cancelIgnore(self):
        """Cancel an ignore by removing the entry"""
        sid, fid, dt = self.tagId, self.tagVal, self.tagTime
        fieldList = self.getList(fid, sid)
        if dt in fieldList:
            # Snapshot BEFORE mutating
            snap = self._snapshotManCell(sid, fid, dt)

            def restore():
                self._restoreManCell(snap, sid, fid, dt)
                self.writeManFile(dt)

            self._pushUndo(f"Cancel ignore {fid} at {sid}", restore)

            del fieldList[dt]
            self.setList(fieldList, fid, sid)
            self.writeManFile(dt)
        self.redisplay()

    def _resetAllEntry(self, description):
        """Reset this cell's all.manual value to '*' (no override), with undo.

        Shared by cancelIgnoreAll and cancelOverrideAll, which had identical
        bodies and differed only in the undo description.
        """
        sid, fid = self.tagId, self.tagVal
        prev = self.getAll(fid, sid)

        def restore():
            self.setAll(prev, fid, sid)
            self.writeAllFile()

        self._pushUndo(description, restore)

        self.setAll("*", fid, sid)
        self.writeAllFile()
        self.redisplay()

    def cancelIgnoreAll(self):
        """Cancel ignoring all by resetting to '*'"""
        self._resetAllEntry(f"Cancel ignore-all {self.tagVal} for {self.tagId}")

    def cancelOverride(self):
        """Cancel an override by removing the entry"""
        sid, fid, dt = self.tagId, self.tagVal, self.tagTime
        fieldList = self.getList(fid, sid)
        if dt in fieldList:
            snap = self._snapshotManCell(sid, fid, dt)

            def restore():
                self._restoreManCell(snap, sid, fid, dt)
                self.writeManFile(dt)

            self._pushUndo(f"Cancel override {fid} at {sid}", restore)

            del fieldList[dt]
            self.setList(fieldList, fid, sid)
            self.writeManFile(dt)
        self.redisplay()

    def cancelOverrideAll(self):
        """Cancel overriding all by resetting to '*'"""
        self._resetAllEntry(f"Cancel override all {self.tagVal} for {self.tagId}")

    # ------------------------------------------------------------------
    # Backups
    # ------------------------------------------------------------------

    def _backupFile(self, path):
        """Copy `path` to `path.bak` if it exists. Best-effort: failure to
        back up should never block writing the new state."""
        try:
            if os.path.exists(path):
                shutil.copy2(path, path + ".bak")
        except Exception as e:
            print(f"Warning: could not back up {path}: {e}")

    # ------------------------------------------------------------------
    # Undo stack
    # ------------------------------------------------------------------

    def _pushUndo(self, description, restore_fn):
        """Push an undo record onto the stack. `restore_fn` is a callable
        with no args that reverts the relevant state and writes any files
        that need to change. The caller is responsible for capturing the
        prior state in the closure BEFORE mutating."""
        self.undoStack.append((description, restore_fn))
        if len(self.undoStack) > self.MAX_UNDO:
            self.undoStack.pop(0)
        if hasattr(self, "undoAction"):
            self.undoAction.setEnabled(True)
            self.undoAction.setText(f"&Undo: {description}")

    def undoLast(self):
        """Pop and apply the most recent undo record."""
        if not self.undoStack:
            return
        description, restore_fn = self.undoStack.pop()
        try:
            restore_fn()
        except Exception as e:
            print(f"Error during undo of '{description}': {e}")
            print(traceback.format_exc())
            QMessageBox.warning(
                self, "Undo Error",
                f"Error undoing '{description}':\n{str(e)}"
            )
        self.redisplay()
        self.statusBar.showMessage(f"Undid: {description}", 4000)
        if hasattr(self, "undoAction"):
            if self.undoStack:
                next_desc = self.undoStack[-1][0]
                self.undoAction.setText(f"&Undo: {next_desc}")
                self.undoAction.setEnabled(True)
            else:
                self.undoAction.setText("&Undo Last QC Action")
                self.undoAction.setEnabled(False)

    def _snapshotManCell(self, stationId, fieldId, datetime):
        """Capture current state of a single self.man cell. Returns a tuple
        (present_bool, value_or_None) suitable for restoration."""
        fieldList = self.getList(fieldId, stationId)
        if datetime in fieldList:
            return (True, fieldList[datetime])
        return (False, None)

    def _restoreManCell(self, snapshot, stationId, fieldId, datetime):
        """Restore a single self.man cell from a snapshot tuple."""
        fieldList = self.getList(fieldId, stationId)
        if snapshot[0]:
            fieldList[datetime] = snapshot[1]
        else:
            fieldList.pop(datetime, None)
        self.setList(fieldList, fieldId, stationId)

    # ------------------------------------------------------------------
    # Auto-detected issue flags (sanity + spike)
    # ------------------------------------------------------------------

    def _effectiveValue(self, stationId, fieldIndex, datetime):
        """Return the current effective string value for a cell, applying
        manual and all-station overrides on top of the raw value. Returns
        None if the cell is empty, ignored, or ignore-all'd."""
        raw_str = ""
        if stationId in self.raw and fieldIndex < len(self.raw[stationId]):
            rawList = self.raw[stationId][fieldIndex]
            if datetime in rawList:
                raw_str = rawList[datetime].strip()

        effective = raw_str if raw_str else None

        # Manual override?
        if stationId in self.man and fieldIndex < len(self.man[stationId]):
            manList = self.man[stationId][fieldIndex]
            if datetime in manList:
                mv = manList[datetime]
                if mv == "":
                    return None  # user-ignored
                elif mv != "*":
                    effective = mv

        # All-station override?
        if stationId in self.all and fieldIndex < len(self.all[stationId]):
            av = self.all[stationId][fieldIndex]
            if av == "":
                return None  # user ignore-all'd
            elif av != "*":
                effective = av

        return effective

    def _computeSanityFlags(self):
        """Compute physically-impossible-value flags across all stations.
        Returns dict {(stationId, datetime): set(field_id_string)}."""
        flags = {}
        field_idx = {f["id"]: i for i, f in enumerate(self.variableFields)}
        limits = self.SANITY_LIMITS

        for stationId in list(self.raw.keys()):
            # Collect every datetime that has any obs for this station
            all_dts = set()
            for fi in range(self.numVariableFields):
                if fi < len(self.raw[stationId]):
                    all_dts.update(self.raw[stationId][fi].keys())

            for dt in all_dts:
                # Pull effective values for each field; numeric only
                vals = {}
                for fid, fi in field_idx.items():
                    eff = self._effectiveValue(stationId, fi, dt)
                    if eff is None:
                        continue
                    try:
                        vals[fid] = float(eff)
                    except ValueError:
                        continue

                cell_flags = set()

                # Range checks
                for fid, v in vals.items():
                    if fid in limits:
                        lim = limits[fid]
                        if v < lim["min"] or v > lim["max"]:
                            cell_flags.add(fid)

                # Cross-field: dewpoint cannot exceed temperature
                if "tmp" in vals and "dew" in vals:
                    if vals["dew"] > vals["tmp"]:
                        cell_flags.add("dew")

                # Cross-field: gust must be >= sustained (when gust present)
                if "spd" in vals and "gst" in vals:
                    if vals["gst"] > 0 and vals["gst"] < vals["spd"]:
                        cell_flags.add("gst")

                if cell_flags:
                    flags[(stationId, dt)] = cell_flags

        return flags

    def _computeSpikeFlags(self):
        """Compute temporal spike flags. A value is flagged when it differs
        from BOTH temporal neighbors by more than the field's threshold —
        this avoids flagging real frontal passages where the new value
        persists. Returns dict {(stationId, datetime, field_id): True}."""
        flags = {}
        if not self.dates or len(self.dates) < 3:
            return flags

        # self.dates is stored newest-first; sort ascending for sequence work
        sorted_dates = sorted(self.dates)
        field_idx = {f["id"]: i for i, f in enumerate(self.variableFields)}

        for stationId in list(self.raw.keys()):
            for fid, fi in field_idx.items():
                if fid not in self.SPIKE_THRESHOLDS:
                    continue
                threshold = self.SPIKE_THRESHOLDS[fid]

                # Build temporal sequence of (dt, float_value) for cells
                # that have an effective value (skipping ignored / missing)
                seq = []
                for dt in sorted_dates:
                    eff = self._effectiveValue(stationId, fi, dt)
                    if eff is None:
                        continue
                    try:
                        seq.append((dt, float(eff)))
                    except ValueError:
                        continue

                if len(seq) < 3:
                    continue

                for i in range(1, len(seq) - 1):
                    _, prev_v = seq[i - 1]
                    cur_dt, cur_v = seq[i]
                    _, next_v = seq[i + 1]

                    if (abs(cur_v - prev_v) > threshold
                            and abs(cur_v - next_v) > threshold):
                        flags[(stationId, cur_dt, fid)] = True

        return flags

    def _computeIssueFlags(self):
        """Recompute self.sanityFlags and self.spikeFlags. Called at the
        start of each redisplay so flags reflect current overrides."""
        try:
            self.sanityFlags = self._computeSanityFlags()
        except Exception as e:
            print(f"Error computing sanity flags: {e}")
            self.sanityFlags = {}
        try:
            self.spikeFlags = self._computeSpikeFlags()
        except Exception as e:
            print(f"Error computing spike flags: {e}")
            self.spikeFlags = {}

        #  Index the flags by time so the "Issues only" filter can ask whether a
        #  station is flagged at the displayed hour in O(1).  It used to scan
        #  both dictionaries once per station, which is quadratic and answered
        #  the wrong question - a 06Z spike kept a station on screen at 18Z.
        self.flaggedStationsByTime = {}
        for (stationId, flagTime) in self.sanityFlags:
            self.flaggedStationsByTime.setdefault(flagTime, set()).add(stationId)
        for (stationId, flagTime, _fieldId) in self.spikeFlags:
            self.flaggedStationsByTime.setdefault(flagTime, set()).add(stationId)

    def _autoStatusForCell(self, stationId, datetime, fieldId, baseStatus):
        """Given a cell's currently-determined base status, upgrade to
        'sanity' or 'spike' if appropriate. Auto-detected statuses only
        apply when the cell would otherwise be 'plain' or 'unchanged'."""
        if baseStatus not in ("plain", "unchanged"):
            return baseStatus
        cell_key = (stationId, datetime)
        if cell_key in self.sanityFlags and fieldId in self.sanityFlags[cell_key]:
            return "sanity"
        if (stationId, datetime, fieldId) in self.spikeFlags:
            return "spike"
        return baseStatus

    # ------------------------------------------------------------------
    # Issues-only filter (Time view)
    # ------------------------------------------------------------------

    def _stationHasIssues(self, stationId, unchangedValuesMap=None, datetime=None):
        """Return True if a station has any of:
        - manual per-time overrides/ignores in self.man
        - all-station overrides/ignores in self.all
        - unchanged-value flags (passed in from display loop)
        - any sanity-flag entries
        - any spike-flag entries

        When datetime is given only that hour counts, because the filter is
        about the row actually on screen.  Ignore-alls in self.all are
        time-independent and always count.
        """
        if stationId in self.man:
            for fieldList in self.man[stationId]:
                if datetime is None:
                    if fieldList:
                        return True
                elif datetime in fieldList:
                    return True

        if stationId in self.all:
            for v in self.all[stationId]:
                if v != "*":
                    return True

        if unchangedValuesMap and stationId in unchangedValuesMap:
            return True

        if datetime is not None:
            return stationId in self.flaggedStationsByTime.get(datetime, ())

        for (sid, _dt) in self.sanityFlags.keys():
            if sid == stationId:
                return True

        for (sid, _dt, _fid) in self.spikeFlags.keys():
            if sid == stationId:
                return True

        return False

    # ------------------------------------------------------------------
    # Multi-cell range actions (shift/ctrl click selection)
    # ------------------------------------------------------------------

    def _selectedVariableCells(self):
        """Return list of (item, stationId, datetime, fieldId) for every
        currently-selected variable cell. Non-variable cells are filtered."""
        result = []
        for item in self.dataTable.selectedItems():
            data = item.data(Qt.UserRole)
            if not data or not isinstance(data, dict):
                continue
            fieldId = data.get("field")
            if not fieldId:
                continue

            row = item.row()
            if self.displayMode == "Time":
                idItem = self.dataTable.item(row, 0)
                if not idItem:
                    continue
                sid = idItem.text()
                dt = self.currentDatetime
            else:
                timeItem = self.dataTable.item(row, 0)
                if not timeItem:
                    continue
                dt = timeItem.data(Qt.UserRole)
                sid = self.currentStation
            if sid and dt:
                result.append((item, sid, dt, fieldId))
        return result

    def ignoreSelectedCells(self):
        """Ignore every currently-selected variable cell in one shot."""
        cells = self._selectedVariableCells()
        if not cells:
            return

        # Snapshot for undo
        snapshots = []  # (sid, fid, dt, snapshot_tuple)
        for _item, sid, dt, fid in cells:
            snapshots.append((sid, fid, dt, self._snapshotManCell(sid, fid, dt)))

        affected_times = set()

        for _item, sid, dt, fid in cells:
            fieldList = self.getList(fid, sid)
            fieldList[dt] = ""
            self.setList(fieldList, fid, sid)
            affected_times.add(dt)

        def restore():
            for sid, fid, dt, snap in snapshots:
                self._restoreManCell(snap, sid, fid, dt)
            for dt in affected_times:
                self.writeManFile(dt)
        self._pushUndo(f"Ignore {len(cells)} selected cell(s)", restore)

        for dt in affected_times:
            self.writeManFile(dt)

        self.redisplay()
        self.statusBar.showMessage(f"Ignored {len(cells)} cell(s)", 4000)

    def cancelIgnoreSelectedCells(self):
        """Remove per-time ignore/override for every selected variable cell
        that has one (only touches self.man entries)."""
        cells = self._selectedVariableCells()
        if not cells:
            return

        snapshots = []
        affected_times = set()
        cleared = 0
        for _item, sid, dt, fid in cells:
            fieldList = self.getList(fid, sid)
            if dt in fieldList:
                snapshots.append((sid, fid, dt, self._snapshotManCell(sid, fid, dt)))
                del fieldList[dt]
                self.setList(fieldList, fid, sid)
                affected_times.add(dt)
                cleared += 1

        if cleared == 0:
            self.statusBar.showMessage(
                "No per-time ignore/override on selected cells", 3000
            )
            return

        def restore():
            for sid, fid, dt, snap in snapshots:
                self._restoreManCell(snap, sid, fid, dt)
            for dt in affected_times:
                self.writeManFile(dt)
        self._pushUndo(
            f"Cancel ignore/override on {cleared} selected cell(s)", restore
        )

        for dt in affected_times:
            self.writeManFile(dt)

        self.redisplay()
        self.statusBar.showMessage(
            f"Cleared {cleared} cell(s) of per-time overrides", 4000
        )

    def hasAnyIgnoreAll(self, stationId):
        """Return True if any variable for this station is currently set to ignore-all ('')."""
        if not stationId or stationId not in self.all:
            return False
        return any(val == "" for val in self.all[stationId])

    def ignoreAllObsForStation(self, stationId=None):
        """Ignore-all every variable for a station in one shot.

        If called from the context menu, uses self.tagId. Can also be called
        programmatically by passing stationId directly (used by the bulk
        non-reporting-station cleanup).
        """
        if stationId is None:
            stationId = self.tagId
        if not stationId:
            return

        # Snapshot prior all-row for undo
        prev_row = list(self.all[stationId]) if stationId in self.all else None

        def restore():
            if prev_row is None:
                self.all.pop(stationId, None)
            else:
                self.all[stationId] = list(prev_row)
            self.writeAllFile()

        self._pushUndo(f"Ignore all obs for {stationId}", restore)

        # Ensure entry exists, then set every field to "" (ignore-all)
        if stationId not in self.all:
            self.all[stationId] = ["*"] * self.numVariableFields
        for fieldIndex in range(self.numVariableFields):
            self.all[stationId][fieldIndex] = ""

        self.writeAllFile()
        self.redisplay()
        self.statusBar.showMessage(
            f"Ignored all observations for station {stationId}", 4000
        )

    def cancelIgnoreAllObsForStation(self, stationId=None):
        """Reset every variable for this station back to '*' (no override)."""
        if stationId is None:
            stationId = self.tagId
        if not stationId or stationId not in self.all:
            return

        # Snapshot for undo
        prev_row = list(self.all[stationId])

        def restore():
            self.all[stationId] = list(prev_row)
            self.writeAllFile()

        self._pushUndo(f"Cancel ignore-all for {stationId}", restore)

        for fieldIndex in range(self.numVariableFields):
            self.all[stationId][fieldIndex] = "*"

        self.writeAllFile()
        self.redisplay()
        self.statusBar.showMessage(
            f"Cleared ignore-all for station {stationId}", 4000
        )

    def findNonReportingStations(self):
        """Return a sorted list of station IDs that have zero raw obs across
        every variable and every loaded time, and which are not already
        fully ignore-alled.

        We only consider stations that are in self.static (i.e. known to the
        app); a station with no raw data AND no override entries is already
        skipped by getIdList, so we focus on stations that would otherwise
        appear in the table.
        """
        nonReporting = []

        for stationId in self.static.keys():
            # Skip if already fully ignore-alled (every variable is "")
            if (stationId in self.all
                    and all(v == "" for v in self.all[stationId])):
                continue

            if self.stationHasNoObs(stationId):
                nonReporting.append(stationId)

        return sorted(nonReporting)

    def stationHasNoObs(self, stationId):
        """True when a station has no non-blank raw observation anywhere in the
        loaded period."""
        if stationId not in self.raw:
            return True

        for fieldIndex in range(self.numVariableFields):
            if fieldIndex < len(self.raw[stationId]):
                for value in self.raw[stationId][fieldIndex].values():
                    if value.strip():
                        return False

        return True

    def findNonReportingIgnoreAlls(self):
        """Return station IDs that carry an all.manual entry but have no raw
        observations anywhere in the loaded period.

        These are ignore-all rows left behind by stations that have since gone
        silent - the entry is no longer suppressing anything.
        """
        return sorted(sid for sid in self.all if self.stationHasNoObs(sid))

    def clearIgnoreAllForNonReporting(self):
        """Remove the all.manual entry for every station with no raw obs in the
        loaded period, after confirmation.  Per-time QC is left alone."""
        try:
            candidates = self.findNonReportingIgnoreAlls()
        except Exception as e:
            print(f"Error scanning for stale ignore-all entries: {e}")
            print(traceback.format_exc())
            QMessageBox.warning(
                self, "Scan Error",
                f"Error scanning for stale ignore-all entries:\n{str(e)}"
            )
            return

        if not candidates:
            QMessageBox.information(
                self, "No Stale Ignore-All Entries",
                "Every station with an all.manual entry has reported at least "
                "once in the loaded period."
            )
            return

        # Build a preview list (cap so the dialog stays sane)
        previewLimit = 40
        previewIds = candidates[:previewLimit]
        more = len(candidates) - len(previewIds)
        previewText = "\n".join(previewIds)
        if more > 0:
            previewText += f"\n... and {more} more"

        msg = (
            f"Found {len(candidates)} station(s) with an all.manual entry but "
            f"NO raw observations across the entire loaded period.\n\n"
            f"Their ignore-all entries will be removed from all.manual.  "
            f"Per-time QC in the hourly .manual files is left alone.\n\n"
            f"Stations:\n{previewText}\n\n"
            f"Proceed?"
        )

        response = QMessageBox.question(
            self, "Clear Ignore-All for Non-Reporting Stations",
            msg,
            QMessageBox.Yes | QMessageBox.Cancel,
            QMessageBox.Cancel
        )

        if response != QMessageBox.Yes:
            return

        # Snapshot prior state for every affected station
        prev_rows = {sid: list(self.all[sid]) for sid in candidates}

        def restore():
            for sid, prev in prev_rows.items():
                self.all[sid] = list(prev)
            self.writeAllFile()

        self._pushUndo(
            f"Clear ignore-all on {len(candidates)} non-reporting station(s)",
            restore
        )

        for stationId in candidates:
            self.all.pop(stationId, None)

        self.writeAllFile()
        self.redisplay()
        self.statusBar.showMessage(
            f"Cleared ignore-all on {len(candidates)} non-reporting station(s)",
            5000
        )

    def ignoreNonReportingStations(self):
        """Find every station with no raw obs in the loaded period and
        ignore-all every variable for them, after user confirmation."""
        try:
            candidates = self.findNonReportingStations()
        except Exception as e:
            print(f"Error scanning for non-reporting stations: {e}")
            print(traceback.format_exc())
            QMessageBox.warning(
                self, "Scan Error",
                f"Error scanning for non-reporting stations:\n{str(e)}"
            )
            return

        if not candidates:
            QMessageBox.information(
                self, "No Non-Reporting Stations",
                "No stations were found with zero observations across the "
                "loaded period (that aren't already fully ignored)."
            )
            return

        # Build a preview list (cap so the dialog stays sane)
        previewLimit = 40
        previewIds = candidates[:previewLimit]
        more = len(candidates) - len(previewIds)
        previewText = "\n".join(previewIds)
        if more > 0:
            previewText += f"\n… and {more} more"

        msg = (
            f"Found {len(candidates)} station(s) with NO raw observations "
            f"across the entire loaded period.\n\n"
            f"Ignoring will set every variable (TMP/DEW/DIR/SPD/GST) to "
            f"ignore-all and persist to all.manual.\n\n"
            f"Stations:\n{previewText}\n\n"
            f"Proceed?"
        )

        response = QMessageBox.question(
            self, "Ignore Non-Reporting Stations",
            msg,
            QMessageBox.Yes | QMessageBox.Cancel,
            QMessageBox.Cancel
        )

        if response != QMessageBox.Yes:
            return

        # Snapshot prior state for every affected station
        prev_rows = {
            sid: (list(self.all[sid]) if sid in self.all else None)
            for sid in candidates
        }

        def restore():
            for sid, prev in prev_rows.items():
                if prev is None:
                    self.all.pop(sid, None)
                else:
                    self.all[sid] = list(prev)
            self.writeAllFile()

        self._pushUndo(
            f"Ignore {len(candidates)} non-reporting station(s)", restore
        )

        # Apply ignore-all to every candidate, then write the file once
        for stationId in candidates:
            if stationId not in self.all:
                self.all[stationId] = ["*"] * self.numVariableFields
            for fieldIndex in range(self.numVariableFields):
                self.all[stationId][fieldIndex] = ""

        self.writeAllFile()
        self.redisplay()
        self.statusBar.showMessage(
            f"Ignored {len(candidates)} non-reporting station(s)", 5000
        )

    def runMoa(self):
        """Run the MatchObsAll process"""
        message = ("This procedure will run in the background and take a few minutes. "
                 "This will create new obs grids from the MatchObsAll procedure. "
                 "Please wait.")

        response = QMessageBox.question(
            self,
            "Run MatchObsAll",
            message,
            QMessageBox.Ok | QMessageBox.Cancel
        )

        if response == QMessageBox.Ok:
            self.runMoaBtn.setEnabled(False)
            self.runMoaBtn.setText("Running...")
            self.runMoaBtn.setStyleSheet("background-color:yellow;")
            self.statusBar.showMessage("Running MatchObsAll...")

            try:
                # Start the process in a separate thread
                self.runnerThread = RunnerThread(f"ssh pv2 {MOADIR}/bin/RunObs.sh")
                self.runnerThread.runFinished.connect(self.moaFinished)
                self.runnerThread.start()
            except Exception as e:
                print(f"Error starting runner thread: {e}")
                self.runMoaBtn.setEnabled(True)
                self.runMoaBtn.setText("Run MatchObsAll")
                self.runMoaBtn.setStyleSheet("background-color:red;")
                QMessageBox.warning(self, "Error", f"Error running MatchObsAll: {str(e)}")
                self.statusBar.showMessage("Error running MatchObsAll", 3000)

    def closeEvent(self, event):
        """Let a running MatchObsAll thread finish before the window goes away.

        Without this, quitting mid-run destroyed a live QThread and Qt aborted
        the process with "QThread: Destroyed while thread is still running".
        """
        thread = getattr(self, "runnerThread", None)
        if thread is not None and thread.isRunning():
            response = QMessageBox.question(
                self, "MatchObsAll Still Running",
                "MatchObsAll is still running.  Wait for it to finish?",
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel,
                QMessageBox.Yes
            )
            if response == QMessageBox.Cancel:
                event.ignore()
                return
            if response == QMessageBox.Yes:
                self.statusBar.showMessage("Waiting for MatchObsAll to finish...")
                thread.wait()
            else:
                thread.terminate()
                thread.wait()

        event.accept()

    def moaFinished(self, success, message):
        """Handle completion of MatchObsAll process"""
        self.runMoaBtn.setEnabled(True)
        self.runMoaBtn.setText("Run MatchObsAll")
        self.runMoaBtn.setStyleSheet("background-color:green;")

        if success:
            QMessageBox.information(self, "Success", "MatchObsAll completed successfully")
            self.refreshData()
        else:
            QMessageBox.warning(self, "Error", f"Error running MatchObsAll: {message}")

        self.statusBar.showMessage("Ready")

    def getGeography(self):
        """Read geography information from geom.dat file"""
        self.curr = {}
        self.currLat = {}
        self.currLon = {}

        fullname = GEOMPATH
        if not os.path.exists(fullname):
            self.statusBar.showMessage(f"Geography file not found: {fullname}")
            return

        #  Each row is parsed on its own.  A single unparseable lat/lon used to
        #  unwind the whole read, so the stations after it were never loaded and
        #  every one of them was silently dropped from the display.
        skipped = 0
        try:
            with open(fullname, "r") as datafile:
                for lineNumber, line in enumerate(datafile, 1):
                    pieces = line.strip().split(",")
                    if len(pieces) != 14:  # Expected format
                        skipped += 1
                        dbg(f"geom.dat line {lineNumber}: expected 14 fields, "
                            f"got {len(pieces)}")
                        continue

                    try:
                        stationId = pieces[0].strip()
                        datalist = [piece.strip() for piece in pieces[1:]]
                        lat = float(pieces[2].strip())
                        lon = float(pieces[3].strip())
                    except ValueError as e:
                        skipped += 1
                        dbg(f"geom.dat line {lineNumber}: {e}")
                        continue

                    if not stationId:
                        skipped += 1
                        continue

                    self.curr[stationId] = datalist
                    self.currLat[stationId] = lat
                    self.currLon[stationId] = lon
        except Exception as e:
            self.statusBar.showMessage(f"Error reading geography file: {e}")
            print(f"Error reading geography file: {e}")
            print(traceback.format_exc())
            return

        if skipped:
            print(f"getGeography: skipped {skipped} unparseable row(s) in {fullname}")

        #  An empty result means every station will be dropped from the table,
        #  and the status bar is overwritten moments later during start-up - so
        #  say so in a way the operator cannot miss.
        if not self.curr:
            QMessageBox.warning(
                self, "Geography File Unusable",
                f"No usable stations were read from:\n{fullname}\n\n"
                f"{skipped} row(s) could not be parsed.  The observation table "
                f"will be empty until this file is fixed."
            )

    def listDataFiles(self):
        """Get list of data files in the data directory"""
        if not os.path.exists(DATADIR):
            print(f"Data directory not found: {DATADIR}")
            return

        try:
            allFiles = os.listdir(DATADIR)
            self.dataFiles = []
            self.manFiles = []
            self.allFile = "all.manual"

            searchpat1 = re.compile(r"\d{12}\.dat$")
            searchpat2 = re.compile(r"\d{12}\.manual$")
            searchpat3 = re.compile(r"all\.manual$")

            for file in allFiles:
                if searchpat1.search(file):
                    self.dataFiles.append(file)
                if searchpat2.search(file):
                    self.manFiles.append(file)
                if searchpat3.search(file):
                    self.allFile = file

            self.dataFiles.sort(reverse=True)
            self.manFiles.sort(reverse=True)

            # Extract dates from filenames
            self.dates = []
            for file in self.dataFiles:
                datetime, ext = file.split(".")
                if datetime not in self.dates:
                    self.dates.append(datetime)

            for file in self.manFiles:
                datetime, ext = file.split(".")
                if datetime not in self.dates:
                    self.dates.append(datetime)

            self.dates.sort(reverse=True)

            # Print debug info
            print(f"Found {len(self.dataFiles)} data files")
            print(f"Found {len(self.manFiles)} manual files")
            print(f"All file: {self.allFile}")
            print(f"Total dates: {len(self.dates)}")

            # Update selectors
            self.updateSelectors()
        except Exception as e:
            print(f"Error listing data files: {e}")
            print(traceback.format_exc())

    def decodeDataLine(self, line):
        """Decode a comma-separated data line, splitting on commas and stripping spaces"""
        try:
            rawvals = line.split(",")
            vals = []
            for rawval in rawvals:
                val = rawval.strip()
                vals.append(val)

            # Ensure the line has enough fields
            while len(vals) < (self.numVariableFields + self.dataFieldStart):
                vals.append("")

            return vals
        except Exception as e:
            print(f"Error decoding data line: {e}")
            print(f"Line: {line}")
            return []

    def checkLatLonLimits(self, vals):
        """Check if lat/lon limits need to be expanded"""
        try:
            lat = float(vals[2])
            lon = float(vals[3])

            if lat < self.minLat:
                self.minLat = lat
                self.changedLatLonLimits = True

            if lat > self.maxLat:
                self.maxLat = lat
                self.changedLatLonLimits = True

            if lon < self.minLon:
                self.minLon = lon
                self.changedLatLonLimits = True

            if lon > self.maxLon:
                self.maxLon = lon
                self.changedLatLonLimits = True
        except (ValueError, IndexError) as e:
            print(f"Error checking lat/lon limits: {e}")

    def readShapes(self):
        """Fill in self.geom for every station that needs it.

        This used to be gated on changedLatLonLimits - and cleared that flag
        inside its own loop - so a station that started reporting inside the
        existing lat/lon box never got a geom entry and was dropped from every
        display until the app was restarted.  A limits change still forces a
        full rebuild; otherwise only the stations we have no entry for are done.
        """
        try:
            if self.changedLatLonLimits:
                stations = list(self.static.keys())
                self.changedLatLonLimits = False
            else:
                stations = [sid for sid in self.static if sid not in self.geom]

            if not stations:
                return

            self.statusBar.showMessage("Updating location information for points")

            notInGeom = []
            for stationId in stations:
                entry = self.curr.get(stationId)
                if entry is None:
                    notInGeom.append(stationId)
                    continue

                self.geom[stationId] = [
                    entry[5],   # zone
                    entry[6],   # zone name
                    entry[4],   # state
                    entry[7],   # county
                    entry[8],   # county name
                    entry[9],   # cwa
                    entry[10],  # fire zone
                    entry[12],  # fire zone name
                ]

            if notInGeom:
                print(f"readShapes: {len(notInGeom)} station(s) not in the "
                      f"geometry file, e.g. {', '.join(notInGeom[:5])}")
        except Exception as e:
            print(f"Error reading shapes: {e}")
            print(traceback.format_exc())

    def readOneDataFile(self, fullname, datetime):
        """Read data from one raw datafile"""
        self.statusBar.showMessage(f"Reading {fullname}")

        try:
            with open(fullname, "r") as datafile:
                for line in datafile:
                    if line.startswith("#"):
                        continue

                    vals = self.decodeDataLine(line)
                    if not vals or len(vals) < 2:
                        continue

                    stationId = vals[0]
                    if not stationId:  # Skip empty IDs
                        continue

                    if stationId not in self.raw:
                        # New station - add static data
                        fielddata = []
                        for field in range(self.numStaticAttrs):
                            fielddata.append(vals[field + 1] if field + 1 < len(vals) else "")
                        self.static[stationId] = fielddata
                        self.checkLatLonLimits(vals)

                        # New station - add room for raw data
                        fielddata = []
                        for field in range(self.numVariableFields):
                            data = {}
                            fielddata.append(data)
                        self.raw[stationId] = fielddata

                    # Store the raw data
                    for field in range(self.numVariableFields):
                        fieldI = field + self.dataFieldStart
                        if fieldI < len(vals):
                            self.raw[stationId][field][datetime] = vals[fieldI]

            # Add stations to the allStations list.  The membership test used
            # to scan the whole list per station, which is quadratic over a
            # few thousand sites; the set keeps the list's order and type.
            known = set(self.allStations)
            for stationId in self.raw:
                if stationId not in known:
                    known.add(stationId)
                    self.allStations.append(stationId)

            # Print number of stations found
            print(f"Processed {fullname}: found {len(self.raw)} stations")
        except Exception as e:
            self.statusBar.showMessage(f"Error reading data file: {e}")
            print(f"Error reading data file: {e}")
            print(traceback.format_exc())

    def readOneManFile(self, fullname, datetime):
        """Read data from one manual override file"""
        self.statusBar.showMessage(f"Reading {fullname}")

        try:
            with open(fullname, "r") as datafile:
                #  Same reason as readAllFile: this file is the complete set of
                #  overrides for its time, so an entry someone else deleted must
                #  not survive in memory and get written back on the next QC action.
                for fieldData in self.man.values():
                    for field in fieldData:
                        field.pop(datetime, None)

                for line in datafile:
                    if line.startswith("#"):
                        continue

                    vals = self.decodeDataLine(line)
                    if not vals or len(vals) < 2:
                        continue

                    stationId = vals[0]
                    if not stationId:  # Skip empty IDs
                        continue

                    if stationId not in self.man:
                        # If new station - add static data if needed
                        if stationId not in self.raw:
                            fielddata = []
                            for field in range(self.numStaticAttrs):
                                fielddata.append(vals[field + 1] if field + 1 < len(vals) else "")
                            self.static[stationId] = fielddata
                            self.checkLatLonLimits(vals)

                        # New man station - add room for manual data
                        fielddata = []
                        for field in range(self.numVariableFields):
                            data = {}
                            fielddata.append(data)
                        self.man[stationId] = fielddata

                    # Store the manual data.  "*" is the no-override placeholder
                    # writeManFile emits for fields a station was never QC'd on;
                    # storing it made an untouched field look like a real override,
                    # so cancelling one ignore left a row of all-"*" values behind
                    # that could never be removed.  Absent and "*" already behave
                    # identically everywhere the value is read.
                    for field in range(self.numVariableFields):
                        fieldI = field + self.dataFieldStart
                        if fieldI < len(vals) and vals[fieldI] != "*":
                            self.man[stationId][field][datetime] = vals[fieldI]

            # Add stations to the allStations list.  The membership test used
            # to scan the whole list per station, which is quadratic over a
            # few thousand sites; the set keeps the list's order and type.
            known = set(self.allStations)
            for stationId in self.man:
                if stationId not in known:
                    known.add(stationId)
                    self.allStations.append(stationId)

            # Print number of stations found
            print(f"Processed {fullname}: found {len(self.man)} stations with manual overrides")
        except Exception as e:
            self.statusBar.showMessage(f"Error reading manual override file: {e}")
            print(f"Error reading manual override file: {e}")
            print(traceback.format_exc())

    def readAllFile(self, fullname):
        """Read data from the all.manual file"""
        self.statusBar.showMessage(f"Reading {fullname}")

        try:
            with open(fullname, "r") as datafile:
                #  The file is the complete set of ignore-alls.  Without dropping
                #  what we already hold, a station another workstation removed
                #  from all.manual survived in memory and writeAllFile put it
                #  straight back, silently undoing their change.  Cleared after
                #  the open succeeds so a read error leaves state untouched.
                self.all = {}

                for line in datafile:
                    if line.startswith("#"):
                        continue

                    vals = self.decodeDataLine(line)
                    if not vals or len(vals) < 2:
                        continue

                    stationId = vals[0]
                    if not stationId:  # Skip empty IDs
                        continue

                    if stationId not in self.all:
                        # If new station - add static data if needed
                        if stationId not in self.raw:
                            fielddata = []
                            for field in range(self.numStaticAttrs):
                                fielddata.append(vals[field + 1] if field + 1 < len(vals) else "")
                            self.static[stationId] = fielddata
                            self.checkLatLonLimits(vals)

                    # Store the all data
                    self.all[stationId] = []
                    for field in range(self.numVariableFields):
                        fieldI = field + self.dataFieldStart
                        if fieldI < len(vals):
                            self.all[stationId].append(vals[fieldI])
                        else:
                            self.all[stationId].append("*")

            # Add stations to the allStations list.  The membership test used
            # to scan the whole list per station, which is quadratic over a
            # few thousand sites; the set keeps the list's order and type.
            known = set(self.allStations)
            for stationId in self.all:
                if stationId not in known:
                    known.add(stationId)
                    self.allStations.append(stationId)

            # Print number of stations found
            print(f"Processed {fullname}: found {len(self.all)} stations with all overrides")
        except Exception as e:
            self.statusBar.showMessage(f"Error reading all.manual file: {e}")
            print(f"Error reading all.manual file: {e}")
            print(traceback.format_exc())

    def readAllDataFiles(self):
        """Read all data files that have changed"""
        # Try to list data files, might fail if directory doesn't exist
        try:
            self.listDataFiles()
        except Exception as e:
            print(f"Error listing data files: {e}")
            print(traceback.format_exc())
            return

        # Process data files
        for filename in self.dataFiles:
            try:
                fullname = os.path.join(DATADIR, filename)
                lastMod = os.path.getmtime(fullname)

                if filename in self.fileTimes and self.fileTimes[filename] == lastMod:
                    continue

                datetime, ext = filename.split(".")
                self.readOneDataFile(fullname, datetime)
                self.fileTimes[filename] = lastMod
            except Exception as e:
                print(f"Error processing data file {filename}: {e}")
                print(traceback.format_exc())

        # Process manual override files
        for filename in self.manFiles:
            try:
                fullname = os.path.join(DATADIR, filename)
                lastMod = os.path.getmtime(fullname)

                if filename in self.manTimes and self.manTimes[filename] == lastMod:
                    continue

                datetime, ext = filename.split(".")
                self.readOneManFile(fullname, datetime)
                self.manTimes[filename] = lastMod
            except Exception as e:
                print(f"Error processing manual file {filename}: {e}")
                print(traceback.format_exc())

        # Process all.manual file
        fullname = os.path.join(DATADIR, self.allFile)
        if os.path.exists(fullname):
            try:
                lastMod = os.path.getmtime(fullname)

                if self.allTimes != lastMod:
                    self.readAllFile(fullname)
                    self.allTimes = lastMod
            except Exception as e:
                print(f"Error processing all.manual file: {e}")
                print(traceback.format_exc())

        # Update station list and shape info
        self.allStations.sort()

        # Store the latest date (chronologically) - which is the first in the dates list
        # since they're sorted in reverse chronological order
        self.latestDataDate = self.dates[0] if self.dates else None

        print(f"Latest date: {self.latestDataDate}")

        # Print stats about the data we loaded
        print(f"Loaded {len(self.raw)} stations with raw data")
        print(f"Loaded {len(self.man)} stations with manual overrides")
        print(f"Loaded {len(self.all)} stations with all overrides")
        print(f"Total stations: {len(self.allStations)}")

        self.readShapes()

        # Update CWA filter with available CWAs
        self.updateCwaFilter()

        # Update selectors
        self.updateSelectors()

    def updateSelectors(self):
        """Repopulate the time and station dropdowns without firing signals.

        clear() followed by the first addItem() moves the current index from -1
        to 0, which emitted currentIndexChanged in the middle of a reload and
        dragged the user back to the newest time before a single file had been
        re-read.  Both selectors are rebuilt with signals blocked and the prior
        selection restored if it still exists; the callers redisplay explicitly.
        """
        previousTime = self.timeSelector.currentData()

        self.timeSelector.blockSignals(True)
        try:
            self.timeSelector.clear()
            for date in self.dates:
                self.timeSelector.addItem(self.getTimeFromDateTime(date), date)

            if previousTime is not None:
                index = self.timeSelector.findData(previousTime)
                if index >= 0:
                    self.timeSelector.setCurrentIndex(index)
        finally:
            self.timeSelector.blockSignals(False)

        allStations = sorted(self.allStations)
        self.stationSelector.setAllStations(allStations)

        dbg(f"Added {len(allStations)} stations to searchable selector")

    def getTimeFromDateTime(self, date):
        """Turn a YYYYMMDDHHMM date into a time string of MM/DD/YY HHZ"""
        try:
            year = int(date[2:4])
            month = int(date[4:6])
            day = int(date[6:8])
            hour = int(date[8:10])

            tstring = f"{month:2d}/{day:2d}/{year:02d} {hour:02d}Z"
            return tstring
        except (ValueError, IndexError, TypeError):
            return date

    def getSplitTime(self, date):
        """Turn a YYYYMMDDHHMM date into a two-line string: MON DD\n DDZ"""
        try:
            year = int(date[2:4])
            mon = int(date[4:6])
            mons = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
            month = mons[mon-1]
            day = int(date[6:8])
            hour = int(date[8:10])
            minute = int(date[10:12])

            ugc = calendar.timegm((2000+year, mon, day, hour, minute, 0, 0, 0, 0))
            timeTuple = time.gmtime(ugc)
            days = ("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
            weekday = timeTuple.tm_wday

            tstring = f"{days[weekday]} {month} {day:2d} {hour:02d}Z"
            return tstring
        except (ValueError, IndexError, TypeError) as e:
            print(f"Error in getSplitTime: {e}")
            return date

    def getList(self, fieldId, stationId):
        """Get list of values for a specific field and station"""
        fieldIndex = -1
        for i, field in enumerate(self.variableFields):
            if field["id"] == fieldId:
                fieldIndex = i
                break

        if fieldIndex == -1 or stationId not in self.man:
            return {}

        return self.man[stationId][fieldIndex]

    def setList(self, fieldList, fieldId, stationId):
        """Set list of values for a specific field and station"""
        fieldIndex = -1
        for i, field in enumerate(self.variableFields):
            if field["id"] == fieldId:
                fieldIndex = i
                break

        if fieldIndex == -1:
            return

        if stationId not in self.man:
            # Initialize man entry for this station
            fielddata = []
            for field in range(self.numVariableFields):
                data = {}
                fielddata.append(data)
            self.man[stationId] = fielddata

        self.man[stationId][fieldIndex] = fieldList

    def getAll(self, fieldId, stationId):
        """Get 'all' value for a specific field and station"""
        fieldIndex = -1
        for i, field in enumerate(self.variableFields):
            if field["id"] == fieldId:
                fieldIndex = i
                break

        if fieldIndex == -1 or stationId not in self.all:
            return "*"

        return self.all[stationId][fieldIndex]

    def setAll(self, allValue, fieldId, stationId):
        """Set 'all' value for a specific field and station"""
        fieldIndex = -1
        for i, field in enumerate(self.variableFields):
            if field["id"] == fieldId:
                fieldIndex = i
                break

        if fieldIndex == -1:
            return

        if stationId not in self.all:
            # Initialize all entry for this station
            self.all[stationId] = ["*"] * self.numVariableFields

        self.all[stationId][fieldIndex] = allValue

    def getIdList(self, datetime):
        """Get list of station IDs in group/elevation order"""
        sortKeys = []

        # Get sort direction (ascending=0, descending=1)
        sortDirection = self.sortDirOptions.checkedId()
        
        sortOpt = self.sortOptions.checkedId()

        # List of stations to include
        includedStations = []

        # Check if any CWAs are selected
        #  Filter whenever the CWA filter exists at all.  Keying off
        #  selectedCwas made an empty selection mean "no filter", so "Select
        #  None" showed every station - identical to "Select All".
        usingCwaFilter = len(self.cwaCheckboxes) > 0

        for stationId in list(self.static.keys()):
            # Skip stations with no data and no all.manual entries
            hasAnyData = False

            # Check if there's any raw data for any field at any time
            for fieldIndex in range(self.numVariableFields):
                if stationId in self.raw:
                    rawList = self.raw[stationId][fieldIndex]
                    for time_key, value in rawList.items():
                        if value.strip():
                            hasAnyData = True
                            break
                if hasAnyData:
                    break

            # If no raw data, check if there's any entries in all.manual
            if not hasAnyData and stationId in self.all:
                for value in self.all[stationId]:
                    if value != "*":  # If there's any non-default value
                        hasAnyData = True
                        break

            # Skip this station if it has no data and no all.manual entries
            if not hasAnyData:
                dbg(f"Skipping station {stationId} - no data or overrides")
                continue

            # Continue with the rest of the processing for stations with data
            if stationId in self.geom:
                try:
                    cwa = self.geom[stationId][5]
                    zone = self.geom[stationId][0]

                    # Sort by temperature or elevation
                    # sortOpt = self.sortOptions.checkedId()

                    # Skip stations from unselected CWAs if filter is active
                    if usingCwaFilter and cwa not in self.selectedCwas:
                        continue

                    if usingCwaFilter:
                        # Only use CWA/Zone for sorting when CWA filter is active
                        cwaSort = cwa
                        zoneSort = zone

                        # Apply CWA/Zone ordering if defined
                        if cwa in CWAORDER:
                            cwaSort = CWAORDER[cwa]
                        if zone in ZONEORDER:
                            zoneSort = ZONEORDER[zone]

                    if sortOpt < 5:  # Field-based sort
                        tmpList, manTmpList, allTmp = self.getVals(stationId, sortOpt)
                        hasData = False
                        tmp = -100
                        if datetime in tmpList:
                            tmpVal = tmpList[datetime]
                            if tmpVal != "":
                                try:
                                    tmp = int(float(tmpVal))
                                    hasData = True
                                except ValueError:
                                    tmp = -100
                            else:
                                tmp = -100

                        if hasData:
                            # Modify the sort key based on direction
                            if sortDirection == 0:  # Ascending
                                tmpSortValue = f"{int(tmp):05d}"
                            else:  # Descending
                                tmpSortValue = f"{400-int(tmp):05d}"
                        else:
                            tmpSortValue = "1_0000"

                        if usingCwaFilter:
                            # When CWA filter is active: sort by CWA, zone, then value
                            sortKey = f"{cwaSort}-{zoneSort}-{tmpSortValue}"
                        else:
                            # When no CWA filter: sort only by value
                            sortKey = tmpSortValue

                    elif sortOpt == 5:  # Elevation sort
                        try:
                            elev = int(float(self.static[stationId][3]))
                        except (ValueError, IndexError):
                            elev = 0

                        # Modify the sort key based on direction
                        if sortDirection == 0:  # Ascending
                            elevSortValue = f"{elev:05d}"
                        else:  # Descending
                            elevSortValue = f"{20000-elev:05d}"

                        if usingCwaFilter:
                            # When CWA filter is active: sort by CWA, zone, then elevation
                            sortKey = f"{cwaSort}-{zoneSort}-{elevSortValue}"
                        else:
                            # When no CWA filter: sort only by elevation
                            sortKey = elevSortValue

                    if sortOpt == 6:  # ID-based sort
                        #  The ID itself is the tiebreaker in the tuple below,
                        #  so the key only has to carry the grouping.
                        sortKey = f"{cwaSort}-{zoneSort}" if usingCwaFilter else ""

                    #  Carry the station ID beside the key rather than encoded
                    #  into it.  Recovering it with split("-")[-1] truncated any
                    #  ID containing a hyphen (e.g. the CoCoRaHS "ID-AD-24"
                    #  form), which then blew up the render with a KeyError.
                    sortKeys.append((sortKey, stationId))
                except (IndexError, ValueError) as e:
                    print(f"Error processing station {stationId} in getIdList: {e}")
            else:
                dbg(f"Skipping site {stationId} not in geometry file")

        sortKeys.sort()

        if sortOpt == 6 and sortDirection == 1:
            #  Descending ID sort: flip the IDs within each group while leaving
            #  the groups themselves in ascending order.  With no CWA filter the
            #  key is "" for everything, so this is one group - which is exactly
            #  the plain reverse sort the old code did in that case.
            grouped = {}
            for key, stationId in sortKeys:
                grouped.setdefault(key, []).append(stationId)

            sortKeys = [
                (key, stationId)
                for key in sorted(grouped)
                for stationId in sorted(grouped[key], reverse=True)
            ]

        return [stationId for _key, stationId in sortKeys]

    def getVals(self, stationId, fieldIndex):
        """Get raw, manual, and all values for a specific station and field index"""
        if stationId in self.raw:
            rawList = self.raw[stationId][fieldIndex]
        else:
            rawList = {}

        if stationId in self.man:
            manList = self.man[stationId][fieldIndex]
        else:
            manList = {}

        if stationId in self.all:
            allVal = self.all[stationId][fieldIndex]
        else:
            allVal = "*"

        return (rawList, manList, allVal)

    def _changedOnDisk(self, fullname, knownMtime):
        """True when the file has moved on since this session last read it.

        Every write here rewrites the whole file from memory, so if another
        workstation saved in the meantime we would silently discard their work.
        The mtimes were already being recorded on every read and write - they
        just were not being consulted before overwriting.
        """
        if not knownMtime or not os.path.exists(fullname):
            return False
        try:
            return os.path.getmtime(fullname) != knownMtime
        except OSError:
            return False

    def _confirmOverwrite(self, fullname):
        """Ask whether to clobber a file someone else has changed."""
        response = QMessageBox.question(
            self, "File Changed On Disk",
            f"{os.path.basename(fullname)} has been changed by someone else "
            f"since this session last read it.\n\n"
            f"Saving now would overwrite their edits with this session's copy."
            f"\n\nOverwrite anyway?  Choosing No reloads the file from disk "
            f"and discards the change you just made.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )
        return response == QMessageBox.Yes

    def _writeAtomic(self, tmpname, fullname):
        """Put the finished temp file in place and return the new mtime.

        open(path, "w") truncates in place, so a crash or a full disk part-way
        through left a corrupt file for every other workstation to read.
        os.replace is atomic, so readers see either the old file or the new one.
        """
        try:
            os.chmod(tmpname, 0o666)
        except Exception as e:
            print(f"Error setting file permissions: {e}")

        os.replace(tmpname, fullname)
        return os.path.getmtime(fullname)

    def _discardTemp(self, tmpname):
        """Remove a part-written temp file after a failed write."""
        try:
            os.remove(tmpname)
        except OSError:
            pass

    def writeManFile(self, datetime):
        """Write data to one manual override file"""
        filename = f"{datetime}.manual"
        fullname = os.path.join(DATADIR, filename)

        if self._changedOnDisk(fullname, self.manTimes.get(filename)):
            if not self._confirmOverwrite(fullname):
                self.readOneManFile(fullname, datetime)
                self.manTimes[filename] = os.path.getmtime(fullname)
                self.redisplay()
                return

        self.statusBar.showMessage(f"Writing {fullname}")

        # Back up the existing file before overwriting
        self._backupFile(fullname)

        tmpname = f"{fullname}.tmp{os.getpid()}"
        try:
            with open(tmpname, "w") as datafile:
                # Debug: Print the data structure
                dbg(f"Writing to {fullname}")
                dbg(f"Field order: ID, name, lat, lon, elev, tmp, dew, dir, spd, gst")

                for stationId in sorted(self.man.keys()):
                    # Check if there's any manual override specifically for this datetime
                    hasOverride = False
                    for fieldIndex in range(self.numVariableFields):
                        fieldData = self.man[stationId][fieldIndex]
                        if datetime in fieldData and fieldData[datetime] != "*":
                            hasOverride = True
                            dbg(f"  Station {stationId} has override for field {fieldIndex} ({self.variableFields[fieldIndex]['id']}): '{fieldData[datetime]}'")
                            break

                    if not hasOverride:
                        continue  # Skip stations with no overrides for this datetime

                    # Starting field list with station ID
                    field_list = [stationId]

                    # IMPORTANT: Add only the first 4 static fields (name, lat, lon, elev)
                    # This ensures we don't accidentally include variable fields in static data
                    if stationId in self.static:
                        # Ensure we only take the first 4 static fields (name, lat, lon, elev)
                        for i in range(min(4, len(self.static[stationId]))):
                            field_list.append(self.static[stationId][i])
                            dbg(f"  Adding static field {i+1}: '{self.static[stationId][i]}'")

                        # If we have fewer than 4 static fields, add placeholders
                        while len(field_list) < 5:  # ID + 4 static fields = 5
                            field_list.append("Unknown")
                            dbg(f"  Adding placeholder static field")
                    else:
                        # Add placeholders for all static fields
                        field_list.extend(["Unknown"] * 4)  # 4 static fields
                        dbg(f"  Adding 4 placeholder static fields")

                    # Add variable data fields
                    for fieldIndex in range(self.numVariableFields):
                        field_id = self.variableFields[fieldIndex]["id"]
                        if stationId in self.man and fieldIndex < len(self.man[stationId]):
                            fieldData = self.man[stationId][fieldIndex]
                            if datetime in fieldData:
                                value = fieldData[datetime]
                                dbg(f"  Adding override for field {field_id}: '{value}'")
                            else:
                                value = "*"
                        else:
                            value = "*"  # Default value if no override
                        field_list.append(value)

                    # Join fields with commas
                    outline = ",".join(field_list) + "\n"

                    # Debug output
                    dbg(f"Writing line: {outline.strip()}")
                    dbg(f"Field count: {len(field_list)}")

                    datafile.write(outline)

            # Swap the finished file into place and record its new time
            self.manTimes[filename] = self._writeAtomic(tmpname, fullname)
            dbg(f"Successfully wrote manual override file: {fullname}")

        except Exception as e:
            self._discardTemp(tmpname)
            self.statusBar.showMessage(f"Error writing manual file: {e}")
            print(f"Error writing manual file: {e}")
            print(traceback.format_exc())

    def writeAllFile(self):
        """Write data to all.manual file"""
        fullname = os.path.join(DATADIR, "all.manual")

        if self._changedOnDisk(fullname, self.allTimes):
            if not self._confirmOverwrite(fullname):
                self.readAllFile(fullname)
                self.allTimes = os.path.getmtime(fullname)
                self.redisplay()
                return

        self.statusBar.showMessage(f"Writing {fullname}")

        # Back up the existing file before overwriting
        self._backupFile(fullname)

        tmpname = f"{fullname}.tmp{os.getpid()}"
        try:
            with open(tmpname, "w") as datafile:
                # Debug: Print the data structure
                dbg(f"Writing to {fullname}")
                dbg(f"Field order: ID, name, lat, lon, elev, tmp, dew, dir, spd, gst")

                for stationId in sorted(self.all.keys()):
                    # Check if there's any non-default data
                    hasNonDefaultData = False
                    for val in self.all[stationId]:
                        if val != "*":
                            hasNonDefaultData = True
                            break

                    if not hasNonDefaultData:
                        continue

                    # Starting field list with station ID
                    field_list = [stationId]

                    # IMPORTANT: Add only the first 4 static fields (name, lat, lon, elev)
                    # This ensures we don't accidentally include variable fields in static data
                    if stationId in self.static:
                        # Ensure we only take the first 4 static fields (name, lat, lon, elev)
                        for i in range(min(4, len(self.static[stationId]))):
                            field_list.append(self.static[stationId][i])
                            dbg(f"  Adding static field {i+1}: '{self.static[stationId][i]}'")

                        # If we have fewer than 4 static fields, add placeholders
                        while len(field_list) < 5:  # ID + 4 static fields = 5
                            field_list.append("Unknown")
                            dbg(f"  Adding placeholder static field")
                    else:
                        # Add placeholders for all static fields
                        field_list.extend(["Unknown"] * 4)  # 4 static fields
                        dbg(f"  Adding 4 placeholder static fields")

                    # Add variable data fields
                    for fieldIndex in range(self.numVariableFields):
                        if fieldIndex < len(self.all[stationId]):
                            value = self.all[stationId][fieldIndex]
                            dbg(f"  Adding override for field {self.variableFields[fieldIndex]['id']}: '{value}'")
                        else:
                            value = "*"  # Default value
                        field_list.append(value)

                    # Join fields with commas
                    outline = ",".join(field_list) + "\n"

                    # Debug output
                    dbg(f"Writing line: {outline.strip()}")
                    dbg(f"Field count: {len(field_list)}")

                    datafile.write(outline)

            # Swap the finished file into place and record its new time
            self.allTimes = self._writeAtomic(tmpname, fullname)
            dbg(f"Successfully wrote all.manual file: {fullname}")

        except Exception as e:
            self._discardTemp(tmpname)
            self.statusBar.showMessage(f"Error writing all.manual file: {e}")
            print(f"Error writing all.manual file: {e}")
            print(traceback.format_exc())

    def displayStation(self, stationId):
        """View data for a specific station"""
        try:
            if not stationId or stationId not in self.static:
                return

            # Recompute auto-detected issue flags for this redisplay.
            self._computeIssueFlags()

            # Update current station
            self.currentStation = stationId

            # Update station selector
            index = self.stationSelector.findText(stationId)
            if index >= 0:
                self.stationSelector.setCurrentIndex(index)

            # Update station info
            try:
                name = self.static[stationId][0]
                lat = float(self.static[stationId][1])
                lon = float(self.static[stationId][2])
                elev = int(float(self.static[stationId][3]))
                wfo = self.geom[stationId][5] if stationId in self.geom else "N/A"

                self.stationInfoLabel.setText(f"{wfo} - {name}\nElev: {elev} ft ({lat:.3f}, {lon:.3f})")
            except (ValueError, IndexError):
                self.stationInfoLabel.setText(f"{stationId}")

            # Clear the table.  clear() does not drop merged cells, so a span
            # set for a zone header on a previous render would otherwise land on
            # an unrelated data row and swallow its columns.
            self.dataTable.clear()
            self.dataTable.clearSpans()
            self.dataTable.setRowCount(0)

            # Set up table headers
            headers = ["Date/Time"]

            # Add field headers if checked
            for fieldId, checkbox in self.fieldCheckboxes.items():
                if checkbox.isChecked():
                    for field in self.variableFields:
                        if field["id"] == fieldId:
                            headers.append(field["shortLabel"])
                            break

            self.dataTable.setColumnCount(len(headers))
            self.dataTable.setHorizontalHeaderLabels(headers)

            # Get field data for the station
            fieldData = []
            for fieldIndex in range(self.numVariableFields):
                rawList, manList, allVal = self.getVals(stationId, fieldIndex)
                fieldData.append((rawList, manList, allVal))

            # Check for unchanged values - NEW CODE
            unchangedFields = {}  # Track which fields have unchanged values
            for fieldIndex, field in enumerate(self.variableFields):
                rawList, manList, allVal = fieldData[fieldIndex]

                # Skip fields with all-station overrides
                if allVal != "*":
                    continue

                # Collect all non-empty values for this field
                values = []
                for datetime in self.dates:
                    # Get the effective value considering raw data and overrides
                    value = ""
                    if datetime in rawList:
                        value = rawList[datetime]

                    # Apply manual overrides
                    if datetime in manList:
                        tmpVal = manList[datetime]
                        if tmpVal != "" and tmpVal != "*":
                            value = tmpVal

                    # Only consider non-empty values
                    if value.strip():
                        values.append(value)

                # Check if all values are the same and we have more than one observation
                if len(values) > 1 and all(v == values[0] for v in values):
                    unchangedFields[fieldIndex] = True

            # Add a row for each date
            row = 0
            for datetime in self.dates:
                # Add the row for every time period, even with no data
                self.dataTable.insertRow(row)

                # Set date/time cell
                timeStr = self.getTimeFromDateTime(datetime)
                timeItem = QTableWidgetItem(timeStr)
                timeItem.setTextAlignment(Qt.AlignLeft | Qt.AlignVCenter)
                timeItem.setData(Qt.UserRole, datetime)  # Store actual datetime for context menu
                self.dataTable.setItem(row, 0, timeItem)

                # Add variable fields if checked
                col = 1
                for fieldIndex, field in enumerate(self.variableFields):
                    if not self.fieldCheckboxes[field["id"]].isChecked():
                        continue

                    rawList, manList, allVal = fieldData[fieldIndex]

                    # Determine value and status
                    value = ""
                    status = "plain"

                    # First set raw value if available
                    if datetime in rawList:
                        value = rawList[datetime]

                    # Check manual overrides
                    if datetime in manList:
                        tmpVal = manList[datetime]
                        if tmpVal == "":
                            # For ignore, keep showing the original value but mark the status
                            status = "ignore"
                        elif tmpVal != "*":
                            value = tmpVal
                            status = "override"

                    # Check all-station overrides (these take precedence)
                    if allVal == "":
                        # For ignore-all, keep showing the original value but mark the status
                        status = "ignoreall"
                    elif allVal != "*":
                        value = allVal
                        status = "overrideall"

                    # Check for unchanged values - NEW CODE
                    # Only apply to fields with normal status (not ignored or overridden)
                    if status == "plain" and fieldIndex in unchangedFields and value.strip():
                        status = "unchanged"

                    # Apply auto-detected sanity / spike flags (only upgrades
                    # 'plain' or 'unchanged' — explicit overrides win visually)
                    status = self._autoStatusForCell(
                        stationId, datetime, field["id"], status
                    )

                    # Create cell item
                    valItem = QTableWidgetItem(f"{value:>3}")
                    valItem.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                    valItem.setBackground(self.dataTable.COLORS[status])

                    # Store field info for context menu
                    valItem.setData(Qt.UserRole, {"status": status, "field": field["id"]})

                    self.dataTable.setItem(row, col, valItem)
                    col += 1

                row += 1

            # Update title and status
            self.setWindowTitle(f"MatchObsAllQC 2.5 - {stationId}")
            self.statusBar.showMessage(f"Viewing data for {stationId}")

            # Set a reasonable width for the table
            self.dataTable.resizeColumnsToContents()
            tableWidth = 0
            for i in range(self.dataTable.columnCount()):
                tableWidth += self.dataTable.columnWidth(i)

            # Add a bit of padding
            tableWidth += 50
            self.dataTable.setMinimumWidth(tableWidth)
            self.dataTable.setMaximumWidth(tableWidth)
        except Exception as e:
            print(f"Error displaying station: {e}")
            print(traceback.format_exc())
            self.statusBar.showMessage(f"Error displaying station: {e}")

    def displayTime(self, datetime):
        """View data for a specific time"""
        try:
            if not datetime:
                return

            # Recompute auto-detected issue flags for this redisplay.
            self._computeIssueFlags()

            # Update current datetime
            self.currentDatetime = datetime

            # Update time selector
            index = self.timeSelector.findData(datetime)
            if index >= 0:
                self.timeSelector.setCurrentIndex(index)

            # Clear the table.  clear() does not drop merged cells, so a span
            # set for a zone header on a previous render would otherwise land on
            # an unrelated data row and swallow its columns.
            self.dataTable.clear()
            self.dataTable.clearSpans()
            self.dataTable.setRowCount(0)

            # Get the stations in order
            stations = self.getIdList(datetime)

            # Set up table headers
            headers = ["ID"]
            if self.showNameCb.isChecked():
                headers.append("Name")
            if self.showCwaCb.isChecked():
                headers.append("CWA")
            if self.showLatCb.isChecked():
                headers.append("Lat")
            if self.showLonCb.isChecked():
                headers.append("Lon")
            if self.showElevCb.isChecked():
                headers.append("Elev")

            # Add field headers if checked
            for fieldId, checkbox in self.fieldCheckboxes.items():
                if checkbox.isChecked():
                    for field in self.variableFields:
                        if field["id"] == fieldId:
                            headers.append(field["shortLabel"])
                            break

            self.dataTable.setColumnCount(len(headers))
            self.dataTable.setHorizontalHeaderLabels(headers)

            # Determine if we should show zone headers
            showZoneHeaders = len(self.availableCwas) > 0

            # Check for unchanged values across all dates - NEW CODE
            unchangedValues = {}  # Dictionary to track stations with unchanged values by field

            for stationId in stations:
                if stationId not in self.raw:
                    continue

                for fieldIndex, field in enumerate(self.variableFields):
                    if not self.fieldCheckboxes[field["id"]].isChecked():
                        continue

                    # Get values for this station and field
                    rawList, manList, allVal = self.getVals(stationId, fieldIndex)

                    # Skip if there's an all-station override
                    if allVal != "*":
                        continue

                    # Collect effective values across all dates
                    values = []
                    for date in self.dates:
                        # Get the effective value
                        value = ""
                        if date in rawList:
                            value = rawList[date]

                        # Apply manual overrides if any
                        if date in manList:
                            tmpVal = manList[date]
                            if tmpVal != "" and tmpVal != "*":
                                value = tmpVal

                        # Only consider non-empty values
                        if value.strip():
                            values.append(value)

                    # Check if all values are identical and we have more than one
                    if len(values) > 1 and all(v == values[0] for v in values):
                        if stationId not in unchangedValues:
                            unchangedValues[stationId] = []
                        unchangedValues[stationId].append(fieldIndex)

            # Apply the "Issues only" filter, if enabled. This trims the
            # station list down to only those with overrides, ignores,
            # unchanged values, sanity flags, or spike flags.
            if self.issuesOnlyCb.isChecked():
                stations = [
                    sid for sid in stations
                    if self._stationHasIssues(sid, unchangedValues, datetime)
                ]

            # Group dictionary to track current group
            currentGroup = None
            row = 0

            # Filter and add data rows
            for stationId in stations:
                if stationId not in self.geom:
                    continue

                # Get CWA for this station
                cwa = self.geom[stationId][5]

                group = self.geom[stationId][0]  # Zone
                wfo = self.geom[stationId][5]    # CWA

                # Add group header if needed and we're showing headers
                if showZoneHeaders and group != currentGroup:
                    groupName = self.geom[stationId][1]  # Zone name
                    self.dataTable.insertRow(row)

                    # Insert a two-column header
                    # First column: CWA
                    cwaCellItem = QTableWidgetItem(wfo)
                    cwaCellItem.setBackground(self.dataTable.COLORS["label"])
                    cwaCellItem.setForeground(Qt.white)
                    cwaCellItem.setTextAlignment(Qt.AlignCenter | Qt.AlignVCenter)
                    font = cwaCellItem.font()
                    font.setBold(True)
                    cwaCellItem.setFont(font)
                    self.dataTable.setItem(row, 0, cwaCellItem)

                    # Second column: Group - GroupName
                    groupCellItem = QTableWidgetItem(f"{group} - {groupName}")
                    groupCellItem.setBackground(self.dataTable.COLORS["label"])
                    groupCellItem.setForeground(Qt.white)
                    groupCellItem.setTextAlignment(Qt.AlignLeft | Qt.AlignVCenter)
                    groupCellItem.setFont(font)

                    # Span all remaining columns
                    numCols = self.dataTable.columnCount()
                    if numCols > 1:
                        self.dataTable.setSpan(row, 1, 1, numCols - 1)
                    self.dataTable.setItem(row, 1, groupCellItem)

                    currentGroup = group
                    row += 1

                # Add station data row
                self.dataTable.insertRow(row)

                # Check if this station has any actual raw data across ALL time periods
                hasAnyRawData = False

                # Check if this station has any overrides or ignores
                hasOverrides = False

                # Check for raw data across all time periods and all fields
                if stationId in self.raw:
                    for fieldIndex in range(self.numVariableFields):
                        if fieldIndex < len(self.raw[stationId]):
                            rawList = self.raw[stationId][fieldIndex]
                            for time_key, value in rawList.items():
                                if value.strip():
                                    hasAnyRawData = True
                                    break
                        if hasAnyRawData:
                            break

                # Check for manual overrides
                if stationId in self.man:
                    for fieldIndex in range(self.numVariableFields):
                        if fieldIndex < len(self.man[stationId]):
                            manList = self.man[stationId][fieldIndex]
                            if manList:  # If there are any manual overrides
                                hasOverrides = True
                                break

                # Check for all-station overrides
                if stationId in self.all:
                    for value in self.all[stationId]:
                        if value != "*":  # If there's any non-default value
                            hasOverrides = True
                            break

                # Highlight if station has overrides but no raw data
                highlightNoData = hasOverrides and not hasAnyRawData

                # Set ID cell
                idItem = QTableWidgetItem(stationId)
                idItem.setTextAlignment(Qt.AlignLeft | Qt.AlignVCenter)

                # Highlight station ID if conditions are met
                if highlightNoData:
                    idItem.setBackground(self.dataTable.COLORS["nodata"])

                self.dataTable.setItem(row, 0, idItem)

                col = 1

                # Add name if checked
                if self.showNameCb.isChecked():
                    name = self.static[stationId][0]
                    if len(name) > self.nameLengthValue:
                        name = name[:self.nameLengthValue]
                    nameItem = QTableWidgetItem(name)
                    nameItem.setTextAlignment(Qt.AlignLeft | Qt.AlignVCenter)

                    # Also highlight name if conditions are met
                    if highlightNoData:
                        nameItem.setBackground(self.dataTable.COLORS["nodata"])

                    self.dataTable.setItem(row, col, nameItem)
                    col += 1

                # Add CWA if checked
                if self.showCwaCb.isChecked():
                    cwa = self.geom[stationId][5] if stationId in self.geom else ""
                    cwaItem = QTableWidgetItem(cwa)
                    cwaItem.setTextAlignment(Qt.AlignCenter | Qt.AlignVCenter)

                    # Highlight CWA if conditions are met
                    if highlightNoData:
                        cwaItem.setBackground(self.dataTable.COLORS["nodata"])

                    self.dataTable.setItem(row, col, cwaItem)
                    col += 1

                # Add lat if checked
                if self.showLatCb.isChecked():
                    try:
                        lat = float(self.static[stationId][1])
                        latItem = QTableWidgetItem(f"{lat:.3f}")
                    except (ValueError, IndexError):
                        latItem = QTableWidgetItem("")
                    latItem.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)

                    # Highlight lat if conditions are met
                    if highlightNoData:
                        latItem.setBackground(self.dataTable.COLORS["nodata"])

                    self.dataTable.setItem(row, col, latItem)
                    col += 1

                # Add lon if checked
                if self.showLonCb.isChecked():
                    try:
                        lon = float(self.static[stationId][2])
                        lonItem = QTableWidgetItem(f"{lon:.3f}")
                    except (ValueError, IndexError):
                        lonItem = QTableWidgetItem("")
                    lonItem.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)

                    # Highlight lon if conditions are met
                    if highlightNoData:
                        lonItem.setBackground(self.dataTable.COLORS["nodata"])

                    self.dataTable.setItem(row, col, lonItem)
                    col += 1

                # Add elevation if checked
                if self.showElevCb.isChecked():
                    try:
                        elev = int(float(self.static[stationId][3]))
                        elevItem = QTableWidgetItem(f"{elev}")
                    except (ValueError, IndexError):
                        elevItem = QTableWidgetItem("")
                    elevItem.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)

                    # Highlight elevation if conditions are met
                    if highlightNoData:
                        elevItem.setBackground(self.dataTable.COLORS["nodata"])

                    self.dataTable.setItem(row, col, elevItem)
                    col += 1

                # Add variable fields if checked
                for fieldIndex, field in enumerate(self.variableFields):
                    if not self.fieldCheckboxes[field["id"]].isChecked():
                        continue

                    rawList, manList, allVal = self.getVals(stationId, fieldIndex)

                    # Determine value and status
                    value = ""
                    status = "plain"

                    # First set raw value if available
                    if datetime in rawList:
                        value = rawList[datetime]

                    # Check manual overrides
                    if datetime in manList:
                        tmpVal = manList[datetime]
                        if tmpVal == "":
                            # For ignore, keep showing the original value but mark the status
                            status = "ignore"
                        elif tmpVal != "*":
                            value = tmpVal
                            status = "override"

                    # Check all-station overrides (these take precedence)
                    if allVal == "":
                        # For ignore-all, keep showing the original value but mark the status
                        status = "ignoreall"
                    elif allVal != "*":
                        value = allVal
                        status = "overrideall"

                    # Check for unchanged values - NEW CODE
                    # Only apply to fields with normal status (not ignored or overridden)
                    if status == "plain" and stationId in unchangedValues and fieldIndex in unchangedValues[stationId] and value.strip():
                        status = "unchanged"

                    # Apply auto-detected sanity / spike flags (only upgrades
                    # 'plain' or 'unchanged' — explicit overrides win visually)
                    status = self._autoStatusForCell(
                        stationId, datetime, field["id"], status
                    )

                    # Create cell item
                    valItem = QTableWidgetItem(f"{value:>3}")
                    valItem.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                    valItem.setBackground(self.dataTable.COLORS[status])

                    # Store field info for context menu
                    valItem.setData(Qt.UserRole, {"status": status, "field": field["id"]})

                    self.dataTable.setItem(row, col, valItem)
                    col += 1

                row += 1

            # Update title and status
            formattedTime = self.getSplitTime(datetime)
            self.setWindowTitle(f"MatchObsAllQC 2.5 - {formattedTime}")
            self.statusBar.showMessage(f"Viewing data for {self.getTimeFromDateTime(datetime)}")

            # Set a reasonable width for the table
            self.dataTable.resizeColumnsToContents()
            tableWidth = 0
            for i in range(self.dataTable.columnCount()):
                tableWidth += self.dataTable.columnWidth(i)

            # Add a bit of padding
            tableWidth += 50
            self.dataTable.setMinimumWidth(tableWidth)
            self.dataTable.setMaximumWidth(tableWidth)
        except Exception as e:
            print(f"Error displaying time: {e}")
            print(traceback.format_exc())
            self.statusBar.showMessage(f"Error displaying time: {e}")

def main():
    """Main application entry point"""
    try:
        # Log diagnostic info
        print(f"Python version: {sys.version}")
        print(f"Current working directory: {os.getcwd()}")

        # Create application
        app = QApplication(sys.argv)

        # Set application style for a modern look
        app.setStyle("Fusion")

        # Create and show the main window
        mainWindow = MatchObsAllQC()
        mainWindow.show()

        # Run the application
        sys.exit(app.exec())
    except Exception as e:
        print(f"Fatal error in main application: {e}")
        print(traceback.format_exc())

        # Show error in GUI if possible
        try:
            if 'app' in locals():
                QMessageBox.critical(None, "Fatal Error",
                                   f"The application has encountered a fatal error and cannot continue:\n\n{str(e)}")
        except:
            pass

        sys.exit(1)

if __name__ == "__main__":
    main()
