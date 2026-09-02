import importlib.util
import sys
import types
import unittest
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_diurnal_module():
    for name in [
        "AbsTime",
        "DatabaseID",
        "LogStream",
        "ProcessVariableList",
        "TimeRange",
    ]:
        sys.modules[name] = types.ModuleType(name)

    smart_script = types.ModuleType("SmartScript")
    smart_script.SmartScript = object
    sys.modules["SmartScript"] = smart_script

    config = types.ModuleType("DiurnalConfig")
    config.configDict = {}
    sys.modules["DiurnalConfig"] = config

    spec = importlib.util.spec_from_file_location(
        "diurnal_under_test", REPO_ROOT / "Diurnal.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DIURNAL = load_diurnal_module()


class FakeAbsTime:
    def __init__(self, value):
        self._value = value

    def unixTime(self):
        return self._value


class FakeTimeRange:
    def __init__(self, start, end):
        self._start = FakeAbsTime(start)
        self._end = FakeAbsTime(end)

    def startTime(self):
        return self._start

    def endTime(self):
        return self._end

    def overlaps(self, other):
        return (
            self._start.unixTime() < other.endTime().unixTime()
            and other.startTime().unixTime() < self._end.unixTime()
        )

    def __hash__(self):
        return hash((self._start.unixTime(), self._end.unixTime()))

    def __eq__(self, other):
        return isinstance(other, FakeTimeRange) and (
            self._start.unixTime(),
            self._end.unixTime(),
        ) == (other.startTime().unixTime(), other.endTime().unixTime())


class FakeDatabase:
    def __init__(self, name, model_time=0, identifier=None):
        self._name = name
        self._model_time = FakeAbsTime(model_time)
        self._identifier = identifier or name

    def modelName(self):
        return self._name

    def modelTime(self):
        return self._model_time

    def modelIdentifier(self):
        return self._identifier

    def __repr__(self):
        return "FakeDatabase({!r}, {!r})".format(
            self._name, self._model_time.unixTime()
        )


class FakeGridInfo:
    def __init__(self, time_range):
        self._time_range = time_range

    def gridTime(self):
        return self._time_range


class FakeParm:
    def __init__(self, name):
        self.name = name


def new_procedure():
    return DIURNAL.Procedure.__new__(DIURNAL.Procedure)


class ObservationAverageTests(unittest.TestCase):
    def test_each_utc_hour_has_an_independent_average(self):
        procedure = new_procedure()
        ranges = [
            FakeTimeRange(0, 3600),
            FakeTimeRange(24 * 3600, 25 * 3600),
            FakeTimeRange(3600, 7200),
        ]
        grids = {
            ranges[0]: np.array([[10.0]]),
            ranges[1]: np.array([[20.0]]),
            ranges[2]: np.array([[100.0]]),
        }
        procedure._excludeTRList = []
        procedure.getWEInventory = lambda *args: ranges
        procedure.getGrids = lambda *args, **kwargs: grids[args[3]]
        procedure.empty = lambda: np.zeros((1, 1), dtype=float)

        averages = procedure.calcObsHourlyAvg("Obs", "T")

        np.testing.assert_array_equal(averages[0], np.array([[15.0]]))
        np.testing.assert_array_equal(averages[1], np.array([[100.0]]))
        self.assertIsNone(averages[2])

    def test_hours_without_observations_are_not_published(self):
        procedure = new_procedure()
        averages = [None] * 24
        averages[0] = np.array([[42.0]])
        procedure.configDict = {"obsModelList": ["Obs"]}
        procedure.calcObsHourlyAvg = lambda *args: averages
        procedure.makeTimeRange = FakeTimeRange
        procedure.statusBarMsg = lambda *args: None

        result = procedure.makeHourlyGrids(
            ["Obs", "GFS"], "T", FakeTimeRange(0, 2 * 3600)
        )

        self.assertEqual(list(result), [FakeTimeRange(0, 3600)])
        np.testing.assert_array_equal(result[FakeTimeRange(0, 3600)], [[42.0]])


class SplineTests(unittest.TestCase):
    def test_first_source_knot_is_reproduced_exactly(self):
        procedure = new_procedure()
        grids = [np.array([[value]], dtype=float) for value in (0, 1, 4, 9)]

        result = procedure._cubicSpline(grids, [0, 1, 2, 3], [0])

        np.testing.assert_allclose(result[0], grids[0])

    def test_one_and_two_source_grids_degrade_to_constant_and_linear(self):
        procedure = new_procedure()
        one = procedure._cubicSpline([np.array([[7.0]])], [0], [0, 1])
        two = procedure._cubicSpline(
            [np.array([[0.0]]), np.array([[6.0]])],
            [0, 3],
            [0, 1, 2, 3],
        )

        np.testing.assert_allclose([grid.item() for grid in one], [7.0, 7.0])
        np.testing.assert_allclose(
            [grid.item() for grid in two], [0.0, 2.0, 4.0, 6.0]
        )

    def test_hourly_output_includes_the_final_source_time(self):
        procedure = new_procedure()
        first = FakeTimeRange(0, 3600)
        last = FakeTimeRange(3 * 3600, 4 * 3600)
        db = object()
        grids = {first: np.array([[0.0]]), last: np.array([[6.0]])}
        procedure.configDict = {"obsModelList": []}
        procedure.getModelInventory = lambda *args: [(first, db), (last, db)]
        procedure.getGrids = lambda *args, **kwargs: grids[args[3]]
        procedure.makeTimeRange = FakeTimeRange

        result = procedure.makeHourlyGrids(
            ["GFS"], "T", FakeTimeRange(0, 4 * 3600)
        )

        self.assertIn(last, result)
        np.testing.assert_allclose(result[last], [[6.0]])

    def test_three_hour_guidance_includes_end_knot_for_complete_period(self):
        procedure = new_procedure()
        database = FakeDatabase("GFS", 24 * 3600)
        inventory = [
            FakeTimeRange((24 + hour) * 3600, (25 + hour) * 3600)
            for hour in (0, 3, 6, 9, 12)
        ]
        grids = {
            tr: np.array([[float(index)]]) for index, tr in enumerate(inventory)
        }
        procedure.configDict = {
            "obsModelList": [],
            "highlightColors": ["yellow"],
        }
        procedure.getAvailableDBIDs = lambda name: [database]
        procedure.getWEInventory = lambda db, *args: inventory
        procedure.getWEFrequency = lambda *args: 3
        procedure.getGrids = lambda db, name, level, tr, **kwargs: grids[tr]
        procedure.makeTimeRange = FakeTimeRange
        procedure.expandTimeRange = lambda tr, step: tr
        procedure.statusBarMsg = lambda *args: None

        result = procedure.makeHourlyGrids(
            ["GFS"], "T", FakeTimeRange(24 * 3600, 36 * 3600)
        )
        procedure._modelHrGridDict = result

        self.assertIn(FakeTimeRange(36 * 3600, 37 * 3600), result)
        self.assertTrue(
            procedure._hourlyCoverageComplete(
                FakeTimeRange(24 * 3600, 36 * 3600)
            )
        )


class ExtremaTests(unittest.TestCase):
    def test_hourly_model_extrema_keep_curve_anchored_to_official_values(self):
        procedure = new_procedure()
        period = FakeTimeRange(0, 2 * 3600)
        procedure._modelHrGridDict = {
            FakeTimeRange(0, 3600): np.array([[10.0]]),
            FakeTimeRange(3600, 7200): np.array([[15.0]]),
        }
        procedure.configDict = {"useHourlyMinMax": True}
        procedure.getWEInventory = lambda model, name: [period]
        procedure.getGrids = lambda model, name, *args, **kwargs: {
            "MinT": np.array([[5.0]]),
            "MaxT": np.array([[20.0]]),
        }[name]

        minimums, maximums = procedure.getMinMaxGrids(
            "GFS", "T", [period], [period]
        )

        np.testing.assert_array_equal(minimums[period], [[10.0]])
        np.testing.assert_array_equal(maximums[period], [[15.0]])

    def test_abutting_hourly_inventory_does_not_return_a_sentinel_grid(self):
        procedure = new_procedure()
        procedure._modelHrGridDict = {
            FakeTimeRange(3600, 7200): np.array([[12.0]])
        }
        procedure.newGrid = lambda value: np.full((1, 1), value, dtype=float)

        result = procedure.calcMinMaxFromHourlies(
            FakeTimeRange(0, 3600), "min"
        )

        self.assertIsNone(result)

    def test_partial_hourly_period_is_not_accepted_as_daily_extremum(self):
        procedure = new_procedure()
        procedure._modelHrGridDict = {
            FakeTimeRange(hour * 3600, (hour + 1) * 3600): np.array([[hour]])
            for hour in range(6, 12)
        }

        result = procedure.calcMinMaxFromHourlies(
            FakeTimeRange(0, 12 * 3600), "min"
        )

        self.assertIsNone(result)

    def test_time_before_first_extreme_never_wraps_to_last_period(self):
        procedure = new_procedure()
        first = FakeTimeRange(10 * 3600, 20 * 3600)
        last = FakeTimeRange(30 * 3600, 40 * 3600)
        procedure._modelMinGridDict = {first: np.array([[1]]), last: np.array([[2]])}

        result = procedure.calcAdjacentTimeRanges(
            "min", FakeTimeRange(5 * 3600, 6 * 3600)
        )

        self.assertEqual(result, [first, first])

    def test_same_day_execution_still_validates_min_and_max_inventory(self):
        procedure = new_procedure()
        messages = []
        procedure.getParm = lambda db, name, level: FakeParm(name)
        procedure.getParmTimeConstraints = lambda parm: (0, 12 * 3600)
        procedure.makeTimeRange = FakeTimeRange
        procedure.statusBarMsg = lambda message, level: messages.append(message)

        result = procedure.checkFcstMinMaxGrids(
            "T", FakeTimeRange(3600, 2 * 3600), {}, {}
        )

        self.assertFalse(result)
        self.assertTrue(any("MinT" in message for message in messages))
        self.assertTrue(any("MaxT" in message for message in messages))

    def test_only_the_extrema_periods_that_can_bracket_output_are_required(self):
        procedure = new_procedure()
        min_periods = [
            FakeTimeRange(0, 12 * 3600),
            FakeTimeRange(24 * 3600, 36 * 3600),
        ]
        max_periods = [
            FakeTimeRange(-12 * 3600, 0),
            FakeTimeRange(12 * 3600, 24 * 3600),
        ]
        procedure.getParm = lambda db, name, level: FakeParm(name)
        procedure.getParmTimeConstraints = lambda parm: (
            (0, 12 * 3600)
            if parm.name == "MinT"
            else (12 * 3600, 12 * 3600)
        )
        procedure.getParmMinMaxLimits = lambda *args: (-100.0, 150.0)
        procedure.makeTimeRange = FakeTimeRange
        procedure.statusBarMsg = lambda *args: None
        min_grids = {tr: np.array([[30.0]]) for tr in min_periods}
        max_grids = {tr: np.array([[80.0]]) for tr in max_periods}

        result = procedure.checkFcstMinMaxGrids(
            "T", FakeTimeRange(0, 24 * 3600), min_grids, max_grids
        )

        self.assertTrue(result)

    def test_missing_leading_bracket_is_rejected(self):
        procedure = new_procedure()
        current_min = FakeTimeRange(0, 12 * 3600)
        current_max = FakeTimeRange(12 * 3600, 24 * 3600)
        messages = []
        procedure.getParm = lambda db, name, level: FakeParm(name)
        procedure.getParmTimeConstraints = lambda parm: (
            (0, 12 * 3600)
            if parm.name == "MinT"
            else (12 * 3600, 12 * 3600)
        )
        procedure.makeTimeRange = FakeTimeRange
        procedure.statusBarMsg = lambda message, level: messages.append(message)

        result = procedure.checkFcstMinMaxGrids(
            "T",
            FakeTimeRange(0, 3600),
            {current_min: np.array([[30.0]])},
            {current_max: np.array([[80.0]])},
        )

        self.assertFalse(result)
        self.assertTrue(any("MaxT" in message for message in messages))

    def test_extrema_qc_detects_selected_sentinel_and_nonfinite_cells(self):
        procedure = new_procedure()
        selected = np.array([[True, True]], dtype=bool)

        self.assertIsNotNone(
            procedure._extremaGridError(
                np.array([[-100.0, -100.0]]), -100.0, selected
            )
        )
        self.assertIsNotNone(
            procedure._extremaGridError(
                np.array([[np.nan, 40.0]]), -100.0, selected
            )
        )
        self.assertIsNone(
            procedure._extremaGridError(
                np.array([[40.0, -100.0]]), -100.0, selected
            )
        )
        self.assertIsNone(
            procedure._extremaGridError(
                np.array([[0.0, 10.0]]), 0.0, selected
            )
        )

    def test_fallback_hourlies_do_not_mix_with_primary_daily_extrema(self):
        procedure = new_procedure()
        period = FakeTimeRange(0, 2 * 3600)
        procedure._modelHrGridDict = {
            FakeTimeRange(0, 3600): np.array([[10.0]]),
            FakeTimeRange(3600, 7200): np.array([[15.0]]),
        }
        procedure._hourlyUsedFallback = True
        procedure.configDict = {"useHourlyMinMax": False}
        procedure.getWEInventory = lambda model, name: [period]
        procedure.getGrids = lambda model, name, *args, **kwargs: {
            "MinT": np.array([[5.0]]),
            "MaxT": np.array([[20.0]]),
        }[name]

        minimums, maximums = procedure.getMinMaxGrids(
            "GFS", "T", [period], [period]
        )

        np.testing.assert_array_equal(minimums[period], [[10.0]])
        np.testing.assert_array_equal(maximums[period], [[15.0]])

    def test_incomplete_fallback_hourlies_reject_primary_daily_extrema(self):
        procedure = new_procedure()
        period = FakeTimeRange(0, 3 * 3600)
        procedure._modelHrGridDict = {
            FakeTimeRange(3600, 2 * 3600): np.array([[12.0]])
        }
        procedure._hourlyUsedFallback = True
        procedure.configDict = {"useHourlyMinMax": False}
        procedure.getWEInventory = lambda model, name: [period]
        procedure.getGrids = lambda model, name, *args, **kwargs: {
            "MinT": np.array([[5.0]]),
            "MaxT": np.array([[20.0]]),
        }[name]

        minimums, maximums = procedure.getMinMaxGrids(
            "GFS", "T", [period], [period]
        )

        self.assertEqual(minimums, {})
        self.assertEqual(maximums, {})

    def test_missing_official_extrema_are_not_synthesized_from_model_hourlies(self):
        procedure = new_procedure()
        period = FakeTimeRange(0, 2 * 3600)
        procedure._modelHrGridDict = {
            FakeTimeRange(0, 3600): np.array([[10.0]]),
            FakeTimeRange(3600, 7200): np.array([[15.0]]),
        }
        procedure.configDict = {"useHourlyMinMax": True}
        procedure.getWEInventory = lambda *args: []

        minimums, maximums = procedure.getMinMaxGrids(
            "Fcst", "T", [period], [period]
        )

        self.assertEqual(minimums, {})
        self.assertEqual(maximums, {})

    def test_stale_optional_daily_inventory_uses_complete_hourly_extrema(self):
        procedure = new_procedure()
        period = FakeTimeRange(0, 2 * 3600)
        procedure._modelHrGridDict = {
            FakeTimeRange(0, 3600): np.array([[10.0]]),
            FakeTimeRange(3600, 7200): np.array([[15.0]]),
        }
        procedure.configDict = {"useHourlyMinMax": True}
        procedure.getWEInventory = lambda *args: [period]
        procedure.getGrids = lambda *args, **kwargs: (_ for _ in ()).throw(
            RuntimeError("stale daily inventory")
        )

        minimums, maximums = procedure.getMinMaxGrids(
            "GFS", "T", [period], [period]
        )

        np.testing.assert_array_equal(minimums[period], [[10.0]])
        np.testing.assert_array_equal(maximums[period], [[15.0]])

    def test_stale_official_inventory_is_returned_as_invalid_not_raised(self):
        procedure = new_procedure()
        period = FakeTimeRange(0, 2 * 3600)
        procedure._modelHrGridDict = {}
        procedure.configDict = {"useHourlyMinMax": True}
        procedure.getWEInventory = lambda *args: [period]
        procedure.getGrids = lambda *args, **kwargs: (_ for _ in ()).throw(
            RuntimeError("stale Official inventory")
        )

        minimums, maximums = procedure.getMinMaxGrids(
            "Fcst", "T", [period], [period]
        )

        self.assertIsNone(minimums[period])
        self.assertIsNone(maximums[period])

    def test_model_curve_is_rescaled_to_reach_official_extrema(self):
        procedure = new_procedure()
        first_hour = FakeTimeRange(0, 3600)
        second_hour = FakeTimeRange(3600, 7200)
        extrema_period = FakeTimeRange(0, 7200)
        procedure._modelHrGridDict = {
            first_hour: np.array([[10.0]]),
            second_hour: np.array([[20.0]]),
        }
        procedure._modelMinGridDict = {extrema_period: np.array([[10.0]])}
        procedure._modelMaxGridDict = {extrema_period: np.array([[20.0]])}
        procedure._fcstMinGridDict = {extrema_period: np.array([[5.0]])}
        procedure._fcstMaxGridDict = {extrema_period: np.array([[25.0]])}
        procedure._fTrendDict = {
            first_hour: np.array([[1]]),
            second_hour: np.array([[-1]]),
        }
        procedure._smoothAreas = []
        procedure.configDict = {"DEBUG": False}
        procedure.getParmMinMaxLimits = lambda *args: (-100.0, 150.0)
        procedure.highlightGrids = lambda *args: None
        procedure.getGrids = lambda *args, **kwargs: (_ for _ in ()).throw(
            KeyError("no prior hourly grid")
        )
        created = []
        procedure.createGrid = lambda *args, **kwargs: created.append(args)

        procedure.calcDiurnalGrids(
            "T",
            FakeTimeRange(0, 7200),
            np.ones((1, 1), dtype=bool),
            np.ones((1, 1), dtype=bool),
        )

        np.testing.assert_array_equal(created[0][3], [[5.0]])
        np.testing.assert_array_equal(created[1][3], [[25.0]])


class ModelInventoryTests(unittest.TestCase):
    def test_missing_element_frequency_with_database_id_returns_zero(self):
        procedure = new_procedure()
        database = FakeDatabase("NBM")
        procedure.getParm = lambda *args: (_ for _ in ()).throw(KeyError("RH"))
        procedure.statusBarMsg = lambda *args: None

        self.assertEqual(procedure.getWEFrequency(database, "RH"), 0)

    def test_database_versions_are_returned_newest_first(self):
        procedure = new_procedure()
        old = FakeDatabase("GFS", 100)
        newest = FakeDatabase("GFS", 300)
        middle = FakeDatabase("GFS", 200)
        d2d = FakeDatabase("GFS", 400, identifier="D2D_GFS")
        procedure.availableDatabases = lambda: [old, newest, d2d, middle]

        result = procedure.getAvailableDBIDs("GFS")

        self.assertEqual(result, [newest, middle, old])

    def test_secondary_supplies_inventory_when_primary_database_is_absent(self):
        procedure = new_procedure()
        secondary = FakeDatabase("GFS")
        inventory = [
            FakeTimeRange(0, 3600),
            FakeTimeRange(3 * 3600, 4 * 3600),
        ]
        procedure.configDict = {"highlightColors": ["yellow", "orange", "red"]}
        procedure.getAvailableDBIDs = lambda name: {
            "NBM": [],
            "GFS": [secondary],
        }.get(name, [])
        procedure.getWEInventory = lambda db, *args: inventory if db is secondary else []
        procedure.getWEFrequency = lambda db, *args: 3 if db is secondary else 0
        procedure.getGrids = lambda *args, **kwargs: np.array([[1.0]])
        procedure.makeTimeRange = FakeTimeRange
        procedure.expandTimeRange = lambda tr, step: tr
        procedure.statusBarMsg = lambda *args: None

        result = procedure.getModelInventory(
            ["NBM", "GFS"], "T", FakeTimeRange(0, 4 * 3600)
        )

        self.assertEqual(result, [(inventory[0], secondary), (inventory[1], secondary)])
        self.assertEqual(procedure._highlightDict[inventory[0]][0], "yellow")

    def test_secondary_inventory_does_not_fill_between_primary_cadence_knots(self):
        procedure = new_procedure()
        primary = FakeDatabase("NBM")
        secondary = FakeDatabase("GFS")
        primary_inventory = [
            FakeTimeRange(0, 3600),
            FakeTimeRange(3 * 3600, 4 * 3600),
        ]
        secondary_inventory = [
            FakeTimeRange(hour * 3600, (hour + 1) * 3600)
            for hour in range(4)
        ]
        procedure.configDict = {"highlightColors": ["yellow", "orange"]}
        procedure.getAvailableDBIDs = lambda name: {
            "NBM": [primary],
            "GFS": [secondary],
        }[name]
        procedure.getWEInventory = lambda db, *args: (
            primary_inventory if db is primary else secondary_inventory
        )
        procedure.getWEFrequency = lambda db, *args: 3 if db is primary else 1
        procedure.getGrids = lambda *args, **kwargs: np.array([[1.0]])
        procedure.makeTimeRange = FakeTimeRange
        procedure.expandTimeRange = lambda tr, step: tr
        procedure.statusBarMsg = lambda *args: None

        result = procedure.getModelInventory(
            ["NBM", "GFS"], "T", FakeTimeRange(0, 4 * 3600)
        )

        self.assertEqual(result, [(tr, primary) for tr in primary_inventory])

    def test_stale_primary_inventory_retries_the_same_time_from_fallback(self):
        procedure = new_procedure()
        primary = FakeDatabase("NBM")
        secondary = FakeDatabase("GFS")
        inventory = [
            FakeTimeRange(0, 3600),
            FakeTimeRange(3 * 3600, 4 * 3600),
        ]
        calls = []
        procedure.configDict = {
            "obsModelList": [],
            "highlightColors": ["yellow", "orange"],
        }
        procedure.getAvailableDBIDs = lambda name: {
            "NBM": [primary],
            "GFS": [secondary],
        }[name]
        procedure.getWEInventory = lambda *args: inventory
        procedure.getWEFrequency = lambda *args: 3

        def get_grid(db, name, level, tr, **kwargs):
            calls.append(db)
            if db is primary:
                return None
            return np.array([[float(tr.startTime().unixTime() / 3600)]])

        procedure.getGrids = get_grid
        procedure.makeTimeRange = FakeTimeRange
        procedure.expandTimeRange = lambda tr, step: tr
        procedure.statusBarMsg = lambda *args: None

        result = procedure.makeHourlyGrids(
            ["NBM", "GFS"], "T", FakeTimeRange(0, 4 * 3600)
        )

        self.assertIn(FakeTimeRange(0, 3600), result)
        self.assertIn(FakeTimeRange(3 * 3600, 4 * 3600), result)
        self.assertEqual(calls, [primary, secondary, primary, secondary])
        self.assertTrue(procedure._hourlyUsedFallback)

    def test_coarser_fallback_knots_bracket_requested_hourly_period(self):
        procedure = new_procedure()
        primary = FakeDatabase("NBM")
        secondary = FakeDatabase("GFS")
        primary_inventory = [FakeTimeRange(0, 3600)]
        fallback_inventory = [
            FakeTimeRange(12 * 3600, 13 * 3600),
            FakeTimeRange(15 * 3600, 16 * 3600),
        ]
        fallback_grids = {
            fallback_inventory[0]: np.array([[0.0]]),
            fallback_inventory[1]: np.array([[3.0]]),
        }
        procedure.configDict = {
            "obsModelList": [],
            "highlightColors": ["yellow", "orange"],
        }
        procedure.getAvailableDBIDs = lambda name: {
            "NBM": [primary],
            "GFS": [secondary],
        }[name]
        procedure.getWEInventory = lambda db, *args: (
            primary_inventory if db is primary else fallback_inventory
        )
        procedure.getWEFrequency = lambda db, *args: 1 if db is primary else 3
        procedure.getGrids = lambda db, name, level, tr, **kwargs: fallback_grids[tr]
        procedure.makeTimeRange = FakeTimeRange
        procedure.expandTimeRange = lambda tr, step: tr
        procedure.statusBarMsg = lambda *args: None

        result = procedure.makeHourlyGrids(
            ["NBM", "GFS"],
            "T",
            FakeTimeRange(13 * 3600, 14 * 3600),
        )

        np.testing.assert_allclose(result[FakeTimeRange(13 * 3600, 14 * 3600)], [[1.0]])
        np.testing.assert_allclose(result[FakeTimeRange(14 * 3600, 15 * 3600)], [[2.0]])

    def test_timeline_aligns_to_model_cadence_instead_of_subtracting_six_hours(self):
        procedure = new_procedure()
        base = 10 * 24 * 3600
        latest = FakeDatabase("GFS", base + 12 * 3600)
        older = FakeDatabase("GFS", base)
        older_inventory = [
            FakeTimeRange(base + hour * 3600, base + (hour + 1) * 3600)
            for hour in (0, 3, 6, 9)
        ]
        latest_inventory = [
            FakeTimeRange(base + 12 * 3600, base + 13 * 3600)
        ]
        procedure.configDict = {"highlightColors": ["yellow", "orange", "red"]}
        procedure.getAvailableDBIDs = lambda name: [latest, older]
        procedure.getWEInventory = lambda db, *args: (
            latest_inventory if db is latest else older_inventory
        )
        procedure.getWEFrequency = lambda *args: 3
        procedure.getGrids = lambda *args, **kwargs: np.array([[1.0]])
        procedure.makeTimeRange = FakeTimeRange
        procedure.expandTimeRange = lambda tr, step: tr
        procedure.statusBarMsg = lambda *args: None

        result = procedure.getModelInventory(
            ["GFS"], "T", FakeTimeRange(base, base + 15 * 3600)
        )

        self.assertEqual(
            result,
            [(tr, older) for tr in older_inventory]
            + [(latest_inventory[0], latest)],
        )
        self.assertTrue(procedure._hourlyUsedFallback)

    def test_missing_expected_model_time_is_marked_for_fallback(self):
        procedure = new_procedure()
        database = FakeDatabase("GFS", 24 * 3600)
        inventory = [
            FakeTimeRange(24 * 3600, 25 * 3600),
            FakeTimeRange(27 * 3600, 28 * 3600),
        ]
        procedure.configDict = {"highlightColors": ["yellow"]}
        procedure.getAvailableDBIDs = lambda name: [database]
        procedure.getWEInventory = lambda *args: inventory
        procedure.getWEFrequency = lambda *args: 3
        procedure.getGrids = lambda *args, **kwargs: np.array([[1.0]])
        procedure.makeTimeRange = FakeTimeRange
        procedure.expandTimeRange = lambda tr, step: tr
        procedure.statusBarMsg = lambda *args: None

        procedure.getModelInventory(
            ["GFS"], "T", FakeTimeRange(24 * 3600, 33 * 3600)
        )

        self.assertEqual(
            procedure._highlightDict[FakeTimeRange(30 * 3600, 31 * 3600)],
            ("red", "Missing guidance"),
        )


class ConfigurationTests(unittest.TestCase):
    def test_local_config_can_add_keys_not_present_in_base_defaults(self):
        procedure = new_procedure()
        procedure.configDict = {"known": 1}
        procedure.statusBarMsg = lambda *args: None
        original = DIURNAL.DiurnalConfig.configDict
        try:
            DIURNAL.DiurnalConfig.configDict = {"fTable": [[1, 2, 3]]}
            procedure.replaceConfigValues()
        finally:
            DIURNAL.DiurnalConfig.configDict = original

        self.assertEqual(procedure.configDict["fTable"], [[1, 2, 3]])

    def test_site_config_anchors_curve_to_official_extrema(self):
        spec = importlib.util.spec_from_file_location(
            "diurnal_config_under_test", REPO_ROOT / "DiurnalConfig.py"
        )
        config_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(config_module)

        self.assertTrue(config_module.configDict["useHourlyMinMax"])


class TimeRangeTests(unittest.TestCase):
    def test_dst_is_evaluated_at_each_named_endpoint(self):
        procedure = new_procedure()
        winter = int(datetime(2026, 1, 1, tzinfo=timezone.utc).timestamp())
        summer = int(datetime(2026, 7, 1, tzinfo=timezone.utc).timestamp())
        procedure.configDict = {
            "adjustDST": True,
            "timeZone": "America/Boise",
        }
        procedure.makeTimeRange = FakeTimeRange
        procedure.getTimeRange = lambda name: {
            "start": FakeTimeRange(winter, winter + 3600),
            "end": FakeTimeRange(summer - 3600, summer),
        }[name]
        procedure.statusBarMsg = lambda *args: None

        result = procedure.makeExecuteTR(
            "start", "end", FakeTimeRange(0, 1)
        )

        self.assertEqual(result, FakeTimeRange(winter, summer + 3600))

    def test_invalid_named_range_is_rejected_before_time_range_construction(self):
        procedure = new_procedure()
        procedure.configDict = {"adjustDST": False}
        procedure.getTimeRange = lambda name: {
            "start": FakeTimeRange(3600, 7200),
            "end": FakeTimeRange(0, 1800),
        }[name]
        procedure.makeTimeRange = lambda *args: self.fail(
            "an invalid range must not be constructed"
        )
        messages = []
        procedure.statusBarMsg = lambda message, level: messages.append(message)

        result = procedure.makeExecuteTR("start", "end", FakeTimeRange(0, 1))

        self.assertIsNone(result)
        self.assertTrue(any("Ending time" in message for message in messages))

    def test_unknown_named_period_returns_a_clear_error(self):
        procedure = new_procedure()
        procedure.configDict = {"adjustDST": False}
        procedure.getTimeRange = lambda name: None
        procedure.makeTimeRange = lambda *args: self.fail(
            "an unresolved named period must not be constructed"
        )
        messages = []
        procedure.statusBarMsg = lambda message, level: messages.append(message)

        result = procedure.makeExecuteTR(
            "Misspelled Period", "Selected TimeRange", FakeTimeRange(0, 3600)
        )

        self.assertIsNone(result)
        self.assertTrue(any("Misspelled Period" in message for message in messages))


class ClimatologyFallbackTests(unittest.TestCase):
    def test_disjoint_fallback_spans_do_not_overwrite_primary_guidance(self):
        procedure = new_procedure()
        fallback_one = FakeTimeRange(0, 3600)
        primary = FakeTimeRange(3600, 7200)
        fallback_two = FakeTimeRange(7200, 10800)
        procedure._highlightDict = {
            fallback_one: ("yellow", "02/00Z GFS"),
            primary: ("gray", "02/00Z NBM"),
            fallback_two: ("yellow", "02/00Z GFS"),
        }
        procedure._modelHrGridDict = {
            fallback_one: np.array([[1.0]]),
            primary: np.array([[1.0]]),
            fallback_two: np.array([[1.0]]),
        }
        procedure.makeTimeRange = FakeTimeRange
        calls = []
        procedure.callProcedure = lambda *args, **kwargs: calls.append(kwargs)

        procedure.overwriteWithClimo(
            "NBM", FakeTimeRange(0, 10800), "Entire Domain"
        )

        self.assertEqual(
            [call["timeRange"] for call in calls],
            [fallback_one, fallback_two],
        )

    def test_uncovered_tail_is_filled_with_climatology(self):
        procedure = new_procedure()
        covered = FakeTimeRange(0, 3600)
        procedure._highlightDict = {}
        procedure._modelHrGridDict = {covered: np.array([[1.0]])}
        procedure.makeTimeRange = FakeTimeRange
        calls = []
        procedure.callProcedure = lambda *args, **kwargs: calls.append(kwargs)

        procedure.overwriteWithClimo(
            "NBM", FakeTimeRange(0, 3 * 3600), "Entire Domain"
        )

        self.assertEqual(
            [call["timeRange"] for call in calls],
            [FakeTimeRange(3600, 3 * 3600)],
        )

    def test_multiword_primary_name_is_not_mistaken_for_fallback(self):
        procedure = new_procedure()
        covered = FakeTimeRange(0, 3600)
        procedure._highlightDict = {covered: ("gray", "02/00Z My Model")}
        procedure._modelHrGridDict = {covered: np.array([[1.0]])}
        procedure.makeTimeRange = FakeTimeRange
        calls = []
        procedure.callProcedure = lambda *args, **kwargs: calls.append(kwargs)

        procedure.overwriteWithClimo(
            "My Model", FakeTimeRange(0, 3600), "Entire Domain"
        )

        self.assertEqual(calls, [])

    def test_climatology_fallback_preserves_caller_options(self):
        procedure = new_procedure()
        fallback = FakeTimeRange(0, 3600)
        procedure._highlightDict = {fallback: ("yellow", "02/00Z GFS")}
        procedure._modelHrGridDict = {fallback: np.array([[1.0]])}
        procedure.makeTimeRange = FakeTimeRange
        calls = []
        procedure.callProcedure = lambda *args, **kwargs: calls.append(kwargs)

        procedure.overwriteWithClimo(
            "NBM",
            fallback,
            "Entire Domain",
            {"Climo-Only Areas": ["Snake Plain"], "Custom Option": 7},
        )

        self.assertEqual(calls[0]["varDict"]["Climo-Only Areas"], ["Snake Plain"])
        self.assertEqual(calls[0]["varDict"]["Custom Option"], 7)
        self.assertEqual(calls[0]["varDict"]["Beginning"], "Selected TimeRange")


class SmoothingTests(unittest.TestCase):
    def test_smoothing_changes_only_selected_cells(self):
        procedure = new_procedure()
        procedure.statusBarMsg = lambda *args: None
        grid = np.zeros((3, 3), dtype=float)
        grid[1, 1] = 16.0
        mask = np.zeros((3, 3), dtype=bool)
        mask[1, 1] = True

        result = procedure.smoothGrid(grid, 1, mask)

        self.assertEqual(result[1, 1], 4.0)
        np.testing.assert_array_equal(result[~mask], grid[~mask])

    def test_nonfinite_and_excessive_weights_are_rejected(self):
        procedure = new_procedure()
        procedure.configDict = {"maxSmoothPasses": 5}
        procedure.statusBarMsg = lambda *args: None
        grid = np.arange(9.0).reshape((3, 3))
        mask = np.ones((3, 3), dtype=bool)

        for weight in (float("nan"), float("inf"), 1000):
            with self.subTest(weight=weight):
                np.testing.assert_array_equal(
                    procedure.smoothGrid(grid, weight, mask), grid
                )


class ExecutionGuardTests(unittest.TestCase):
    def _base_procedure(self):
        procedure = new_procedure()
        procedure.encodeEditArea = lambda area: np.ones((1, 1), dtype=bool)
        procedure.loadedParms = lambda: []
        procedure.makeTimeRange = FakeTimeRange
        procedure.statusBarMsg = lambda *args: None
        return procedure

    def test_empty_element_selection_returns_before_derived_products(self):
        procedure = self._base_procedure()
        procedure.makeTd = lambda *args: self.fail("derived products must not run")

        procedure.execute(
            FakeTimeRange(0, 3600),
            "Entire Domain",
            {
                "Beginning": "Selected TimeRange",
                "Ending": "Selected TimeRange",
                "Element": [],
                "Primary Model for Diurnal:": "NBM",
                "Secondary Model for Diurnal:": "GFS",
            },
        )

    def test_none_and_noniterable_element_values_return_cleanly(self):
        for element in (None, 7):
            with self.subTest(element=element):
                procedure = self._base_procedure()
                procedure.loadedParms = lambda: self.fail(
                    "invalid element input must stop before loading"
                )
                procedure.execute(
                    FakeTimeRange(0, 3600),
                    "Entire Domain",
                    {
                        "Beginning": "Selected TimeRange",
                        "Ending": "Selected TimeRange",
                        "Element": element,
                        "Primary Model for Diurnal:": "NBM",
                        "Secondary Model for Diurnal:": "GFS",
                    },
                )

    def test_reversed_execution_range_returns_before_loading(self):
        procedure = self._base_procedure()
        messages = []
        procedure.statusBarMsg = lambda message, level: messages.append(message)
        procedure.makeExecuteTR = lambda *args: FakeTimeRange(3600, 0)
        procedure.loadedParms = lambda: self.fail(
            "invalid execution range must stop before loading"
        )

        procedure.execute(
            FakeTimeRange(0, 3600),
            "Entire Domain",
            {
                "Beginning": "Selected TimeRange",
                "Ending": "Selected TimeRange",
                "Element": ["T"],
                "Primary Model for Diurnal:": "NBM",
                "Secondary Model for Diurnal:": "GFS",
            },
        )

        self.assertTrue(any("Ending time" in message for message in messages))

    def test_programmatic_smooth_areas_accept_none_and_a_single_name(self):
        for raw_areas, weights, expected in (
            (None, {}, []),
            ("ridge", {"ridge": 1}, ["ridge"]),
        ):
            with self.subTest(raw_areas=raw_areas):
                procedure = self._base_procedure()
                procedure.loadParm = lambda *args: None
                procedure.mutableID = lambda: "Fcst"
                procedure.getParm = lambda db, name, level: FakeParm(name)
                procedure.getParmTimeConstraints = lambda parm: (0, 12 * 3600)
                procedure.makeHourlyGrids = lambda *args: {}
                procedure.execute(
                    FakeTimeRange(0, 3600),
                    "Entire Domain",
                    {
                        "Beginning": "Selected TimeRange",
                        "Ending": "Selected TimeRange",
                        "Element": ["T"],
                        "Primary Model for Diurnal:": "NBM",
                        "Secondary Model for Diurnal:": "GFS",
                        "Custom Smooth Areas": raw_areas,
                        "Custom Smooth Weights": weights,
                    },
                )
                self.assertEqual(procedure._smoothAreas, expected)

    def test_programmatic_obs_exclusion_builds_its_own_date_mapping(self):
        procedure = self._base_procedure()
        procedure._gmtime = lambda: FakeAbsTime(0)
        original = DIURNAL.DiurnalConfig.configDict
        try:
            DIURNAL.DiurnalConfig.configDict = {
                "displayDaysToExclude": True,
                "daysToExclude": 1,
            }
            procedure.execute(
                FakeTimeRange(0, 3600),
                "Entire Domain",
                {
                    "Beginning": "Selected TimeRange",
                    "Ending": "Selected TimeRange",
                    "Element": [],
                    "Primary Model for Diurnal:": "Obs",
                    "Secondary Model for Diurnal:": "GFS",
                    "Exclude Obs:": ["1970/01/01"],
                },
            )
        finally:
            DIURNAL.DiurnalConfig.configDict = original

        self.assertEqual(procedure._excludeTRList, [FakeTimeRange(0, 24 * 3600)])

    def test_unsupported_programmatic_element_returns_before_loading(self):
        procedure = self._base_procedure()
        messages = []
        procedure.statusBarMsg = lambda message, level: messages.append(message)
        procedure.loadedParms = lambda: self.fail("unsupported element must not load")

        procedure.execute(
            FakeTimeRange(0, 3600),
            "Entire Domain",
            {
                "Beginning": "Selected TimeRange",
                "Ending": "Selected TimeRange",
                "Element": ["Td"],
                "Primary Model for Diurnal:": "NBM",
                "Secondary Model for Diurnal:": "GFS",
            },
        )

        self.assertTrue(any("Unsupported" in message for message in messages))

    def test_no_hourly_guidance_returns_before_derived_products(self):
        procedure = self._base_procedure()
        messages = []
        procedure.statusBarMsg = lambda message, level: messages.append(message)
        procedure.loadParm = lambda *args: None
        procedure.mutableID = lambda: "Fcst"
        procedure.makeExecuteTR = lambda *args: FakeTimeRange(0, 3600)
        procedure.getMinMaxTRLists = lambda *args: ([], [])
        procedure.getParm = lambda db, name, level: FakeParm(name)
        procedure.getParmTimeConstraints = lambda parm: (0, 12 * 3600)
        procedure.makeHourlyGrids = lambda *args: {}
        procedure.makeTd = lambda *args: self.fail("derived products must not run")

        procedure.execute(
            FakeTimeRange(0, 3600),
            "Entire Domain",
            {
                "Beginning": "Selected TimeRange",
                "Ending": "Selected TimeRange",
                "Element": ["T"],
                "Primary Model for Diurnal:": "NBM",
                "Secondary Model for Diurnal:": "GFS",
            },
        )

        self.assertTrue(any("No T guidance" in message for message in messages))

    def test_partial_hourly_guidance_aborts_instead_of_leaving_stale_hours(self):
        procedure = self._base_procedure()
        messages = []
        procedure.statusBarMsg = lambda message, level: messages.append(message)
        procedure.loadParm = lambda *args: None
        procedure.mutableID = lambda: "Fcst"
        procedure.getParm = lambda db, name, level: FakeParm(name)
        procedure.getParmTimeConstraints = lambda parm: (0, 12 * 3600)
        procedure.makeHourlyGrids = lambda *args: {
            FakeTimeRange(3600, 2 * 3600): np.array([[10.0]])
        }
        procedure.getMinMaxGrids = lambda *args: self.fail(
            "partial hourly coverage must stop before extrema processing"
        )
        procedure.makeTd = lambda *args: self.fail(
            "partial hourly coverage must stop before derived products"
        )

        procedure.execute(
            FakeTimeRange(0, 3 * 3600),
            "Entire Domain",
            {
                "Beginning": "Selected TimeRange",
                "Ending": "Selected TimeRange",
                "Element": ["T"],
                "Primary Model for Diurnal:": "NBM",
                "Secondary Model for Diurnal:": "GFS",
            },
        )

        self.assertTrue(any("continuously cover" in message for message in messages))

    def test_primary_climatology_rejects_rh_before_writing_temperature(self):
        procedure = self._base_procedure()
        messages = []
        procedure.statusBarMsg = lambda message, level: messages.append(message)
        procedure.loadedParms = lambda: self.fail("invalid selection must stop early")
        procedure.callProcedure = lambda *args, **kwargs: self.fail(
            "climatology must not write a partial run"
        )

        procedure.execute(
            FakeTimeRange(0, 3600),
            "Entire Domain",
            {
                "Beginning": "Selected TimeRange",
                "Ending": "Selected TimeRange",
                "Element": ["T", "RH"],
                "Primary Model for Diurnal:": "F-table (Climo)",
                "Secondary Model for Diurnal:": "GFS",
            },
        )

        self.assertTrue(any("only be the primary source for T" in m for m in messages))


class ApparentTemperatureTests(unittest.TestCase):
    def test_missing_wind_skips_the_period_instead_of_crashing(self):
        procedure = new_procedure()
        time_range = FakeTimeRange(0, 3600)
        procedure.getGridInfo = lambda *args, **kwargs: [FakeGridInfo(time_range)]
        procedure.getGrids = lambda db, name, *args, **kwargs: {
            "T": np.array([[90.0]]),
            "RH": np.array([[30.0]]),
            "Wind": None,
        }[name]
        created = []
        procedure.createGrid = lambda *args, **kwargs: created.append(args)
        deleted = []
        procedure.deleteGrid = lambda *args, **kwargs: deleted.append(args)

        procedure.makeAT("Fcst", time_range)

        self.assertEqual(created, [])
        self.assertEqual(deleted[0][1], "ApparentT")

    def test_missing_apparent_temperature_input_removes_stale_outputs(self):
        procedure = new_procedure()
        hourly = FakeTimeRange(0, 3600)
        extrema = FakeTimeRange(0, 12 * 3600)
        procedure.getGridInfo = lambda db, name, *args, **kwargs: {
            "T": [FakeGridInfo(hourly)],
            "MinT": [FakeGridInfo(extrema)],
            "MaxT": [FakeGridInfo(extrema)],
        }.get(name, [])
        procedure.getGrids = lambda db, name, *args, **kwargs: {
            "T": np.array([[90.0]]),
            "RH": None,
            "Wind": (np.array([[5.0]]), np.array([[0.0]])),
            "ApparentT": np.array([[90.0]]),
        }[name]
        deleted = []
        procedure.deleteGrid = lambda db, name, level, tr: deleted.append((name, tr))
        procedure.createGrid = lambda *args, **kwargs: None

        procedure.makeAT("Fcst", extrema)
        procedure.maxminAT("Fcst", extrema)

        self.assertEqual(
            deleted,
            [("ApparentT", hourly), ("MinApT", extrema), ("MaxApT", extrema)],
        )

    def test_dry_heat_index_uses_the_actual_temperature(self):
        procedure = new_procedure()
        time_range = FakeTimeRange(0, 3600)
        temperature = np.array([[95.0]])
        humidity = np.array([[5.0]])
        procedure.getGridInfo = lambda *args, **kwargs: [FakeGridInfo(time_range)]
        procedure.getGrids = lambda db, name, *args, **kwargs: {
            "T": temperature,
            "RH": humidity,
            "Wind": (np.array([[0.0]]), np.array([[0.0]])),
        }[name]
        created = []
        procedure.createGrid = lambda *args, **kwargs: created.append(args)

        procedure.makeAT("Fcst", time_range)

        base = (
            -42.379
            + 2.04901523 * 95.0
            + 10.14333127 * 5.0
            - 0.22475541 * 95.0 * 5.0
            - 0.00683783 * 95.0 ** 2
            - 0.05481717 * 5.0 ** 2
            + 0.00122874 * 95.0 ** 2 * 5.0
            + 0.00085282 * 95.0 * 5.0 ** 2
            - 0.00000199 * 95.0 ** 2 * 5.0 ** 2
        )
        expected = base - ((13.0 - 5.0) / 4.0)
        np.testing.assert_allclose(created[0][3], [[expected]])

    def test_min_and_max_apparent_temperature_bootstrap_from_temperature_periods(self):
        procedure = new_procedure()
        min_range = FakeTimeRange(0, 12 * 3600)
        max_range = FakeTimeRange(12 * 3600, 24 * 3600)

        def grid_info(db, name, *args, **kwargs):
            return {
                "MinT": [FakeGridInfo(min_range)],
                "MaxT": [FakeGridInfo(max_range)],
                "MinApT": [],
                "MaxApT": [],
            }.get(name, [])

        procedure.getGridInfo = grid_info
        procedure.getGrids = lambda db, name, level, tr, mode=None, **kwargs: {
            "Min": np.array([[10.0]]),
            "Max": np.array([[90.0]]),
        }[mode]
        created = []
        procedure.createGrid = lambda *args, **kwargs: created.append(args)

        procedure.maxminAT("Fcst", FakeTimeRange(0, 24 * 3600))

        self.assertEqual([call[1] for call in created], ["MinApT", "MaxApT"])

    def test_legacy_heat_index_helper_uses_numpy_array_operations(self):
        procedure = new_procedure()
        time_range = FakeTimeRange(0, 3600)
        procedure.getGridInfo = lambda *args, **kwargs: [FakeGridInfo(time_range)]
        procedure.getGrids = lambda db, name, *args, **kwargs: {
            "T": np.array([[95.0]]),
            "Td": np.array([[55.0]]),
        }[name]
        created = []
        procedure.createGrid = lambda *args, **kwargs: created.append(args)

        procedure.makeHI("Fcst", time_range)

        self.assertEqual(created[0][1], "HeatIndex")
        self.assertTrue(np.isfinite(created[0][3]).all())

    def test_legacy_wind_chill_helper_uses_numpy_array_operations(self):
        procedure = new_procedure()
        time_range = FakeTimeRange(0, 3600)
        procedure.getGridInfo = lambda *args, **kwargs: [FakeGridInfo(time_range)]
        procedure.getGrids = lambda db, name, *args, **kwargs: {
            "T": np.array([[20.0]]),
            "Wind": (np.array([[10.0]]), np.array([[0.0]])),
        }[name]
        created = []
        procedure.createGrid = lambda *args, **kwargs: created.append(args)

        procedure.makeWC("Fcst", time_range)

        self.assertEqual(created[0][1], "WindChill")
        self.assertTrue(np.isfinite(created[0][3]).all())


class SolarCalculationTests(unittest.TestCase):
    def test_negative_solar_flux_is_initialized_to_zero(self):
        procedure = new_procedure()
        procedure.empty = lambda dtype=float: np.full((1, 1), -123.0, dtype=dtype)
        procedure.eqnOfTime = lambda *args: (
            0.0,
            1.0,
            0.0,
            np.array([[-1.0]]),
            np.array([[0.0]]),
            1,
        )

        flux, *_ = procedure.solarFlux(np.array([[90.0]]), 2026, 12, 21)

        np.testing.assert_array_equal(flux, [[0.0]])

    def test_polar_night_sunrise_and_sunset_are_finite(self):
        procedure = new_procedure()
        procedure.empty = lambda dtype=float: np.full((1, 1), 123.0, dtype=dtype)

        with np.errstate(all="ignore"):
            sunrise, sunset = procedure.sunrSuns(
                2026,
                12,
                21,
                np.array([[0.0]]),
                np.array([[90.0]]),
            )

        self.assertTrue(np.isfinite(sunrise).all())
        self.assertTrue(np.isfinite(sunset).all())
        np.testing.assert_array_equal(sunset - sunrise, [[0.0]])

    def test_sunrise_day_number_preserves_integer_calendar_formula(self):
        procedure = new_procedure()
        procedure.empty = lambda dtype=float: np.empty((1, 1), dtype=dtype)

        sunrise, sunset = procedure.sunrSuns(
            2026,
            9,
            2,
            np.array([[-116.2]]),
            np.array([[43.6]]),
        )

        np.testing.assert_allclose(sunrise, [[13.16814560894]], atol=1e-9)
        np.testing.assert_allclose(sunset, [[26.31118058239]], atol=1e-9)

    def test_zenith_roundoff_is_clipped_with_a_full_grid_mask(self):
        procedure = new_procedure()
        shape = (2, 2)
        procedure.empty = lambda dtype=float: np.empty(shape, dtype=dtype)
        procedure._Procedure__pointDiagnostics = False
        matching_angle = -88.99287992879928

        zenith, cosine = procedure.getZenAng(
            np.deg2rad(matching_angle),
            0.0,
            np.array([[matching_angle, 0.0], [10.0, 20.0]]),
            np.zeros(shape),
            12,
            np.full(shape, 1013.0),
            np.full(shape, 25.0),
        )

        self.assertEqual(zenith.shape, shape)
        self.assertTrue(np.isfinite(zenith).all())
        self.assertTrue(np.logical_and(cosine >= -1.0, cosine <= 1.0).all())


class WBGTTests(unittest.TestCase):
    def test_diffuse_radiation_scales_sky_and_reflected_components(self):
        procedure = new_procedure()
        stefan_boltzmann = 5.67e-8

        result = procedure._wbgtRadiationFactor(
            np.array([[0.0]]),
            np.array([[1.0]]),
            np.array([[0.5]]),
            np.array([[0.2]]),
            stefan_boltzmann,
        )

        np.testing.assert_allclose(result, [[1.2 / stefan_boltzmann]])

    def test_subfreezing_night_does_not_mirror_globe_temperature_above_freezing(self):
        procedure = new_procedure()
        procedure._Procedure__pointDiagnostics = False
        procedure.empty = lambda dtype=float: np.empty((1, 1), dtype=dtype)
        procedure.newGrid = lambda value: np.full((1, 1), value, dtype=float)
        procedure.getZenAng = lambda *args: (
            np.array([[90.0]]),
            np.array([[1.0]]),
        )
        temperature_f = np.array([[-20.0]])
        temperature_c = (temperature_f - 32.0) * 5.0 / 9.0
        dewpoint_c = np.array(temperature_c, copy=True)

        result = procedure.wbgtMainCalc(
            np.array([[0.0]]),
            np.array([[0.0]]),
            1,
            0,
            np.array([[0.0]]),
            np.array([[45.0]]),
            dewpoint_c,
            temperature_c,
            np.array([[1013.0]]),
            np.array([[2.5 * 3600.0]]),
            np.array([[2.5]]),
            temperature_f,
            FakeTimeRange(0, 3600),
            ["Low", "Elevated", "Moderate", "High", "Extreme"],
            np.array([[0.2]]),
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
        )

        self.assertTrue(np.isfinite(result).all())
        np.testing.assert_allclose(result, temperature_f, atol=1e-9)

    def test_missing_pressure_database_returns_none(self):
        procedure = new_procedure()
        procedure.availableDatabases = lambda: []

        result = procedure._findLongestModelDB(
            "D2D_GFS_", "p", "SFC", FakeTimeRange(0, 3600)
        )

        self.assertIsNone(result)

    def test_failed_wbgt_preflight_does_not_delete_existing_grids(self):
        procedure = new_procedure()
        procedure._prepareWBGTInputs = lambda: None
        deleted = []
        procedure.deleteCmd = lambda *args, **kwargs: deleted.append(args)

        result = procedure._updateWBGT()

        self.assertFalse(result)
        self.assertEqual(deleted, [])

    def test_failed_wbgt_calculation_validation_does_not_delete_existing_grids(self):
        procedure = new_procedure()
        procedure._prepareWBGTInputs = lambda: {}
        procedure._validateWBGTCalculations = lambda inputs: False
        procedure.statusBarMsg = lambda *args: None
        deleted = []
        procedure.deleteCmd = lambda *args, **kwargs: deleted.append(args)

        result = procedure._updateWBGT()

        self.assertFalse(result)
        self.assertEqual(deleted, [])

    def test_wbgt_threat_keys_must_have_the_expected_order(self):
        procedure = new_procedure()

        self.assertTrue(
            procedure._validWBGTThreatKeys(
                ["Low", "Elevated", "Moderate", "High", "Extreme"]
            )
        )
        self.assertTrue(
            procedure._validWBGTThreatKeys(
                ["Extreme", "High", "Moderate", "Elevated", "Low"]
            )
        )
        self.assertFalse(procedure._validWBGTThreatKeys(["Low", "High"]))

    def test_wbgt_output_grid_types_are_checked_before_mutation(self):
        procedure = new_procedure()

        class GridInfo:
            def __init__(self, grid_type):
                self._grid_type = grid_type

            def getGridType(self):
                return self._grid_type

        class Parm:
            def __init__(self, grid_type):
                self._info = GridInfo(grid_type)

            def getGridInfo(self):
                return self._info

        self.assertTrue(procedure._validWBGTParmType(Parm("SCALAR"), "SCALAR"))
        self.assertFalse(procedure._validWBGTParmType(Parm("DISCRETE"), "SCALAR"))

    def test_wbgt_write_failure_propagates_abort_after_warning(self):
        procedure = new_procedure()
        procedure._prepareWBGTInputs = lambda: {}
        procedure._validateWBGTCalculations = lambda inputs: True
        procedure._writeWBGTGrids = lambda inputs: (_ for _ in ()).throw(
            RuntimeError("write failed")
        )
        messages = []
        procedure.statusBarMsg = lambda message, level: messages.append((message, level))

        with self.assertRaisesRegex(RuntimeError, "write failed"):
            procedure._updateWBGT()

        self.assertEqual(messages[-1][1], "U")
        self.assertIn("Reload WBGT", messages[-1][0])

    def test_wbgt_outputs_are_recomputed_after_memory_bounded_validation(self):
        procedure = new_procedure()
        periods = [FakeTimeRange(0, 3600), FakeTimeRange(3600, 7200)]
        all_time = FakeTimeRange(-1000, 10000)
        inputs = {
            "timeRanges": periods,
            "threatKeys": ["Low", "Elevated", "Moderate", "High", "Extreme"],
            "allTR": all_time,
            "fullTR": FakeTimeRange(0, 7200),
        }
        procedure._prepareWBGTInputs = lambda: inputs
        calculations = []
        procedure._calculateWBGTGrid = lambda data, tr, cache: (
            calculations.append(tr) or np.array([[75.0]])
        )
        procedure._makeWBGTRiskGrid = lambda grid: np.array([[1.0]])
        deleted = []
        procedure.deleteCmd = lambda elements, tr: deleted.append((elements, tr))
        procedure.createGrid = lambda *args, **kwargs: None
        procedure._createMaxGrids = lambda *args: None
        procedure.statusBarMsg = lambda *args: None

        result = procedure._updateWBGT()

        self.assertTrue(result)
        self.assertEqual(calculations, periods + periods)
        self.assertEqual(deleted[0][1], all_time)

    def test_max_wbgt_risk_remaps_reordered_keys_without_double_mapping(self):
        procedure = new_procedure()
        hourly = FakeTimeRange(0, 3600)
        maximum = FakeTimeRange(0, 12 * 3600)
        source_keys = ["Extreme", "High", "Moderate", "Elevated", "Low"]
        source_grid = np.array([[4, 3, 2, 1, 0]])
        procedure.createFromScratchCmd = lambda *args: None
        procedure.fragmentCmd = lambda *args: None
        procedure.getGrids = lambda *args, **kwargs: {
            hourly: (source_grid, source_keys)
        }
        procedure.getGridInfo = lambda *args, **kwargs: [FakeGridInfo(maximum)]
        procedure.newGrid = lambda value: np.full((1, 5), value, dtype=np.int8)
        created = []
        procedure.createGrid = lambda *args, **kwargs: created.append(args)
        procedure.deleteGrid = lambda *args, **kwargs: None

        procedure._createMaxGrids(
            maximum,
            [hourly],
            ["Low", "Elevated", "Moderate", "High", "Extreme"],
        )

        np.testing.assert_array_equal(created[0][3][0], [[0, 1, 2, 3, 4]])
        np.testing.assert_array_equal(source_grid, [[4, 3, 2, 1, 0]])

    def test_risk_grid_is_zero_outside_configured_regions(self):
        procedure = new_procedure()
        procedure.empty = lambda dtype=float: np.empty((2, 2), dtype=dtype)
        procedure.newGrid = lambda value: np.full((2, 2), value, dtype=float)
        procedure.regionsDict = {
            "Region1": {
                "low": 70.0,
                "elevated": 75.0,
                "moderate": 80.0,
                "high": 85.0,
            }
        }
        procedure.validEAs = [("Region1", "inside")]
        procedure.encodeEditArea = lambda name: np.array(
            [[True, False], [False, False]], dtype=bool
        )

        result = procedure._makeWBGTRiskGrid(
            np.array([[90.0, 90.0], [90.0, 90.0]])
        )

        np.testing.assert_array_equal(result, [[4.0, 0.0], [0.0, 0.0]])


if __name__ == "__main__":
    unittest.main()
