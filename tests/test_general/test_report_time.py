from datetime import timedelta

from source.general.report_time import parse_time_delta


def test_0_day_0_hour_0_min_0_5_sec():
    delta = timedelta(seconds=0.5)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 0
    assert hours == 0
    assert minutes == 0
    assert seconds == 0


def test_0_day_0_hour_0_min_0_sec():
    delta = timedelta(seconds=0)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 0
    assert hours == 0
    assert minutes == 0
    assert seconds == 0


def test_0_day_0_hour_0_min_10_sec():
    delta = timedelta(seconds=10)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 0
    assert hours == 0
    assert minutes == 0
    assert seconds == 10


def test_0_day_0_hour_1_min_0_sec():
    delta = timedelta(seconds=60)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 0
    assert hours == 0
    assert minutes == 1
    assert seconds == 0


def test_0_day_2_hour_0_min_0_sec():
    delta = timedelta(seconds=2 * 3600)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 0
    assert hours == 2
    assert minutes == 0
    assert seconds == 0


def test_3_day_0_hour_0_min_0_sec():
    delta = timedelta(seconds=3 * 3600 * 24)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 3
    assert hours == 0
    assert minutes == 0
    assert seconds == 0


def test_0_day_0_hour_1_min_2_sec():
    delta = timedelta(seconds=62)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 0
    assert hours == 0
    assert minutes == 1
    assert seconds == 2


def test_0_day_3_hour_0_min_2_sec():
    delta = timedelta(seconds=3 * 3600 + 2)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 0
    assert hours == 3
    assert minutes == 0
    assert seconds == 2


def test_4_day_0_hour_0_min_2_sec():
    delta = timedelta(seconds=4 * 3600 * 24 + 2)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 4
    assert hours == 0
    assert minutes == 0
    assert seconds == 2


def test_0_day_3_hour_1_min_2_sec():
    delta = timedelta(seconds=3 * 3600 + 62)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 0
    assert hours == 3
    assert minutes == 1
    assert seconds == 2


def test_1_day_0_hour_1_min_2_sec():
    delta = timedelta(seconds=24 * 3600 + 62)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 1
    assert hours == 0
    assert minutes == 1
    assert seconds == 2


def test_1_day_3_hour_1_min_2_sec():
    delta = timedelta(seconds=24 * 3600 + 3 * 3600 + 62)
    days, hours, minutes, seconds = parse_time_delta(delta)
    assert days == 1
    assert hours == 3
    assert minutes == 1
    assert seconds == 2
