from datetime import timedelta

from source.general.report_time import parse_time_delta
from source.general.report_time import get_string_for_time_unit
from source.general.report_time import print_runtime


def test_print_runtime_0_day_0_hour_0_min_0_sec():
    output = print_runtime(0, "title")
    answer = "title:\n0\tseconds\n\n"
    assert output == answer


def test_print_runtime_0_day_0_hour_0_min_1_sec():
    output = print_runtime(1, "title")
    answer = "title:\n1\tsecond\n\n"
    assert output == answer


def test_print_runtime_0_day_0_hour_0_min_2_sec():
    output = print_runtime(2, "title")
    answer = "title:\n2\tseconds\n\n"
    assert output == answer


def test_print_runtime_0_day_0_hour_1_min_0_sec():
    output = print_runtime(60, "title")
    answer = "title:\n1\tminute\n\n"
    assert output == answer


def test_print_runtime_0_day_0_hour_1_min_2_sec():
    output = print_runtime(62, "title")
    answer = "title:\n1\tminute\n2\tseconds\n\n"
    assert output == answer


def test_get_string_for_time_unit_0():
    output = get_string_for_time_unit(0, "unit")
    assert output == "0\tunits\n"


def test_get_string_for_time_unit_1():
    output = get_string_for_time_unit(1, "unit")
    assert output == "1\tunit\n"


def test_get_string_for_time_unit_2():
    output = get_string_for_time_unit(2, "unit")
    assert output == "2\tunits\n"


def test_0_day_0_hour_0_min_0_5_sec():
    delta = timedelta(seconds=0.5)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 0
    assert delta_t["hour"] == 0
    assert delta_t["minute"] == 0
    assert delta_t["second"] == 0


def test_0_day_0_hour_0_min_0_sec():
    delta = timedelta(seconds=0)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 0
    assert delta_t["hour"] == 0
    assert delta_t["minute"] == 0
    assert delta_t["second"] == 0


def test_0_day_0_hour_0_min_10_sec():
    delta = timedelta(seconds=10)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 0
    assert delta_t["hour"] == 0
    assert delta_t["minute"] == 0
    assert delta_t["second"] == 10


def test_0_day_0_hour_1_min_0_sec():
    delta = timedelta(seconds=60)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 0
    assert delta_t["hour"] == 0
    assert delta_t["minute"] == 1
    assert delta_t["second"] == 0


def test_0_day_2_hour_0_min_0_sec():
    delta = timedelta(seconds=2 * 3600)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 0
    assert delta_t["hour"] == 2
    assert delta_t["minute"] == 0
    assert delta_t["second"] == 0


def test_3_day_0_hour_0_min_0_sec():
    delta = timedelta(seconds=3 * 3600 * 24)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 3
    assert delta_t["hour"] == 0
    assert delta_t["minute"] == 0
    assert delta_t["second"] == 0


def test_0_day_0_hour_1_min_2_sec():
    delta = timedelta(seconds=62)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 0
    assert delta_t["hour"] == 0
    assert delta_t["minute"] == 1
    assert delta_t["second"] == 2


def test_0_day_3_hour_0_min_2_sec():
    delta = timedelta(seconds=3 * 3600 + 2)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 0
    assert delta_t["hour"] == 3
    assert delta_t["minute"] == 0
    assert delta_t["second"] == 2


def test_4_day_0_hour_0_min_2_sec():
    delta = timedelta(seconds=4 * 3600 * 24 + 2)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 4
    assert delta_t["hour"] == 0
    assert delta_t["minute"] == 0
    assert delta_t["second"] == 2


def test_0_day_3_hour_1_min_2_sec():
    delta = timedelta(seconds=3 * 3600 + 62)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 0
    assert delta_t["hour"] == 3
    assert delta_t["minute"] == 1
    assert delta_t["second"] == 2


def test_1_day_0_hour_1_min_2_sec():
    delta = timedelta(seconds=24 * 3600 + 62)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 1
    assert delta_t["hour"] == 0
    assert delta_t["minute"] == 1
    assert delta_t["second"] == 2


def test_1_day_3_hour_1_min_2_sec():
    delta = timedelta(seconds=24 * 3600 + 3 * 3600 + 62)
    delta_t = parse_time_delta(delta)
    assert delta_t["day"] == 1
    assert delta_t["hour"] == 3
    assert delta_t["minute"] == 1
    assert delta_t["second"] == 2
