from datetime import timedelta

def get_string_for_time_unit(delta_time, time_unit):
    delta_time_str = f"{delta_time}\t {time_unit}"
    if delta_time > 1:
        delta_time_str + "s"
    delta_time_str += "\n"
    return delta_time_str

def parse_time_delta(delta_time):
    seconds_per_hour = 3600
    seconds_per_min = 60
    days = delta_time.days
    total_seconds = delta_time.seconds
    hours, remainder = divmod(total_seconds, seconds_per_hour)
    minutes, seconds = divmod(remainder, seconds_per_min)
    return days, hours, minutes, seconds

def print_runtime(delta_clock, title):
    

    delta = timedelta(seconds=delta_clock)
    days, hours, minutes, seconds = parse_time_delta(delta)

    print(title)
    for delta_time, time_unit in zip([days, hours, minutes], ["day", "hour", "minute"]):
        if delta_time != 0:
            output_str = get_string_for_time_unit(delta_time, time_unit)
            print(output_str)
    seconds_str = get_string_for_time_unit(seconds, "second")
    print(seconds_str)