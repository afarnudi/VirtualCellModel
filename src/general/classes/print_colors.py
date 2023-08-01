class TerminalColors:
    TRESET = "\033[0m"
    TBLINK = "\033[5m"
    TBOLD = "\033[1m"
    TGRAY = "\033[38;5;249m"
    TRED = "\033[38;5;196m"
    TGREEN = "\033[38;5;118m"
    TYELLOW = "\033[38;5;226m"
    TBLUE = "\033[38;5;33m"
    TPINK = "\033[38;5;200m"
    TCYAN = "\033[38;5;123m"
    TWHITE = "\033[38;5;255m"
    TORANGE = "\033[38;5;214m"
    TGB = "\033[38;5;50m"  # Greenish Blue
    TDARKGB = "\033[38;5;6m"  # Dark Greenish Blue
    TPURPLE = "\033[38;5;93m"
    TLGREEN = "\033[38;5;155m"  # Light Green
    TBORANGE = "\033[38;5;166m"  # Blood Orange

    TSUCCESS = TGREEN
    TON = TGREEN
    TFAILED = TRED
    TOFF = TRED
    TFILE = TDARKGB
    TWARN = TORANGE
    TWWARN = TBLINK + TRED

    TMEM = TBLUE
    TACT = TCYAN
    TECM = TDARKGB
    TCHR = TPINK
    TOMM = TGB

    TOCL = TPURPLE
    TCUD = TGREEN
    TCPU = TBORANGE

    tc_color = {}
    tc_color["OpenCL"] = TOCL
    tc_color["CUDA"] = TCUD
    tc_color["CPU"] = TCPU
