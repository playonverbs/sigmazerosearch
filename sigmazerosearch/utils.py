def file_ok(filename: str, mode: str = "read") -> bool:
    """
    Test file for read/write availability and convert Exceptions to boolean
    Defaults to checking for read availability
    """
    try:
        fp = open(filename) if mode == "read" else open(filename, "w")
    except FileNotFoundError:
        return False
    else:
        fp.close()

    return True
