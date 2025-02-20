def print_info(str, banner="="):
    """
    Print information with a banner to help
    with visibility during debug printing
    """
    print(banner * 80)
    print(f"{str}".center(80))
    print(banner * 80)