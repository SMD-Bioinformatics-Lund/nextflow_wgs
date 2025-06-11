#!/usr/bin/env python3

import argparse


description = """
Parse ROH output and UPD output into a single bed file
This bed will then be loaded into Gens as a sample annotation track
"""

def main():
    pass



def parse_roh():
    pass


def parse_upd():
    pass


def write_bed():
    pass


def parse_arguments():
    parser = argparse.ArgumentParser();

    # FIXME: Understand the difference roh / upd

    parser.add_argument("--roh", required=True)
    parser.add_argument("--upd", required=True)

    parser.add_argument("--roh_quality_threshold", default=85)

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
