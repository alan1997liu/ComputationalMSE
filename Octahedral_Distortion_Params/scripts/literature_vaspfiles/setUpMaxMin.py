#!/usr/bin/env python3

import csv
import sys

with open('data.csv', 'w', newline = '') as csv_file:
    crystaldata_writer = csv.writer(csv_file, delimiter = ",", quotechar = '"', \
            quoting = csv.QUOTE_MINIMAL)
    crystaldata_writer.writerow(['CIFName', 'MaxInPlane', 'MinInPlane', 'MaxOutPlane', 'MinOutPlane', 'MaxTilt', 'MinTilt', 'MaxOctAngleVariance', 'MinOctAngleVariance'])


