#!/usr/bin/env python3

import csv
import sys

with open('data.csv', 'w', newline = '') as csv_file:
    crystaldata_writer = csv.writer(csv_file, delimiter = ",", quotechar = '"', \
            quoting = csv.QUOTE_MINIMAL)
    crystaldata_writer.writerow(['CIFName', "Octahedral Number", 'AvgIn', 'AvgOut', 'AvgTilt', 'MaxIn', 'MaxOut', 'MaxTilt', 'MinIn', 'MinOut', 'MinTilt', 'Bond_Distortion'])

