#!/usr/bin/env python3

import csv
import sys

with open('data.csv', 'w', newline = '') as csv_file:
    crystaldata_writer = csv.writer(csv_file, delimiter = ",", quotechar = '"', \
            quoting = csv.QUOTE_MINIMAL)
    crystaldata_writer.writerow(['CIFName', "Octahedral Number", 'InPlaneDistortion', 'Out_Plane_Distortion', 'Tilt_Angle_Distortion', 'Bond_Distortion'])

