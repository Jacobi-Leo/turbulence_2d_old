#!/usr/bin/python

fspec = open("spectrum.d")
tmp = []
time = []
energy = []
spectrum = []
wavenum = []
counter = 0

while True:
    dataline = fspec.readline().split()

    if len(dataline) == 0:
        break

    if dataline[0] == "#k":
        time.append(float(dataline[6]))
        energy.append(float(dataline[8]))
        if counter > 0:
            spectrum.append(tmp)
            tmp = []
        counter += 1
    else:
        if counter == 1:
            wavenum.append(int(dataline[0]))
        tmp.append(float(dataline[1]))

fspec.close()


