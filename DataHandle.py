#!/usr/bin/python

import numpy as np
import scipy.interpolate as sip
import matplotlib.pyplot as plt

plt.close('all')

def plot_vort(numstring):
    ntime = 1000
    startpoint = 100

    filename = prefix + 'vorstat.d'
    filename1 = prefix + "vortfield" + numstring

    raw_data = np.loadtxt(filename)
    time = raw_data.transpose()[0]
    num_of_vort = raw_data.transpose()[1]

    vortfield_raw = np.loadtxt(filename1)
    vortfield = vortfield_raw[0:n, 0:n]

    f, ax = plt.subplots(1, 2)
    ax[0].loglog(time[startpoint:ntime], num_of_vort[startpoint:ntime])
    #ax[0].loglog(time, num_of_vort)
    ax[0].set_title("Number of vortex vs time")
    ax[1].matshow(vortfield, cmap=plt.cm.Spectral)
    #ax[1].set_title("Vorticity distribution")
    plt.show()

#L = scatter_plot(zip(time[startpoint:ntime],num_of_vort[startpoint:ntime]))
#L.show(scale="loglog")

def plot_energy(step):
    fspec = open("spectrum.d")
    tmp = []
    time = []
    energy = []
    spectrum = []
    wavenum = []
    counter = 0
    istep = int(step) / ieout

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

    f, ax = plt.subplots(1, 2)
    ax[0].loglog(time,energy)
    ax[0].set_title("Energy decay with time")
    ax[1].plot(wavenum, spectrum[istep])
    ax[1].set_title("spectrum")
    plt.show()

prefix = '/home/liu/dat_inverse_512/'
n = 512
step = '009000'
ieout = 100

plot_vort(step)
plot_energy(step)

