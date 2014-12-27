#!/usr/bin/python

import numpy as np
import scipy.interpolate as sip
import matplotlib.pyplot as plt
import os
nstep = 9000

if 'info.txt' in [f for f in os.listdir('.') if os.path.isfile(f)]:
    prefix = os.getcwd() + '/'
else:
    prefix = os.getcwd() + '/dat/'

n = 512
ieout = 100
def readindata():
    global n
    global ieout
    finfo = open(prefix+'info.txt')
    info = finfo.readlines()
    n = int(info[1].split()[1])
    ieout = int(info[4].split()[5])
readindata()

plt.close('all')

def plot_vort():
    '''numstring is the string form of the number of steps, which
    should be controlled by the variable step.'''
    ntime = 1000
    startpoint = 100

    step = '{0:0>#6}'.format(nstep)
    filename = prefix + 'vorstat.d'
    filename1 = prefix + "vortfield" + step

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
    f.colorbar(ax[1])
    #ax[1].set_title("Vorticity distribution")
    plt.show()

def plot_energy():
    '''step is the string form of the number of the step you want
    to plot, which is recommended to be controlled by the variable
    DataHandle.step.'''
    fspec = open("spectrum.d")
    tmp = []
    time = []
    energy = []
    spectrum = []
    wavenum = []
    counter = 0
    istep = nstep / ieout

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

#plot_vort()
#plot_energy()

