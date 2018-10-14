#!/usr/bin/env python
#-*- coding: utf-8 -*-

import pylab as py
import numpy as np

#--get data
path = './'
file1 = path + 'G4dataCheck.txt'
#file1 = path + 'G4data-src-1.txt'

data1 = np.loadtxt(file1)
nTotEvents = 100000
print "Total energy in plastic: ", (data1[:,:28]).sum()
print "Total energy in BGO: ", (data1[:,28:]).sum()

#--process data
trigP = (data1[:,:28] > 30).sum(1)
trigB = (data1[:,28:] > 50).sum(1)
ud = (data1 > 300).sum(1)
threshP = data1[:,:28] > 2
threshB = data1[:,28:] > 50
thresh = np.concatenate((threshP, threshB), axis=1)
dataProc = data1*thresh
dataProc = np.compress(((trigP+trigB > 0) & (ud == 0)), dataProc, axis=0)


events = (dataProc != 0).sum(1)
eventsP = (dataProc[:,:28] != 0).sum(1)
eventsB = (dataProc[:,28:] != 0).sum(1)
P_etot = (dataProc[:,:28]).sum(1)
B_etot = (dataProc[:,28:]).sum(1)

surface_per_event = (np.arccos(-1)*20**2)/nTotEvents

singles = np.compress((events == 1), dataProc, axis=0)
doubles = np.compress((events == 2), dataProc, axis=0)
doublesPB = np.compress((eventsP == 1) & (eventsB == 1), dataProc, axis=0)
doublesPcBp = np.compress((eventsP == 1) & (eventsB == 1) & (P_etot < B_etot), dataProc, axis=0)

singles_etot = np.sum(singles,axis=1)
doubles_etot = np.sum(doubles,axis=1)


#--visualization
py.figure()

emax = 1000
binning = 10

#py.subplot(2,1,1)
##py.loglog()
#py.semilogy()
#py.grid()
#bins=np.arange(0,emax+1,binning)
#singles_spectrum, bins, p = py.hist(singles_etot, bins=bins, histtype='step')
#py.ylabel('Singles')

#py.subplot(2,1,2)
##py.loglog()
#py.semilogy()
#py.grid()
#bins=np.arange(0,emax+1,binning)
#doubles_spectrum, bins, p = py.hist(doubles_etot, bins=bins, histtype='step')
#py.ylabel('Doubles')
#py.xlabel('keV')

#--prints
print
print "singles:"
#print "photopeak =", np.amax(singles_spectrum), "ph @ ", bins[np.argmax(singles_spectrum)], "keV"
#print "photopeak % = ", float(np.amax(singles_spectrum))/np.sum(singles_spectrum)
print "effective area = ", singles.shape[0]*surface_per_event, "cm²"
print
print "doubles:"
#print "photopeak =", np.amax(doubles_spectrum), "ph @ ", bins[np.argmax(doubles_spectrum)], "keV"
#print "photopeak % = ", float(np.amax(doubles_spectrum))/np.sum(doubles_spectrum)
print "effective area all = ", doubles.shape[0]*surface_per_event, "cm²"
print "effective area P<->B = ", doublesPB.shape[0]*surface_per_event, "cm²"
print "effective area P->B = ", doublesPcBp.shape[0]*surface_per_event, "cm²"

#py.show()


#l = []
#for i in range(0,10,2):
#    print i,
#    l.append(i)
#print l
#l = l + l + [13]
#print l, type(l)
#l = np.array(l)
#print l, type(l)
#
#for i in np.arange(0,10.5,0.5):
#    print "coucou", i, "machin %s %d ca va"%(i, i+1)
#
#if i == 10:
#    print "genial"
#elif i > 3:
#    print "azerty"
#else:
#    print "fin"
#
#vecteur = np.array([])
#vect = np.zeros(1000)
#matrice = np.zeros((2,4))
#print vect
# Python
# list => l = []  = string, nombre, vecteur, matrice
# tuple => t = ()
# dictionnaire => d = {}
# Numpy
# ndarray => fonction numpy.machin ou np.array(tuple, list)

