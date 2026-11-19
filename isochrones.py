#!/usr/bin/env python

import numpy as np
from astropy.table import Table,hstack,vstack
#from dlnpyutils import utils as dln,plotting as pl

class Isochrone:

    def __init__(self,isolist):
        metals = np.zeros(len(isolist),float)
        ages = np.zeros(len(isolist),float)
        data = []
        for i in range(len(isolist)):
            iso1 = isolist[i]
            metals[i] = iso1['mh']
            ages[i] = iso1['agegyr']
            data.append(iso1['data'])
        self._metals = metals
        self._ages = ages
        self._data = data

    def __call__(self,age,metal):
        ind, = np.where((np.abs(self._ages-age)<0.01) & (np.abs(self._metals-metal)<0.01))
        if len(ind)==0:
            print('Age,Metal=({:.2f},{:.2f}) not found'.format(age,metal))
        return self._data[ind[0]]

    def __repr__(self):
        """ Print out the string representation of the Isochrone object."""
        s = "Isochrone \n"
        s += "Ages: "+','.join(np.char.array(self.ages).astype(str))+"\n"
        s += "Metals: "+','.join(np.char.array(self.metals).astype(str))+"\n"
        return s

    @property
    def metals(self):
        return np.unique(self._metals)

    @property
    def ages(self):
        return np.unique(self._ages)
