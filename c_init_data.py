# -*- coding: utf-8 -*-
"""
Created on 30.03.23

@author: shv
@description:
    Meant for reading any needed type of files
"""

class initData(object):
    """This is a doc string.
    Class ist meant for FE data of 1D-SFE solver"""
    def __init__(self):
        self.E = 0.0
        self.nu = 0.0
        self.ElementType = 0
        self.degree = 1
        self.numElements = 0
        self.x = []
        self.A = []
        self.forceX = []
        
    def read(self, filename):
        file = open(filename, 'r')
        lines = file.readlines()
        for line in lines:
            if len(line.rsplit()) > 0:
                indx = line.rsplit()[0]
                if (indx == '**'):
                    continue
                if (indx == '*0001'):
                    self.E = float(line.rsplit()[1])
                if (indx == '*0002'):
                    self.nu = float(line.rsplit()[1])
                if (indx == '*0003'):
                    self.ElementType = int(line.rsplit()[1])
                if (indx == '*0004'):
                    self.degree = int(line.rsplit()[1])
                if (indx == '*1000'):
                    self.x.append(float(line.rsplit()[1]))
                    height = float(line.rsplit()[2])
                    width = float(line.rsplit()[3])
                    self.A.append(height*width)
                    self.forceX.append(float(line.rsplit()[4]))
        self.numElements = len(self.x) - 1
        file.close()
                