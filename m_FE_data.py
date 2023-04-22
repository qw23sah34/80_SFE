# -*- coding: utf-8 -*-
"""
Created on 30.03.23

@author: shv
@description:
    Stores FE Datatype and function for FE caculations
"""
import numpy as np

class FEData(object):
    """This is a doc string.
    Class ist meant for FE data of 1D-SFE solver"""
    def __init__(self):
        self.e = []
        self.p = FE_Point()
        self.E = 0.0
        self.nu = 0.0
        self.nPoints = 0
        self.nElements = 0
        self.ElementType = 0
        
    def mesh(self, initData):
        if (initData.ElementType == 1):
            #
            # Save general data
            #
            self.E  = initData.E
            self.nu  = initData.nu
            self.ElementType = initData.ElementType
            self.nPoints = len(initData.x)
            # 
            # So far only main points are known. Dependent on approximation
            # degree additional points can be added. Everything should be 
            # first defined at element level. Define all elements
            #
            self.nElements = self.nPoints - 1
            if (initData.ElementType == 1):
                for iElement in range(self.nElements):
                    element = FE_Element()
                    element.type = 1 # Beam element. 2 main points.
                    element.degree = 1
                    element.nPoints = 2
                    element.pointMainIX =  np.zeros(element.nPoints, int)
                    if (iElement == 0):
                        ix_i = iElement
                        ix_j = ix_i + 1
                    else:
                        ix_i = self.e[iElement-1].pointMainIX[-1]
                        ix_j = ix_i + 1
                        
                    ix_i_o = iElement
                    ix_j_o = ix_i_o + 1
                    element.pointMainIX[0] = ix_i
                    element.pointMainIX[1] = ix_j
                    element.p = FE_Point()
                    #
                    element.p.x = np.array([initData.x[ix_i_o], initData.x[ix_j_o]])
                    element.p.A = np.array([initData.A[ix_i_o], initData.A[ix_j_o]])
                    element.p.F = np.array([initData.forceX[ix_i_o], initData.forceX[ix_j_o]])
                    element.p.u = np.array([0.0, 0.0])
                    if (initData.degree == 2):
                        #
                        # Second order form function (quadratic). One 
                        # additional point in the middle of element
                        #
                        element.degree = 2
                        element.nPoints = element.nPoints + 1
                        element.p.x = np.insert(element.p.x, 1, sum(element.p.x)/2.0)
                        element.p.A = np.insert(element.p.A, 1, sum(element.p.A)/2.0)
                        element.p.F = np.insert(element.p.F, 1, 0.0)
                        element.p.u = np.insert(element.p.u, 1, 0.0)
                        element.pointMainIX = \
                            np.insert(element.pointMainIX, 2, element.pointMainIX[-1]+1)
                    #
                    element.A = (initData.A[ix_i_o] + initData.A[ix_j_o])/2.0
                    element.length = initData.x[ix_j_o] - initData.x[ix_i_o]
                    #
                    self.e.append(element)
            
    def calc_element_stiffness_matrix(self):
        """The stiffness matrix for each element is calculated here.
        So far the integration of the stiffness matrix is hardcoded"""
        #
        def get_gauss_integral_B(i,j,a,b,degree):
            """Calculates the integral. 
            Borders: from a to b.
            Integral: B_i, B_j where B are the B-functions of the degree"""
            def get_B(ix, x):
                if (ix == 0):
                    B = -3.0/L + 4.0*x/pow(L,2)
                elif (ix == 1):
                    B =  4.0/L - 8.0*x/pow(L,2)
                elif (ix == 2):
                    B = -1.0/L + 4.0*x/pow(L,2)
                return B
            
            # Define X_i and alpha_i (from Table)
            if (degree == 1):
                X_i = np.array([0.0], float)
                alpha_i = np.array([2.0], float)
            elif (degree == 2):
                X_i = np.array([-np.sqrt(1.0/3.0),np.sqrt(1.0/3.0)], float)
                alpha_i = np.array([1.0,1.0], float)
            #
            GI = 0.0
            for iD in range(degree):
                x = (b-a)/2.0*X_i[iD] + (a+b)/2.0
                GI = GI + get_B(i,x)*get_B(j,x)*alpha_i[iD]
            #
            GI = GI*(b-a)/2.0
            return GI
            #
        E = self.E
        for iElement in range(self.nElements):
            nPoints = self.e[iElement].nPoints
            #!shv
            # Define element stifness matrix. It depend on element order
            # Is it defined ONLY by the element's order? Can it be hard coded
            # only in dependence of the order?
            #
            # I assume that stiffness matrix can be hard coded for every type
            # of element defined, and the is no reason to execute the following
            # steps. But for the sake of comprehence i'll do it.
            #
            # Everything that is needed to calculate the stiffness matrix for 
            # any type of element is:
            #   B[1:nElementPoints] - function, which is the first derivative of 
            #                         the form function N in respect to x.
            #   E, L, Am            - constants
            #   gedree              - degree of the form function is necessary 
            #                         for numerical Gauss-Legendre-Integration 
            #
            L = self.e[iElement].length
            A = self.e[iElement].A
            degree = self.e[iElement].degree
            if (degree == 1):
                #
                # C = E*A/L
                #
                C = self.E*self.e[iElement].A/self.e[iElement].length
                self.e[iElement].K = np.array([[ C, -C], [-C,  C]], float)
            elif (degree == 2):
                self.e[iElement].K = np.zeros((nPoints, nPoints), float)
                for i in range(nPoints):
                    for j in range(nPoints):
                        self.e[iElement].K[i,j] = E*A*get_gauss_integral_B(i,j,0.0,L,degree)
                        #
                #
            #
         
class FE_Element(object):
    def __init__(self):
        # 1 - beam
        # 2 - beam :-)
        # 3 - beam :-)
        self.type = 0
        self.degree = 1
        self.nPoints = 0
        self.A = []
        self.pointMainIX = []
        self.length = 0.0
        self.K = []
        self.p = []
        
class FE_Point(object):
    def __init__(self):
        self.x = []
        self.A = []
        self.u = []
        self.F = []
        
