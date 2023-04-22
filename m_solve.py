# -*- coding: utf-8 -*-
"""
Created on 30.03.23

@author: shv
@description:
    Main module to solve the problems that have to be solved.
"""
import numpy as np

class sol(object):
    """ All FE solution data and routines"""
    def __init__(self):
        self.K_tot = []
        self.u = []
        self.F = []
    #
    def create_total_system (self, FE_data):
        #
        # Calculate number of points
        #
        nPoints = 1
        for iElement in range(FE_data.nElements):
            #
            # -1 because it is assumed (for ease purposes) that 
            # elements share exactly one single point.
            #
            nPoints = nPoints + FE_data.e[iElement].nPoints - 1
        #
        # Allocate all data structures
        #
        self.K_tot = np.zeros((nPoints, nPoints), float)
        self.u = np.zeros(nPoints, float)
        self.F = np.zeros(nPoints, float)
        #
        # Calculate the total striffness matrix
        #
        for iElement in range(FE_data.nElements):
            n = FE_data.e[iElement].pointMainIX[0]
            delta = FE_data.e[iElement].pointMainIX[-1]+1 - FE_data.e[iElement].pointMainIX[0]
            for i in range(delta):
                for j in range(delta):
                    self.K_tot[n+i,n+j] = self.K_tot[n+i,n+j] + FE_data.e[iElement].K[i,j]
        #
        # Compile force vector. The force can be either defined,
        # undefined (in other words equal to zero), or can be -999.
        # In the case of -999 the clamped boundary condition of the 
        # point is meant. Force unknown, and displacement is zero.
        #
        for iElement in range(FE_data.nElements):
            for iPoint in range(FE_data.e[iElement].nPoints):
                ix_glob = FE_data.e[iElement].pointMainIX[iPoint]
                self.u[ix_glob] = -999.0
                self.F[ix_glob] = -999.0
                #
                if (FE_data.e[iElement].p.F[iPoint] >= 0.0):
                    self.F[ix_glob] = FE_data.e[iElement].p.F[iPoint]
                elif (np.absolute(FE_data.e[iElement].p.F[iPoint]-(-999.0)) < 1.0e-6):
                    self.u[ix_glob] = 0.0
        
    #
    def solve(self):
        """ Solve the linear equation problem """
        #
        # The equation system is not uniform. At some places only displacement 
        # is known. At other - only force. Decompose the system into two 
        # subsystems. 
        # First - where forces are known.
        # Second - where displacements are known. 
        #
        # self.K_tot*self.u = self.F
        #
        displ_unkn = \
            np.logical_not((np.absolute(self.F - (-999.0)) < 1.0e-6))
        K1 = self.K_tot[displ_unkn,:]
        #
        K1 = K1[:,displ_unkn]
        u1 = np.zeros(np.count_nonzero(displ_unkn))
        F1 = self.F[displ_unkn]
        #
        # Now we have K1*u1 = F1
        # where K1 and F1 are known, and u1 must be calculated
        #
        u1 = np.linalg.solve(K1, F1)
        #
        # After u1 displacements are known the second system should be solved
        #
        K2 = self.K_tot[~displ_unkn]
        u2 = self.u
        counter = 0
        for i in range(len(self.u)):
            if (displ_unkn[i]):
                u2[i] = u1[counter]
                counter = counter + 1
        #
        F2 = np.matmul(K2,u2)
        self.u = u2
        counter = 0
        for i in range(len(self.F)):
            if (displ_unkn[i]):
                #
                # Where displacement is unknown, force is known
                #
                self.F[i] = self.F[i]
            elif (~displ_unkn[i]):
                self.F[i] = F2[counter]
                counter = counter + 1
        #
        pass
    #
    def evaluate_results(self, FE_data):
        """ Calculate displacements and stresses along whole part """
        import matplotlib.pyplot as plt
        from SFE_main import BeamType
        #
        # Displacement 
        #
        length = FE_data.e[-1].p.x[-1]
        displacement = np.zeros(int(length), float)
        stresses = np.zeros(int(length), float)
        #
        # Displacement function has the form: u(x) = (1-x/L)*u_i + (x/L)*u_j
        # for each element, if the form function was of first order.
        #
        for iElement in range(FE_data.nElements):
            u = np.zeros(FE_data.e[iElement].nPoints, float) # nodes displacements
            x_n = np.zeros(FE_data.e[iElement].nPoints, int) # nodes x-coordinates
            for iPoint in range(FE_data.e[iElement].nPoints):
                u[iPoint] = self.u[FE_data.e[iElement].pointMainIX[iPoint]]
                x_n[iPoint] = int(FE_data.e[iElement].p.x[iPoint])
            #
            L = FE_data.e[iElement].length
            for i in range(x_n[0], x_n[-1]):
                x = float(i-x_n[0])
                if (FE_data.e[iElement].degree == 1):
                    displacement[i] = (1.0-x/L)*u[0] + (x/L)*u[1]
                    stresses[i] = (FE_data.E/L)*(-u[0] + u[1])
                    pass
                elif (FE_data.e[iElement].degree == 2):
                    displacement[i] = (1.0 - 3.0*(x/L) + 2.0*pow(x/L,2))*u[0] + \
                        (0.0 + 4.0*(x/L)*(1.0-x/L))*u[1] + (-x/L + 2.0*pow(x/L,2))*u[2]
                    stresses[i] = (-3.0/L + 4.0*x/pow(L,2))*u[0] + \
                        (4.0/L - 8.0*x/pow(L,2))*u[1] + (-1.0/L + 4.0*x/pow(L,2))*u[2]
                    stresses[i] = stresses[i]*FE_data.E
                    pass
            #
        #
        plt.figure(1)
        plt.plot(range(int(length)), displacement)
        plt.grid(True)
        plt.figure(2)
        #plt.show()
        plt.plot(range(int(length)), stresses)
        plt.grid(True)
        

                
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        