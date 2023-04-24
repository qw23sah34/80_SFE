# -*- coding: utf-8 -*-
"""
Created on 30.03.23

@author: shv

@description:
    Super simple FE-Solver for learning practice and better understanding
    of FE fundamentals.
"""

import os
os.chdir(os.path.dirname(os.path.realpath(__file__)))
import c_init_data
import m_FE_data
import m_solve as FE_SOLUTION

global BeamType
BeamType = 1
#
# Read control data. There is defined:
#    Material constants: E, nu
#    Element type: beam element type
#    FE-Point information:  - x Position
#                           - heigth of the part at position x
#                           - width of the part at position x
#                           - applied external force,
#                             (-999 means clamped boudary condition)
#
initData = c_init_data.initData()
initData.read('.\control.ste')
#
# From FE-points data calculate element data
#
FE_data = m_FE_data.FEData()
FE_data.mesh(initData)
#
# Create element stiffness matrix
#
FE_data.calc_element_stiffness_matrix()
#
# Solution initialization
#
sol = FE_SOLUTION.sol()
#
# Composing the system to be solved
#
sol.create_total_system(FE_data)
#
# Solve the problem.
#
sol.solve()
sol.evaluate_results(FE_data)
print("The displacements are:")
print(sol.u)
print("The forces are:")
print(sol.F)
