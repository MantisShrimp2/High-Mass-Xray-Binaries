#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 15:04:32 2024

@author: karan
"""

import numpy as np
import pandas as pd
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import matplotlib.pyplot as plt
from astropy.table import Column, join, Table, vstack, hstack
from astropy.io import ascii
class GalacticTraceback:
    def __init__(self,table):
        
        # self.l = long
        # self.b = lat
        # self.mu_l = mu_l_cosb
        # self.mu_b = mu_b
        self.table = table
        
    def traceback_time(self):
        '''
        Calculate the traceback time for a star
        The time is takes to return to the galactic midplane in years
        
        input:
        self
        b - galactic longitude degrees
        mu_b- proper motion in b mas/yr

        Returns
        -------
        None.

        '''

        lat = self.table['b']
        mu_b = self.table['pm_b_poleski']
        
        #if lat == 0:
           #raise ValueError('Proper motion in b cannot be zero ')
        #convert mu_b to degree/yr 1deg = 3.6 million milliarcseconds
        mu_b_deg = mu_b /(3.6e6) 
        trace_time = np.abs(lat/mu_b_deg)
        
        self.table['Trace Time'] = np.array(trace_time)/float(1e6)
        self.table['Trace Time'].unit = 'Million years'
        
        
        return self.table
        
    def trace_linear_path(self, time_step=1000):
        """
        Trace the path of the star in Galactic coordinates until b = 0, using the Euler method.

        Parameters:
        - time_step (float): Step size in years for tracing the path.
        - max_steps (int): Maximum number of steps for tracing.

        Returns:
        - path (list of tuples): List of (l, b) pairs representing the path.
        """
        # Convert proper motions to degrees per year
        # long = self.table['l']
        # lat = self.table['b']
        # mu_l = self.table['pm_l_poleski']
        # mu_b = self.table['pm_b_poleski']
        # mu_l_deg_per_year = mu_l/ 3.6e6
        # mu_b_deg_per_year = mu_b / 3.6e6

        # Initialize path
        traced_paths = []
        for row in self.table:
            l = row['l']
            b = row['b']
            mu_l = row['pm_l_poleski']
            mu_b = row['pm_b_poleski']
            mu_l_deg_per_year = mu_l/ 3.6e6
            mu_b_deg_per_year = mu_b / 3.6e6
            path = [(l,b)]

            max_steps = int(1e6)
        # Trace back in time using Euler's method
            current_l = l
            current_b = b
            for _ in range(max_steps):
                # Update coordinates using Euler method
                #fixed timestep
                current_l -= mu_l_deg_per_year * time_step  # Adjust longitude
                current_b -= mu_b_deg_per_year * time_step  # Adjust latitude
    
                # Wrap longitude to [0, 360)
                current_l %= 360
                current_b %= 90
            

            # Add to path
                path.append((current_l, current_b))
                
            # Stop if b is close to zero (or negative)
            if current_b <1e-3:
                break
            traced_paths.append(path)
        
        #add to table
        self.table['Galactic Path'] = traced_paths
            
            
        return self.table

# galactic_longitude = 150.0  # degrees
# galactic_latitude = 10.0  # degrees
# proper_motion_l_cos_b = 5.0  # mas/yr
# proper_motion_b = -2.0  # mas/yr (negative to simulate motion below the Galactic plane)



