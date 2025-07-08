#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 11:57:18 2025

@author: caiazzo
"""

import numpy as np
import matplotlib.pyplot as plt
import os,sys

def plot_mean_phenotype(filename, yline_value, output_filename):
    # Load data
    dt_write = 1000
    mean_pheno = np.loadtxt(filename)
    endval = len(mean_pheno) - 1
    x = np.arange(0, endval + 1) * dt_write  # timestep

    # Plotting
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x/1e6, mean_pheno, color='black', linewidth=3, linestyle='-',label="Results")
    ax.axhline(y=yline_value, linewidth=3, color='red', linestyle='--',label="Asymptotic")
    
    ax.set_xlim(0, 10)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel('Timestep ($\cdot 10^6$)', fontsize=20)
    ax.set_ylabel('Average phenotype', fontsize=18)
    ax.tick_params(labelsize=20)
    ax.legend(framealpha=1,frameon=False,bbox_to_anchor=(.85,1.0),
                    loc='upper center').set_draggable(True)
    #plt.box(False)
    ax.yaxis.grid(color='grey', linestyle='--', linewidth=0.5) # vertical lines
    ax.xaxis.grid(color='grey', linestyle='--', linewidth=0.5) # vertical lines
    # Save figure
    if not os.path.exists(output_directory):
        print(" ** cannot print picture:", filename, ". Path does not exist.")
    else:
        fig.savefig(output_filename, bbox_inches='tight')
    fig.show()
    plt.close(fig)

# Run plots
output_directory = '../tests/paper-1/figures/'
input_directory = '../tests/paper-1/output/'

plot_mean_phenotype(input_directory + 'constant_02_100_mean_phenotype.txt', 20/23, output_directory + 'average_phenotype_o2_100.png')
plot_mean_phenotype(input_directory + 'constant_02_20_mean_phenotype.txt', 20/27, output_directory + 'average_phenotype_o2_20.png')
plot_mean_phenotype(input_directory + 'constant_02_6_mean_phenotype.txt', 15/29, output_directory + 'average_phenotype_o2_6.png')
plot_mean_phenotype(input_directory + 'constant_02_1_mean_phenotype.txt', 10/61, output_directory + 'average_phenotype_o2_1.png')
