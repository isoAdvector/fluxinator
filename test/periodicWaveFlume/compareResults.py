#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import subprocess
import os
plt.rcParams.update({
#    "text.usetex": True,
    "font.family": "monospace",#"Times New Roman",
    "font.monospace": "DejaVu Sans Mono",
    "font.size":20,
    'lines.linewidth':1
})
import argparse

import re

def get_parameter_value(text, parameter_name):
    # Create a pattern to match the parameter name followed by its value
    pattern = r'{}:\s*(\S+)'.format(re.escape(parameter_name))
    
    # Search for the pattern in the text
    match = re.search(pattern, text)
    
    # Return the matched value if found, None otherwise
    if match:
        return match.group(1)
    else:
        return None

def count_time_lines(log_file_path):
    count = 0
    
    with open(log_file_path, 'r') as file:
        for line in file:
            if line.strip().startswith("Time = "):
                count += 1
    
    return int(count)

def get_last_clocktime(log_file_path):
    # Regular expression to extract ClockTime value
    pattern = r'ClockTime = (\d+\.?\d*) s'
    last_clocktime = None
    
    with open(log_file_path, 'r') as file:
        for line in file:
            if "ExecutionTime" in line and "ClockTime" in line:
                match = re.search(pattern, line)
                if match:
                    last_clocktime = match.group(1)
    
    return float(last_clocktime)

# Set up argument parser
parser = argparse.ArgumentParser(description='Process some case names.')

# Add an argument for caseNames
parser.add_argument(
    '--scanDir',
#    nargs='+',  # '+' means one or more arguments
    required=True,  # Make it a required argument
    help='Directory where cases are found'
)

# Add an argument for caseNames
parser.add_argument(
    '--caseNames',
    nargs='+',  # '+' means one or more arguments
    required=True,  # Make it a required argument
    help='List of case names'
)

parser.add_argument(
    '--figName',
#    type=str,  # '+' means one or more arguments
    required=True,  # Make it a required argument
    help='List of case names'
)

# Parse the arguments
args = parser.parse_args()

# Access the caseNames argument
caseNames = args.caseNames

# Access the figName argument
figName = args.figName

#Plotting position
fig, ax = plt.subplots(2)
fig.set_size_inches(20, 14)

scanName = args.scanDir

zorders1 = [0, 1, 2, -1, 4, 5, 6, 7, 8]
zorders2 = [0, -1, 2, 3, 4, 5, 6, 7, 8]
title_string = ""
Umax = 0
Umin = 1e6
maxT = 0
refTime = 1
refSteps = 1
tRight = 200
tLeft = -1
for count, caseName in enumerate(caseNames):
    caseDir=scanName+'/'+caseNames[count]
    print(caseDir)
    os.chdir(caseDir)
    t, WG1 = np.loadtxt('postProcessing/interfaceHeight1/0/height.dat', comments="#", unpack=True, usecols = (0, 2))
    mask = (t <= tRight) & (t >= tLeft)
    WG1 = WG1[mask]
    t = t[mask]
    maxT = max(max(t), maxT)
    with open("parmScannerParms", 'r') as file:
        parmScannerParms = file.read().strip()
    application = get_parameter_value(parmScannerParms, "application")
#    application = "interIsoGoam"
    legend="Id: "+caseNames[count]
    legend+=", Solver: "+application.ljust(12)
    if count < len(caseNames)-1: #Skip last analytic case here
        nTimeSteps = count_time_lines("log."+application)
        ClockTime = get_last_clocktime("log."+application)
        if count == 0: #Use 1st case as reference
            refTime = ClockTime
            refSteps = nTimeSteps
        relTime = ClockTime/refTime
        relSteps = nTimeSteps/refSteps
        legend+=" Time: "+str(format(ClockTime,".0f"))+"s ("+str(format(relTime,".2f"))+"),"
        legend+=" Steps: "+str(nTimeSteps)+" ("+str(format(relSteps,".2f"))+")"
        nz = get_parameter_value(parmScannerParms, "nz")
        maxCo = get_parameter_value(parmScannerParms, "maxCo")
        deltaT = get_parameter_value(parmScannerParms, "deltaT")
    #Plotting position
    line, = ax[0].plot(t, WG1, '-', label=legend, zorder=zorders1[count])

    t, maxU = np.loadtxt('postProcessing/fieldMinMax1/0/fieldMinMax.dat', unpack=True, skiprows=2, usecols = (0, 2))
    mask = (t <= tRight) & (t >= tLeft)
    maxU = maxU[mask]
    t = t[mask]
    Umax = max(Umax,max(maxU))
    Umin = min(Umin,min(maxU))
    line, = ax[1].plot(t, maxU, '.', label=legend, zorder=zorders2[count])

    os.chdir("../..")

u_crest = 0.44453480326804;
Umax = min(Umax,10*u_crest)
Umax = 2.4
hmax = 0.081866967464134;
hmin = -0.043133032535866;
ax[0].plot([0, maxT],[hmax, hmax],'k')
ax[0].plot([0, maxT],[hmin, hmin],'k')
ax[0].set_ylim(1.1*hmin,1.2*hmax)
#ax[0].set_title(f"Stream Function Wave (nz = {nz}, maxCo = {maxCo})")
ax[0].set_title(f"Stream Function Wave (nz = {nz}, deltaT = {deltaT})")
#plt.rcParams['text.usetex'] = False
#plt.rcParams['text.usetex'] = True
ax[0].set(xlabel=r'$t$ [s]', ylabel=r'$\eta$ [m]')
ax[1].plot([0, maxT],[u_crest, u_crest],'k')
ax[1].set(xlabel=r'$t$ [s]', ylabel=r'max$|U|$ [m/s]')
ax[1].legend(loc='center left',fontsize=12)
##ax[1].set_ylim(Umin, Umax)
#ax[0].xlabel(r'$t$ [s]')
#ax[0].ylabel(r'$\eta$ [m]')
ax[0].grid(True)
ax[1].grid(True)
plt.show()
#plt.savefig(figName+'.png', bbox_inches='tight', dpi=300)
