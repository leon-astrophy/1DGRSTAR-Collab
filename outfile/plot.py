#import#
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#search for all files #
filename = []
for root, dirs, files in os.walk('.'):
  for file in files:
    if file.endswith("dat"):
      filename.append(os.path.join(file))

# plot #
for i in range (0, len(filename)):

  #read#
  mr = pd.read_csv(filename[i], delim_whitespace=True, header=None)

  # look for how many fluids #
  fluid = filename[i].split('-')
  
  #case by case#
  if(fluid[0] == '1F'):

    #assign#
    nmlogrho = mr[0].to_numpy()
    nmmass = mr[1].to_numpy()
    nmrad = mr[2].to_numpy()

    #plot#
    plt.plot(nmlogrho, nmmass, label='NM')
    plt.xlabel('Log10 NM Central Density (g cm$^{-3}$)', size=15)
    plt.ylabel(r'NM Mass ($M_{\odot}$)', size=15)
    plt.tick_params(axis="both", labelsize=15)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('rhovsmass.png')
    plt.clf()

    #plot#
    plt.plot(nmrad, nmmass, label='NM')
    plt.xlabel(r'NM Radius ($R_{\odot}$)', size=15)
    plt.ylabel(r'NM Mass ($M_{\odot}$)', size=15)
    plt.tick_params(axis="both", labelsize=15)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('radvsmass.png')
    plt.clf()

  # 2f case #
  elif(fluid[0] == '2F'):

    #assign#
    nmlogrho = mr[0].to_numpy()
    dmmass = mr[1].to_numpy()
    nmmass = mr[2].to_numpy()
    dmrad = mr[3].to_numpy()
    nmrad = mr[4].to_numpy()

    #plot#
    plt.plot(nmlogrho, nmmass+dmmass, label='NM')
    plt.xlabel('Log10 NM Central Density (g cm$^{-3}$)', size=15)
    plt.ylabel(r'NM Mass ($M_{\odot}$)', size=15)
    plt.tick_params(axis="both", labelsize=15)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('rhovsmass.png')
    plt.clf()

    #plot#
    plt.plot(nmrad, nmmass+dmmass, label='NM')
    plt.xlabel(r'NM Radius ($R_{\odot}$)', size=15)
    plt.ylabel(r'NM Mass ($M_{\odot}$)', size=15)
    plt.tick_params(axis="both", labelsize=15)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('radvsmass.png')
    plt.clf()
