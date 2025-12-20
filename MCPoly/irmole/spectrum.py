import re
import os
import matplotlib.pyplot as plt
import numpy as np

def spectrum(files, loc='./', legend_rename=[]):
    fig, ax = plt.subplots()

    for file in files:
        f0 = open(loc+file, 'r')
        
        freq = []
        intensity = []
        
        for line in f0:
            a = re.findall(r'[0-9]+\.[0-9]+', line)
            if len(a) == 2:
                freq.append(eval(a[0]))
                intensity.append(eval(a[1])-1000)
        
        ax.plot(freq, intensity)
        
        f0.close()
    
    ax.set_xlim(4000, 300)
    ax.set_xlabel('wavenumber (cm$^{-1}$)')
    ax.set_ylabel('Relative Intensity')
    if legend_rename == []:
        ax.legend([files], bbox_to_anchor=[1.22,0.9])
    else:
        ax.legend(legend_rename, bbox_to_anchor=[1.22,0.9])
    #plt.title('')
    plt.show()