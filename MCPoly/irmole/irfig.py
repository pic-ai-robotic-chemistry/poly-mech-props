from scipy.stats import rankdata
import os
import re

class irfig:
    def __init__(self, file, loc='./'):
        f = open(loc+file, 'r')
        freq = []
        intensity = []
        peaks = []
        peaks_freq = []

        for line in f:
            a = re.findall(r'[0-9]+\.[0-9]+', line)
            if len(a) == 2:
                freq.append(eval(a[0]))
                intensity.append(eval(a[1])-1000)

        self.freq = freq
        self.intensity = intensity
        self.I = intensity
        self.intensity_R = [-i for i in intensity]
        self.I_R = [-i for i in intensity]

        for i, fq in enumerate(self.intensity):
            if i == 0 or i == len(self.intensity) - 1:
                continue
            try:
                if self.intensity[i] < self.intensity[i-1] and self.intensity[i] < self.intensity[i+1]:
                    peaks.append(self.intensity[i])
                    peaks_freq.append(self.freq[i])
            except:
                continue

        self.peaks = peaks
        self.peaks_R = [-i for i in peaks]
        self.peaks_freq = peaks_freq
        
        f.close()

    def toprank(self, num):
        rankorder = rankdata(self.peaks)
        #print(rankorder)
        
        tops = []
        tops_freq = []
        
        for j in range(1, num+1):
            for i,rank in enumerate(rankorder):
                if rank == j:
                    tops.append(self.peaks[i])
                    tops_freq.append(self.peaks_freq[i])

        self.tops = tops
        self.tops_R = [-i for i in tops]
        self.tops_freq = tops_freq
        return tops, tops_freq