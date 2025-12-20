from scipy.stats import rankdata
import re
import os

class irout:

    """
    The method to get the information about infrared spectrum in output file.
    irout(file, loc='./'))
    file: The name of output file.
    loc: Your location of output file.

    Functions:
    freq: All frequencies in the output file.
    intensity/I: All intensities of relevant frequencies in the output file.

    Example:
        Input:
            from MCPoly.irmole import irfig
            a = irfig('et.out')
            print(a.freq)
            print(a.intensity) # You can use a.I instead.
        Output:
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 303.02, 827.41, 827.6, 996.26, 1221.94, 1222.0, 1411.48, 1421.8, 1501.2, 1501.27, 1504.87,
             1504.91, 3026.99, 3028.07, 3070.5, 3070.74, 3096.03, 3096.23]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.55, 3.55, 0.0, 0.0, 0.0, 1.23, 0.0, 0.0,
             0.0, 10.08, 10.14, 0.0, 61.64, 0.0, 0.0, 67.1, 67.14]

    Tip: You can use toprank to get the highest peak in the spectrum.
    
    """

    def __init__(self, file, loc='./'):
        f = open(loc+file+'.out', 'r')
        n0 = 0
        n2 = 0
        n3 = 0
        freq = []
        intensity = []
        
        for line in f:
            if n0 == 1 and n2 == 1:
                if n3 == 0:
                    c1 = re.search(r'[0-9]+\:', line)
                    #print(eval(c1.group(0)[:-1]))
                    freq = eval(c1.group(0)[:-1])*[0.0]
                    intensity = eval(c1.group(0)[:-1])*[0.0]
                    n3 = 1
                    
                if n3 == 1:
                    c2 = re.findall(r'\-?[0-9]+\.[0-9]+', line)
                    if len(c2) > 6:
                        freq.append(eval(c2[0]))
                        intensity.append(eval(c2[2]))
                    n3 = 1
        
                x3 = re.search('The epsilon (eps)', line)
                if x3:
                    break
            
            x0 = re.search('IR SPECTRUM', line)
            if x0:
                n0 = 1
            x2 = re.search('-------------------------------------', line)
            if n0 == 1 and x2:
                n2 = 1

        self.freq = freq
        self.intensity = intensity
        self.I = intensity
        
        f.close()

    def toprank(self, num):
        """
        The method to get the highest peak in the spectrum.
        toprank(num)
        num: numbers of peaks. e.g. If num = 3, it will output the first three highest peak in the spectrum.

        After using it, two new data will be created.
        tops_freq: The frequencies of highest peaks in the Hessian data file.
        tops: The intensities of highest peaks in the Hessian data file.

        Example:
            Input:
                from MCPoly.irmole import irfig
                a = irfig('et.out')
                a.toprank(3)
                print(a.tops)
                print(a.tops_freq)
            Output:
                [67.14, 67.1, 61.64]
                [3096.23, 3096.03, 3028.07]

        Tip: If you use code 'a.toprank(3)' above, it will directly output a.tops and a.tops_freq.
        """
        
        rankorder = rankdata([-i for i in self.intensity])
        #print(rankorder)
        
        tops = []
        tops_freq = []
        
        for j in range(1, num+1):
            for i,rank in enumerate(rankorder):
                if rank == j:
                    tops.append(self.intensity[i])
                    tops_freq.append(self.freq[i])

        self.tops = tops
        self.tops_freq = tops_freq
        return tops, tops_freq