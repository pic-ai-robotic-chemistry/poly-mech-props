def relat_peaks(selected_freq, intensity, freq, tolerance=10.0, neednum=False):
    real_tops = []
    real_tops_freq = []
    js = []
    
    for t_fq in selected_freq:
        base = 0
        base_freq = 0
        j = 0
        for i,fq in enumerate(freq):
            if abs(t_fq-fq) < tolerance:
                #print(t_fq, fq, I[i])
                if intensity[i] > base:
                    base = intensity[i]
                    base_freq = freq[i]
                    j = i
        real_tops.append(base)
        real_tops_freq.append(base_freq)
        js.append(j)

    if neednum == True:
        return real_tops, real_tops_freq, js
    else:
        return real_tops, real_tops_freq