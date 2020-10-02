import numpy as np
class Peak:
    def __init__(self, startidx):
        self.born = self.left = self.right = startidx
        self.died = None

    def get_persistence(self, seq):
        return float("inf") if self.died is None else seq[self.born] - seq[self.died]

def get_persistent_homology(seq):
    peaks = [] #all peaks as defined in the above class
    # Maps indices to peaks
    idxtopeak = [None for s in seq] # each index has a main peak it is attributed to i.e. the index of the top of the island
    # Sequence indices sorted by values bottom up, ala points emerging if the water level is reduced
    indices = range(len(seq))
    indices = sorted(indices, key = lambda i: seq[i], reverse=True)

    # Each point is investigated with respect to its neighbours in descending order
    for idx in indices:
        lftdone = (idx > 0 and idxtopeak[idx-1] is not None) #idx is the index of the highest wo lowest point, if the index was 0 there would be no idx-1
        rgtdone = (idx < len(seq)-1 and idxtopeak[idx+1] is not None) # if the adjacent index left or right already has a peak attributed lftdone/rightdone is True
        il = idxtopeak[idx-1] if lftdone else None #il and ir give the peak index of the left or right neighbour respectively
        ir = idxtopeak[idx+1] if rgtdone else None

        # New peak born
        if not lftdone and not rgtdone: #if neither the left nor the right neighbor has a peak attributed the point must be the top of an island
            peaks.append(Peak(idx)) #a new peak is added to peaks with born left and right all equal to the index
            idxtopeak[idx] = len(peaks)-1 # new peaks are numbered from 0 up, and the index is attributed the peaknumber

        # Directly merge to next peak left
        if lftdone and not rgtdone: #only a left neighbor, the point is added to the peak on the left
            peaks[il].right += 1 #the dominant peak gets on more left neighbour, which is kinda useless if i look at it noch
            idxtopeak[idx] = il #the index gets attributed whatever peak the left neighbour is attributed to

        # Directly merge to next peak right, same theing as above only to the right
        if not lftdone and rgtdone:
            peaks[ir].left -= 1
            idxtopeak[idx] = ir

        # Merge left and right peaks
        if lftdone and rgtdone:
            # Left was born earlier: merge right to left
            if seq[peaks[il].born] > seq[peaks[ir].born]:
                peaks[ir].died = idx
                peaks[il].right = peaks[ir].right
                idxtopeak[peaks[il].right] = idxtopeak[idx] = il
            else:
                peaks[il].died = idx
                peaks[ir].left = peaks[il].left
                idxtopeak[peaks[ir].left] = idxtopeak[idx] = ir

    # This is optional convenience
    return sorted(peaks, key=lambda p: p.get_persistence(seq), reverse=True)

def find_peaks(number,x,data):
    x = np.array(x)
    data = np.array(data)
    peak_ind = []
    peak_x = []
    for i in range(number):
        #print(get_persistent_homology(data)[i].born)
        peak_ind += [get_persistent_homology(data)[i].born]
        
       # print(x)
        peak_x += [x[get_persistent_homology(data)[i].born] ]
    return np.array(peak_x),np.array(data[peak_ind])
