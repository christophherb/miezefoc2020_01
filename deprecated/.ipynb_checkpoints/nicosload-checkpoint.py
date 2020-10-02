import time
import numpy as np
from numpy import loadtxt


def read_nicos_file(filename):
    '''
    Load a NICOS dat file (Georg Brandl style)
    returns: tuple of:
    colnames    Names of colums of the scan
    colunits    Units of colums of the scan
    cols        actual scan data
    meta        metadata from NICOS header as dictonary
    '''
    fp = open(filename, 'r')
    meta = {}
    dtline = fp.readline()
    if not dtline.startswith('### NICOS data file'):
        raise RuntimeError('%r does not appear to be a NICOS data file' %
                           filename)
    ctime = time.mktime(time.strptime(
        dtline[len('### NICOS data file, created at '):].strip(),
        '%Y-%m-%d %H:%M:%S'))
    meta['created'] = ctime
    remark = ''
    for line in iter(fp.readline, ''):
        line = line.decode('ascii', 'ignore').encode('ascii')
        if line.startswith('### Scan data'):
            break
        if line.startswith('# '):
            items = line.strip().split(None, 3)
            try:
                oval, unit = items[3].split(None, 1)
                val = float(oval)
            except (IndexError, ValueError):
                try:
                    oval = items[3]
                    val = float(oval)
                except ValueError:
                    val = items[3]
                except IndexError:
                    continue
                unit = None
            key = items[1]
            if key.endswith(('_offset', '_precision')):
                # we don't need these
                continue
            if key.endswith('_value'):
                key = key[:-6]
            meta[key] = val
    if remark and 'title' in meta:
        meta['title'] += ', ' + remark
    colnames = fp.readline()[1:].split()
    colunits = fp.readline()[1:].split()
    def convert_value(s):
        try:
            return float(s)
        except ValueError:
            return 0.0
    cvdict = dict((i, convert_value) for i in range(len(colnames))
                  if colnames[i] != ';')
    colnames = [name for name in colnames if name != ';']
    colunits = [unit for unit in colunits if unit != ';']
    usecols = cvdict.keys()
    coldata = loadtxt(fp, converters=cvdict, usecols=usecols, ndmin=2,
                      comments='#')
    fp.close()
    if not coldata.size:
        raise RuntimeError('empty data file')
    cols = dict((name, coldata[:,i]) for (i, name) in enumerate(colnames))
    # return values:
    # - column names
    # - column units
    # - data by columns
    # - metadata (file header) as a dictionary
    return colnames, colunits, cols, meta

#--------------------------------------------------------------------------------------------------------------------

def read_tof_file(filename, pixels=(128,128), timechannels=16, foils=8):
    '''Read a tof file and return data and metadata
    filename: Name of the file
    pixels: tupel of detector pixel size
    timechannels: number of timechannels
    foils: number of detector foils
    
    Returns:
    data: Raw detector data reshaped (e.g. 128x128x16x18)
    sum_img: summed over all time channels and foils
    meta: metadata collected from file'''
    
    with open(filename, 'rb') as data:
        data = np.fromfile(data, dtype=np.int32)[0:128*128*16*8]
        data = data.reshape((8, 16, 128,128))
    
    #Sum up over the image
    sum_image = np.zeros((128,128))
    for i in list(range(6)):
        for j in list(range(15)):
            sum_image = sum_image + data[i][j]
    
    with open(filename,'r') as fh:
        meta = {}
        startnum =1e9
        for num, curline in enumerate(fh, 1):
            # check if the current line starts with "#"
            if curline.startswith("### NICOS"):
                startnum = num
            if num >= startnum:
                if curline.startswith('#'):
                    pass
                else:
                    #print(curline.split())
                    items = curline.split()
                    try:
                        key = items[0]
                        if key.endswith(('_offset', '_precision', 'status', '_userlimits', '_alias')) or key.startswith(('Exp_', 'Sample', 'reseda_')):
                            # we don't need these
                            continue
                        try:
                            val = float(items[2])
                        except:
                            val = items[2]
                    except:
                        pass
                    meta[key] = val
                    
    return data, sum_image, meta

#--------------------------------------------------------------------------------------------------------------------

def read_pad_file(filename, pixels=(128,128)):
    '''Read a pad file and return data and metadata
    filename: Name of the file
    pixels: tupel of detector pixel size
    
    Returns:
    data: Raw detector data reshaped (e.g. 128x128px)
    meta: metadata collected from file'''
    
    with open(filename, 'rb') as data:
        data = np.fromfile(data, dtype=np.int32)[0:pixels[0]*pixels[1]]
        data = data.reshape(pixels)
    
  
    with open(filename,'r', errors='ignore') as fh:
        meta = {}
        startnum =1e9
        for num, curline in enumerate(fh, 1):
            # check if the current line starts with "#"
            if curline.startswith("### NICOS"):
                startnum = num
            if num >= startnum:
                if curline.startswith('#'):
                    pass
                else:
                    #print(curline.split())
                    items = curline.split()
                    try:
                        key = items[0]
                        if key.endswith(('_offset', '_precision', 'status', '_userlimits', '_alias')) or key.startswith(('Exp_', 'Sample', 'reseda_')):
                            # we don't need these
                            continue
                        try:
                            val = float(items[2])
                        except:
                            val = items[2]
                    except:
                        pass
                    meta[key] = val
                    
    return data, meta
