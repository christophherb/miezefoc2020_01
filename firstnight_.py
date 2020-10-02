#angles sth 
import numpy as np
sth_vals = np.linspace(-0.2,-1.8,81)

dtx_vals = [k for k in range(10,260,30)]+[400,600]
for dtx_val in dtx_vals:
    maw(dtx,dtx_val)
    for sth_val in sth_vals:
        maw(sth,sth_val)
        count(50)
Shutter.close()