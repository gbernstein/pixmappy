#!/usr/bin/env python

# Example code for using DELVE WCS solutions for DECAM data.

import sys
import numpy as np
from astropy.table import Table
import argparse
import pixmappy as pm


def run(dmc, expnum, catfile):
    tab = Table.read(catfile)
    ccdnums = np.unique(tab['CCDNUM'])
    x = tab['XWIN_IMAGE']
    y = tab['YWIN_IMAGE']
    try:
        c = tab['GI']
    except:
        c = np.ones_like(x)*pm.REF_COLOR

    # New array for ra,dec
    radec = np.zeros((len(x),2),dtype=float)
    # Also record the size of trap for each point
    trap = np.zeros_like(x)
    
    for ccdnum in ccdnums:
        use = tab['CCDNUM']==ccdnum
        # Acquire and use the WCS
        wcs = dmc.getDelveWCS(expnum,ccdnum)
        r,d = wcs.toSky(x[use],y[use],c[use])
        radec[use] = np.stack([r,d],axis=-1)
        
        # Now record the trap size for each source
        trapName = dmc.trapMapFor(expnum,ccdnum)
        trapMap = dmc.getMap(trapName)
        trap[use] = trapMap.mas(x[use],y[use])

    # Save the information as new columns in the table
    tab['radec'] = radec
    tab['trap'] = trap
    return tab

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=\
        'Apply astrometric mapping to a DELVE catalog')
    parser.add_argument("expnum", help='exposure number', type=int)
    parser.add_argument('-i','--infile', help='path to input catalog',type=str)
    parser.add_argument('-o','--outfile', help='path to save output catalog',type=str)
    args = parser.parse_args()

    if args.infile is None or args.outfile is None:
        print('ERROR: -i and -o file specs required')
        sys.exit(1)

    dmc = pm.DelveMaps()

    tab = run(dmc,args.expnum, args.infile)
    tab.write(args.outfile)
    sys.exit(0)
