# Copyright (c) 2016 by Mike Jarvis and the other collaborators on GitHub at
# https://github.com/rmjarvis/Piff  All rights reserved.
#
# Piff is free software: Redistribution and use in source and binary forms
# with or without modification, are permitted provided that the following
# conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the disclaimer given in the accompanying LICENSE
#    file.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the disclaimer given in the documentation
#    and/or other materials provided with the distribution.

"""
.. module:: decaminfo
"""

import numpy as np

# Translation from CCDNUM to DETPOS for DECam:
ccdnum2detpos = {1:'S29',  2:'S30',  3:'S31',  4:'S25',  5:'S26',  6:'S27',
                 7:'S28',  8:'S20',  9:'S21', 10:'S22', 11:'S23', 12:'S24',
                13:'S14', 14:'S15', 15:'S16', 16:'S17', 17:'S18', 18:'S19',
                19:'S8',  20:'S9',  21:'S10', 22:'S11', 23:'S12', 24:'S13',
                25:'S1',  26:'S2',  27:'S3',  28:'S4',  29:'S5',  30:'S6',
                31:'S7',  32:'N1',  33:'N2',  34:'N3',  35:'N4',  36:'N5',
                37:'N6',  38:'N7',  39:'N8',  40:'N9', 41:'N10', 42:'N11',
                43:'N12', 44:'N13', 45:'N14', 46:'N15', 47:'N16', 48:'N17',
                49:'N18', 50:'N19', 51:'N20', 52:'N21', 53:'N22', 54:'N23',
                55:'N24', 56:'N25', 57:'N26', 58:'N27', 59:'N28', 60:'N29',
                61:'N30', 62:'N31'}
#...and inverse
detpos2ccdnum = {v:k for k,v in ccdnum2detpos.items()}

# Dictionary of valid pixels, with 0-based indexing and
# "end" pixel not included.
xyBounds = {'xstart':25,'xend':2023, 'ystart':15, 'yend':4081}

uvBounds =  {'N1': (-1.0811, -0.782681, -0.157306, -0.00750506),
             'N2': (-0.771362, -0.472493, -0.157385, -0.00749848), 
             'N3': (-0.461205, -0.161464, -0.157448, -0.00749265), 
             'N4': (-0.150127, 0.149894, -0.15747, -0.00749085), 
             'N5': (0.161033, 0.460796, -0.157638, -0.0074294), 
             'N6': (0.472171, 0.771045, -0.157286, -0.00740563), 
             'N7': (0.782398, 1.08083, -0.157141, -0.0074798), 
             'N8': (-0.92615, -0.627492, -0.321782, -0.172004), 
             'N9': (-0.616455, -0.317043, -0.322077, -0.172189), 
             'N10': (-0.305679, -0.00571999, -0.322071, -0.17217), 
             'N11': (0.00565427, 0.305554, -0.322243, -0.172254), 
             'N12': (0.31684, 0.616183, -0.322099, -0.172063), 
             'N13': (0.627264, 0.925858, -0.321792, -0.171887), 
             'N14': (-0.926057, -0.62726, -0.485961, -0.336213), 
             'N15': (-0.616498, -0.317089, -0.486444, -0.336606), 
             'N16': (-0.30558, -0.00578257, -0.486753, -0.336864), 
             'N17': (0.00532179, 0.305123, -0.486814, -0.33687), 
             'N18': (0.316662, 0.616018, -0.486495, -0.336537), 
             'N19': (0.62708, 0.92578, -0.485992, -0.336061), 
             'N20': (-0.770814, -0.471826, -0.650617, -0.500679), 
             'N21': (-0.460777, -0.161224, -0.650817, -0.501097), 
             'N22': (-0.149847, 0.149886, -0.650816, -0.501308), 
             'N23': (0.161001, 0.460566, -0.650946, -0.501263), 
             'N24': (0.47163, 0.770632, -0.650495, -0.500592), 
             'N25': (-0.615548, -0.316352, -0.814774, -0.665052), 
             'N26': (-0.305399, -0.00591217, -0.814862, -0.665489), 
             'N27': (0.00550714, 0.304979, -0.815022, -0.665418), 
             'N28': (0.316126, 0.615276, -0.814707, -0.664908), 
             'N29': (-0.46018, -0.16101, -0.97887, -0.829315), 
             'N31': (0.160884, 0.460147, -0.978775, -0.829426), 
             'S1': (-1.08096, -0.782554, 0.00715956, 0.15689), 
             'S2': (-0.7713, -0.47242, 0.0074194, 0.157269), 
             'S3': (-0.4611, -0.161377, 0.00723009, 0.157192), 
             'S4': (-0.149836, 0.150222, 0.00737069, 0.157441), 
             'S5': (0.161297, 0.461031, 0.0072399, 0.1572), 
             'S6': (0.472537, 0.771441, 0.00728934, 0.157137), 
             'S7': (0.782516, 1.08097, 0.00742809, 0.15709), 
             'S8': (-0.92583, -0.627259, 0.171786, 0.32173), 
             'S9': (-0.616329, -0.31694, 0.171889, 0.321823), 
             'S10': (-0.305695, -0.00579187, 0.172216, 0.322179), 
             'S11': (0.00556739, 0.305472, 0.172237, 0.322278), 
             'S12': (0.316973, 0.61631, 0.172015, 0.322057), 
             'S13': (0.627389, 0.925972, 0.171749, 0.321672), 
             'S14': (-0.925847, -0.627123, 0.335898, 0.48578), 
             'S15': (-0.616201, -0.316839, 0.336498, 0.486438), 
             'S16': (-0.305558, -0.00574858, 0.336904, 0.486749), 
             'S17': (0.00557115, 0.305423, 0.33675, 0.486491), 
             'S18': (0.316635, 0.615931, 0.33649, 0.486573), 
             'S19': (0.627207, 0.925969, 0.336118, 0.485923), 
             'S20': (-0.770675, -0.471718, 0.500411, 0.65042), 
             'S21': (-0.46072, -0.161101, 0.501198, 0.650786), 
             'S22': (-0.149915, 0.14982, 0.501334, 0.650856), 
             'S23': (0.160973, 0.460482, 0.501075, 0.650896), 
             'S24': (0.47167, 0.770647, 0.50045, 0.650441), 
             'S25': (-0.615564, -0.316325, 0.66501, 0.814674), 
             'S26': (-0.30512, -0.0056517, 0.665531, 0.81505), 
             'S27': (0.00560886, 0.305082, 0.665509, 0.815022), 
             'S28': (0.316158, 0.615391, 0.665058, 0.814732), 
             'S29': (-0.46021, -0.160988, 0.829248, 0.978699), 
             'S30': (-0.150043, 0.149464, 0.829007, 0.978648), 
             'S31': (0.160898, 0.460111, 0.82932, 0.978804) }

''' Class for determining the calibration epoch that should be used for
an exposure taken at a given MJD.
'''
from astropy.time import Time
import numpy as np

def mjdOfEpoch(epoch):
    # Return mjd of date specified by 8-character epoch
    return Time(epoch[:4]+'-'+epoch[4:6]+'-'+epoch[6:8], 
                format='fits',scale='utc').mjd
class EpochFinder:
    '''Function class which returns the star flat epoch nearest to specifed input MJD
    that does not have an intervening camera event.  Returns '00000000' if no star flats
    occur in the same interval between events.
    '''
    # Epochs of star flats and of camera "events" when calibration changes.
    sfEpochs = ['20121120','20121223','20130221','20130829','20131115','20140118',
                '20140807','20141105','20150204','20150926','20160209',
                '20160223','20160816','20161117','20170111','20170214',
                '20170411','20170814','20170906','20171129','20180103',
                '20180327','20180829','20181123','20181218','20190116']
        # Skipping  20181025, bad registration
        # Also note only ugri are usable 20180103, no zY.
    warmups = ['20121230','20130512','20130722','20131015','20140512',
               '20141201',
               # remove, see below: '20150625','20150725',
               '20150809',
               # remove '20150825',
               '20160219','20161013',
               # remove '20161214',
               '20161226','20170714',
               # remove '20170803', # This was changing r filter positions
               '20170903','20171103','20171215','20180318', #rizY filters moved
               '20180619', #Y filter moved 20180718, g on 20180814
               '20181118']
    # Missing a SF set between 20121226 and 20121230; -> remove former cooldown
    # 20150625,0725,0809,0825;  -> remove first 2, last one?
    # 20161214 and 1226;  -> omit first one
    # 20170714 and 0803;  -> omit latter
    # 20171215 and 20180314 (missing zY only in former) and 20180318 ->drop 0314 as last is
    # optics work; will need to kludge a 20180103 solution for zY from 20171129,
    # which is preferable to going to 20180327 star flat because filter/shutter service
    # just before the latter (which also is missing Y band)
    
    cooldowns=['20151126']  # Omitting '20121226','20180314'
    nogood = '00000000'
    def __init__(self):
        # Place the epochs at ~midday Chile time of their stated date.
        self.sfMjds = np.array([mjdOfEpoch(e) for e in self.sfEpochs]) + 0.7
        self.eventMjds = np.array([mjdOfEpoch(e) for e in self.warmups + self.cooldowns]) + 0.7
        self.eventMjds.sort() 
        return
    def __call__(self, mjd):
        if mjd is None:
            return self.nogood
        # which events are before, after our mjd?
        before = mjd >= self.eventMjds
        # Mark which star flat MJDs are in same interval between events
        if not np.any(before):
            # Our mjd is before any events
            same = self.sfMjds < self.eventMjds[0]
        elif np.all(before):
            # Our mjd is after all events
            same = self.sfMjds >= self.eventMjds[-1]
        else:
            # Our mjd is between two events, get index of preceding one
            precede = np.where(before)[0][-1]
            same = np.logical_and(self.sfMjds>=self.eventMjds[precede],
                                  self.sfMjds <self.eventMjds[precede+1])
        
        if not same.any():
            # No star flats in the same event interval.
            return self.nogood
        sameIndices = np.where(same)[0]
        closest = np.argmin(np.abs(self.sfMjds[same]-mjd))
        return self.sfEpochs[sameIndices[closest]]
