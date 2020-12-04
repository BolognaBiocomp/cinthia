'''
Created on 31/gen/2011

@author: Castrense Savojardo
'''

import numpy

def logistic(x):
    # return math.tanh(x)
    return 1 / (1 + numpy.exp(-x))

class InvalidCheckpointFileError(Exception):
    def __init__(self):
        pass

def BlastCheckPointProfile(checkpointFile, newFormat = True):
    try:
        checkpointFile = open(checkpointFile).readlines()
    except IOError:
        print("Error while open/reading checkpoint file.")
        raise
    profile = None
    if newFormat:
        try:
            profile = _profileParseNew(checkpointFile)
        except:
            raise
    return profile

def _profileParseNew(checkpoint):
    headerSize = 3
    footerSize = -6
    shift = 22
    aaOrder = checkpoint[2].split()[:20]
    try:
        _check(checkpoint[1])
    except:
        raise

    profile = []

    for line in checkpoint[headerSize:footerSize]:
        line = line.split()
        pos = numpy.zeros(20)
        for j in range(20):
            val = float(line[j + shift]) / 100.0
            if not val == 0.0:
                pos[j] = val
        if numpy.sum(pos) == 0.0:
            aa = line[1]
            try:
                pos[aaOrder.index(aa)] = 0.001
            except ValueError:
                pass
        profile.append(pos)
    return numpy.array(profile)

def _check(line):
    import re
    if not re.search('Last position-specific scoring matrix computed', line):
        raise InvalidCheckpointFileError


def BlastCheckPointPSSM(checkpointFile, newFormat = True, transform = True):
    try:
        checkpointFile = open(checkpointFile).readlines()
    except IOError:
        print("Error while open/reading checkpoint file.")
        raise
    pssm = None
    if newFormat:
        try:
            pssm = _pssmParseNew(checkpointFile, transform)
        except:
            raise
    return pssm

def _pssmParseNew(checkpoint, transform):
    headerSize = 3
    footerSize = -6

    try:
        _check(checkpoint[1])
    except:
        raise

    pssm = []

    for line in checkpoint[headerSize:footerSize]:
        line = line.split()[2:22]
        pos = numpy.zeros(20)
        for j in range(20):
            if transform:
                pos[j] = logistic(float(line[j]))  # 1.0 / (1.0 + math.exp(-1 * float(line[j + shift])))
            else:
                pos[j] = float(line[j])
        pssm.append(pos)
    return numpy.array(pssm)
