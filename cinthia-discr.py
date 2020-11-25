#!/usr/bin/env python

import sys
import os
import tempfile
from time import localtime, strftime
import argparse
import modules.cinthiaalgo as cinthia
import modules.config as cinthiaconfig
import shutil
import re
import numpy
from Bio import SeqIO


def printDate(msg):
  print "[%s] %s" % (strftime("%a, %d %b %Y %H:%M:%S", localtime()), msg)


def SetUpTemporaryEnvironment():
  tempfile.tempdir = os.path.abspath(tempfile.mkdtemp(prefix="job.tmpd.",
                                                      dir="."))
  printDate("Setting up job temporary enviroment [%s]" % tempfile.tempdir)


def DestroyTemporaryEnvironment():
  if not tempfile.tempdir == None:
    shutil.rmtree(tempfile.tempdir)
    printDate("Destroying job temporary enviroment [%s]" % tempfile.tempdir)

def getNewTmpFile(prefix, suffix):
  outTmpFile = tempfile.NamedTemporaryFile(mode   = 'write',
                                           prefix = prefix,
                                           suffix = suffix,
                                           delete = False)
  outTmpFileName = outTmpFile.name
  outTmpFile.close()
  return outTmpFileName

def writeNTerminus(profile, outfile, n, newline=True):
  try:
    assert(n < len(profile) and n > 0)
  except AssertionError:
    n = len(profile)
  of = open(outfile, 'w')
  for i in range(n):
    of.write(profile[i])
  if newline:
    of.write("\n")
  of.close()

def writeCTerminus(profile, outfile, n, newline=True):
  try:
    assert(n < len(profile) and n > 0)
  except AssertionError:
    n = 0
  of = open(outfile, 'w')
  for i in range(n, len(profile)):
    of.write(profile[i])
  if newline:
    of.write("\n")
  of.close()

def writeResFile(observation, prediction, outfile):
  ofs = open(outfile, 'w')
  ofs.write('# obs hmm_ma\n')
  for i in range(len(observation)):
    ofs.write("\t".join(observation[i], prediction[i]) + '\n')
  ofs.close()

## Models and parameters
CRFmodel=os.path.join(cinthiaconfig.CINTHIA_HOME, 'models','CRF.modnew')
HMMUmodel=os.path.join(cinthiaconfig.CINTHIA_HOME, 'models','HMM.modDR')
HMMWmodel=os.path.join(cinthiaconfig.CINTHIA_HOME, 'models','PES.modDR')
decoding='posterior-viterbi-sum'

## Parsing input arguments
DESC="Cinthia: Combined predictor of signal peptide, transit peptide and transmembrane topology"
parser = argparse.ArgumentParser(description=DESC)

parser.add_argument("-f", "--fasta",
                    help = "The input protein sequence in FASTA format",
                    dest = "fasta", required = True)
parser.add_argument("-s", "--swissprot",
                    help = "The input sequence profile computed using the SwissProt database",
                    dest = "swissprot", required = True)
parser.add_argument("-u", "--uniref",
                    help = "The input sequence profile computed using the Uniref90 database",
                    dest = "uniref", required = True)
parser.add_argument("-o", "--out",
                    help = "The output prediction file",
                    dest = "output", required = True)
parser.add_argument("-t", "--forcetopo",
                    help = "Force topology to contain at least one TM segment",
                    dest = "forcetopo", action="store_true")

ns = parser.parse_args()

acc = None
sequence = None
for record in SeqIO.parse(ns.fasta, 'fasta'):
  acc = record.id
  sequence = str(record.seq)
  # Only the first sequence is processed
  break

swissprotProf = open(ns.swissprot).readlines()
unirefProf = open(ns.uniref).readlines()

assert(len(swissprotProf) == len(sequence))
assert(len(unirefProf) == len(sequence))

SetUpTemporaryEnvironment()

CRFprofile  = getNewTmpFile("crf.", ".prof")
HMMSprofile = getNewTmpFile("hmm.", ".prof")
CRFraw      = getNewTmpFile("crf.raw.", ".res")
CRFpost     = getNewTmpFile("crf.posterior.", ".res.")
CINTHIAInp  = getNewTmpFile("cinthia.", ".input")
CINTHIAOut  = getNewTmpFile("cinthia.", ".output")

cleavageSite = 0
peptide = 'none'

writeCTerminus(unirefProf, HMMSprofile, cleavageSite, newline=False)
writeCTerminus(unirefProf, CRFprofile, cleavageSite)

printDate("Predicting trans-membrane helices on mature protein sequence.")
_, CRFprediction, CRFprobs = cinthia.runCRF(CRFmodel, CRFprofile, CRFraw, CRFpost, decoding, 8)
_, HMMUprediction = cinthia.runHMM(HMMUmodel, HMMSprofile)
_, HMMWprediction = cinthia.runHMM(HMMWmodel, HMMSprofile)

pred = None
if ('T' in CRFprediction) or ('T' in HMMUprediction) or ('T' in HMMWprediction):
  ofs = open(CINTHIAInp, 'w')
  ofs.write("# M1 M2 M3 M4\n")
  if 'T' in CRFprediction:
    pred = CRFprediction
  elif 'T' in HMMUprediction:
    pred = HMMUprediction
  else:
    pred = HMMWprediction
  score = str(round(numpy.mean([CRFprobs[i] for i in range(len(CRFprobs)) if pred[i] == "T"]), 2))

  for i in range(len(CRFprediction)):
    ofs.write("\t".join([pred[i], pred[i], pred[i], pred[i]]) + '\n')
  ofs.close()

  DP,names=cinthia.readPreds(CINTHIAInp)
  tsymb=cinthia.tmsymbols()
  tmseg,topSeg,topSum,mVote=cinthia.topology0(DP,names,12,35,0.0,tsymb)
  pLen=len(DP[names[0]])

  cinthia.writeConsensus(tmseg,pLen,topSeg,topSum,mVote,tsymb,CINTHIAOut)
  topology = "".join([x.strip() for x in open(CINTHIAOut).readlines()])
  if not re.match("-GLOBULAR", topology):
    topology = "".join([x.strip() for x in open(CINTHIAOut).readlines()]).replace("l", "i").replace("L", "o")
else:
  topology = "-GLOBULAR"
  score = str(round(numpy.mean([CRFprobs[i] for i in range(len(CRFprobs)) if CRFprediction[i] != "T"]), 2))

peptideStr = ""
sigStr = "SP=NO"
trStr = "TR=NO"
ofs = open(ns.output, 'w')
ofs.write("\t".join([acc, sigStr, trStr, peptideStr+topology, score]) + '\n')

DestroyTemporaryEnvironment()

sys.exit(0)

