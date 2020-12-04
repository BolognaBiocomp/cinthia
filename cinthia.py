#!/usr/bin/env python
import sys
import os
import logging
if 'CINTHIA_ROOT' in os.environ:
    sys.path.append(os.environ['CINTHIA_ROOT'])
else:
    logging.error("CINTHIA_ROOT environment varible is not set")
    logging.error("Please, set and export CINTHIA_ROOT to point to cinthia root folder")
    sys.exit(1)
import argparse
from modules import cinthiaalgo as cinthia
from modules import config as cfg
from modules import cpparser as bcp
from modules import workenv
from modules import utils
import re
import numpy
from Bio import SeqIO


def main():
    ## Parsing input arguments
    DESC="Cinthia: Combined predictor of signal peptide, transit peptide and transmembrane topology"
    parser = argparse.ArgumentParser(description=DESC)
    parser.add_argument("-f", "--fasta",
                        help = "The input protein sequence in FASTA format",
                        dest = "fasta", required = True)
    parser.add_argument("-p", "--pssm",
                        help = "The input PSSM/Profile file (PSIBLAST)",
                        dest = "pssm", required = True)
    parser.add_argument("-o", "--out",
                        help = "The output GFF3 prediction file",
                        dest = "output", required = True)
    parser.add_argument("-t", "--forcetopo",
                        help = "Force topology to contain at least one TM segment",
                        dest = "forcetopo", action="store_true")
    ns = parser.parse_args()

    we = workenv.TemporaryEnv()

    try:
        record = SeqIO.read(ns.fasta, "fasta")
    except:
        logging.exception("Error reading FASTA: file is not FASTA or more than one sequence is present")
        we.destroy()
        sys.exit(1)
    else:
        acc = record.id
        sequence = str(record.seq)
        try:
            utils.check_sequence_pssm_match(sequence, ns.pssm)
        except:
            logging.exception("Error in PSSM: sequence and provided PSSM do not match.")
            we.destroy()
            sys.exit(1)
        else:
            profile = bcp.BlastCheckPointProfile(ns.pssm)
            if ns.forcetopo:
                CRFprediction, CRFprobs = cinthia.runCRF(cfg.CRFFORCEDMODEL, profile, we)
                

SPprofile   = getNewTmpFile("sp.", ".prof")
SPresults   = getNewTmpFile("sp.", ".res")
TRprofile   = getNewTmpFile("tr.", ".prof")
TRresults   = getNewTmpFile("tr.", ".res")
HMMSprofile = getNewTmpFile("hmm.", ".prof")
CRFprofile  = getNewTmpFile("crf.", ".prof")
CRFraw      = getNewTmpFile("crf.raw.", ".res")
CRFpost     = getNewTmpFile("crf.posterior.", ".res.")
CINTHIAInp  = getNewTmpFile("cinthia.", ".input")
CINTHIAOut  = getNewTmpFile("cinthia.", ".output")

cleavageSite = 0
peptide = 'none'

writeCTerminus(unirefProf, HMMSprofile, cleavageSite, newline=False)
writeCTerminus(unirefProf, CRFprofile, cleavageSite)


if ns.forcetopo:
  _, CRFprediction, CRFprobs = cinthia.runCRF(CRFForcedModel, CRFprofile, CRFraw, CRFpost, decoding, 8, we)
  ofs = open(CINTHIAInp, 'w')
  ofs.write("# M1 M2 M3 M4\n")
  for i in range(len(CRFprediction)):
    ofs.write("\t".join([CRFprediction[i], CRFprediction[i], CRFprediction[i], CRFprediction[i]]) + '\n')
  ofs.close()
  DP,names=cinthia.readPreds(CINTHIAInp)
  tsymb=cinthia.tmsymbols()
  tmseg,topSeg,topSum,mVote=cinthia.topology0(DP,names,12,35,0.0,tsymb)
  pLen=len(DP[names[0]])

  cinthia.writeConsensus(tmseg,pLen,topSeg,topSum,mVote,tsymb,CINTHIAOut)
  topology = "".join([x.strip() for x in open(CINTHIAOut).readlines()]).replace("l", "i").replace("L", "o")
else:

  printDate("Predicting trans-membrane helices on mature protein sequence.")
  _, CRFprediction, CRFprobs = cinthia.runCRF(CRFmodel, CRFprofile, CRFraw, CRFpost, decoding, 8)
  _, HMMUprediction = cinthia.runHMM(HMMUmodel, HMMSprofile)
  _, HMMWprediction = cinthia.runHMM(HMMWmodel, HMMSprofile)

  if 'T' in CRFprediction and 'T' in HMMUprediction and 'T' in HMMWprediction:
    ofs = open(CINTHIAInp, 'w')
    ofs.write("# M1 M2 M3 M4\n")
    for i in range(len(CRFprediction)):
      ofs.write("\t".join([HMMUprediction[i], HMMWprediction[i], CRFprediction[i], CRFprediction[i]]) + '\n')
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

peptideStr = ""
sigStr = "SP=NO"

trStr = "TR=NO"

score = str(round(numpy.mean([CRFprobs[i] for i in range(len(CRFprobs)) if CRFprediction[i] == "T"]), 2))
ofs = open(ns.output, 'w')
ofs.write( "\t".join([acc, sigStr, trStr, peptideStr+topology, score]) + '\n' )

DestroyTemporaryEnvironment()

sys.exit(0)
