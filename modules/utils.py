import logging
from time import localtime, strftime
from . import cpparser as bcp

def printDate(msg):
  print("[%s] %s" % (strftime("%a, %d %b %Y %H:%M:%S", localtime()), msg))

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

def check_sequence_pssm_match(sequence, psiblast_pssm):
  try:
    pssm_mat = bcp.BlastCheckPointPSSM(psiblast_pssm)
  except:
    logging.error("Failed reading/parsing PSSM file")
    raise
  else:
    try:
      assert(len(sequence) == pssm_mat.shape[0])
    except:
      logging.error("Sequence and PSSM have different lengths")
      raise
  return True
