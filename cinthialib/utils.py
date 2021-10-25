import logging
from time import localtime, strftime
import re
import numpy as np
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

def one_hot_encoding(sequence, alph="VLIMFWYGAPSTCHRKQEND"):
    profile = numpy.zeros((len(sequence), 20))
    for (i, aa) in enumerate(sequence):
        try:
            j = alph.index(aa)
            profile[i,j] = 1.0
        except:
            pass
    return profile

def rearrange_profile(profile, in_alph, out_alph):
    new_profile = np.zeros(profile.shape)
    for (i, a) in enumerate(out_alph):
        new_profile[:,i] = profile[:,in_alph.index(a)]
    return new_profile

def get_data_cache(cache_dir):
    import os
    from . import datacache
    ret = None
    if cache_dir is not None:
        if os.path.isdir(cache_dir):
            ret = datacache.DataCache(cache_dir)
    return ret

def write_gff_output(acc, sequence, output_file, topology, scores):
    if topology != "":
        for mo in re.finditer("T+", topology):
            s = mo.start()
            e = mo.end()
            score = np.mean(scores[s:e])
            print(acc, "CINTHIA/ENSEMBLE3.0", "Transmembrane", s+1, e, round(score,2), ".", ".", "Note=Helical;evidence=ECO:0000256")


def get_json_output(acc, sequence, topology, scores):
    acc_json = {'accession': acc, 'features': []}
    acc_json['sequence'] = {
                              "length": len(sequence),
                              "sequence": sequence
                           }
    if topology != "":
        for mo in re.finditer("T+", topology):
            s = mo.start()
            e = mo.end()
            score = np.mean(scores[s:e])
            acc_json['features'].append({
                  "type": "TRANSMEM",
                  "category": "TOPOLOGY",
                  "description": "Helical",
                  "begin": s+1,
                  "end": e,
                  "score": round(float(score),2),
                  "evidences": [
                    {
                      "code": "ECO:0000256",
                      "source": {
                        "name": "SAM",
                        "id": "CINTHIA/ENSEMBLE3.0",
                        "url": "https://busca.biocomp.unibo.it"
                      }
                    }
                  ]
            })
    return acc_json
