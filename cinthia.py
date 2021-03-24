#!/usr/bin/env python
import sys
import os
import logging
import json
if 'CINTHIA_ROOT' in os.environ:
    sys.path.append(os.environ['CINTHIA_ROOT'])
else:
    logging.error("CINTHIA_ROOT environment varible is not set")
    logging.error("Please, set and export CINTHIA_ROOT to point to cinthia root folder")
    sys.exit(1)
import argparse
from cinthialib import cinthiaalgo as cinthia
from cinthialib import config as cfg
from cinthialib import cpparser as bcp
from cinthialib import workenv
from cinthialib import utils
from cinthialib import blast
import re
import numpy
from Bio import SeqIO

def run_json(ns):
    we = workenv.TemporaryEnv()
    ifs = open(ns.i_json)
    input_json = json.load(ifs)
    ifs.close()
    i = 0
    protein_jsons = []
    #for record in SeqIO.parse(ns.fasta, 'fasta'):
    for i_json in input_json:
        acc = i_json['accession']
        if not ns.mature:
            sequence, cleavage = utils.cut_peptide(i_json)
        else:
            sequence = i_json['sequence']['sequence']
            cleavage = 0
        prefix = "seq%d" % i
        fastaSeq  = we.createFile(prefix+".", ".fasta")
        #SeqIO.write([record], fastaSeq, 'fasta')
        fsofs=open(fastaSeq,'w')
        #SeqIO.write([fasta], fsofs, 'fasta')
        print(">%s" % acc, file=fsofs)
        print(sequence, file=fsofs)
        fsofs.close()
        pssm = blast.runPsiBlast(prefix, cfg.BLASTDB, fastaSeq, we)
        profile = bcp.BlastCheckPointProfile(pssm)
        profile = utils.rearrange_profile(profile, cfg.BLASTALPH, cfg.HSSPALPH)
        topology = ""
        if ns.forcetopo:
            CRFprediction, CRFprobs = cinthia.runCRF(cfg.CRFFORCEDMODEL, profile, we)
            cinthia_input_tmp_file = we.createFile("cinthia.", ".input.dat")
            cinthia_output_tmp_file = we.createFile("cinthia.", ".output.dat")
            ofs = open(cinthia_input_tmp_file, 'w')
            ofs.write("# M1 M2 M3 M4\n")
            for i in range(len(CRFprediction)):
                ofs.write("\t".join([CRFprediction[i], CRFprediction[i], CRFprediction[i], CRFprediction[i]]) + '\n')
            ofs.close()
            DP,names=cinthia.readPreds(cinthia_input_tmp_file)
            tsymb=cinthia.tmsymbols()
            tmseg,topSeg,topSum,mVote=cinthia.topology0(DP,names,12,35,0.0,tsymb)
            pLen=len(DP[names[0]])
            cinthia.writeConsensus(tmseg,pLen,topSeg,topSum,mVote,tsymb,cinthia_output_tmp_file)
            topology = "".join([x.strip() for x in open(cinthia_output_tmp_file).readlines()]).replace("l", "i").replace("L", "o")
        else:
            CRFprediction, CRFprobs = cinthia.runCRF(cfg.CRFMODEL, profile, we)
            HMMUprediction = cinthia.runHMM(cfg.HMMUMODEL, profile, we)
            HMMWprediction = cinthia.runHMM(cfg.HMMWMODEL, profile, we)
            if 'T' in CRFprediction and 'T' in HMMUprediction and 'T' in HMMWprediction:
                cinthia_input_tmp_file = we.createFile("cinthia.", ".input.dat")
                cinthia_output_tmp_file = we.createFile("cinthia.", ".output.dat")
                ofs = open(cinthia_input_tmp_file, 'w')
                ofs.write("# M1 M2 M3 M4\n")
                for i in range(len(CRFprediction)):
                    ofs.write("\t".join([HMMUprediction[i], HMMWprediction[i], CRFprediction[i], CRFprediction[i]]) + '\n')
                ofs.close()
                DP,names=cinthia.readPreds(cinthia_input_tmp_file)
                tsymb=cinthia.tmsymbols()
                tmseg,topSeg,topSum,mVote=cinthia.topology0(DP,names,12,35,0.0,tsymb)
                pLen=len(DP[names[0]])

                cinthia.writeConsensus(tmseg,pLen,topSeg,topSum,mVote,tsymb,cinthia_output_tmp_file)
                topology = "".join([x.strip() for x in open(cinthia_output_tmp_file).readlines()])
                if not re.match("-GLOBULAR", topology):
                    topology = "".join([x.strip() for x in open(cinthia_output_tmp_file).readlines()]).replace("l", "i").replace("L", "o")
                else:
                    topology = ""
            else:
                topology = ""
        if cleavage > 0:
            topology = "P" * cleavage + topology
            CRFprobs = [0.0] * cleavage + CRFprobs
        acc_json = utils.get_json_output(i_json, topology, CRFprobs)
        protein_jsons.append(acc_json)
        i = i + 1
    ofsout = open(ns.outf, 'w')
    json.dump(protein_jsons, ofsout, indent=5)
    ofsout.close()
    we.destroy()
    sys.exit(0)

def main():
    ## Parsing input arguments
    DESC="Cinthia: Predictor of helical transmembrane topology"
    parser = argparse.ArgumentParser(description=DESC)
    parser.add_argument("-i", "--i-json", help = "The input JSON file name", dest = "i_json", required = True)
    parser.add_argument("-o", "--outf", help = "The output file", dest = "outf", required = True)
    parser.add_argument("-t", "--forcetopo", help = "Force topology to contain at least one TM segment", dest = "forcetopo", action="store_true")
    parser.add_argument("-m", "--mature", help = "Sequences are mature. Do not try to cut SIGNAL or TRANSIT peptides before prediction", dest= "mature", action="store_true")
    ns = parser.parse_args()
    run_json(ns)

if __name__ == "__main__":
    main()
