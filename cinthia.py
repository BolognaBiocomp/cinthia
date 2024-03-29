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

def run_multifasta(ns):
    we = workenv.TemporaryEnv()
    data_cache = utils.get_data_cache(ns.cache_dir)
    i = 0
    ofsout = open(ns.outf, 'w')
    protein_jsons = []
    for record in SeqIO.parse(ns.fasta, 'fasta'):
        acc = record.id
        sequence = str(record.seq)
        prefix = "seq%d" % i
        fastaSeq  = we.createFile(prefix+".", ".fasta")
        SeqIO.write([record], fastaSeq, 'fasta')
        pssm = blast.runPsiBlast(prefix, ns.dbfile, fastaSeq, we, data_cache=data_cache)
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
        if ns.outfmt == "json":
            acc_json = utils.get_json_output(acc, sequence, topology, CRFprobs)
            #json.dump([acc_json], ofsout, indent=5)
            protein_jsons.append(acc_json)
        else:
            utils.write_gff_output(acc, sequence, ofsout, topology, CRFprobs)
        i = i + 1
    if ns.outfmt == "json":
        json.dump(protein_jsons, ofsout, indent=5)
    ofsout.close()
    we.destroy()
    sys.exit(0)

def run_pssm(ns):
    we = workenv.TemporaryEnv()
    record = SeqIO.read(ns.fasta, "fasta")
    acc = record.id
    sequence = str(record.seq)

    profile = bcp.BlastCheckPointProfile(ns.pssm)
    profile = utils.rearrange_profile(profile, cfg.BLASTALPH, cfg.HSSPALPH)
    topology = ""
    ofsout = open(ns.outf, 'w')
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
        if ns.outfmt == "json":
            acc_json = utils.get_json_output(acc, sequence, topology, CRFprobs)
            json.dump([acc_json], ofsout, indent=5)
        else:
            utils.write_gff_output(acc, sequence, ofsout, topology, CRFprobs)
    ofsout.close()
    we.destroy()
    sys.exit(0)

def main():
    ## Parsing input arguments
    DESC="Cinthia: Predictor of helical transmembrane topology"
    parser = argparse.ArgumentParser(description=DESC)
    subparsers   = parser.add_subparsers(title = "subcommands", description = "valid subcommands", required = True)
    multifasta  = subparsers.add_parser("multi-fasta", help = "Multi-FASTA input module", description = "DeepMito: Multi-FASTA input module.")
    pssm  = subparsers.add_parser("pssm", help = "PSSM input module (one sequence at a time)", description = "DeepMito: PSSM input module.")
    multifasta.add_argument("-f", "--fasta", help = "The input multi-FASTA file name", dest = "fasta", required = True)
    multifasta.add_argument("-d", "--dbfile", help = "The PSIBLAST DB file", dest = "dbfile", required= True)
    multifasta.add_argument("-o", "--outf", help = "The output file", dest = "outf", required = True)
    multifasta.add_argument("-m", "--outfmt", help = "The output format: json or gff3 (default)", choices=['json', 'gff3'], required = False, default = "gff3")
    multifasta.add_argument("-t", "--forcetopo", help = "Force topology to contain at least one TM segment", dest = "forcetopo", action="store_true")
    multifasta.add_argument("-c", "--cache-dir", help="Cache dir for alignemnts", dest="cache_dir", required=False, default=None)
    multifasta.set_defaults(func=run_multifasta)
    pssm.add_argument("-f", "--fasta", help = "The input FASTA file name (one sequence)", dest = "fasta", required = True)
    pssm.add_argument("-p", "--pssm", help = "The PSIBLAST PSSM file", dest = "pssm", required= True)
    pssm.add_argument("-o", "--outf", help = "The output file", dest = "outf", required = True)
    pssm.add_argument("-m", "--outfmt", help = "The output format: json or gff3 (default)", choices=['json', 'gff3'], required = False, default = "gff3")
    pssm.add_argument("-t", "--forcetopo", help = "Force topology to contain at least one TM segment", dest = "forcetopo", action="store_true")
    pssm.set_defaults(func=run_pssm)
    if len(sys.argv) == 1:
      parser.print_help()
    else:
      ns = parser.parse_args()
      ns.func(ns)

if __name__ == "__main__":
    main()
