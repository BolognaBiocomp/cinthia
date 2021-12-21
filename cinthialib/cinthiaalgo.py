#!/usr/bin/python
import os
import sys
import subprocess
import string
import math
import re
from . import maxss # WARNING ! This must be in also in other places
#sys.path.append('/home/savojard/LocPipeline/CINTHIA-MOD/hmm')
from .hmm import tr_obj as tr_obj
from .hmm import HMM_IO as HMM_IO
from .hmm import algo_HMM as algo_HMM
from . import config as cfg

def runCRF_multi(model, profiles, we, num_threads = 1):
  crfdat = we.createFile("crf.", ".dat")
  cdofs=open(crfdat,'w')
  for profile in profiles:
      for i in range(profile.shape[0]):
          for j in range(profile.shape[1]):
              cdofs.write("%f " % profile[i][j])
          cdofs.write("0.0 l\n")
      cdofs.write("\n")
  cdofs.close()
  crfpred = we.createFile("crf.", ".pred")
  crfplabel = we.createFile("crf.",".plabel")

  subprocess.call([cfg.BIOCRF, '-test',
                   '-m', model,
                   '-o', crfpred,
                   '-d', cfg.CRF_DECONDING,
                   '-q', crfplabel,
                   '-w', str(cfg.CRF_WINDOW),
                   '-a', str(num_threads), crfdat],
                   stdout=open('/dev/null', 'w'),
                   stderr=open('/dev/null', 'w'))
  probs = []
  for i in range(len(profiles)):
      probs.append([float(line.split()[1]) for line in open(crfplabel+"_%d" % i)])
  prediction = []
  p = []
  for line in open(crfpred):
      line = line.split()
      if len(line) > 0:
          p.append(line[1])
      else:
          p = "".join(p)
          prediction.append(p)
          p = []
  return prediction, probs

def runCRF(model, profile, we):
  crfdat = we.createFile("crf.", ".dat")
  cdofs=open(crfdat,'w')
  for i in range(profile.shape[0]):
      for j in range(profile.shape[1]):
          cdofs.write("%f " % profile[i][j])
      cdofs.write("0.0 l\n")
  cdofs.write("\n")
  cdofs.close()
  crfpred = we.createFile("crf.", ".pred")
  crfplabel = we.createFile("crf.",".plabel")

  subprocess.call([cfg.BIOCRF, '-test',
                   '-m', model,
                   '-o', crfpred,
                   '-d', cfg.CRF_DECONDING,
                   '-q', crfplabel,
                   '-w', str(cfg.CRF_WINDOW),
                   '-a', '1', crfdat],
                   stdout=open('/dev/null', 'w'),
                   stderr=open('/dev/null', 'w'))
  prediction  = ''.join([x.split()[1] for x in open(crfpred).readlines()[:-1]])
  probs       = [float(line.split()[1]) for line in open(crfplabel+"_0").readlines()]
  return prediction, probs

def runHMM(model, profile, we):
  #profile = open(inpf).readlines()
  hmmdat = we.createFile("hmm.", ".dat")
  hdofs=open(hmmdat,'w')
  for i in range(profile.shape[0]):
      for j in range(profile.shape[1]):
          hdofs.write("%f " % profile[i][j])
      hdofs.write("0.0 l\n")
  hdofs.close()
  dp = open(hmmdat).readlines()
  labels  = ["l" for i in range(len(dp))]
  seq     = [[float(x) for x in dp[i].split()[:-1]] for i in range(len(dp))]
  obj     = tr_obj.TR_OBJ(seq, labels, name=dp)
  hmm     = HMM_IO.get_hmm(model)
  bestpath, bestval = algo_HMM.maxAcc_decoder(hmm, obj.seq)
  return ''.join(list(bestpath))

class tmsymbols:
      def __init__(self,tm='T',inside='l',outside='L',globular='g'):
          '''__init__(self,tm='T',inside='l',outside='L',globular='g') '''
          self.tm=tm
          self.inside=inside
          self.outside=outside
          self.globular=globular

def readPreds(fname):
    ''' readPreds(fname) '''
    try:
        lines=open(fname,'r').readlines()
    except:
        sys.stderr.write('cannot open file '+fname+'\n')
        sys.exit()
    v=lines[0].split() # assuming first line # seq obs m1 m2 ...mn
    P={}
    names=[]
    for e in v[1:]:  # '# m1 m2 ...mn'
        P[e]=''
        names.append(e)
    for line in lines[1:]:
        v=line.split()
        for i in range(len(names)):
            P[names[i]]+=v[i]
    return P,names

def topology0(DP,names,minl,maxl,th,tSymb):
    ''' - topology(DP,minl,maxl,th,tSymb):
          computes the average helix score and uses MaxSubseq to get the tm-segments
          evaluates the best agrement on the loops to assign topology
    '''
    score=[]
    lp=len(DP[names[0]])
    loops=[0.0]*lp
    # conpute the TM-Helix score
    majorVote=0
    for n in names:
        if DP[n][0] == tSymb.inside:
            majorVote+=1
        elif DP[n][0] == tSymb.outside:
            majorVote-=1
    ds={tSymb.tm:1}
    for i in range(lp):
        stmp=0.0
        for n in names:
            stmp+=ds.get(DP[n][i],-1)
            if DP[n][i] == tSymb.inside:
               loops[i]+=1
            elif DP[n][i] == tSymb.outside:
               loops[i]-=1
        score.append(stmp-th)
    #print score
    # optimize the segment position with maxsubseq
    m=maxss.MaxSubSeq(score,minl,maxl)
    seglist=m.max_seq()
    #print seglist
#    pt=[-1]*lp
    hseg=[]
    ind=0
    isLoop=0
    topScore=0.0
    topSegScore=0
    evenloop=1
    for i in range(len(seglist)):
        isLoop = isLoop + 1
        isLoop = isLoop%2
        if isLoop : # compute loop score
           tmpsum=evenloop*sum(loops[ind:ind+seglist[i]])
           topScore+=tmpsum
           if tmpsum > 0:
              topSegScore+=1
           elif tmpsum < 0:
              topSegScore-=1
           evenloop*=-1 # change phase
        else:
           hseg.append((ind,ind+seglist[i]))
        ind=ind+seglist[i]
#        for j in range(seglist[i]):
#            pt[ind]=isLoop
#            ind=ind +1
    if hseg == []:
        topSegScore=topScore=None
    return hseg,topSegScore,topScore,majorVote

def writeConsensus(tmSeg,pLen,topSeg,topSum,mVote,tSymb,outf):
    '''printConsensus(tmSeg,pLen,topSeg,topSum,mVote,tSymb)'''
    ofs = open(outf, 'w')
    if tmSeg == []:
        #for i in range(pLen):
            #print tSymb.globular
        ofs.write("-GLOBULAR\n")
        return
    topLab=[tSymb.inside,tSymb.outside]
    if topSeg > 0:
        it=0
    elif topSeg < 0:
        it=1
    else:
        if topSum > 0:
            it=0
        elif topSum < 0:
            it=1
        else:
            if mVote > 0:
               it=0
            elif mVote < 0:
               it=1
            else:
               it=0
               sys.stderr.write("ACTUNG set in but undefined topology\n")
    bOld=0
    for b,e in tmSeg:
        for i in range(bOld,b):
            ofs.write( str(topLab[it]) + '\n')
        it=(it+1)%2
        for i in range(b,e):
            ofs.write( str(tSymb.tm) + '\n')
        bOld=e
    for i in range(bOld,pLen):
        ofs.write( str(topLab[it]) + '\n')
    return 0
