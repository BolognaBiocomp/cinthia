#!/usr/bin/python
import os
import sys
import subprocess
import string
import sets
import math
import re
import maxss # WARNING ! This must be in also in other places 
#sys.path.append('/home/savojard/LocPipeline/CINTHIA-MOD/hmm')
import hmm.tr_obj as tr_obj
import hmm.HMM_IO as HMM_IO
import hmm.algo_HMM as algo_HMM
import config as cinthiaconfig

def runCRF(model, inpf, outf, postf, decoding, window):

  """
  crfplabel = tempfile.NamedTemporaryFile()
    crfplabelname = crfplabel.name
    crfplabel.close
    crfpstate = tempfile.NamedTemporaryFile()
    crfpstatename = crfpstate.name
    crfpstate.close
    # Run prediction to get P(State)
    sp.check_call(['/home/savojard/BUSCA/tools/tppred3/tppred3/tools/biocrf',
                   '-test', '-w', '5', '-m',
                   crfmodel, '-a', '1', '-d',
                   'posterior-viterbi-max', '-q',
                   crfpstatename, '-o', crfpredname, crfdatname], stderr = open("/dev/null", 'w'), stdout = open("/dev/null", 'w'))
    state_prob = [float(line.split()[45]) for line in open(crfpstatename+"_0").readlines()]
  """
  subprocess.call([cinthiaconfig.BIOCRF, '-test',
                   '-m', model,
                   '-o', outf,
                   '-d', decoding,
                   '-q', postf,
                   '-w', str(window),
                   '-a', '1', inpf],
                   stdout=open('/dev/null', 'w'),
                   stderr=open('/dev/null', 'w'))
  observation = ''.join([x.split()[0] for x in open(outf).readlines()[:-1]])
  prediction  = ''.join([x.split()[1] for x in open(outf).readlines()[:-1]])
  probs       = [float(line.split()[1]) for line in open(postf+"_0").readlines()]
  return observation, prediction, probs

def runHMM(model, inpf):
  profile = open(inpf).readlines()
  labels  = [profile[i].split()[-1] for i in range(len(profile))]
  seq     = [map(float, profile[i].split()[:-1]) for i in range(len(profile))]
  obj     = tr_obj.TR_OBJ(seq, labels, name=profile)
  hmm     = HMM_IO.get_hmm(model)
  bestpath, bestval = algo_HMM.maxAcc_decoder(hmm, obj.seq)
  return ''.join(labels), ''.join(list(bestpath))

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

