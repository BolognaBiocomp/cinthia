import os

CINTHIA_ROOT = os.environ['CINTHIA_ROOT']

BIOCRF = os.path.join(CINTHIA_ROOT, 'tools' 'biocrf-static')


## Models and parameters
SPMODEL=os.path.join(CINTHIA_ROOT, 'models','CRF-SP.mod')
TRMODEL=os.path.join(CINTHIA_ROOT, 'models','CRF-TR.mod')
CRFMODEL=os.path.join(CINTHIA_ROOT, 'models','CRF.modnew')
CRFFORCEDMODEL=os.path.join(CINTHIA_ROOT, 'models','CRF_force-tm.modnew')
HMMUMODEL=os.path.join(CINTHIA_ROOT, 'models','HMM.modDR')
HMMWMODEL=os.path.join(CINTHIA_ROOT, 'models','PES.modDR')

CRF_DECONDING='posterior-viterbi-sum'
CRF_WINDOW = 8
