'''maxss

This module defines an object that computes *maximal subsequences*
in a sense made precise in "Calcolo della sottosequenza Ottima":maxss.dvi.

Objects:

  MaxSubSeq -- gives information about maximal subsequence of a given
  _reward sequence_
'''

# Ivan Fri, 20 Jul 2001 11:58:44 +0200
# rinominato in MaxSubSeq
class MaxSubSeq :
    '''MaxSubSeq Object -- Compute and store information about optimal
    subsequence.

    Provide list interface to access score values, does not support slicing.
    '''

    def __init__(self, R, lmin, lmax) :
        '''\
        Compute optimal susquences and store results in newly created
        object.

        Arguments:
        o R -- a reward sequence of len N, \'R[i]\' is reward for 1
               in position i.
        o lmin -- minimum length of a subsequence of 1.
        o lmax -- maximum length of a subsequence of 0.'''

        N = len(R)
        assert((0 < lmin <= lmax) and (N >= lmin))
        # ntries is the number of subseqs to be considered for each index
        ntries = lmax - lmin + 1
        # r_sum[n] is sum_{j=0}^{n-1} R[j]
        r_sum = [0]*(N+1)
        bottom = -1
        for n in range(N) :
            r_sum[n+1] = r_sum[n] + R[n]
            bottom = bottom - abs(R[n])
        bottom = 2 * bottom
        # max_k is the maximum k s.t. S(N, k) > bottom
        max_k = int((N+1) / (lmin+1))
        # init score matrix: s[k][n] = S(n, k+1)
        self.s = [None] * (max_k+1)
        sk1 = self.s[0] = [0] * (N+1)
        # DIRTY: we use s[0][-1] = s[0][N] = 0, by Python list semantics

        # INV: sk1[n] = S(n, k-1) for each n,
        # INV: self.s[h,n] = S(n, h) for each n, for all h < k - 1
        # max_model contains the column with the maximum score
        self.max_model=0
        for k in range(1, max_k + 1) :
            sk = [bottom] * (N+1)
            lastinv = 0
            v = [bottom] * ntries

            #INV: v[lastinv ... lastinv + lmax % nchoices] = V(n-1,k)
            #INV: sk[i] = S(k,i) for i < n
            for n in range((lmin+1)*k-1, N+1) :
                v[lastinv] = sk1[n-lmin-1] - r_sum[n-lmin]
                lastinv = (lastinv + 1) % ntries
                sk[n] = max(sk[n-1], r_sum[n] + max(v))
            self.s[k] = sk
            sk1 = sk
            if(self.s[k][N-1] > self.s[self.max_model][N-1]):
               self.max_model=k

        # Store results in the newly created object
        self.len = N
        self.lmin = lmin
        self.lmax = lmax
        self.max_k = max_k
        self.r_sum = r_sum

    def max_seq(self, k = None) :
        '''\
        Returns a sequence L of len 2*k+1 such that the sequence
        0^\'L[0] + 1^\'L[1]\' + ... + 0^\'L[2*k]\' has score S(n,k).

        Arguments:
        o k -- number of 1 subsequences (default to \'len(self) - 1\').

        Raises IndexError for invalid indexes k.
        '''

        if (k == None) : k = self.max_model
        n = self.len

        seq = [None] * (2*k+1)
        while (k > 0) :
            maxscore = self.s[k][n]
            # Count number of zeros
            m = 0
            while (self.s[k][n-m-1] == maxscore) : m = m + 1
            seq[2*k] = m + 1
            # Count number of ones
            n = n - m
            m = self.lmin
            maxscore = maxscore - self.r_sum[n]
            sk1 = self.s[k-1]
            while (abs(maxscore - sk1[n-m-1] + self.r_sum[n-m]) > 10e-9) : m = m + 1
            seq[2*k-1] = m
            k = k - 1
            n = n - m -1
        seq[0] = n + 1
        seq[-1] = seq[-1] - 1
        return seq

    def __len__(self) :
        '''Return the number of k (starting from 0) for which Seq(N, k) is
        non empty.'''
        return self.max_k + 1

    def __getitem__(self, k) :
        '''Returns Score(N,k).

        Raises IndexError for invalid indexes k.
        '''
        return self.s[k][self.len]
