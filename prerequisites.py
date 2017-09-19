"""
This file contains functions for probabilistic association rules
mining for prerequisite relationships discovery

Step 1: To find all the possible pairs of KCs
Step 2: Each KC pair as a pattern, calculate the support probability distribution
        function fx(k) for each pattern
        (Dynamic-Programming Algorithm, see (Sun et al. KDD 2010))
Step 3: To compute association rule probability
        Given minsup and minconf, calculate P(X=>Y)
        (This algorithm can also be found in (Sun et al. KDD 2010))
        Given minprob, if P(X=>Y) >= minprob, then return prAR

Copyright reserved by Yang CHEN (yangchen312@gmail.com)
"""

def getPatterns(df):

    """get patterns from dataset, here, patterns are pairs of KCs
    Input: dataframe
    Output: patterns list -- palist
    """

    palist = []
    KCs = df.columns
    for kc1 in KCs:
        for kc2 in KCs:
            if kc1 != kc2 and sorted((kc1, kc2)) not in palist:
                palist.append(sorted((kc1, kc2)))
    return palist


def DP(df, pat):

    """Dynamic-Programming algorithm to calculate the pmf fx(k) for a pattern
    Input: dataframe, a pattern
    Output:support pmf -- fx(k)
    """

    fx = [1] + [0 for i in range(len(df))]
    fx2 = [0 for i in range(len(df) + 1)]
    for i in df.index:
        pix = 1
        for e in pat:
            pix = pix * df.ix[i][e]
        fx2[0] = (1-pix) * fx[0]
        for k in range(1, len(fx)):
            fx2[k] = pix * fx[k-1] + (1 - pix) * fx[k]
        fx = [a for a in fx2]
    return fx

def DP2(df, pat):

    """Dynamic-Programming algorithm to calculate the pmf fx(k) for a pattern
    Input: dataframe, a pattern
    Output:support pmf -- fx(k)
    """

    fx = [1] + [0 for i in range(len(df))]
    fx2 = [0 for i in range(len(df) + 1)]
    for i in df.index:
        pix = df.ix[i][pat[0]] * (1 - df.ix[i][pat[1]])
        fx2[0] = (1-pix) * fx[0]
        for k in range(1, len(fx)):
            fx2[k] = pix * fx[k-1] + (1 - pix) * fx[k]
        fx = [a for a in fx2]
    return fx


def ARP(fxy, fxny, minsup, minconf):

    """Computing the probability of an association rule
    Input: support pmfs fxy & fxny, thresholds minsup & minconf
    Output: probability of rule x=>y
    """

    prAR = 0
    for i in range(minsup, len(fxy)):
        j = 0
        prCum = 0
        while j < len(fxny):
            if j > ((1 - minconf) / minconf) * i:
                break
            else:
                prCum = prCum + fxny[j]
                j += 1
        prAR = prAR + prCum * fxy[i]
    return prAR


def ARM(df, minsup, minconf, minprob):

    """Mining probabilistic association rules
    Input: dataframe, three thresholds minsup, minconf, minprob
    Output: pairs (x, y), x=>y
    """

    ars = []
    pairs = getPatterns(df)
    revpairs = [list(reversed(pair)) for pair in pairs]
    pairs = pairs + revpairs
    for xy in pairs:
        fxy = DP(df, xy)
        fxny = DP2(df, xy)
        arp = ARP(fxy, fxny, minsup, minconf)
        if arp >= minprob:
            ars.append(xy)
    return ars

def ARO(df, minsup, minconf):

    """Rank association rules according to their probabilities
    Input: data, two thresholds minsup, minconf
    Output: pairs (x, y), x=>y and their probabilities
    """

    import bisect
    ars = []
    arps = []
    pairs = getPatterns(df)
    revpairs = [list(reversed(pair)) for pair in pairs]
    pairs = pairs+revpairs
    for xy in pairs:
        fxy = DP(df, xy)
        fxny = DP2(df, xy)
        arp = ARP(fxy, fxny, minsup, minconf)
        bisect.insort(arps, arp)
        ars.insert(arps.index(arp), xy)
    return list(reversed(ars)), list(reversed(arps))



if __name__ == '__main__':
    #Given threshold values
    minsupct = 150   #minimum support count
    minconf = 0.8     #minimum confidence threshold
    minprob = 0.9     #minimum probability threshold

    import pandas as pd

    df = pd.read_csv('data_demo.csv')
    ars = ARM(df, minsupct, minconf, minprob)
    print "The rules in the form of Sj=1=>Si=1 are found."

    df1 = 1 - df
    ars1 = ARM(df1, minsupct, minconf, minprob)
    print "The rules in the form of Si=0=>Sj=0 are found."

    print "\nDiscovered prerequisite links:"
    for ar in ars:
        if list(reversed(ar)) in ars1:
            print ar[1] + " is a prerequisite of " + ar[0]


    #Only given minsupct & minconf values, without minprob value,
    #the rules are ranked according to their probabilities
    print "\nThe ranking of rules [x, y] x=1=>y=1:"
    ars_rank, arps = ARO(df, minsupct, minconf)
    for i in range(len(ars_rank)):
        print str(ars_rank[i]) + ' probability: ' + str(arps[i])

    print "\nThe ranking of rules [x, y] x=0=>y=0:"
    ars1_rank, arps1 = ARO(df1, minsupct, minconf)
    for i in range(len(ars1_rank)):
        print str(ars1_rank[i]) + ' probability: ' + str(arps1[i])
