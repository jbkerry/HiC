#!/usr/bin/env python

import re
import numpy as np
import pandas as pd
from Bio import SeqIO

rs_dict = {'DpnII': 'GATC',
           'NlaIII': 'CATG',
           'HindIII': 'AAGCTT'}

def digest_genome(fa, chromosome, enzyme='DpnII', window=10, step=5):
    '''in silico digests the provided fasta file using the recognition
    sequence of the specifed enzyme
    
    Parameters
    ----------
    fa: reference genome/chromosome fasta file
    chromosome: the number/letter of the chromsome to digest
    enzyme: restriction enzyme to use in the in-silico digest, default=DpnII
    window: int, the size (in kb) of the window within which to perform
        ligations, default = 10kb
    step: int, the size (in kb) to move the window for each ligation set,
        default = 5kb
    
    Returns
    -------
    frag_df: a dictionary containing the fragement coordinates as the key
        and the fragment sequence as the value
    start_boundaries: a list containing all start coordinates of the windows
        within which to perform ligations
    stop_boundaries: a list containing all stop coordinates of the windows
        within which to perform ligations
    
    '''
    
    chr_name = 'chr'+str(chromosome)
    seq_dict = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    seq = seq_dict[chr_name].seq.upper()
    
    start = 0; stop = len(seq)
    p = re.compile(rs_dict[enzyme])
    pos_list = [m.start() for m in p.finditer(str(seq[start:stop]))]
    pos_list = np.array(pos_list)
    
    win_size = 1000*window
    step_size = 1000*step
    start_boundaries = []
    stop_boundaries = []
    this_start = pos_list[0]
    final_stop = pos_list[-1]
    while this_start <= final_stop:
        start_boundary = pos_list[pos_list>=this_start][0]
        if start_boundary <= final_stop-win_size:
            start_boundaries.append(start_boundary)
            this_stop = start_boundary+win_size
            stop_boundary = pos_list[pos_list>=this_stop][0]
            stop_boundaries.append(stop_boundary)
        this_start=start_boundary+step_size
    
    fragments = []
    for x, y in enumerate(pos_list[:-1]):
        z = pos_list[x+1]
        fragments.append([chr_name, y, z, str(seq[y:z])])
    frag_df = pd.DataFrame(fragments, columns = ['chr', 'start',
                                                 'stop', 'seq'])
    
    return frag_df, start_boundaries, stop_boundaries

def ligate(df, starts, stops):
    '''Performs self-ligations of fragments within the same window based on
    coordinates supplied from start_boundaries and stop_boundaries lists
    '''
    #p = re.compile(rs_dict[enzyme])
    p = re.compile('GATC')
    fa_seqs = {}
    for x, y in zip(starts, stops):
        seqs = df[(df['start']>=x) & (df['stop']<=y)]['seq']
        lig_frags = []
        for i, j in enumerate(seqs):
            start = i+1
            for k in seqs[start:]:
                lig_frags.append(j + k)
                lig_frags.append(k + j)
        for frag in lig_frags:
            if len(frag)>=200:
                junction = [m.start() for m in p.finditer(frag)][1]
                read_start = junction-100 if junction-100>=0 else 0
                read_stop = junction+100 if junction<=len(frag) else len(frag)
                key = '{}-{}'.format(read_start, read_stop)
                fa_seqs[key] = frag[read_start:read_stop]
    return fa_seqs
    
    