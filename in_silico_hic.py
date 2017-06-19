#!/usr/bin/env python

import re
from Bio import SeqIO

rs_dict = {'DpnII': 'GATC',
           'NlaIII': 'CATG',
           'HindIII': 'AAGCTT'}

def digest_genome(fa, chromosome, enzyme='DpnII'):
    '''in silico digests the provided fasta file using the recognition
    sequence of the specifed enzyme
    
    Parameters
    ----------
    fa = reference genome/chromosome fasta file
    chromosome = the number/letter of the chromsome to digest
    enzyme = restriction enzyme to use in the in-silico digest, default=DpnII
    
    Returns
    -------
    fragments = a dictionary containing the fragement coordinates as the key
        and the fragment sequence as the value
    
    '''
    
    chr_name = 'chr'+str(chromosome)
    seq_dict = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    seq = seq_dict[chr_name].seq.upper()
    
    start = 0; stop = len(seq)
    p = re.compile(rs_dict[enzyme])
    pos_list = []
    for m in p.finditer(str(seq[start:stop])):
        pos_list.append(m.start())
    pos_list = np.array(pos_list)
    
    this_start = pos_list[0]
    start_boundaries = []
    stop_boundaries = []
    while this_start<= pos_list[-1]:
        start_boundary = pos_list[pos_list>=this_start][0]
        if start_boundary<=pos_list[-1]-10000:
            start_boundaries.append(start_boundary)
            this_stop = start_boundary+10000
            stop_boundary = pos_list[pos_list>=this_stop][0]
            stop_boundaries.append(stop_boundary)
        this_start=start_boundary+5000
    
    if len(start_boundaries)!=len(stop_boundaries):
        start_boundaries = start_boundaries[:-1]
    seq_list = []
    cut_size = len(rs_dict[enzyme])
    fragments = {}
    for i, j in enumerate(pos_list[:-1]):
        coor = '{}:{}-{}'.format(chr_name, j, pos_list[i+1])
        fragments[coor] = str(seq[j:pos_list[i+1]])
        seq_list.append(str(seq[j:pos_list[i+1]]))
    #seq_list = np.array(seq_list)
    
    return seq_list, fragments, start_boundaries, stop_boundaries