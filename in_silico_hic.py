#!/usr/bin/env python

import re
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import math
import subprocess

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
    
    #start = 0; stop = len(seq)
    p = re.compile(rs_dict[enzyme])
    pos_list = [m.start() for m in p.finditer(str(seq))]
    bin_list = pos_list[::10]
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
    
    return frag_df, start_boundaries, stop_boundaries, bin_list

def ligate(df, starts, stops):
    '''Performs ligations of fragments within the same window based on
    coordinates supplied from start_boundaries and stop_boundaries lists
    
    Parameters
    ----------
    df: pandas data-frame containing coordinates and sequences of all digested
        fragments (returned from digest_genome)
    starts: a list containing all start coordinates of the windows
        within which to perform ligations (returned from digest_genome)
    stops: a list containing all stop coordinates of the windows
        within which to perform ligations (returned from digest_genome)
        
    Returns
    -------
    final_frags: a list containing sonicated, pulled-down and size-selected
        fragments (all fragments will be in the range 240-400bp and contain
        one restricion site)
    
    '''
    
    #p = re.compile(rs_dict[enzyme])
    p = re.compile('GATC')
    fa_seqs = {}
    final_frags = []
    for x, y in zip(starts, stops):
        seqs = df[(df['start']>=x) & (df['stop']<=y)]['seq']
        lig_frags = []
        for i, j in enumerate(seqs):
            start = i+1
            for k in seqs[start:]:
                lig_frags.append(j + k)
                lig_frags.append(k + j)
        for frag in lig_frags:
            # remove leading restriction site
            frag = frag[len('GATC'):]
            # sonicate
            sonic_frags = [frag[i:i+400] for i in range(0, len(frag), 400)]
            #pulldown
            junc_frags = [i for i in sonic_frags if bool(
                                                re.search('GATC', i))==True] 
            #size-select    
            for i in junc_frags:
                if len(i)>=240: final_frags.append(i)

    return final_frags

def sequence(seqs):
    '''Paired end sequences a list of supplied fragments. Outputs forward and
    reverse FASTA files
    
    Parameters
    ----------
    seqs: a list of fragment sequences to sequence
    
    Returns:
    --------
    for_seqs:
    rev_seqs:
    
    '''
    for_fasta = open('seqs_1.fa', 'w')
    rev_fasta = open('seqs_2.fa', 'w')
    for i, seq in enumerate(seqs):
        seq_name = 'seq.'+str(i)
        for_seq = seq[0:120]
        rev_seq = Seq(seq[-120:]).reverse_complement()
        for_fasta.write('>{}\n{}\n'.format(seq_name, for_seq))
        rev_fasta.write('>{}\n{}\n'.format(seq_name, rev_seq))
    
    return

def fill_matrix(bins):
    bound_coor = ['{}-{}'.format(pos, bins[i+1]) for i, pos in enumerate(bins[:-1])]
    hic_matrix = pd.DataFrame(columns=bound_coor[:500], index=bound_coor[:500])
    
    file1 = '/t1-data1/WTSA_Dev/jkerry/HiC/HiCUP/jon_scripts/Dev/seqs_1.sorted.bam'
    file2 = '/t1-data1/WTSA_Dev/jkerry/HiC/HiCUP/jon_scripts/Dev/seqs_2.sorted.bam'
    
    for i, coor in enumerate(hic_matrix.index.values):
        coor1 = 'chr22:'+coor
        in_first = subprocess.Popen('samtools view {} {} | cut -f 1'.format(file1, coor1),
                                        shell=True,
                                        stdout=subprocess.PIPE).communicate()[0]
        first_list = in_first.decode('utf-8').split('\n')[:-1]
        for j in hic_matrix.columns[i:]:
            coor2 = 'chr22:'+j
            in_second = subprocess.Popen('samtools view {} {} | cut -f 1'.format(file2, coor2),
                                        shell=True,
                                        stdout=subprocess.PIPE).communicate()[0]
            second_list = in_second.decode('utf-8').split('\n')[:-1]
            
            hic_matrix.loc[coor, j] = len([i for i in first_list if i in second_list])
            
    return hic_matrix
    
    
    