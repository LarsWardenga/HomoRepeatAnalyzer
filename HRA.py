#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 12:56:10 2018

@author: lars

"""

from __future__ import division
#from Bio.Align.Applications import MuscleCommandline
#from Bio.Emboss.Applications import NeedleCommandline
from collections import defaultdict
from os import listdir
from os.path import isfile, join
from scipy import stats

import copy
import filecmp
import itertools
import math
import matplotlib.pyplot as plt; plt.close('all')
import numpy as np
import os
import pandas as pd
import pickle
import pybam
import re
import resource
import seaborn as sns
import shutil
import subprocess
import sys
import time








t0 = time.time()

#order: A T C G




###############################################################################
#start of house keeping, helper, and clean-up functions, also variables
###############################################################################

# this function returns the peak memory usage of the main modules if called
#@profile
def get_mem(name):
    #process = psutil.Process(os.getpid())
    MaxMem1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    #MaxMem2 = process.memory_info()[0]
    #return name + ': 1) ' + str(MaxMem1) + 'kb | ' + '2) ' + str(MaxMem2/1000) + 'kb'
    return name + ' peak memory usage [kb]: ' + str(MaxMem1)

#taken from https://github.com/bosswissam/pysize
def get_size(obj, seen=None):
        """Recursively finds size of objects"""
        size = sys.getsizeof(obj)
        if seen is None:
            seen = set()
        obj_id = id(obj)
        if obj_id in seen:
            return 0
        seen.add(obj_id)
        if isinstance(obj, dict):
            size += sum([get_size(v, seen) for v in obj.values()])
            size += sum([get_size(k, seen) for k in obj.keys()])
        elif hasattr(obj, '__dict__'):
            size += get_size(obj.__dict__, seen)
        elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
            size += sum([get_size(i, seen) for i in obj])
        return size

# this function takes a dictionary as the ones returned by the main modules and
# merges all other dictionaries of its kind into one inclusive dictionary
#@profile    
def merge_dicts(dict_list):
    first_dict = dict_list[0].copy()
    for any_further_dict in xrange(1, len(dict_list), 1):
        first_dict.update(dict_list[any_further_dict])
    return first_dict

# this function checks if the "pre decimal" digit of an ENST (the last figure 
# of the ENST name, since the figure after the period-mark denotes the version 
# number) and if it matches it returns True 
#@profile
def predecimal_ones_check(ones_pre_dec, target):
    #print('target: ' + target + ' | ' + 'REF ones_pre_dec: ' + ones_pre_dec)
    if re.match(r"ENS.*" + ones_pre_dec + "\..+", target):
        return True

def predecimal_tens_and_ones_check(tens_pre_dec, ones_pre_dec, target):
    #print('target: ' + target + ' | ' + 'REF ones_pre_dec: ' + ones_pre_dec)
    if re.match(r"ENS.*" + tens_pre_dec + ones_pre_dec + "\..+", target):
        return True

#OUTDATED: this function is a back-up function that catches sequences that were
# missing for the reference upon start of compare_subject_vs_reference
#@profile
def fetch_seq_from_consensus_pickle(working_directory, ChrOI, bed_name, bam_name, unique_ENST):
    #proc_ID = 'P_' + ChrOI + '_' + unique_ENST.split('.')[0][-1]
    csv_path = working_directory + 'raw_consensus_matrices/' + unique_ENST + '.csv'
    if os.path.exists(working_directory + bed + bam + '_HRA_mapping_analyzer_merged_dict_pickle.txt') == True:
        pickle_opened_for_reading = open(working_directory + bed + bam + '_HRA_mapping_analyzer_merged_dict_pickle.txt')
        HRA_mapping_analyzer_dict = pickle.load(pickle_opened_for_reading)
        consensus_seq = HRA_mapping_analyzer_dict[unique_ENST + '|' + ChrOI + '|consensus_seq']
        pickle_opened_for_reading.close()
        return consensus_seq
    elif os.path.exists(csv_path) == True:
        csv_opened_for_reading = open(csv_path, 'r')
        consensus_seq = ''.join([line.split(',')[0] for line in csv_opened_for_reading])
        return consensus_seq

#TODO: Check if this table is correct
AAS_table = { 
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
            } 

def greater_func(first_repeat_length, first_repeat_start, first_repeat_stop, first_repeat_seq, second_repeat_length, second_repeat_start, second_repeat_stop, second_repeat_seq):
    if first_repeat_length > second_repeat_length:
        return True
    else:
        return False
    
def lesser_func(first_repeat_length, first_repeat_start, first_repeat_stop, first_repeat_seq, second_repeat_length, second_repeat_start, second_repeat_stop, second_repeat_seq):
    if first_repeat_length < second_repeat_length:
        return True
    else:
        return False
    
def equal_func(first_repeat_length, first_repeat_start, first_repeat_stop, first_repeat_seq, second_repeat_length, second_repeat_start, second_repeat_stop, second_repeat_seq):
    if first_repeat_length == second_repeat_length:
        return True
    else:
        return False
        
def any_seq_changes_func(first_repeat_length, first_repeat_start, first_repeat_stop, first_repeat_seq, second_repeat_length, second_repeat_start, second_repeat_stop, second_repeat_seq):
    if first_repeat_seq != second_repeat_seq:
        return True
    else:
        return False
        
def greater_length_and_any_seq_changes_func(first_repeat_length, first_repeat_start, first_repeat_stop, first_repeat_seq, second_repeat_length, second_repeat_start, second_repeat_stop, second_repeat_seq):
    if (first_repeat_length > second_repeat_length) and (first_repeat_seq != second_repeat_seq):
        return True
    else:
        return False
    
def lesser_length_and_any_seq_changes_func(first_repeat_length, first_repeat_start, first_repeat_stop, first_repeat_seq, second_repeat_length, second_repeat_start, second_repeat_stop, second_repeat_seq):
    if (first_repeat_length < second_repeat_length) and (first_repeat_seq != second_repeat_seq):
        return True
    else:
        return False
    
def equal_length_and_any_seq_changes_func(first_repeat_length, first_repeat_start, first_repeat_stop, first_repeat_seq, second_repeat_length, second_repeat_start, second_repeat_stop, second_repeat_seq):
    if (first_repeat_length == second_repeat_length) and (first_repeat_seq != second_repeat_seq):
        return True
    else:
        return False

#returns True if new_element matches an element in a given list, disregarding both their [1] placed subcontents
#TODO: verify if needed, else remove this function 
def check_if_in_list_regardless_of_second_element(new_element, list_to_be_checked):
    new_element_minus_ENST = new_element
    list_to_be_checked_minus_ENST = list_to_be_checked
    del new_element_minus_ENST[1]
    for old_element_minus_ENST in list_to_be_checked_minus_ENST:
        del old_element_minus_ENST[1]
        if new_element_minus_ENST == old_element_minus_ENST:
            return True
    return False

def check_if_list_elements_not_in_string(string, checklist):
    for character in checklist:
        if character in string:
            return False
    return True

###############################################################################
#end of house keeping, helper, and clean-up functions, also variables
###############################################################################


###############################################################################
#start of parsing functions
###############################################################################


def simple_fasta_parser(fasta_file_object_opened_for_reading):    
    header = fasta_file_object_opened_for_reading.readline().rstrip('\n')
    seq = ''.join((line.rstrip('\n') for line in fasta_file_object_opened_for_reading))
    fasta_file_object_opened_for_reading.close()
    return [(header, seq)]


# this function parses the fasta files as reference (ensembl format should be 
#used) - it is called by HRA_ref_analyzer()
#@profile
def fasta_parser(directory, bed, fasta_file_object_opened_for_reading, ChrOI, ones_pre_dec, HRA_bed_parser_dict, GRCh37_cds = False, GRCh37_dna = False, specific_enst = False):
    return_list = []
    if GRCh37_cds == True:
        if specific_enst == False:
            line = fasta_file_object_opened_for_reading.readline()
            while line != '':
                if line.startswith('>'):
                    #example of a header line ">ENST00000361390.2 mt_genbank_import:known chromosome:GRCh37:MT:3307:4262:1 gene:ENSG00000198888.2 [...]"
                    #as explained here: http://may2012.archive.ensembl.org/info/website/tutorials/module3_feb2009_ensembl.pdf
                    header = line.rstrip('\n')
                    header_elements = header.split(' ')
                    enst = header_elements[0].lstrip('>')
                    chromosome_element = header_elements[2].split(':')
                    chromosome = chromosome_element[2]
                    ensg_element = header_elements[3].split(':')
                    ensg = ensg_element[1]
                    enst_start = chromosome_element[3]
                    enst_stop = chromosome_element[4]
                    if chromosome_element[5] == '1':
                        strand = '+'
                    elif chromosome_element[5] == '-1':
                        strand = '-'
                    else:
                        strand = 'NA'
                        raise Exception('strand of' + str(enst) + 'couldnt be determined')
                    seq = ''
                    line = fasta_file_object_opened_for_reading.readline()
                    
                    while (not line.startswith('>') and line != ''):
                        seq += line.rstrip('\n')
                        line = fasta_file_object_opened_for_reading.readline()
                    return_list.append((chromosome, ensg, enst, seq, strand, enst_start, enst_stop))
                    
            return return_list

        else:            
            line = fasta_file_object_opened_for_reading.readline()
            while line != '':
                if line.startswith('>'):
                    #example of a header line ">ENST00000361390.2 mt_genbank_import:known chromosome:GRCh37:MT:3307:4262:1 gene:ENSG00000198888.2 [...]"
                    #as explained here: http://may2012.archive.ensembl.org/info/website/tutorials/module3_feb2009_ensembl.pdf
                    header = line.rstrip('\n')
                    header_elements = header.split(' ')
                    enst = header_elements[0].lstrip('>')
                    chromosome_element = header_elements[2].split(':')
                    chromosome = chromosome_element[2]
                    ensg_element = header_elements[3].split(':')
                    ensg = ensg_element[1]
                    enst_start = chromosome_element[3]
                    enst_stop = chromosome_element[4]
                    if chromosome_element[5] == '1':
                        strand = '+'
                    elif chromosome_element[5] == '-1':
                        strand = '-'
                    else:
                        strand = 'NA'
                        print('strand of' + str(enst) + 'couldnt be determined')
                    seq = ''
                    line = fasta_file_object_opened_for_reading.readline()
                    
                    while (not line.startswith('>') and line != ''):
                        seq += line.rstrip('\n')
                        line = fasta_file_object_opened_for_reading.readline()
                    if specific_enst == enst:
                        return [(chromosome, ensg, enst, seq, strand, enst_start, enst_stop)]
    
    elif GRCh37_dna == True:
        if specific_enst == False:
            header = fasta_file_object_opened_for_reading.readline().rstrip('\n')
            seq = ''.join((line.rstrip('\n') for line in fasta_file_object_opened_for_reading))
            fasta_file_object_opened_for_reading.close()
            return [(header, seq)]

        else:    
            header = fasta_file_object_opened_for_reading.readline().rstrip('\n')
            seq = ''.join((line.rstrip('\n') for line in fasta_file_object_opened_for_reading))
            proc_ID = 'P_' + ChrOI + '_' + ones_pre_dec
            if os.path.exists(directory + 'bed_pickles/' + bed + '_' + proc_ID + '_ENST_split.txt') == True:
                HRA_bed_parser_dict = pickle.load(open(directory + 'bed_pickles/' + bed + '_' + proc_ID + '_ENST_split.txt'))
            Chromosome, Start, Stop, strand = HRA_bed_parser_dict[specific_enst + '|bed_data']
            ENST_start = int(Start)
            ENST_stop = int(Stop) + 1
            #fasta_file_object_opened_for_reading.close()
            return [(header, seq[ENST_start:ENST_stop])]

# this function parses the gff3 files to get information on where coding exons
#(CDS) are located, and it returns a dict with the coordinates as values to 
# ENST-keys
#@profile
def GFF3_parser(GFF3_path):
    #sub_feature_phase allows 0,1,2 and is given in relation to the feature's coding strand
    #ensembl_start_phase and ensembl_stop_phase are absolute in the sense that they refer to the reference genomes plus-strand regardless of the feature's orientation
    #the manual reads: Column 8: "phase"
    #                            For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon. In other words, a phase of "0" indicates that the next codon begins at the first base of the region described by the current line, a phase of "1" indicates that the next codon begins at the second base of this region, and a phase of "2" indicates that the codon begins at the third base of this region. This is NOT to be confused with the frame, which is simply start modulo 3. If there is no phase, put a "." (a period) in this field.
    #                            For forward strand features, phase is counted from the start field. For reverse strand features, phase is counted from the end field.
    #                            The phase is required for all CDS features.    
    #common errors: for calculating the length of a feature the start and end of the feature are substracted, but a +1 has to be added to account for the actual number of nucleotides; the phase in CDS entries in gff3 files is the number of bases one needs to skip to reach the first codon starting base in the feature - a better use for phase would be to annotate the first bases codon position ... this can be derived from the gff3 "phase" by calculating: (3 - gff3_phase) % 3
    #The following formulas were tested extensively and work for both strands:
    #next_exons_first_bases_phase = ((((feature_stop - feature_start + 1) % 3) + real_phase) % 3)
    #last_bases_phase = ((((feature_stop - feature_start + 1) % 3) + real_phase + 2) % 3)
    #as shown here: ("P" = phase, len = real length, "X" = intermediate, "B-phase" = base-of-interess's phase)
            
    #P	     len % 3 = X | (X + P) % 3 = B-phase
    #012012012 0  
    #agtgggttt B    9 % 3 = 0 | (0 + 0) % 3 = 0    y
    #120120120 1
    #agtgggttt B    9 % 3 = 0 | (0 + 1) % 3 = 1    y
    #201201201 2
    #agtgggttt B    9 % 3 = 0 | (0 + 2) % 3 = 2    y
            
    #2012012012 0
    #xagtgggttt B    10 % 3 = 1 | (1 + 2) % 3 = 0    y
    #0120120120 1
    #xagtgggttt B    10 % 3 = 1 | (1 + 0) % 3 = 1    y
    #1201201201 2
    #xagtgggttt B    10 % 3 = 1 | (1 + 1) % 3 = 2    y
            
    #12012012012 0
    #xxagtgggttt B    11 % 3 = 2 | (2 + 1) % 3 = 0    y
    #20120120120 1
    #xxagtgggttt B    11 % 3 = 2 | (2 + 2) % 3 = 1    y
    #01201201201 2
    #xxagtgggttt B    11 % 3 = 2 | (2 + 0) % 3 = 2    y
    
    #P	     len % 3 = X | (X + P + 2) % 3 = B-phase
    #012012012 0  
    #agtgggttB t    9 % 3 = 0 | (0 + 0 + 2) % 3 = 2    y
    #120120120 1
    #agtgggttB t    9 % 3 = 0 | (0 + 1 + 2) % 3 = 0    y
    #201201201 2
    #agtgggttB t    9 % 3 = 0 | (0 + 2 + 2) % 3 = 1    y
    
    #2012012012 0
    #xagtgggttB t    10 % 3 = 1 | (1 + 2 + 2) % 3 = 2    y
    #0120120120 1
    #xagtgggttB t    10 % 3 = 1 | (1 + 0 + 2) % 3 = 0    y
    #1201201201 2
    #xagtgggttB t    10 % 3 = 1 | (1 + 1 + 2) % 3 = 1    y
    
    #12012012012 0
    #xxagtgggttB t    11 % 3 = 2 | (2 + 1 + 2) % 3 = 2    y
    #20120120120 1
    #xxagtgggttB t    11 % 3 = 2 | (2 + 2 + 2) % 3 = 0    y
    #01201201201 2
    #xxagtgggttB t    11 % 3 = 2 | (2 + 0 + 2) % 3 = 1    y 

    #          P    len % 3 = X | (X + P) % 3 =      B-phase
    #0 210210210  
    #B agtgggttt      9 % 3 = 0 | (0 + 0) % 3 = 0    y
    #2 102102102  
    #B agtgggttt      9 % 3 = 0 | (0 + 2) % 3 = 2    y
    #1 021021021  
    #B agtgggttt      9 % 3 = 0 | (0 + 1) % 3 = 1    y
    
    #0 2102102102  
    #B xagtgggttt      10 % 3 = 1 | (1 + 2) % 3 = 0    y
    #1 0210210210    
    #B xagtgggttt      10 % 3 = 1 | (1 + 0) % 3 = 1    y  
    #2 1021021021    
    #B xagtgggttt      10 % 3 = 1 | (1 + 1) % 3 = 2    y  
    
    #2 10210210210    
    #B xxagtgggttt      11 % 3 = 2 | (2 + 0) % 3 = 2    y  
    #0 21021021021    
    #B xxagtgggttt      11 % 3 = 2 | (2 + 1) % 3 = 0    y  
    #1 02102102102    
    #B xxagtgggttt      11 % 3 = 2 | (2 + 2) % 3 = 1    y  
    
    #          P    len % 3 = X | (X + P + 2) % 3 =      B-phase
    #0 210210210  
    #a Bgtgggttt      9 % 3 = 0 | (0 + 0 + 2) % 3 = 2    x
    #2 102102102  
    #a Bgtgggttt      9 % 3 = 0 | (0 + 2 + 2) % 3 = 1    x
    #1 021021021  
    #a Bgtgggttt      9 % 3 = 0 | (0 + 1 + 2) % 3 = 0    x
    
    #0 2102102102  
    #a Bagtgggttt      10 % 3 = 1 | (1 + 2 + 2) % 3 = 2    y
    #1 0210210210    
    #a Bagtgggttt      10 % 3 = 1 | (1 + 0 + 2) % 3 = 0    y  
    #2 1021021021    
    #a Bagtgggttt      10 % 3 = 1 | (1 + 1 + 2) % 3 = 1    y  
    
    #2 10210210210    
    #a Bxagtgggttt      11 % 3 = 2 | (2 + 0 + 2) % 3 = 1    y  
    #0 21021021021    
    #a Bxagtgggttt      11 % 3 = 2 | (2 + 1 + 2) % 3 = 2    y  
    #1 02102102102    
    #a Bxagtgggttt      11 % 3 = 2 | (2 + 2 + 2) % 3 = 0    y
    
    GFF3_opened = open(GFF3_path.rstrip('\r'), 'r')
    exon_coordinate_dict = defaultdict(list)
    line = GFF3_opened.readline()
    if re.match('##gff-version\s*3', line):
        while line != '':
            if 'ID=transcript:ENST' in line:
                
                columns = line.rstrip('\n').split('\t')
                if len(line.rstrip('\n').split('\t')) == 9:
                    #feature_seqid = columns[0]
                    #feature_source = columns[1]
                    #feature_type = columns[2]
                    #feature_start = columns[3]
                    #feature_stop = columns[4]
                    #feature_score = columns[5]
                    #feature_strand = columns[6]
                    #feature_phase = columns[7]
                    feature_attributes = columns[-1]
                    
                    attributes = feature_attributes.split(';')
                    general_line_ENST = attributes[0].split(':')[1]
                    version_number_line_ENST = attributes[-1].split('=')[1]
                    line_ENST = general_line_ENST + '.' + version_number_line_ENST
                    
                    line = GFF3_opened.readline()
                    while (not 'ID=transcript:ENST' in line) and (not line.startswith('###')) and (line != ''):
                        
                        if '\tCDS\t' in line:
                            sub_columns = line.rstrip('\n').split('\t')
                            #sub_feature_seqid = sub_columns[0]
                            #sub_feature_source = sub_columns[1]
                            sub_feature_type = sub_columns[2]
                            CDS_range_start = sub_columns[3]
                            CDS_range_stop = sub_columns[4]
                            #sub_feature_score = sub_columns[5]
                            #sub_feature_strand = sub_columns[6]
                            sub_feature_phase = sub_columns[7]
                            sub_feature_attributes = sub_columns[8]
                            sub_attributes = sub_feature_attributes.split(';')
                            sub_feature_parent = sub_attributes[1].lstrip('Parent=transcript:')
                                                        
                            if sub_feature_type == 'CDS' and (general_line_ENST == sub_feature_parent):
                                #confusingly enough, out of the information in ensembl's gff3 "version" only the exon option contains phase start AND end in its attributes, not the CDS entries though which contain the phase start as column[7]...
                                
                                #omitted because I decidedc against including Exons, since we have no use for UTRs
                                #exon_coordinate_dict[line_ENST].append((int(CDS_range_start), int(CDS_range_stop), sub_feature_type, sub_feature_phase, ensembl_start_phase, ensembl_stop_phase))
                                exon_coordinate_dict[line_ENST].append((int(CDS_range_start), int(CDS_range_stop), ((3 - int(sub_feature_phase)) % 3)))
                                line = GFF3_opened.readline()
                            else:
                                raise Exception('GFF3_parser(ERROR): transcript (' + general_line_ENST + ') name and CDS-parent name (' + sub_feature_parent + ') do not match ...')
                        else:
                            line = GFF3_opened.readline()
                else:
                    raise Exception('GFF3_parser(ERROR): snip-bit "ID=transcript:ENST" detected in line and the line was expected to have 9 tab-delimited columns, instead ' + str(len(line.rstrip('\n').split('\t'))) + ' columns found ... \nline reads: ' + str(line))
            else:
                line = GFF3_opened.readline()
    
    else:
        raise Exception('GFF3_parser(ERROR): the provided file is not of the format "gff3"; first line does not start with "##gff-version 3\n"')    
    GFF3_opened.close()
    return exon_coordinate_dict

# this function parses the results of polyQ_finder() which are stored as csv
#and prepares them for use by compare_suject_vs_reference()
#@profile
def csv_to_exonic_lot_and_sot_parser(csv_path, ChrOI):
    
    #the following replaces the next hashed line:
    csv_genomic_pre_lol = [(line.rstrip('\n').split(','))[:-1] for line in open(csv_path)]
    
    #because of the repeats having been identified by regex matches, there should be no overlap of the repeats
    #meaning sorting them by starting position and ENST will provide the definite linear order of their physical location
    csv_genomic_sorted_for_position_lol = (sorted(csv_genomic_pre_lol, key=return6th_as_integer))
    csv_genomic_sorted_for_ENST_lol = (sorted(csv_genomic_sorted_for_position_lol, key=return3rd))
    former_ENST = ''
    new_ENST = ''
    position_count = 1
    csv_genomic_lol = []
    for line in csv_genomic_sorted_for_ENST_lol:
        new_ENST = line[2]
        if new_ENST != former_ENST:
            position_count = 1
        elif new_ENST == former_ENST:
            position_count += 1

            
        #as a reminder, here is what a line in csv_genomic_sorted_for_ENST_lol looks like:
        #'4', 'ENSG00000145242.9', 'ENST00000511294.1', '-', '5.0', '66521831', '66521846', 'TTGCTGCTGCTGCTG',
        # so the new_line will look like this:
        #   '4'(chromosome) = line[0], 
        #   'ENSG00000145242.9' (ENSG) = line[1], 
        #   'ENST00000511294.1' (ENST) = line[2], 
        #   '-' (strand) = line[3], 
        #   '5.0'(number of repeats) = line[4], 
        #   '1.0' ('rank') == (placement in all repeats identified in this ENST) = str(position_count),
        #   'TTGCTGCTGCTGCTG' (repeat seq) = line[7],
        #   '66521831' (repeat start) = line[5], 
        #   '66521846' (repeat end) = line[6],

        new_line = (line[0], line[1], line[2], line[3], line[4], str(position_count), line[7].upper(), line[5], line[6])
        csv_genomic_lol.append(new_line)
        former_ENST = new_ENST
        
    #to find ENST with changes in repeat structure the following will be needed: chromosome, ENSG, ENST, strand, repeat length, number of repeats in this ENST, repeat-seq
    csv_exonic_lot = [tuple(csv_genomic_lol_entry) for csv_genomic_lol_entry in csv_genomic_lol if csv_genomic_lol_entry[0] == ChrOI]# and (check_if_cds(ENST_of_interest = csv_genomic_lol_entry[2], left_read_coordinate = csv_genomic_lol_entry[7], right_read_coordinate = csv_genomic_lol_entry[8], GFF3_dict = reference_exon_ranges_dict)]
    #csv_exonic_lot = [tuple(csv_genomic_lol_entry[:-2]) for csv_genomic_lol_entry in csv_genomic_lol if (check_if_cds(ENST_of_interest = csv_genomic_lol_entry[2], left_read_coordinate = csv_genomic_lol_entry[7], right_read_coordinate = csv_genomic_lol_entry[8], GFF3_dict = reference_exon_ranges_dict) and csv_genomic_lol_entry[0] == ChromosomeOfInterest)]
    
    csv_exonic_ENST_Seq_sot = set([tuple([(line.rstrip('\n').split(','))[2], (line.rstrip('\n').split(','))[8]]) for line in open(csv_path) if (line.split(',')[0] == ChrOI)])#(check_if_cds(ENST_of_interest = (line.rstrip('\n').split(','))[2], left_read_coordinate = (line.rstrip('\n').split(','))[5], right_read_coordinate = (line.rstrip('\n').split(','))[6], GFF3_dict = reference_exon_ranges_dict) and  (line.rstrip('\n').split(',')[0][0] == ChromosomeOfInterest))])

    return (csv_exonic_lot, csv_exonic_ENST_Seq_sot)

###############################################################################
#end of parsing functions
###############################################################################



###############################################################################
#start of functions needed as 'key' for sorted()
###############################################################################

#@profile
def returnself(anything):
    return anything
#@profile
def return1st(anything):
    return anything[0]
#@profile
def return2nd(anything):
    return anything[1]
#@profile
def return3rd(anything):
    return anything[2]
#@profile
def return4th(anything):
    return anything[3]
#@profile
def return5th(anything):
    return anything[4]
#@profile
def return6th(anything):
    return anything[5]
#@profile
def return7th(anything):
    return anything[6]
#@profile
def return8th(anything):
    return anything[7]
#@profile
def return9th(anything):
    return anything[8]
#@profile
def return10th(anything):
    return anything[9]
#@profile
def return1st_case_insensitive(anything):
    return anything[0].lower()

#@profile
def returnself_as_absolute_value(anything):
    if anything < 0:
        return (anything * (-1))
    else:
        return anything
    #@profile
def return6th_as_integer(anything):
    return int(anything[5])    

#@profile
def return_read_start_as_absolute_value(anything):
    #the first base satrts at position [1], and its genomic location is at position [1][2]
    if anything[1][2] < 0:
        return (anything[1][2] * (-1))
    else:
        return anything[1][2]
    
###############################################################################
#end of functions needed as 'key' for sorted()
###############################################################################



###############################################################################
#start of read-matrix analyzing functions
###############################################################################
        
# this function calls the bases from all aligned reads to build the consensus
# matrix - it will pick the heighest weighed base and will call iupac-bases for
# ambigious cases with two or more equally often aligned bases per position
#@profile
def call_base(weighed_position):
    weighed_position_labels = [[weighed_position[0], 'A'], [weighed_position[1], 'T'], [weighed_position[2], 'C'], [weighed_position[3], 'G'], [weighed_position[4], '-'], [weighed_position[5], 'N']]
    highest_score = (sorted(weighed_position_labels, key=return1st, reverse=True))
    #test if ambiguity code is needed but maybe I should get the length of all cahracters that match the highest score first
    
    if highest_score[0][0] != highest_score[1][0]:
        return highest_score[0][1]
    
    elif highest_score[0][0] == highest_score[1][0]!= highest_score[2][0]:
        duo_one = highest_score[0][1]
        duo_two = highest_score[1][1]
        if set((duo_one, duo_two)) == set(('C', 'T')):
            return 'R'
        elif set((duo_one, duo_two)) == set(('A', 'G')):
            return 'Y'
        elif set((duo_one, duo_two)) == set(('A', 'T')):
            return 'W'
        elif set((duo_one, duo_two)) == set(('G', 'C')):
            return 'S'
        elif set((duo_one, duo_two)) == set(('T', 'G')):
            return 'M'
        elif set((duo_one, duo_two)) == set(('C', 'A')):
            return 'K'
        elif (set((duo_one, duo_two)) == set(('N', 'A'))) or (set((duo_one, duo_two)) == set(('N', 'T'))) or (set((duo_one, duo_two)) == set(('N', 'C'))) or (set((duo_one, duo_two)) == set(('N', 'G'))):
            return 'N'
        elif (set((duo_one, duo_two)) == set(('-', 'A'))) or (set((duo_one, duo_two)) == set(('-', 'T'))) or (set((duo_one, duo_two)) == set(('-', 'C'))) or (set((duo_one, duo_two)) == set(('-', 'G'))) or (set((duo_one, duo_two)) == set(('-', 'N'))):
            return '-'
        
    elif highest_score[0][0] == highest_score[1][0] == highest_score[2][0] != highest_score[3][0]:
        triplet_one = highest_score[0][1]
        triplet_two = highest_score[1][1]
        triplet_three = highest_score[2][1]
        if set((triplet_one, triplet_two, triplet_three)) == set(('A','G','T')):
            return 'H'
        elif set((triplet_one, triplet_two, triplet_three)) == set(('A','C','G')):
            return 'B'
        elif set((triplet_one, triplet_two, triplet_three)) == set(('A','C','T')):
            return 'D'
        elif set((triplet_one, triplet_two, triplet_three)) == set(('C','G','T')):
            return 'V'
        else:
            return 'N'
    
    elif (highest_score[0][0] == highest_score[1][0] == highest_score[2][0] == highest_score[3][0] != highest_score[4][0]) or (highest_score[0][0] == highest_score[1][0] == highest_score[2][0] == highest_score[3][0] == highest_score[4][0] != highest_score[5][0]) or (highest_score[0][0] == highest_score[1][0] == highest_score[2][0] == highest_score[3][0] == highest_score[4][0] == highest_score[5][0]):
        return 'N'
    
    else:
        return 'N'
    #'(A|T|C|G|-|N|R|Y|W|S|M|K|X|H|B|D|V)'
    #from https://www.gendx.com/SBTengine/Help_220/hs310.htm:
    # >These Rules are as close as possible to the published version [see Biochem. J., 1985, 229, 281-286; Eur. J. Biochem., 1985, 150, 1-5; J. Biol. Chem., 1986, 261, 13-17; Mol. Biol. Evol., 1986, 3, 99-108; Nucl. Acids Res., 1985, 13, 3021-3030; Proc. Nat. Acad. Sci. (U. S.), 1986, 83, 4-8<

# this function builds the matrix to be used by call_base() and directly calls
#call_base() to action 
#@profile
def consensus_writer(matrix_of_read_lots):
    new_matrix = []
    matrix_of_read_lols_sorted_by_pos1_per_read = sorted(matrix_of_read_lots, key=return_read_start_as_absolute_value)
    first_position = matrix_of_read_lols_sorted_by_pos1_per_read[0][1][2]
    
    last_position_per_read_list = []
    for single_read_lot in matrix_of_read_lots:
        #add counted Insertions and the position of the last base for true last alignment position
        last_position_per_read_list.append(single_read_lot[-1][2])# + single_read_lot[0])
    last_position_first = sorted(last_position_per_read_list, key=returnself, reverse=True)
    last_position = last_position_first[0]
    
    for i in xrange(first_position, last_position + 1, 1):
        # A, T, C, G, -, X
        new_matrix.append([0,0,0,0,0,0,])
     
    for single_read_lot in matrix_of_read_lols_sorted_by_pos1_per_read:
        
        insertion_log = 0
        for base, reference_position, read_position, InDel_counter in single_read_lot[1:]:
            
            if read_position > 0:
                column = read_position - first_position
            elif read_position < 0:
                insertion_log += 1
                column = (read_position * -1) - first_position
    
            if base == 'A':
                new_matrix[column][0] += 1
            elif base == 'T':
                new_matrix[column][1] += 1
            elif base == 'C':
                new_matrix[column][2] += 1
            elif base == 'G':
                new_matrix[column][3] += 1
            elif base == '-':
                new_matrix[column][4] += 1
            elif base == 'X' or 'N':
                new_matrix[column][5] += 1
            else:
                new_matrix[column][5] += 1
                
            new_matrix[column].append(reference_position)
            new_matrix[column].append(read_position)
            new_matrix[column].append(InDel_counter)
    
    consensus = ''
    for position in new_matrix:
        consensus += call_base(position[0:6])
    return [consensus, first_position, new_matrix]
    
###############################################################################
#end of read-matrix analyzing functions
###############################################################################



###############################################################################
#start of polyQ-detecting functions
###############################################################################

# this function finds repeat sections via regex re.finditer search and chekcs
# the hit's validity with check_if_in_phase (which checks for phase AND coding 
# regions)
#@profile
def polyQ_finder(seq_in, chromosome_in, ENSG_in, ENST_in, strand_in, starting_coordinate, exon_dict):
    #this finds clean repeats of "CAA" and "CAG" that could code for Q on the "+" strand_in, 
    #or their reverse compliments "TTG" and "CTG" on the "-" strand_in 
    #and returns them as a csv-row
    return_string = ''
    if strand_in == '+':
        pattern = re.finditer(r'((CA[AG])+.{3})*(CA[AG]){4,}(.{3}(CA[AG])+)*', seq_in, flags = re.I)
    elif strand_in == '-':
        pattern = re.finditer(r'(([TC]TG)+.{3})*([TC]TG){4,}(.{3}([TC]TG)+)*', seq_in, flags = re.I)

    for repeat in pattern:    
        repeat_seq = repeat.group()
        
        repeat_start = repeat.start()
        #this too is to counter pythons 0-based ranges
        repeat_stop = (repeat.end() - 1)
        #TODO: change the [0][0], changes in HRA_bed_parser() and a rerun of it will be required
                
        len_seq = (len(seq_in))
        left_fig = 50
        right_fig = 50
        if repeat_start <= left_fig:
            left_fig = repeat_start
        if (len_seq - repeat_stop) <= right_fig:
            right_fig = (len_seq - repeat_stop)
        seq_plus_minus_50 = (seq_in[(repeat_start - left_fig):(repeat_stop + right_fig)])
    
        if check_if_in_phase(ENST_of_interest = ENST_in, left_query_coordinate = (int(repeat_start) + int(starting_coordinate)), right_query_coordinate = (int(repeat_stop) + int(starting_coordinate)), strand = strand_in, GFF3_dict = exon_dict) == True:
            return_string += str(chromosome_in) + ',' + str(ENSG_in) + ',' + str(ENST_in) + ',' + str(strand_in) + ',' + str(len(repeat_seq)/3) + ',' + str(int(repeat_start) + int(starting_coordinate)) + ',' + str(int(repeat_stop) + int(starting_coordinate)) + ',' + str(repeat_seq) + ',' + seq_in[(int(repeat_start)):(int(repeat_stop))] + ',' + str(starting_coordinate) + ',' + str(starting_coordinate + len(seq_in)) + ',' + str(seq_plus_minus_50) + ',' + str(seq_in) +'\n'
    return return_string

# this function is called by check_if_in_phase() and checks if a repeat hit is
# in a CDS region defined by the gff3_dict provided by gff3_parser()
#@profile
def check_if_in_cds(ENST_of_interest, left_query_coordinate, right_query_coordinate, CDS_range_start, CDS_range_stop):
    print('\nchecking for ' + str(ENST_of_interest) + ' if repeat ' +  str(left_query_coordinate) + '-' + str(right_query_coordinate) + ' matches coding exon ' + str(CDS_range_start) + '-' + str(CDS_range_stop) + ' and its phase')
    query_start = int(left_query_coordinate)
    query_stop = int(right_query_coordinate)
    print('check_if_in_cds(): checking for ' + str(ENST_of_interest) + ' if the section in question ' +  str(left_query_coordinate) + '-' + str(right_query_coordinate) + ' is in exon ' + str(CDS_range_start) + ' - ' + str(CDS_range_stop))

    if (not ((query_start < CDS_range_start) and (query_stop < CDS_range_start)) and not ((query_start > CDS_range_stop) and (query_stop > CDS_range_stop))):
        print('Yes!')
        return True
    else:
        print('No!')
        return False
    
# this function checks if a repeat hit is in phase as defined by the gff3_dict 
# provided by gff3_parser()
def check_if_in_phase(ENST_of_interest, left_query_coordinate, right_query_coordinate, strand, GFF3_dict):
    query_start = int(left_query_coordinate)
    query_stop = int(right_query_coordinate)
    print('repeat range : ' + str(left_query_coordinate) + '|' + str(right_query_coordinate))
    gff3_CDS_list = GFF3_dict[ENST_of_interest]
    print('gff3_cds_list: ' + str(gff3_CDS_list))
    for CDS in gff3_CDS_list:
        #sub_feature == coding exon == CDS
        CDS_range_start = CDS[0] 
        CDS_range_stop = CDS[1]
        sub_feature_phase = CDS[2]

        if check_if_in_cds(ENST_of_interest = ENST_of_interest, left_query_coordinate = left_query_coordinate, right_query_coordinate = right_query_coordinate, CDS_range_start = CDS_range_start, CDS_range_stop = CDS_range_stop):

            print('check_if_in_phase(): checking for ' + str(ENST_of_interest) + ' if the repeat ' +  str(left_query_coordinate) + '-' + str(right_query_coordinate) + ' is in phase of the following exon ' + str(CDS_range_start) + ',' + str(strand) + ',' + str(sub_feature_phase))

            if strand == '+':
                print('((' + str(query_start) + ' - (' + str(CDS_range_start) + ' - ' + str(sub_feature_phase) + ')) % 3) = ' + str(((query_start - (CDS_range_start - sub_feature_phase)) % 3)) + ')')
                print('((query_start - (CDS_range_start - sub_feature_phase)) % 3) == 0)')
                if ((query_start - (CDS_range_start - sub_feature_phase)) % 3) == 0:
                    print('Yes!')
                    return True
                else:
                    print('No!')
            if strand == '-':
                print('(((' + str(CDS_range_stop) + ' + ' + str(sub_feature_phase) + ') - ' + str(query_stop) + ') % 3) = ' + str((((CDS_range_stop + sub_feature_phase) - query_stop) % 3)) + ')')
                print('(((CDS_range_stop + sub_feature_phase) - query_stop) % 3) == 0')
                if (((CDS_range_stop + sub_feature_phase) - query_stop) % 3) == 0:
                    print('Yes!')
                    return True
                else:
                    print('No!')
    return False

#this function translates a reference seq to AAS, not usefull yet, as the subject can't be translated ...
def translate_seq_to_AA(ENST, strand, full_seq, GFF3_dict):
    nucleotide_string = ''
    AA_Seq = ''
    gff3_CDS_list = GFF3_dict[ENST]
    for CDS in gff3_CDS_list:
        CDS_range_start = CDS[0] 
        CDS_range_stop = CDS[1]
        #sub_feature_phase = int(CDS[2])
        nucleotide_string += full_seq[CDS_range_start:CDS_range_stop]
    
    nucleotide_string = nucleotide_string.upper()
    
    if strand == '-':
        nucleotide_string == nucleotide_string[::-1]
    
    for step in xrange(0,len(nucleotide_string),3):
        triplet_start = step
        triplet_stop = step + 3
        codon = nucleotide_string[triplet_start:triplet_stop]
        if re.compile('(A|T|C|G){3}').match(codon):
            AA_Seq += AAS_table[codon]
        else:
            AA_Seq += 'X'

    return (nucleotide_string, AA_Seq)

def simple_translate_seq_to_AA(strand, string_triplet):
    AA_Seq = '|'
    string_triplet = [string_triplet[0].upper(), string_triplet[1].upper(), string_triplet[2].upper(), ]
    
    if string_triplet[1] != '':
        print('string triplet: ' + str(AAS_table[string_triplet[1]]))
    
    if strand == '-':
        string_triplet = [string_triplet[2][::-1].replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper(), string_triplet[1][::-1].replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper(), string_triplet[0][::-1].replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper(),]
    for nucleotide_string in string_triplet:
            
        for step in xrange(0,len(nucleotide_string),3):
            triplet_start = step
            triplet_stop = step + 3
            codon = nucleotide_string[triplet_start:triplet_stop]
            if re.compile('(A|T|C|G){3}').match(codon):
                AA_Seq += AAS_table[codon]
            else:
                AA_Seq += 'X'
        AA_Seq += '|'
    
    return (AA_Seq)

###############################################################################
#end of polyQ-detecting functions
###############################################################################



###############################################################################
#start of comparison processing functions that return statistical data on the results
###############################################################################
    
def cluster_hit_if_overlapping(list_of_clusters, submitted_hit):
    submitted_start = submitted_hit[6]
    submitted_stop = submitted_hit[7]
    if list_of_clusters == ['$']:
        return [[submitted_hit]]
    for cluster in list_of_clusters:
        for hit in cluster:
            if hit[0] == 'REF':
                hit_start = hit[6]
                hit_stop = hit[7]
                if not ((submitted_stop < hit_start) and (submitted_start < hit_start)) and not ((submitted_stop > hit_stop) and (submitted_start > hit_stop)):
                    old_cluster = cluster
                    cluster.append(submitted_hit)                    
                    new_list_of_clusters = [element for element in list_of_clusters if element != old_cluster]
                    new_list_of_clusters.append(cluster)
                    return new_list_of_clusters
    list_of_clusters.append([submitted_hit])
    return list_of_clusters
    
def get_statistics(working_directory, data_table, experiment_list, desired_information, bed_name, ref_name, gff3, list_of_passing_first_lengths = 'NA', list_of_passing_second_lengths = 'NA', aligner = 'NA', ):
    
    if list_of_passing_first_lengths == 'NA' or list_of_passing_first_lengths == 'NA':
        raise Exception('get_statistics(ERROR): please enter a list of integers or "any" as list_of_passing_first_lengths and list_of_passing_first_lengths parameters')
    
    if desired_information == 'greater':
        desired_operator = greater_func
    elif desired_information == 'lesser':
        desired_operator = lesser_func 
    elif desired_information == 'equal':
        desired_operator = equal_func
    elif desired_information == 'any_seq_changes':
        desired_operator = any_seq_changes_func
    if desired_information == 'greater_length_and_any_seq_changes':
        desired_operator = greater_length_and_any_seq_changes_func
    elif desired_information == 'lesser_length_and_any_seq_changes':
        desired_operator = lesser_length_and_any_seq_changes_func 
    elif desired_information == 'equal_length_and_any_seq_changes':
        desired_operator = equal_length_and_any_seq_changes_func
    
    list_of_ROI = []
    
    df = pd.DataFrame.from_dict(data_table)
    
    # the data frame will have the following structure (example)
    # number example value                                          variable name
    #	0	6   												    Chromosome
    #	1	ENSG00000010017.9									    ENSG
    #	2 	ENST00000011619.3  									ENST
    #	3 	AATGTAGGTAGTCTTCCACTGTGGCAAATGCGCAGGATCCAATTCC...    	ENST_seq_plus_introns
    #	4 	13622596   											ENST_start
    #	5 	13711783                                                ENST_stop
    #	6 	0                                                       by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed
    #	7 	0                          							by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed
    #	8 	ZVE  												experiment_label
    #	9 	AGGACGACTCCGGAGACNNNNNNNNNNNNNNNNNNNGGTGGCGGCG...       lflank_seq                     
    #	10	CTGCTGCTGCTGTTGCTGCTGCTG              					middle_repeat_seq
    #	11	8                           							repeat_length
    #	12	CTGCTGCTGCTGTTGCTGCTGCTG  							repeat_seq
    #	13	AGGACGACTCCGGAGACNNNNNNNNNNNNNNNNNNNGGTGGCGGCG...       repeat_seq_plus_50_bases_left_and_right 
    #	14	13711684                                  				repeat_start
    #	15	89088     											repeat_start_in_ENST_seq_plus_introns
    #	16 	13711707                                 				repeat_stop
    #	17 	89112  												repeat_stop_in_ENST_seq_plus_introns
    #	18 	CGGCGGCGGCGGCGGCGGCTGCCCGGACATCCCGGCCGCGACTCAG...       rflank_seq   
    #	19 	89038         										start_minus_50
    #	20 	89162      											stop_plus_50
    #	21 	-													strand
    
    print('df:')
    print(df.to_string())
    all_ENSGs = (set(df['ENSG'].tolist()))
    for single_ENSG in all_ENSGs:
        all_ENSTs_per_ENSG_df = df[df.ENSG == single_ENSG]
        iterable_rows = all_ENSTs_per_ENSG_df.iterrows()
        for first_index, first_row in iterable_rows:
            first_Chromosome = str(first_row[0])
            first_ENSG = str(first_row[1])
            first_ENST = str(first_row[2])
            #first_ENST_seq_plus_introns = str(first_row[3])
            #first_ENST_start = str(first_row[4])
            #first_ENST_stop = str(first_row[5])
            #first_by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed = str(first_row[6])
            #first_by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed = str(first_row[7])
            first_experiment_label = str(first_row[8])
            #first_lflank_seq = str(first_row[9])
            #first_middle_repeat_seq = str(first_row[10])
            first_repeat_length = int(first_row[11])
            first_repeat_seq = str(first_row[12])
            #first_repeat_seq_plus_50_bases_left_and_right = str(first_row[13])
            first_repeat_start = int(first_row[14])
            #first_repeat_start_in_ENST_seq_plus_introns = str(first_row[15])
            first_repeat_stop = int(first_row[16])
            #first_repeat_stop_in_ENST_seq_plus_introns = str(first_row[17])
            #first_rflank_seq = str(first_row[18])
            #first_start_minus_50 = str(first_row[19])
            #first_stop_plus_50 = str(first_row[20])
            #first_strand = str(first_row[21])
            for second_index, second_row in iterable_rows:
                second_Chromosome = str(second_row[0])
                second_ENSG = str(second_row[1])
                second_ENST = str(second_row[2])
                #second_ENST_seq_plus_introns = str(second_row[3])
                #second_ENST_start = str(second_row[4])
                #second_ENST_stop = str(second_row[5])
                #second_by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed = str(second_row[6])
                #second_by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed = str(second_row[7])
                second_experiment_label = str(second_row[8])
                #second_lflank_seq = str(second_row[9])
                #second_middle_repeat_seq = str(second_row[10])
                second_repeat_length = int(second_row[11])
                second_repeat_seq = str(second_row[12])
                #second_repeat_seq_plus_50_bases_left_and_right = str(second_row[13])
                second_repeat_start = int(second_row[14])
                #second_repeat_start_in_ENST_seq_plus_introns = str(second_row[15])
                second_repeat_stop = int(second_row[16])
                #second_repeat_stop_in_ENST_seq_plus_introns = str(second_row[17])
                #second_rflank_seq = str(second_row[18])
                #second_start_minus_50 = str(second_row[19])
                #second_stop_plus_50 = str(second_row[20])
                #second_strand = str(second_row[21])
                if (first_experiment_label != second_experiment_label) and (first_ENST == second_ENST):
                    if ((list_of_passing_first_lengths == 'any') or (int(first_repeat_length) in list_of_passing_first_lengths)) and ((list_of_passing_second_lengths == 'any') or (int(second_repeat_length) in list_of_passing_second_lengths)):
                        if desired_operator(first_repeat_length = first_repeat_length, first_repeat_start = first_repeat_start, first_repeat_stop = first_repeat_stop, first_repeat_seq = first_repeat_seq, second_repeat_length = second_repeat_length, second_repeat_start = second_repeat_start, second_repeat_stop = second_repeat_stop, second_repeat_seq = second_repeat_seq):
                            if first_repeat_seq != 'NA':
                                list_of_ROI.append((first_ENST, first_repeat_start, first_repeat_stop, first_ENSG, first_Chromosome))
                            if second_repeat_seq != 'NA':
                                list_of_ROI.append((second_ENST, second_repeat_start, second_repeat_stop, second_ENSG, second_Chromosome))
    
    #elminiate duplicates
    set_of_ROI = set(list_of_ROI)
    #turn into list for indexing and sorting
    unique_list_of_ROI = list (set_of_ROI)
    #sort
    sorted_stop_list_of_ROI = sorted(unique_list_of_ROI, key=return3rd)
    sorted_start_list_of_ROI = sorted(sorted_stop_list_of_ROI, key=return2nd)
    sorted_ENST_list_of_ROI = sorted(sorted_start_list_of_ROI, key=return1st)
    
    print('\nsorted_ENST_list_of_ROI: ' + str(sorted_ENST_list_of_ROI))
    
    #identify clusters by ovelap within ENSTs
    lol_of_clusters = []
    first_member_of_cluster = sorted_ENST_list_of_ROI[0]
    last_member_of_cluster = (0, 0, 0)
    farthest_reach_of_cluster = 0
    
    for entry in sorted_ENST_list_of_ROI:
        #if this is the last entry wrap up the final cluster
        if sorted_ENST_list_of_ROI.index(entry) == (len(sorted_ENST_list_of_ROI) - 1):
            if farthest_reach_of_cluster == 0:
                lol_of_clusters.append(list(sorted_ENST_list_of_ROI[sorted_ENST_list_of_ROI.index(first_member_of_cluster):sorted_ENST_list_of_ROI.index(entry)]))
                lol_of_clusters.append([(entry)])
            elif entry[0] == last_member_of_cluster[0] and (entry[1] <= last_member_of_cluster[2] or entry[1] <= farthest_reach_of_cluster):
                lol_of_clusters.append(list(sorted_ENST_list_of_ROI[sorted_ENST_list_of_ROI.index(first_member_of_cluster):(sorted_ENST_list_of_ROI.index(entry) + 1)]))
            elif (entry[0] != last_member_of_cluster[0]) or (entry[0] == last_member_of_cluster[0] and entry[1] > last_member_of_cluster[2]):
                #for the end coordinate we take the entry index, since it is exclusive
                lol_of_clusters.append(list(sorted_ENST_list_of_ROI[sorted_ENST_list_of_ROI.index(first_member_of_cluster):sorted_ENST_list_of_ROI.index(entry)]))
                lol_of_clusters.append([(entry)])
            else:
                raise Exception('HRA_compare_all_QRYs_and_REF(ERROR): unexpected clustering possibility')
        #continue an established cluster if ...
        elif (entry[0] == last_member_of_cluster[0] and (entry[1] <= last_member_of_cluster[2] or entry[1] <= farthest_reach_of_cluster)):
            last_member_of_cluster = entry
            if entry[2] > farthest_reach_of_cluster:
                farthest_reach_of_cluster = entry[2]
        #wrap up the cluster if ...
        elif (entry[0] != last_member_of_cluster[0]) or (entry[0] == last_member_of_cluster[0] and entry[1] > last_member_of_cluster[2]):
            lol_of_clusters.append(list(sorted_ENST_list_of_ROI[sorted_ENST_list_of_ROI.index(first_member_of_cluster):sorted_ENST_list_of_ROI.index(entry)]))
            first_member_of_cluster = entry
            last_member_of_cluster = entry
            farthest_reach_of_cluster = 0
        else:
            raise Exception('HRA_compare_all_QRYs_and_REF(ERROR): unexpected clustering possibility')
    
    #remove the initial empty one
    if lol_of_clusters[0] == []:
        del lol_of_clusters[0]
    #double check if everything went well
    number_of_hits = 0
    for cluster in lol_of_clusters:
        for hit in cluster:
            number_of_hits += 1
    if number_of_hits != len(sorted_ENST_list_of_ROI):
        raise Exception('HRA_compare_all_QRYs_and_REF(ERROR): it seems that some hits were lost during clustering')
    #gather the sequences for MS-alignments
    cluster_counter = 0
    for cluster in lol_of_clusters:
        cluster_counter += 1
        fa_aln_text = ''
        experiment_log = []
        farthes_left_position = int(cluster[0][1])
        farthes_right_position = int(cluster[0][2])
        cluster_ENST = cluster[0][0]
        cluster_ENSG = cluster[0][3]
        cluster_chromosome = cluster[0][4]
        cluster_pre_dec = cluster_ENSG.split('.')[0][-1]
        align_pickle_name = ''
        if (cluster_chromosome == 'Y') or (cluster_chromosome == 'MT'):
            if cluster_pre_dec == '1' or cluster_pre_dec == '2':
                align_pickle_name = 'P_' + cluster_chromosome + '_(1|2)_align_pickle.txt'
            if cluster_pre_dec == '3' or cluster_pre_dec == '4':
                align_pickle_name = 'P_' + cluster_chromosome + '_(3|4)_align_pickle.txt'
            if cluster_pre_dec == '5' or cluster_pre_dec == '6':
                align_pickle_name = 'P_' + cluster_chromosome + '_(5|6)_align_pickle.txt'
            if cluster_pre_dec == '7' or cluster_pre_dec == '8':
                align_pickle_name = 'P_' + cluster_chromosome + '_(7|8)_align_pickle.txt'
            if cluster_pre_dec == '9' or cluster_pre_dec == '0':
                align_pickle_name = 'P_' + cluster_chromosome + '_(9|0)_align_pickle.txt'
        else:
            align_pickle_name = 'P_' + cluster_chromosome + '_' + cluster_pre_dec + '_align_pickle.txt'
            
        print('cluster_ENST: ' + cluster_ENST)
        print('cluster_ENSG: ' + cluster_ENSG)
        for hit_ENST, hit_start, hit_stop, hit_ENSG, hit_Chr in cluster:
            print('\n\nnon-redundant cluster entry ->   hit_ENST: ' + str(hit_ENST) + ' |hit_start: ' + str(hit_start) + ' |hit_stop: ' + str(hit_stop) + ' |hit_ENSG: ' + str(hit_ENSG) + ' |hit_Chr: ' + str(hit_Chr))
            #keep track of the farthest outreaches among all samples
            if int(hit_start) < farthes_left_position:
                farthes_left_position = int(hit_start)
            if int(hit_stop) > farthes_right_position:
                farthes_right_position = int(hit_stop)
                
            #the cluster is non-redundant, therefor all matching entries can be extracted in the upcoming for-loop
            cluster_ENST_df = df[(df.ENST == cluster_ENST) & (df.repeat_start == hit_start) & (df.repeat_stop == hit_stop)]
            print('cluster_ENST_df.to_string(): ' + cluster_ENST_df.to_string())
            iterable_rows = cluster_ENST_df.iterrows()
                
            for index, row in iterable_rows:
                    
                header = '>'
                aln_Chromosome = str(row[0])
                aln_ENSG = str(row[1])
                aln_ENST = str(row[2])
                #aln_ENST_seq_plus_introns = str(row[3])
                aln_ENST_start = str(row[4])
                aln_ENST_stop = str(row[5])
                #aln_by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed = str(row[6])
                #aln_by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed = str(row[7])
                aln_experiment_label = str(row[8])
                #aln_lflank_seq = str(row[9])
                #aln_middle_repeat_seq = str(row[10])
                aln_repeat_length = str(row[11])
                #aln_repeat_seq = str(row[12])
                aln_repeat_seq_plus_50_bases_left_and_right = str(row[13])
                aln_repeat_start = str(row[14])
                #aln_repeat_start_in_ENST_seq_plus_introns = str(row[15])
                aln_repeat_stop = str(row[16])
                #aln_repeat_stop_in_ENST_seq_plus_introns = str(row[17])
                #aln_rflank_seq = str(row[18])
                #aln_start_minus_50 = str(row[19])
                #aln_stop_plus_50 = str(row[20])
                aln_strand = str(row[21])
                #print('aln_strand: ' + str(aln_strand))
                header += 'experiment_label:' + aln_experiment_label
                header += '_repeat_start:' + str(aln_repeat_start)
                header += '_repeat_stop:' + str(aln_repeat_stop)
                header += '_Chromosome:' + aln_Chromosome
                header += '_strand:' + aln_strand
                header += '_ENSG:' + aln_ENSG
                header += '_ENST:' + aln_ENST
                #header += '_ENST_seq_plus_introns:' + aln_ENST_seq_plus_introns
                header += '_ENST_start:' + aln_ENST_start
                header += '_ENST_stop:' + aln_ENST_stop
                header += '_repeat_length:' + aln_repeat_length
                #header += '_repeat_seq:' + aln_repeat_seq
                #header += '_repeat_seq_plus_50_bases_left_and_right:' + aln_repeat_seq_plus_50_bases_left_and_right
                header += '_'
                fa_aln_text += header + '\n'
                fa_aln_text += aln_repeat_seq_plus_50_bases_left_and_right + '\n'
                print('fa_aln_text: ' + fa_aln_text)
                experiment_log.append(aln_experiment_label)
        list_in_case_of_REF_not_being_in_the_df = [] 
        for experiment in experiment_list:
            print('experiment: ' + experiment)
            if experiment not in experiment_log and experiment != 'REF':
                print('experiment is neither in experiment_log nor "REF"')
                HRA_mapping_analyzer_dict_path = working_directory + experiment + '/analyze_pickles/' + align_pickle_name
                print('HRA_mapping_analyzer_dict_path: ' + str(HRA_mapping_analyzer_dict_path))
                HRA_mapping_analyzer_dict_path_open_for_reading = open(HRA_mapping_analyzer_dict_path, 'rb')
                HRA_mapping_analyzer_dict = pickle.load(HRA_mapping_analyzer_dict_path_open_for_reading)
                    
                try:
                    print('trying to retrieve consensus_ENST_path')
                    consensus_ENSTs_path = HRA_mapping_analyzer_dict[cluster_ENST + '|consensus_ENSTs_path']
                    print('HRA_mapping_analyzer_dict[' + cluster_ENST + "|consensus_ENSTs_path'] was retrieved ... consensus_ENSTs_path: " + str(consensus_ENSTs_path))
                    consensus_leftmost_coordinate = HRA_mapping_analyzer_dict[cluster_ENST + '|consensus_ENST_leftmost_coordinate']                    
                    print('HRA_mapping_analyzer_dict[' + cluster_ENST + "|consensus_ENST_leftmost_coordinate'] was retrieved ... consensus_leftmost_coordinate: " + str(consensus_leftmost_coordinate))
                    
                    if consensus_ENSTs_path != [] and consensus_ENSTs_path != '':
                        print('consensus_ENSTs_path != [] and consensus_ENSTs_path != ""')
                        consensus_file_open_for_reading = open(consensus_ENSTs_path, 'r')
                        consensus_file_read = consensus_file_open_for_reading.read()
                        consensus_seq = consensus_file_read[0]
                        print('first 50 bases of consensus_seq: ' + consensus_seq[0:50])
                        print('consensus_leftmost_coordinate: ' + str(consensus_leftmost_coordinate))
                        print('farthes_left_position: ' + str(farthes_left_position))
                        print('consensus_leftmost_coordinate: ' + str(consensus_leftmost_coordinate))
                        print('farthes_right_position: ' + str(farthes_right_position))
                        print("fa_aln_text += '>' + experiment + ':seq_retrieved_from_consensus_seq\n' + consensus_seq[(farthes_left_position - consensus_leftmost_coordinate):(farthes_right_position - consensus_leftmost_coordinate)] + '\n'")
                        print("the following are the  seq's coordinates: " +  str(farthes_left_position) + ' - ' + str(consensus_leftmost_coordinate) + ' : ' + str(farthes_right_position) + ' - ' + str(consensus_leftmost_coordinate))
                        header = '>experiment_label:' + experiment.lower() + '_repeat_start:' + str(farthes_left_position) + '_repeat_stop:' + str(farthes_right_position) + '_Chromosome:' + cluster_chromosome + '_strand:/_ENSG:' + cluster_ENSG + '_ENST:' + cluster_ENST + '_consensus_leftmostcoordinate:' + str(consensus_leftmost_coordinate) + '_repeat_length:/_\n'

                        if ((int(farthes_left_position) - 50) - int(consensus_leftmost_coordinate)) < 0 and ((int(farthes_right_position) + 50) - int(consensus_leftmost_coordinate)) > len(consensus_seq):
                            fa_aln_text += header + consensus_seq + '\n'
                        elif ((int(farthes_left_position) - 50) - int(consensus_leftmost_coordinate)) < 0:
                            fa_aln_text += header + consensus_seq[:((int(farthes_right_position) + 50) - int(consensus_leftmost_coordinate))] + '\n'
                        elif ((int(farthes_right_position) + 50) - int(consensus_leftmost_coordinate)) > len(consensus_seq):
                            fa_aln_text += header + consensus_seq[((int(farthes_left_position) - 50) - int(consensus_leftmost_coordinate)):] + '\n'
                        else:
                            fa_aln_text += header + consensus_seq[((int(farthes_left_position) - 50) - int(consensus_leftmost_coordinate)):((int(farthes_right_position) + 50) - int(consensus_leftmost_coordinate))] + '\n'
                        print('fa_aln_text: ' + fa_aln_text)
                        consensus_file_open_for_reading.close()
                    else:
                        print('the following was not fulfilled (so the program went with the "else" option): ... consensus_ENSTs_path != [] and consensus_ENSTs_path != ""')
                        fa_aln_text += '>' + experiment.lower() + ':no_consensus_seq\n' + ((farthes_right_position - farthes_left_position) * '-') + '\n'
                    
                except KeyError:
                    print('encountered KeyError ... check above whether HRA_mapping_analyzer_dict[' + cluster_ENST + "|consensus_ENSTs_path'] was sucessfuly retrieved ... if yes, then check if HRA_mapping_analyzer_dict[" + cluster_ENST + "|consensus_ENST_leftmost_coordinate'] was successully retrieved")
                    #TODO: is this obsolete? -> fa_aln_text += '>' + experiment + ':no_consensus_seq\n' + ((farthes_right_position - farthes_left_position) * '-') + '\n'
                    fa_aln_text += '>' + experiment.lower() + ':no_consensus_seq\n' + ((farthes_right_position - farthes_left_position) * '-') + '\n'
                    
                del HRA_mapping_analyzer_dict
                HRA_mapping_analyzer_dict_path_open_for_reading.close()
                
            elif experiment not in experiment_log and experiment == 'REF':
                print('experiment: ' + experiment)
                print('experiment_log: ' + str(experiment_log))
                print('experiment is not in experiment_log but it is indeed "REF"')
                print('experiment folder for ' + experiment_list[1] + ' will be used as random source for reference retrieval')
                fasta_REF_path = working_directory + experiment_list[1] + '/' + ref_name + cluster_chromosome + '.fa'
                print('the following fasta REF path will be opened: ' + fasta_REF_path)
                ref_open_for_reading = open(fasta_REF_path, 'r')
                    
                bed_pickle_path = ''
                    
                if (cluster_chromosome != 'Y') and (cluster_chromosome != 'MT'):
                    bed_pickle_path = working_directory + experiment_list[1] + '/bed_pickles/' + bed_name + '_P_' + cluster_chromosome + '_' + cluster_pre_dec + '_ENST_split.txt'
                elif (cluster_chromosome == 'Y') or (cluster_chromosome == 'MT'):
                    for ones_pre_dec in ['(1|2)','(3|4)','(5|6)','(7|8)','(9|0)']:
                        if cluster_pre_dec in ones_pre_dec:
                            bed_pickle_path = working_directory + experiment_list[1] + '/bed_pickles/' + bed_name + '_P_' + cluster_chromosome + '_' + ones_pre_dec + '_ENST_split.txt'       
                print('the following bed_pickle_path will be opened: ' + bed_pickle_path)
                bed_pickle_path_open_for_loading = open(bed_pickle_path, 'rb')
                HRA_bed_parser_dict = pickle.load(bed_pickle_path_open_for_loading)
                print('HRA_bed_parser_dict succesfully loaded from pickle')
                retrieved_Chromosome, retrieved_Start, retrieved_Stop, retrieved_strand = HRA_bed_parser_dict[cluster_ENST + '|bed_data']
                print('experiment_list: ' + str(experiment_list))
                fasta_parser_output = fasta_parser(directory = working_directory + experiment_list[1] + '/', bed = bed_name, fasta_file_object_opened_for_reading = ref_open_for_reading, ChrOI = cluster_chromosome, ones_pre_dec = cluster_pre_dec, HRA_bed_parser_dict = HRA_bed_parser_dict, GRCh37_cds = False, GRCh37_dna = True, specific_enst = cluster_ENST)
                print('fasta_parser_output: ' + str(fasta_parser_output))
                REF_ENST_seq = fasta_parser_output[0][1]
                print("the following are the ref seq's coordinates: ((" +  str(farthes_left_position) + ' - ' + str(retrieved_Start) + ') - (50 - ' + str() + ')) : (' + str(farthes_right_position) + ' - ' + str(retrieved_Start)) + ') + (50 - ' + str() + '))'
                header = '>experiment_label:ref_repeat_start:' + str(farthes_left_position) + '_repeat_stop:' + str(farthes_right_position) + '_\n'
                if ((int(farthes_left_position) - int(retrieved_Start)) - 50) < 0 and ((int(farthes_right_position) - int(retrieved_Start)) + 50) > len(REF_ENST_seq):
                    fa_aln_text += header + REF_ENST_seq + '\n'
                    list_in_case_of_REF_not_being_in_the_df = [(50 - (int(farthes_left_position) - int(retrieved_Start))), (((int(farthes_right_position) - int(retrieved_Start)) + 50) - len(REF_ENST_seq)), int(farthes_left_position), int(farthes_right_position)]
                elif ((int(farthes_left_position) - int(retrieved_Start)) - 50) < 0:
                    fa_aln_text += header + REF_ENST_seq[:((int(farthes_right_position) - int(retrieved_Start)) + 50)] + '\n'
                    list_in_case_of_REF_not_being_in_the_df = [(50 - (int(farthes_left_position) - int(retrieved_Start))), 0, int(farthes_left_position), int(farthes_right_position)]
                elif ((int(farthes_right_position) - int(retrieved_Start)) + 50) > len(REF_ENST_seq):
                    fa_aln_text += header + REF_ENST_seq[((int(farthes_left_position) - int(retrieved_Start)) - 50):] + '\n'
                    list_in_case_of_REF_not_being_in_the_df = [0, (((int(farthes_right_position) - int(retrieved_Start)) + 50) - len(REF_ENST_seq)), int(farthes_left_position), int(farthes_right_position)]
                else:
                    fa_aln_text += header + REF_ENST_seq[((int(farthes_left_position) - int(retrieved_Start)) - 50):((int(farthes_right_position) - int(retrieved_Start)) + 50)] + '\n'
                    list_in_case_of_REF_not_being_in_the_df = [0, 0, int(farthes_left_position), int(farthes_right_position)]

                bed_pickle_path_open_for_loading.close()
                print('fa_aln_text: ' + fa_aln_text)
                                
        fa_path = working_directory + cluster_ENST + '.cc' + str(cluster_counter) + '.fa'
        opened_path = open(fa_path, 'w')
        opened_path.write(fa_aln_text)
        opened_path.close()
        fa_aln_path = ''
            
        #create the alignment
        if aligner == 'mafft':
            #https://mafft.cbrc.jp/alignment/software/anysymbol.html
            fa_aln_path = working_directory + cluster_ENST + '.cc' + str(cluster_counter) + '.mafft.nc.fa.aln'
            try:
                cmd_line = 'mafft --nuc ' + fa_path + ' > ' + fa_aln_path
                return_code = subprocess.call(cmd_line, shell=True)
                if return_code != 0:
                    exception_string = 'Child was terminated by signal; return_code = ' + str(return_code)
                    raise Exception(exception_string)
                else:
                    return_zero_string = 'Child returned; return_code = ' + str(return_code)
                    print(return_zero_string)
            except OSError as e:
                exception_string = 'Execution failed; Error message = ' + str(e)
                raise Exception(exception_string)
    
        elif aligner == 'muscle':
            #https://www.drive5.com/muscle/muscle.html#_Toc81224840
            fa_aln_path = working_directory + cluster_ENST + '.cc' + str(cluster_counter) + '.muslce.nc.fa.aln'
            muscle_cline = MuscleCommandline(input=fa_path, out=fa_aln_path, seqtype='nucleo')
            print(muscle_cline)
            stdout, stderr = muscle_cline()
                
        alignment_open = open(fa_aln_path, 'r')
        
        #this creates a file that contains the alignment without linebreaks throughout the sequence, but maintaining the general FASTA format
        read_aln_txt = ''
        line = alignment_open.readline()
        while line != '':
            if line.startswith('>'):
                #example of a header line ">ENST00000361390.2 mt_genbank_import:known chromosome:GRCh37:MT:3307:4262:1 gene:ENSG00000198888.2 [...]"
                #as explained here: http://may2012.archive.ensembl.org/info/website/tutorials/module3_feb2009_ensembl.pdf
                read_aln_txt += '\n' + line
                line = alignment_open.readline()
                    
                while (not line.startswith('>') and line != ''):
                    read_aln_txt += line.rstrip('\n')
                    line = alignment_open.readline()
            
        read_aln_txt = read_aln_txt.lstrip('\n')
        no_line_breaks_path = fa_aln_path + '.fa.aln.no_linebreaks'
        no_line_breaks_opened = open(no_line_breaks_path, 'w')
        no_line_breaks_opened.write(read_aln_txt)
        no_line_breaks_opened.close()
        alignment_open.close

        #Reminder-example
        #>experiment_label:REF_Chromosome:12_strand:+_ENSG:ENSG00000065970.4_ENST:ENST00000162391.3_ENST_start:8185299_ENST_stop:8208100_repeat_length:4_repeat_start:8200531_repeat_stop:8200542_coordinate_helper_string:lflank=TTTCTTCTCTCCTGGGGGACATCCCACCCTCGAACAAC
        #these margins will be needed to verify the correct repeat among several hits from the same ENST
        REF_repeat_start_in_chr = 'NA'
        REF_repeat_stop_in_chr = 'NA'
        
        #this creates a list that contains tuples of the alignment without linebreaks throughout the sequence, but maintaining the general FASTA format
        alignment_open = open(fa_aln_path, 'r')
        re_read_aln_txt = ''
        line = alignment_open.readline()
        while line != '':
            if line.startswith('>'):
                print(line)
                if ':no_consensus_seq' in line or ':seq_retrieved_from_consensus_seq' in line:
                    line = alignment_open.readline()
                    while (not line.startswith('>') and line != ''):
                        line = alignment_open.readline()
                elif line.startswith('>experiment_label:REF') or line.startswith('>experiment_label:ref'):
                    REF_repeat_start_in_chr = int(line.split('_')[3].split(':')[1])
                    #REF_repeat_stop_in_chr = int(line.split('_')[5].split(':')[1])
                    label = line.split('_')[1].split(':')[1]
                    re_read_aln_txt += '\n' + label + ':'
                    line = alignment_open.readline()
                    while (not line.startswith('>') and line != ''):
                        re_read_aln_txt += line.rstrip('\n')
                        line = alignment_open.readline()
                else:
                    label = line.split('_')[1].split(':')[1]
                    re_read_aln_txt += '\n' + label + ':'
                    line = alignment_open.readline()
                    while (not line.startswith('>') and line != ''):
                        re_read_aln_txt += line.rstrip('\n')
                        line = alignment_open.readline()
        
        re_read_list = re_read_aln_txt.lstrip('\n').rstrip('\n').split('\n')
        alignment_lot = []
        for line in re_read_list:
            label = line[:4]
            original_seq = line[4:]
            gap_free_seq = original_seq.replace('-', '')
            alignment_lot.append((label, original_seq, gap_free_seq))
                
        print('alignment_lot: ' + str(alignment_lot))
        
        #this parsed exon_dict or "gff3_dict" is needed for checking once more the whole aligned sequences' and the homorepeats' phase and exon-location
        gff3_path = working_directory + experiment_list[1] + '/' + gff3 + cluster_chromosome + '.gff3'
        print('the following gff3_path will be opened: ' + gff3_path)
        gff3_dict = GFF3_parser(gff3_path)        
        gff3_CDS_list = gff3_dict[cluster_ENST]
        
        #TODO: there must be a better way to extract the strand here, correct?
        cluster_strand = ''
        if '_strand:+_' in fa_aln_text:
            cluster_strand = '+'
        elif '_strand:-_' in fa_aln_text:
            cluster_strand = '-'
        
        #initiate these variables
        by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed = None
        by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed = None
        left_query_coordinate = None
        right_query_coordinate = None
        
        #It isn't pretty, but the df can be used to retrieve some data that will help in the following analysis
        start_curtail_containing_list = [string for string in df['by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed'][df['ENST']==cluster_ENST][df['experiment_label']=='REF'][df['repeat_start']==REF_repeat_start_in_chr]]
        print('start_curtail_containing_list: ' + str(start_curtail_containing_list))
        if len(start_curtail_containing_list) == 1:
            by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed = start_curtail_containing_list[0]
            print('by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed: ' + str(by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed))
        elif len(start_curtail_containing_list) > 1:
            raise Exception('get_statistics(ERROR): start_containing_list contains more than one helper strings')
        elif len(start_curtail_containing_list) == 0:
            start_curtail_containing_list = False

        stop_curtail_containing_list = [string for string in df['by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed'][df['ENST']==cluster_ENST][df['experiment_label']=='REF'][df['repeat_start']==REF_repeat_start_in_chr]]
        print('stop_curtail_containing_list: ' + str(stop_curtail_containing_list))        
        if len(stop_curtail_containing_list) == 1:
            by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed = stop_curtail_containing_list[0]
            print('by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed: ' + str(by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed))            
        elif len(stop_curtail_containing_list) > 1:
            raise Exception('get_statistics(ERROR): start_containing_list contains more than one helper strings')
        elif len(stop_curtail_containing_list) == 0:
            stop_curtail_containing_list = False
        
        start_containing_list = [string for string in df['repeat_start'][df['ENST']==cluster_ENST][df['experiment_label']=='REF'][df['repeat_start']==REF_repeat_start_in_chr]]
        print('start_containing_list: ' + str(start_containing_list))
        if len(start_containing_list) == 1:
            left_query_coordinate = start_containing_list[0] - (50 - by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed)
            print('left_query_coordinate: ' + str(left_query_coordinate))
        elif len(start_containing_list) > 1:
            raise Exception('get_statistics(ERROR): start_containing_list contains more than one helper strings')
        elif len(start_containing_list) == 0:
            start_containing_list = False

        stop_containing_list = [string for string in df['repeat_stop'][df['ENST']==cluster_ENST][df['experiment_label']=='REF'][df['repeat_start']==REF_repeat_start_in_chr]]
        print('stop_containing_list: ' + str(stop_containing_list))
        if len(stop_containing_list) == 1:
            right_query_coordinate = stop_containing_list[0] + (50 - by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed)
            print('right_query_coordinate: ' + str(right_query_coordinate))
        elif len(stop_containing_list) > 1:
            raise Exception('get_statistics(ERROR): stop_containing_list contains more than one helper strings')
        elif len(stop_containing_list) == 0:
            stop_containing_list = False
        
        if start_curtail_containing_list == False:
            print('start_curtail_containing_list == False')
            by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed = list_in_case_of_REF_not_being_in_the_df[0]
        if stop_curtail_containing_list == False:
            print('stop_curtail_containing_list == False')
            by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed = list_in_case_of_REF_not_being_in_the_df[1]
        if start_containing_list == False:
            print('start_containing_list == False')
            left_query_coordinate = list_in_case_of_REF_not_being_in_the_df[2] - (50 - by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed)
        if stop_containing_list == False:
            print('stop_containing_list == False')
            right_query_coordinate = list_in_case_of_REF_not_being_in_the_df[3] + (50 - by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed)


        phase_string = ''
        #phase_line = ''
        REF_seq = ''
        
        #build a string that helps to understand the alignment's suggested reading frame and that will help with the following step of identifying the coding repeats
        for label, original_seq, gap_free_seq in alignment_lot:
            if label.upper() == 'REF:':
                print('REF found!')
                REF_seq = gap_free_seq
                for CDS in gff3_CDS_list:
                    #sub_feature == coding exon == CDS
                    
                    CDS_range_start = CDS[0]
                    CDS_range_stop = CDS[1]
                    sub_feature_phase = int(CDS[2])
                    print('CDS_range_start: ' + str(CDS_range_start))
                    print('CDS_range_stop: ' + str(CDS_range_stop))
                    print('sub_feature_phase: ' + str(sub_feature_phase))
            
                    if check_if_in_cds(ENST_of_interest = cluster_ENST, left_query_coordinate = left_query_coordinate, right_query_coordinate = right_query_coordinate, CDS_range_start = CDS_range_start, CDS_range_stop = CDS_range_stop) == True:
                        print('coding sequence found!')
                        initial_left_query_coordinate = left_query_coordinate
                        initial_right_query_coordinate = right_query_coordinate
                        print('left_query_coordinate: ' + str(left_query_coordinate))
                        print('CDS_range_start: ' + str(CDS_range_start))
                        print('right_query_coordinate: ' + str(right_query_coordinate))
                        print('CDS_range_stop: ' + str(CDS_range_stop))
                        if (left_query_coordinate < CDS_range_start) and (right_query_coordinate > CDS_range_stop):
                            coded_seq = REF_seq[(CDS_range_start - left_query_coordinate):(CDS_range_stop - left_query_coordinate)]
                            left_query_coordinate = CDS_range_start
                            right_query_coordinate = CDS_range_stop
                            print('outcome: left_query_coordinate < CDS_range_start) and (right_query_coordinate > CDS_range_stop')
                        elif left_query_coordinate < CDS_range_start:
                            coded_seq = REF_seq[(CDS_range_start - left_query_coordinate):]
                            left_query_coordinate = CDS_range_start
                            print('outcome: left_query_coordinate < CDS_range_start')
                        elif right_query_coordinate > CDS_range_stop:
                            coded_seq = REF_seq[:(CDS_range_stop - left_query_coordinate)]
                            right_query_coordinate = CDS_range_stop
                            print('outcome: right_query_coordinate > CDS_range_stop')
                        
                        if cluster_strand == '+':
                            phase_string = ''
                            phase_of_left_query_coordinate = ((((left_query_coordinate - CDS_range_start) % 3) + sub_feature_phase) % 3)
                            for step in range(0,((right_query_coordinate - left_query_coordinate) + 1),1):
                                phase_string += str((step + phase_of_left_query_coordinate) % 3)
                            phase_string = (' ' * (left_query_coordinate - initial_left_query_coordinate)) + phase_string + (' ' * (initial_right_query_coordinate - right_query_coordinate))
                        elif cluster_strand == '-':
                            phase_string = ''
                            phase_of_right_query_coordinate = ((((CDS_range_stop - right_query_coordinate) % 3) + sub_feature_phase) % 3)
                            print('creating the phase string')
                            print('left_query_coordinate: ' + str(left_query_coordinate))
                            print('initial_left_query_coordinate: ' + str(initial_left_query_coordinate))
                            print('(left_query_coordinate - initial_left_query_coordinate)): ' + str((left_query_coordinate - initial_left_query_coordinate)))
                            print('right_query_coordinate: ' + str(right_query_coordinate))
                            print('initial_right_query_coordinate: ' + str(initial_right_query_coordinate))
                            print('(initial_right_query_coordinate - right_query_coordinate): ' + str((initial_right_query_coordinate - right_query_coordinate)))
                            print('(right_query_coordinate - left_query_coordinate): ' + str((right_query_coordinate - left_query_coordinate)))
                            print('len(REF_seq): ' + str(len(REF_seq)))
                            print('REF_seq: ' + str(REF_seq))
                            print('len(coded_seq): ' + str(len(coded_seq)))
                            print('coded_seq: ' + str(coded_seq))
                            print('(' ' * (left_query_coordinate - initial_left_query_coordinate)): ' + str((' ' * (left_query_coordinate - initial_left_query_coordinate))))
                            print('phase_string[::-1]: ' + str(phase_string[::-1]))
                            print('(' ' * (initial_right_query_coordinate - right_query_coordinate)): ' + str((' ' * (initial_right_query_coordinate - right_query_coordinate))))
                            for step in range(0,((right_query_coordinate - left_query_coordinate) + 1),1):
                                phase_string += str((step + phase_of_right_query_coordinate) % 3)
                            phase_string = (' ' * (left_query_coordinate - initial_left_query_coordinate)) + phase_string[::-1] + (' ' * (initial_right_query_coordinate - right_query_coordinate))
                        print('phase_string: ' + str(phase_string))
                        break
        #next correct the cases of the alignment with uppercase indicating coding repeat sequences and lowercase indicating 
        new_txt = ''
                    
        for label, original_seq, gap_free_seq in alignment_lot:
            if label == 'REF:':
                if cluster_strand == '+':
                    pattern = re.finditer(r'((CA[AG])+.{3})*(CA[AG]){4,}(.{3}(CA[AG])+)*', gap_free_seq, flags = re.I)
                elif cluster_strand == '-':
                    pattern = re.finditer(r'(([TC]TG)+.{3})*([TC]TG){4,}(.{3}([TC]TG)+)*', gap_free_seq, flags = re.I)
                former_stop = 0
                case_corrected_txt = ''
                for repeat in pattern:
                    repeat_start = repeat.start()
                    #this too is to counter pythons 0-based ranges
                    repeat_stop = (repeat.end() - 1)
                    #this checks for coding repeats
                    print('cluster_strand: ' + str(cluster_strand))
                    print('phase_string: ' + str(phase_string))
                    print('repeat_start: ' + str(repeat_start))
                    print('repeat_stop: ' + str(repeat_stop))
                    if (cluster_strand == '+' and phase_string[repeat_start] == '0') or (cluster_strand == '-' and phase_string[repeat_stop] == '0'):
                        case_corrected_txt += gap_free_seq[former_stop:repeat_start].lower() + gap_free_seq[repeat_start:(repeat_stop + 1)].upper()
                        former_stop = (repeat_stop + 1)
                case_corrected_txt += gap_free_seq[former_stop:].lower()
                print('case_corrected_txt REF: ' + case_corrected_txt)
                
                list_of_gaps = []
                gap_pattern = re.finditer(r'-+', original_seq, flags = re.I)
                last_stop = 0
                for gap in gap_pattern:
                    prior_length = (gap.start() - last_stop)
                    gap_seq = gap.group()
                    last_stop = gap.end()
                    list_of_gaps.append((prior_length, gap_seq))
                    
                REF_reconstruction = ''
                last_stop = 0
                for prior_length, gap_seq in list_of_gaps:
                    REF_reconstruction += case_corrected_txt[:(prior_length)] + gap_seq
                    case_corrected_txt = case_corrected_txt[prior_length:]
                REF_reconstruction += case_corrected_txt
                print('REF_reconstruction: ' + REF_reconstruction)
                
                ORF_string_copy = phase_string
                print('ORF_string_copy: ' + ORF_string_copy)
                ORF_gap_construction = ''
                last_stop = 0
                for prior_length, gap_seq in list_of_gaps:
                    ORF_gap_construction += ORF_string_copy[:(prior_length)] + gap_seq
                    ORF_string_copy = ORF_string_copy[prior_length:]
                ORF_gap_construction += ORF_string_copy
                print('ORF_gap_construction: ' + ORF_gap_construction)
                
                new_txt += 'ORF:' + ORF_gap_construction + '\nREF:' + REF_reconstruction + '\n'
            
        for label, original_seq, gap_free_seq in alignment_lot:
            if label != 'REF:':
                print(label)
                if cluster_strand == '+':
                    pattern = re.finditer(r'((CA[AG])+.{3})*(CA[AG]){4,}(.{3}(CA[AG])+)*', gap_free_seq, flags = re.I)
                elif cluster_strand == '-':
                    pattern = re.finditer(r'(([TC]TG)+.{3})*([TC]TG){4,}(.{3}([TC]TG)+)*', gap_free_seq, flags = re.I)
                former_stop = 0
                case_corrected_txt = ''
                for repeat in pattern:
                    repeat_start = repeat.start()
                    #this too is to counter pythons 0-based ranges
                    repeat_stop = (repeat.end() - 1)
                    #this checks for coding repeats
                    print(repeat_start)
                    print(repeat_stop)
                    print(phase_string)
                    if (cluster_strand == '+' and phase_string[repeat_start] == '0') or (cluster_strand == '-' and phase_string[repeat_stop] == '0'):
                        case_corrected_txt += gap_free_seq[former_stop:repeat_start].lower() + gap_free_seq[repeat_start:(repeat_stop + 1)].upper()
                        former_stop = (repeat_stop + 1)
                    
                case_corrected_txt += gap_free_seq[former_stop:].lower()
                print('case_corrected_txt: ' + case_corrected_txt)
                
                list_of_gaps = []
                gap_pattern = re.finditer(r'-+', original_seq, flags = re.I)
                last_stop = 0
                for gap in gap_pattern:
                    prior_length = (gap.start() - last_stop)
                    gap_seq = gap.group()
                    last_stop = gap.end()
                    list_of_gaps.append((prior_length, gap_seq))
                    
                reconstruction = ''
                last_stop = 0
                for prior_length, gap_seq in list_of_gaps:
                    reconstruction += case_corrected_txt[:(prior_length)] + gap_seq
                    case_corrected_txt = case_corrected_txt[prior_length:]
                reconstruction += case_corrected_txt
                new_txt += label + reconstruction + '\n'
                print('reconstruction: ' + case_corrected_txt)
        
        print('new_txt: ' + new_txt)
        new_path = fa_aln_path + '.fa.aln.phased_case_change'
        x = open(new_path, 'w')
        x.write(new_txt)
        x.close()
        alignment_open.close()
        translation_dict = {}
        translation_lot = []
        list_of_case_changed_lines = new_txt.split('\n')
        for line in list_of_case_changed_lines:
            name, seq = line.split(':')
            translation_dict[name] = ''    
            if cluster_strand == '+':
                if name != 'ORF':
                    translation_lot.append((name, seq.upper()))
                else: 
                    orf = seq
            elif cluster_strand == '-':
                if name != 'ORF':
                    translation_lot.append((name, seq[::-1].upper()))
                else: 
                    orf = seq
                    
        character_counter = -1
        for phase_counter in orf:
            character_counter += 1
            if (character_counter - phase_counter) < 0:
                for name, seq in translation_lot:
                    translation_dict[str(phase_counter) + '|' + name] = 'X'
            else:
                    
                if phase_counter == '0':
                    for name, seq in translation_lot:
                        translation_dict[str(phase_counter) + '|' + name] = seq[character_counter]
                if phase_counter == '1':
                    for name, seq in translation_lot:
                        translation_dict[str(phase_counter) + '|' + name] = seq[character_counter - 1] + seq[character_counter]
                if phase_counter == '2':
                    for name, seq in translation_lot:
                        translation_dict[str(phase_counter) + '|' + name] = seq[character_counter - 2] + seq[character_counter - 1] + seq[character_counter]
                        if check_if_list_elements_not_in_string(string = translation_dict[str(phase_counter) + '|' + name], checklist = ['-','N','R','Y','W','S','M','K','X','H','B','D','V']) == True:
                            new_AA = AAS_table[translation_dict[str(phase_counter) + '|' + name]]
                        else:
                            new_AA = 'X'
                        translation_dict[name] = translation_dict[name] + new_AA
        if orf[-1] == '2':
            pass
        else:
            for name, seq in translation_lot:
                translation_dict[name] = translation_dict[name] + 'X'
        AA_txt = ''
        for name, seq in translation_lot:
            AA_txt += name + ':' + translation_dict[name] + '\n'

        print('AA_txt: ' + AA_txt)
        new_path = fa_aln_path + '.translated_AA.aln'
        x = open(new_path, 'w')
        x.write(AA_txt)
        x.close()
                            

    print('aborting')    
    return True


###############################################################################
#end of comparison processing functions that return statistical data on the results
###############################################################################



###############################################################################
#start of major module-functions callable via cmd-agrs
###############################################################################

#in my thesis referred to as Function 1
#@profile
def HRA_bed_parser(working_directory, bed_name, ChrOI_list, overwrite_old_files):
    
    if (os.path.exists(working_directory + 'bed_pickles/') == False):
        os.mkdir(working_directory + 'bed_pickles/')
    else:
        if overwrite_old_files == 'True':
            shutil.rmtree(working_directory + 'bed_pickles/')
            os.mkdir(working_directory + 'bed_pickles/')
        else:
            raise Exception('HRA_bed_parser(ERROR): directory ' + working_directory + 'bed_pickles/ already exists and "overwrite" is not set to "True" ... ')
    
    if os.path.exists(working_directory + 'bam_pickles/') == False:
        os.mkdir(working_directory + 'bam_pickles/')
    else:
        if overwrite_old_files == 'True':
            shutil.rmtree(working_directory + 'bam_pickles/')
            os.mkdir(working_directory + 'bam_pickles/')
        else:
            raise Exception('HRA_bed_parser(ERROR): directory ' + working_directory + 'bam_pickles/ already exists and "overwrite" is not set to "True" ... ')
    
    if (os.path.exists(working_directory + 'indexes/') == False):
        os.mkdir(working_directory + 'indexes/')
    else:
        if overwrite_old_files == 'True':
            shutil.rmtree(working_directory + 'indexes/')
            os.mkdir(working_directory + 'indexes/')
        else:
            raise Exception('HRA_bed_parser(ERROR): directory ' + working_directory + 'indexes/ already exists and "overwrite" is not set to "True" ... ')
    
    if (os.path.exists(working_directory + 'copied_cds_fastas/') == False):
        os.mkdir(working_directory + 'copied_cds_fastas/')
    else:
        if overwrite_old_files == 'True':
            shutil.rmtree(working_directory + 'copied_cds_fastas/')
            os.mkdir(working_directory + 'copied_cds_fastas/')
        else:
            raise Exception('HRA_bed_parser(ERROR): directory ' + working_directory + 'copied_cds_fastas/ already exists and "overwrite" is not set to "True" ... ')

    if (os.path.exists(working_directory + 'analyze_pickles/') == False):
        os.mkdir(working_directory + 'analyze_pickles/')
    else:
        if overwrite_old_files == 'True':
            shutil.rmtree(working_directory + 'analyze_pickles/')
            os.mkdir(working_directory + 'analyze_pickles/')
        else:
            raise Exception('HRA_bed_parser(ERROR): directory ' + working_directory + 'analyze_pickles/ already exists and "overwrite" is not set to "True" ... ')

    if (os.path.exists(working_directory + 'consensus_ENSTs/') == False):
        os.mkdir(working_directory + 'consensus_ENSTs/')
    else:
        if overwrite_old_files == 'True':
            shutil.rmtree(working_directory + 'consensus_ENSTs/')
            os.mkdir(working_directory + 'consensus_ENSTs/')
        else:
            raise Exception('HRA_bed_parser(ERROR): directory ' + working_directory + 'consensus_ENSTs/ already exists and "overwrite" is not set to "True" ... ')

    if (os.path.exists(working_directory + 'raw_consensus_matrices/') == False):
        os.mkdir(working_directory + 'raw_consensus_matrices/')
    else:
        if overwrite_old_files == 'True':
            shutil.rmtree(working_directory + 'raw_consensus_matrices/')
            os.mkdir(working_directory + 'raw_consensus_matrices/')
        else:
            raise Exception('HRA_bed_parser(ERROR): directory ' + working_directory + 'raw_consensus_matrices/ already exists and "overwrite" is not set to "True" ... ')

    if (os.path.exists(working_directory + 'csv_results/') == False):
        os.mkdir(working_directory + 'csv_results/')
    else:
        if overwrite_old_files == 'True':
            shutil.rmtree(working_directory + 'csv_results/')
            os.mkdir(working_directory + 'csv_results/')
        else:
            raise Exception('HRA_bed_parser(ERROR): directory ' + working_directory + 'csv_results/ already exists and "overwrite" is not set to "True" ... ')


    #ChrOI = ChromosomeOfInterest
    #ones_pre_dec = last pre-decimal digit
    for ChrOI in ChrOI_list:
        pre_dec_list = []
        
        if (ChrOI != 'Y' and ChrOI != 'MT'):
            pre_dec_list = ['1','2','3','4','5','6','7','8','9','0']
        else:
            pre_dec_list = ['(1|2)','(3|4)','(5|6)','(7|8)','(9|0)']
        
        
        for ones_pre_dec in pre_dec_list:
            print('input chromosome: ' + str(ChrOI) + ' | input ones_pre_dec: ' + str(ones_pre_dec))
                
            proc_ID = 'P_' + ChrOI + '_' + ones_pre_dec
            
            if (os.path.exists(working_directory + bed_name + '_' + proc_ID + '_ENST_split.txt') == True):
                HRA_bed_parser_dict_pickle = open(working_directory + bed_name + '_' + proc_ID + '_ENST_split.txt', 'r')
                print('Pre-existing pickle.dump() was found under the following path: ' + working_directory + bed_name + '_' + proc_ID + '_ENST_split.txt')
                HRA_bed_parser_dict = pickle.load(HRA_bed_parser_dict_pickle)
                return HRA_bed_parser_dict   
            
            elif os.path.exists(working_directory + bed_name + '_' + proc_ID + '_ENST_split.txt') == False:
                print('No pre-existing pickle.dump() could be found under the following path: ' + working_directory + bed_name + '_' + proc_ID + '_ENST_split.txt ... new pickle.dump() will be created')
        
            #### BED FILE PROCESSING ####
            bed_file_location = working_directory + bed_name
        
            if (os.path.exists(bed_file_location)):
                print('working on bed file: ')
                print(bed_file_location)
            else:
                print('cannot find bed file')
                print(bed_file_location)
            
            bed_lol = (line.rstrip('\n').split('\t') for line in open(bed_file_location))
            #if gene entry exists, append row information, if gene entry 
            #doesn't exist "defaultdict(list) + append" will initiate a new 
            #instance and append the first row
            HRA_bed_parser_dict = defaultdict(list)
            ENSGs = []
            ENSTs = []
            
            #creating the following key series: 
            #HRA_bed_parser_dict[ENST|fasta_path]
            #HRA_bed_parser_dict[ENST|bed_data]
            #HRA_bed_parser_dict[Chromosome|ENSG|bed_row_for_ranges]
            #HRA_bed_parser_dict[Chromosome|ENSGxyz]
            #HRA_bed_parser_dict[Chromosome]
            
            for bed_row in bed_lol:
                Chromosome = bed_row[0]
                Start = int(bed_row[1])
                Stop = int(bed_row[2])
                Strand = bed_row[3]
                ENST = bed_row[4]
                ENSG = bed_row[5]
                
                if Chromosome == ChrOI:
                    if predecimal_ones_check(ones_pre_dec, ENSG) == True:
                        #TODO: this causes the dict to return a needless lol such as [['ensg']], no list at all would be favorable, but if that fails at least the outer brackets can easily be removed
                        HRA_bed_parser_dict[ENST + '|find_ENSG'].append([ENSG])
                        HRA_bed_parser_dict[ENST + '|bed_data'] = (str(Chromosome), str(Start), str(Stop), str(Strand))
                        #The list below could contain each bed row as follows: 
                        #[[Chromosome, start, end, strand, ENSTxyz.1, ENSGxyz,], [ etc,]]
                        #but currently its only use is for creating the range keys via iteration 
                        #which is why the full row's content in form of a dict is hashed out:
                        #HRA_bed_parser_dict[ENSG + '|' 'bed_row'].append([Chromosome, Start, Stop, Strand, ENST, ENSG])
                        #check: check done with random ENSGs and all looks good
                        HRA_bed_parser_dict[Chromosome + '|' + ENSG + '|bed_rows_for_ranges'].append([ENST, Start, Stop])
                        
                        #check: checked for length with bed_file_location = '/home/lw/Desktop/GRCh37_protcodgenes/GRCh37_protcodgenes_info3.bed'
                        #all looks good: content is a long list of ENSGs, and len(ENSGs) = 20151
                        ensg_line = [Chromosome, ENSG]
                        enst_line = [Chromosome, ENST]
                        if ensg_line not in ENSGs:
                            ENSGs.append(ensg_line)
                        if enst_line not in ENSTs:
                            ENSTs.append(enst_line)
            HRA_bed_parser_dict[ChrOI + '|all_ENSTs'] = [element[1] for element in ENSTs]
            #write the ENSG, the most downstream TS-Start and most upstream TS-Stop to a list
            #creating HRA_bed_parser_dict[chromosome|MinMax]
            #this step con only be performed when ENSG|bed_rows_for_ranges has collected all relevant rows
            for chromosome, ENSG in ENSGs:
                bed_rows_for_ranges = HRA_bed_parser_dict[chromosome + '|' + ENSG + '|bed_rows_for_ranges']
                temp_start_list = []
                temp_stop_list = []
                for ENST_coordinates in bed_rows_for_ranges:
                    temp_start_list.append(ENST_coordinates[1])
                    temp_stop_list.append(ENST_coordinates[2])
                HRA_bed_parser_dict[chromosome + '|MinMax'].append([ENSG, int(sorted(temp_start_list)[0]), int(sorted(temp_stop_list, reverse=True)[0])])
    
            HRA_bed_parser_dict_pickle = open(working_directory + 'bed_pickles/' + bed_name + '_' + proc_ID + '_ENST_split.txt', 'wb')
            pickle.dump(HRA_bed_parser_dict, HRA_bed_parser_dict_pickle, protocol=pickle.HIGHEST_PROTOCOL)
            HRA_bed_parser_dict_pickle.close()
        print('HRA_bed_parser() completed')
        print(get_mem('HRA_bed_parser'))
    
    return HRA_bed_parser_dict

#in my thesis referred to as Function 2
#@profile    
def HRA_bam_parser(working_directory, bam_name, ChrOI, overwrite_indexes = False):

    #can run as 25 processes
    
    proc_ID = 'P_' + ChrOI
    HRA_bam_parser_dict = defaultdict(list)
    
    if os.path.exists(working_directory + 'bam_pickles/' + bam_name + '_' + proc_ID + '_HRA_bam_parser.txt') == True:
        HRA_bam_parser_dict_pickle = open(working_directory  + 'bam_pickles/' +  bam_name + '_' + proc_ID + '_HRA_bam_parser.txt', 'r')
        print('Pre-existing pickle.dump() was found ...')
        HRA_bam_parser_dict = pickle.load(HRA_bam_parser_dict_pickle)
        if os.path.exists(working_directory + 'indexes/' + proc_ID + '.txt') == True:
            print('index locations found ...')
            if overwrite_indexes == False:
                print('finishing HRA_bam_parser without parsing the bam file anew')
                return HRA_bam_parser_dict
            elif overwrite_indexes == True:
                print('indexes will be overwritten at given location ... bam file is parsed anew')
        elif os.path.exists(working_directory + 'indexes/' + proc_ID + '.txt') == False:
            print('indexes cannot be found at given location ... bam file is parsed anew')
        
    else:
        print('No pre-existing pickle.dump() was found ... bam has to be parsed anew ... ')
    
    #### BAM FILE PROCESSING ####
    
    bam_file_location = working_directory + bam_name
    print('working on bam file: ')
    print(bam_file_location)
    
    #parse the bam file with the pybam module
    parsed_bam = pybam.read(bam_file_location)
    
    #create a dictionary that containins the paths to text files for all bam reads per chromosome 
    chromosomal_data_path = (working_directory + 'indexes/' + proc_ID + '.txt')
    #prepare path to each chromosome text file in dict for fast retrieval and create txt files
    HRA_bam_parser_dict[ChrOI + '|location'] = chromosomal_data_path
    create_file = open(chromosomal_data_path, 'w')
    create_file.close()
    #write the all relevant alignment information to text file in tab-delim manner 
    #because storing a complete sequencing experiment in memory would be too excessive
    total_count = 0
    chromosome_count = 0
    passed_count = 0
    for alignment in parsed_bam:
        total_count += 1
        chr_rname = str(alignment.sam_rname).replace('chr','').replace('Chr','').replace('CHR','')
        if (chr_rname == ChrOI):
            chromosome_count += 1
            if (int(alignment.sam_mapq) >= 35) and check_if_list_elements_not_in_string(string = str(alignment.sam_qual), checklist = ['!','"','#','$','%','&',"'",'(',')','*','+',',','-','.','/','0','1','2','3','4','5']) == True:
                passed_count += 1
                current_text_file = open(HRA_bam_parser_dict[chr_rname + '|location'], 'a')
                current_text_file.write(chr_rname + '\t' + str(alignment.file_alignments_read) + '\t' + str(alignment.sam_pos1) + '\t' + str(alignment.sam_tlen) + '\t' + str(alignment.sam_mapq) + '\t' + alignment.sam_seq + '\t' + alignment.sam_qual + '\t' + alignment.sam_cigar_string + '\n')
                current_text_file.close()

    HRA_bam_parser_dict_pickle = open(working_directory + 'bam_pickles/' + bam_name + '_' + proc_ID + '_HRA_bam_parser.txt', 'w')
    pickle.dump(HRA_bam_parser_dict, HRA_bam_parser_dict_pickle, protocol=pickle.HIGHEST_PROTOCOL)
    HRA_bam_parser_dict_pickle.close()    
    
    print('Checked a total of ' + str(total_count) + ' mapped reads.\nOf these ' + str(chromosome_count) + ' mapped to the chromosome in question (' + str(chromosome_count / total_count) + '%).\n Of these mapped reads ' + str(passed_count) + ' passed the quality threshold set (' + str(passed_count / chromosome_count) + '%).')
    
    print(get_mem('HRA_bam_parser'))
    
    return HRA_bam_parser_dict

#in my thesis referred to as Function 3
#@profile
def HRA_mapping_analyzer(working_directory, bed_name, bam_name, ChrOI, tens_pre_dec, ones_pre_dec, GFF3_path, introns):
    
    HRA_mapping_analyzer0 = time.time()
    csv_dir = working_directory + 'csv_results/'
    
    if introns != 'exclude' and introns != 'include':
        raise Exception('HRA_mapping_analyzer(ERROR): use the "intron" flag and set it either to "exclude" or "include" to specify how to proceed with intronic regions')
    
    if re.match(r'\d', tens_pre_dec):
        proc_ID = 'P_' + ChrOI + '_' + tens_pre_dec + '_' + ones_pre_dec
    
    elif re.match(r'NA', tens_pre_dec):
        proc_ID = 'P_' + ChrOI + '_' + ones_pre_dec
        tens_pre_dec = r'\d'
        
    else:
        raise Exception('HRA_mapping_analyzer(ERROR): cannot work with tens_pre_dec submitted as ' + tens_pre_dec + ' ... has to be either "NA" or single digit.')
    
    print('starting alignment analyzer')
    
    HRA_mapping_analyzer_dict = defaultdict(list)
    
    bed_pickle_path = working_directory + 'bed_pickles/' + bed_name + '_P_' + ChrOI + '_' + ones_pre_dec + '_ENST_split.txt'
    if os.path.exists(bed_pickle_path) == True:
        HRA_bed_parser_dict_original = pickle.load(open(bed_pickle_path, 'rb'))
        HRA_bed_parser_dict = copy.deepcopy(HRA_bed_parser_dict_original)
        print(proc_ID + ': ' + bed_pickle_path + ' will be used as HRA_bed_parser_dict')
    else:
        print(proc_ID + ': ' + bed_pickle_path + ' is missing ...')
    if os.path.exists(working_directory + 'bam_pickles/' + bam_name + '_P_' + ChrOI + '_HRA_bam_parser.txt') == True:
        HRA_bam_parser_dict_original = pickle.load(open(working_directory + 'bam_pickles/' + bam_name + '_P_' + ChrOI + '_HRA_bam_parser.txt', 'rb'))
        HRA_bam_parser_dict = copy.deepcopy(HRA_bam_parser_dict_original)
    else:
        print(proc_ID + ': ' + working_directory + 'bam_pickles/' + bam_name + '_P_' + ChrOI + '_HRA_bam_parser.txt is missing ...')
    
    csv_repeat_analysis_path = csv_dir + 'csv_analysis_for_chromosome_' + proc_ID + '_polyQs_query.csv'
    csv_repeat_analysis_path_initiate = open(csv_repeat_analysis_path, 'w')
    csv_repeat_analysis_path_initiate.close()
    
    #### Sort reads from BED and BAM according to mapping position ####
        
    print(proc_ID + ': assigning all reads to matching regions (ENSTs) on chromosome ' + ChrOI)
        
        
    unfiltered_current_ranges_lol = HRA_bed_parser_dict[ChrOI + '|MinMax']
    filtered_current_ranges_lol = [range_row for range_row in unfiltered_current_ranges_lol if predecimal_tens_and_ones_check(tens_pre_dec, ones_pre_dec, range_row[0]) == True]
    bam_text_file_path = HRA_bam_parser_dict[ChrOI + '|location']
    bam_lol = (line.rstrip('\n').split('\t') for line in open(bam_text_file_path))
    exon_dict = GFF3_parser(working_directory + GFF3_path)

    if introns == 'exclude':
    
        for bam_row in bam_lol:
                
            ### setting variables
            ### RR = "Read Right", RL = "Read Left"
            bam_running_read_counter = bam_row[1]
            RL = int(bam_row[2])
            bam_read_template_length = int(bam_row[3])
            RR = RL + bam_read_template_length - 1
            #bam_read_mapq = bam_row[4]
            bam_read_seq = bam_row[5]
            #bam_read_alignment_qual = bam_row[6]
            bam_read_cigar = bam_row[7]
                    
            ### find all genes (ENSG ranges) on the cromosome that the current read will map into
            ### GRMin = "Gene Range Minimum", GRMax = "Gene Range Maximum"
            matching_ENSGs = (range_row[0] for range_row in filtered_current_ranges_lol if not (((RR < range_row[1]) and (RL < range_row[1])) or ((RR > range_row[2]) and (RL > range_row[2]))))
            for ENSG in matching_ENSGs:
    
                current_ENST_lol = HRA_bed_parser_dict[ChrOI + '|' + ENSG + '|bed_rows_for_ranges']
                    
                for ENST_row in current_ENST_lol:
                    ENST = ENST_row[0]
                    TL = ENST_row[1]
                    TR = ENST_row[2]
                        
                    if not ((RR < TL) and (RL < TL)) and not ((RR > TR) and (RL > TR)):
                        gff3_CDS_list = exon_dict[ENST]
                            
                        for CDS in gff3_CDS_list:
                            CDS_range_start = CDS[0] 
                            CDS_range_stop = CDS[1]
                                
                            if (not ((RL < CDS_range_start) and (RR < CDS_range_start)) and not ((RL > CDS_range_stop) and (RR > CDS_range_stop))):
                                #####add to ENST dictionary entry for containing all reads to reconstruct the alignment
                                HRA_mapping_analyzer_dict[ENST + '|' + ChrOI].append([bam_running_read_counter, bam_row[2], bam_read_seq, bam_read_cigar])
                                strand = (HRA_bed_parser_dict[ENST + '|bed_data'])[3]
                                    
                                if [ENST, strand] not in HRA_mapping_analyzer_dict[ChrOI + '|all_aligned_ENSTs']:
                                    HRA_mapping_analyzer_dict[ChrOI + '|all_aligned_ENSTs'].append([ENST, strand])
                                    
        HRA_mapping_analyzer1 = time.time()
        print(proc_ID + ': asigning all reads to their relevant ENSTs took ' + str(HRA_mapping_analyzer1 - HRA_mapping_analyzer0) + 's')

    elif introns == 'include':
        
        current_ranges_lol = HRA_bed_parser_dict[ChrOI + '|MinMax']
        for bam_row in bam_lol:
                
            ### setting variables
            ### RR = "Read Right", RL = "Read Left"
            bam_running_read_counter = bam_row[1]
            RL = int(bam_row[2])
            bam_read_template_length = int(bam_row[3])
            RR = RL + bam_read_template_length - 1
            #bam_read_mapq = bam_row[4]
            bam_read_seq = bam_row[5]
            #bam_read_alignment_qual = bam_row[6]
            bam_read_cigar = bam_row[7]
                    
            ### find all genes (ENSG ranges) on the cromosome that the current read will map into
            ### GRMin = "Gene Range Minimum", GRMax = "Gene Range Maximum"
            for range_row in current_ranges_lol:
                ENSG = range_row[0]
                if predecimal_tens_and_ones_check(tens_pre_dec, ones_pre_dec, ENSG) == True:
                    GRMin = range_row[1]
                    GRMax = range_row[2]
                    if not (((RR < GRMin) and (RL < GRMin)) or ((RR > GRMax) and (RL > GRMax))):    
                                    
                        ####find all transcript variants (ENSTs) in the ENSGs' range that the current read will map to
                        ####TL = "Transcript Left" or 'ENST_start', TR = "Transcript Right" or 'ENST_stop'
                        #ENSG = range_row[0]
                        current_ENST_lol = HRA_bed_parser_dict[ChrOI + '|' + ENSG + '|bed_rows_for_ranges']
                        
                        for ENST_row in current_ENST_lol:
                                                
                            ENST = ENST_row[0]
                            TL = ENST_row[1]
                            TR = ENST_row[2]
                            if not ((RR < TL) and (RL < TL)) and not ((RR > TR) and (RL > TR)):
                                                                    
                                #####add to ENST dictionary entry for containing all reads to reconstruct the alignment
                                HRA_mapping_analyzer_dict[ENST + '|' + ChrOI].append([bam_running_read_counter, bam_row[2], bam_read_seq, bam_read_cigar])
                                strand = (HRA_bed_parser_dict[ENST + '|bed_data'])[3]
                                if [ENST, strand] not in HRA_mapping_analyzer_dict[ChrOI + '|all_aligned_ENSTs']:
                                    HRA_mapping_analyzer_dict[ChrOI + '|all_aligned_ENSTs'].append([ENST, strand])


    #### retrieve the alignment ####
    HRA_mapping_analyzer2 = time.time()
    print(proc_ID + ': retrieving alignment information per region (ENST) on chromosome: ' + ChrOI)
    all_passed_ENSTs = HRA_mapping_analyzer_dict[ChrOI + '|all_aligned_ENSTs']
    print('all passed ENSTs: ' + str(all_passed_ENSTs))
    
    for ENST, strand in all_passed_ENSTs:
        #HRA_mapping_analyzer1 = time.time()
        print(proc_ID + ': working on ' + ENST)
            
        #reads are sorted so that the genomic start position of the rightmost read and stop of the leftmost read can be determined as start and stop of the alignment
        reads_sorted_by_pos1 = sorted(HRA_mapping_analyzer_dict[ENST + '|' + ChrOI], key=return2nd)
        #the rightmost position equals the last/rightmost read's starting position + its read length
            
        list_of_retrieved_alignment_matrices_for_parsing = []
        ENSG = HRA_bed_parser_dict[ENST + '|' + 'find_ENSG'][0][0]
            
        for bam_running_read_counter, pos1_str, bam_read_seq, bam_read_cigar  in reads_sorted_by_pos1:
            
            pos1 = int(pos1_str)
            list_of_cigar_elements = re.findall(r'[\d]*[MIDNSHP=X]', bam_read_cigar, flags=0)
            #the first position [0] holds a value that will count the total amount of insertions to the reference
            
            alignment_matrix = [0]                
            reference_position_log = pos1
            read_position_log = pos1
                
            #only the leftmost M/D/N/=/X operator accounts for pos1; see this thread: https://github.com/samtools/hts-specs/issues/80
            #therefor the read position log will have to be adjusted: it will start X bases earlier than pos, with X being the sum of all non-M|D|N|=|X operations
            for element in list_of_cigar_elements:
                count = element.rstrip('MIDNSHP=X')
                operator = element[len(count):]   
                if not re.match(r"(M|D|N|=|X)", operator):
                    read_position_log -= (int(count))
                else:
                    break
                
            read_pos1_including_clipping = read_position_log
            
            for element in list_of_cigar_elements:
                InDel_counter = 0
                count = element.rstrip('MIDNSHP=X')
                operator = element[len(count):]
                
                #turn information stored in the CIGAR string into the aligned read matrix, storing Base, reference_positions, read_position
                iterstart = (read_position_log - read_pos1_including_clipping)
                iterstop = ((read_position_log - read_pos1_including_clipping) + int(count))
                base_iterator = itertools.islice(bam_read_seq, iterstart, iterstop, 1)
                for base in base_iterator:
                    if operator == 'M':
                        alignment_column = (base, reference_position_log, read_position_log, InDel_counter)
                        alignment_matrix.append(alignment_column)
                        reference_position_log += 1
                        read_position_log += 1
                    elif operator == 'I':
                        #insertion to the reference
                        InDel_counter += (1 * -1)
                        alignment_matrix[0] += 1
                        #instead of printing the base, '-' indicated by count, and since reference_position isn't updated the loop doesn't skip any bases from the read
                        alignment_column = ('-', (reference_position_log * -1), read_position_log, InDel_counter)
                        alignment_matrix.append(alignment_column)
                        #reference_position_log += 0
                        read_position_log += 1
                    elif operator == 'D':
                        #deletion to the reference
                        InDel_counter += 1
                        alignment_column = (base, reference_position_log, (read_position_log * -1), InDel_counter)
                        alignment_matrix.append(alignment_column)
                        reference_position_log += 1
                        #read_position_log += 0
                    elif operator == 'S':
                        InDel_counter += (1 * -1)
                        alignment_matrix[0] += 1
                        #instead of printing the base, '-' indicated by count, and since reference_position isn't updated the loop doesn't skip any bases from the read
                        alignment_column = ('-', (reference_position_log * -1), read_position_log, InDel_counter)
                        alignment_matrix.append(alignment_column)
                        #reference_position_log += 0
                        read_position_log += 1
                    elif operator == 'H':
                        read_position_log += 1
                        reference_position_log += 1
                    elif operator == '=' or operator == 'X':
                        alignment_column = (base, reference_position_log, read_position_log, InDel_counter)
                        alignment_matrix.append(alignment_column)
                        reference_position_log += 1
                        read_position_log += 1
                    elif operator == 'N':
                        #an insert that is an intron to the purely exonic reference - how we deal to them depends on wether we use a reference and QRY including introns or not
                        InDel_counter += (1 * -1)
                        alignment_matrix[0] += 1
                        alignment_column = (base, reference_position_log, read_position_log, InDel_counter)
                        alignment_matrix.append(alignment_column)
                        reference_position_log += 1
                        read_position_log += 1
                    elif operator == 'P':
                        #a padded insert to the reference that is not consistent amongst the reads
                        InDel_counter += (1 * -1)
                        alignment_matrix[0] += 1
                        alignment_column = (base, (reference_position_log * -1), read_position_log, InDel_counter)
                        alignment_matrix.append(alignment_column)
                        #reference_position_log += 0
                        read_position_log += 1
            
            alignment = '0' + ''.join([base[0] for base in alignment_matrix[1:]])
            if len(alignment) > 1:
                list_of_retrieved_alignment_matrices_for_parsing.append(alignment_matrix)
                
                    
        if len(list_of_retrieved_alignment_matrices_for_parsing) > 0:
            consensus_matrix = consensus_writer(list_of_retrieved_alignment_matrices_for_parsing)
            HRA_mapping_analyzer_dict[ENST + '|consensus_matrix'] = consensus_matrix
            consensus_seq = consensus_matrix[0]
            raw_consensus_matrix = consensus_matrix[2]
            consensus_leftmost_coordinate = consensus_matrix[1]
            consensus_ENSTs_path = working_directory + 'consensus_ENSTs/' + ENST + '.txt'
            consensus_ENSTs_txt = open(consensus_ENSTs_path, 'w')
            consensus_ENSTs_txt.write(consensus_seq)
            consensus_ENSTs_txt.close()
            HRA_mapping_analyzer_dict[ENST + '|consensus_ENSTs_path'] = consensus_ENSTs_path
            HRA_mapping_analyzer_dict[ENST + '|consensus_ENST_leftmost_coordinate'] = consensus_leftmost_coordinate

            if len(consensus_seq) > 0:
                raw_consensus_matrix_csv = open(working_directory + 'raw_consensus_matrices/' + ENST + '.csv', 'w')
                for position in xrange(0, len(consensus_seq), 1):
                    raw_consensus_csv_line = consensus_seq[position] + ',' + ','.join((str(figure) for figure in raw_consensus_matrix[position])) + '\n'
                    raw_consensus_matrix_csv.write(raw_consensus_csv_line)
                raw_consensus_matrix_csv.close()                    
                repeat_analysis_csv_rows = polyQ_finder(seq_in = consensus_seq, chromosome_in = ChrOI, ENSG_in = ENSG, ENST_in = ENST, strand_in = strand, starting_coordinate = consensus_leftmost_coordinate, exon_dict = exon_dict)
                processed_cds_fasta_path = working_directory + 'copied_cds_fastas/' + proc_ID + '_cds_copy.fasta'
                processed_cds_fasta = open(processed_cds_fasta_path, 'w')
                csv_repeat_analysis = open(csv_repeat_analysis_path, 'a')
                #I wasn't sure if you can get the length of a given regex by using len(), so I just hardcoded the 3 ... can't imagine a reason to step out of the triplet code anyway
                csv_repeat_analysis.write(repeat_analysis_csv_rows)
                csv_repeat_analysis.close()
                processed_cds_fasta.close()
                    
        else:
            print('no alignment could be retrieved for ' + ENST)
            
    HRA_mapping_analyzer3 = time.time()
    print(proc_ID + ': finishing repeat extraction from chromosome ' + proc_ID + ' took ' + str(HRA_mapping_analyzer3 - HRA_mapping_analyzer2) + 's')
                    
    new_pickle = open(working_directory + 'analyze_pickles/' + proc_ID + '_align_pickle.txt', 'w')
    pickle.dump(HRA_mapping_analyzer_dict, new_pickle, protocol=pickle.HIGHEST_PROTOCOL)
    new_pickle.close()

    if re.match(r'\d', tens_pre_dec):
        print(proc_ID + ': done with HRA_mapping_analyzer() for tchromosome ' + ChrOI + ' and all ENSTs of the pattern "ENST.*' + tens_pre_dec + ones_pre_dec + '\..*"')
    
    elif re.match(r'NA', tens_pre_dec):
        print(proc_ID + ': done with HRA_mapping_analyzer() for tchromosome ' + ChrOI + ' and all ENSTs of the pattern "ENST.*\d' + ones_pre_dec + '\..*"')
    
    print(get_mem('HRA_mapping_analyzer'))
    
    return HRA_mapping_analyzer_dict

#in my thesis referred to as Function 4
# this function parses the reference taht comes as fasta file and immideately 
# runs a polyQ_finder() analysis on it
#@profile
def HRA_ref_analyzer(working_directory, bed_name, ChrOI, ones_pre_dec, ref_name, GFF3_path):

    csv_dir = working_directory + 'csv_results/'
    
    proc_ID = 'P_' + ChrOI + '_' + ones_pre_dec
    
    bed_pickle_path = working_directory + 'bed_pickles/' + bed_name + '_P_' + ChrOI + '_' + ones_pre_dec + '_ENST_split.txt'
    if os.path.exists(bed_pickle_path) == True:
        HRA_bed_parser_dict_original = pickle.load(open(bed_pickle_path, 'rb'))
        HRA_bed_parser_dict = copy.deepcopy(HRA_bed_parser_dict_original)
        print(proc_ID + ': ' + bed_pickle_path + ' will be used as HRA_bed_parser_dict')
    else:
        print(proc_ID + ': ' + bed_pickle_path + ' is missing ...')     
    
    csv_repeat_analysis_path = csv_dir + ref_name + '_' + proc_ID + '_polyQs_reference.csv'
    csv_repeat_analysis = open(csv_repeat_analysis_path, 'w')
    csv_repeat_analysis.close()
        
    ref_open_for_reading = open(working_directory + ref_name, 'r')
    
    parsed_exon_dict = GFF3_parser(GFF3_path)

    ### launching fasta_parser() ###  
    if 'Homo_sapiens.GRCh37.dna' in ref_name:
        parsed_ref = fasta_parser(directory = working_directory, bed = bed_name, fasta_file_object_opened_for_reading = ref_open_for_reading, ChrOI = str(ChrOI), ones_pre_dec = ones_pre_dec, GRCh37_cds = False, GRCh37_dna = True, specific_enst = False, HRA_bed_parser_dict = HRA_bed_parser_dict)
    elif 'Homo_sapiens.GRCh37.cds' in ref_name:
        parsed_ref = fasta_parser(directory = working_directory, bed = bed_name, fasta_file_object_opened_for_reading = ref_open_for_reading, ChrOI = str(ChrOI), ones_pre_dec = ones_pre_dec, GRCh37_cds = True, GRCh37_dna = False, specific_enst = False, HRA_bed_parser_dict = HRA_bed_parser_dict)
    else:
        print(proc_ID + ' : the reference is not recognized as either a full dna emble fasta file, or coding sequence exclusive emble fasta file. It will be assumed that the reference file is a common full dna fasta file.')
        parsed_ref = fasta_parser(directory = working_directory, bed = bed_name, fasta_file_object_opened_for_reading = ref_open_for_reading, ChrOI = str(ChrOI), ones_pre_dec = ones_pre_dec, GRCh37_cds = False, GRCh37_dna = True, specific_enst = False, HRA_bed_parser_dict = HRA_bed_parser_dict)
    
    if len(parsed_ref[0]) == 7:
        print(proc_ID + ' : coding Sequences fasta file was recognized as reference')
        for chromosome, ensg, enst, seq, strand, enst_start, enst_stop in parsed_ref:
            ### launching polyQ_finder() ###
            new_csv_rows = polyQ_finder(seq_in = seq, chromosome_in = ChrOI, ENSG_in = ensg, ENST_in = enst, strand_in = strand, starting_coordinate = enst_start, exon_dict = parsed_exon_dict)
            csv_repeat_analysis = open(csv_repeat_analysis_path, 'a')
            csv_repeat_analysis.write(new_csv_rows)
            csv_repeat_analysis.close()
        return 'cds'
                
    elif len(parsed_ref[0]) == 2:

        print(proc_ID + ' : full DNA fasta file was recognized as reference')
        
        #the addition of the '0' (could have been any other character) at the beginning of the nucleotide string is convinient because of python being 0-based, but ensembl annotation and sam files being 1-based
        whole_chromosome_seq = '0' + parsed_ref[0][1]
        
        if whole_chromosome_seq.startswith('>'):
            raise Exception(proc_ID + ' : HRA_ref_analyzer: header was not removed from submitted sequence')
        print(proc_ID + ' : these ENSTs are to be parsed from the reference: ' + str(HRA_bed_parser_dict[ChrOI + '|all_ENSTs']))
        for ENST in HRA_bed_parser_dict[ChrOI + '|all_ENSTs']:
            Chromosome, Start, Stop, strand = HRA_bed_parser_dict[ENST + '|bed_data']
            ENST_start = int(Start)
            ENST_stop = int(Stop) + 1
            print(proc_ID + ' : working on ' + ENST)
            print(proc_ID + ' : searching for polyQ on the following ENST: (' + strand + ', ' + str(ENST_start) + ', ' + str(ENST_stop) + ')')# + whole_chromosome_seq[ENST_start:ENST_stop])

            new_csv_rows = polyQ_finder(seq_in = whole_chromosome_seq[ENST_start:ENST_stop], chromosome_in = ChrOI, ENSG_in = HRA_bed_parser_dict[ENST + '|' + 'find_ENSG'][0][0], ENST_in = ENST, strand_in = strand, starting_coordinate = ENST_start, exon_dict = parsed_exon_dict)


            if len(new_csv_rows) > 0:
                print(proc_ID + ' : new csv row found')
            csv_repeat_analysis = open(csv_repeat_analysis_path, 'a')
            csv_repeat_analysis.write(new_csv_rows)
            csv_repeat_analysis.close()
        return 'dna'
    elif os.path.getsize(working_directory + ref_name) == 0:
        raise Exception(proc_ID + ' : reference cannot be parsed since ' + working_directory + ref_name + ' is an empty file.')
    else:
        raise Exception(proc_ID + ' : HRA_ref_analyzer()Error')

#in my thesis referred to as Function 5.1
#@profile
def HRA_merge_outfiles(working_directory, bam_name, ref_name, ChrOI_list, overwrite):
    
    ### merge query csv ###
    for ChrOI in ChrOI_list:
        csv_files = []
        csv_dir = working_directory + 'csv_results/'
        print('working on the merged query csv for chromosome: ' + ChrOI)
        new_query_csv = csv_dir + bam_name + '_Chromosome_' + ChrOI + '_merged_query.csv'
        if os.path.exists(csv_dir) == False:
            raise Exception('HRA_merge_outfiles(Error) : "csv_results/" does not exist, but the target files would be expected to be in this directory ...')
        elif os.path.getsize(csv_dir) == 0:
            raise Exception('HRA_merge_outfiles(Error) : "csv_results/" exists, but appears to be empty since its size is 0')
            if str(listdir(csv_dir)) == '[]':
                print('the process will continue since the folder returns an empty list of files')
            else:
                print('since the folder does not return an empty list of files it will be checked if the contained file appears to be a merged csv_file of chromosome ' + ChrOI + '...')
        if os.path.exists(new_query_csv) == False:
            print('It appears that no file called "' + new_query_csv + ' exists, therefore csv files of chromosome ' + ChrOI + ' are written to a single file: ' + new_query_csv)
            initiating = open(new_query_csv, 'w')
            initiating.close()
            csv_files = [(csv_dir + any_file) for any_file in listdir(csv_dir) if (isfile(join(working_directory  + 'csv_results/', any_file)) and re.match(r'csv_analysis_for_chromosome_P_' + ChrOI + r'(_\d){,2}_polyQs_query\.csv', any_file))]
            print('list of detected QRY csv files to be merged: ' + str(csv_files))
            if csv_files == []:
                print('no macthing query csv files were found at the given location ... ')
            for csv_file in csv_files:
                read_file = open(csv_file, 'r')
                target = open(new_query_csv, 'a')
                for line in read_file:
                    target.write(line)
                read_file.close()
                target.close()
        elif os.path.exists(new_query_csv) == True and overwrite != 'True':
            print('It appears that the QRY csv files of chromosome ' + ChrOI + ' have already been written to a single file called "' + new_query_csv + '" ... if you want to merge csv files anyway remove or rename "' + new_query_csv + '"')
        elif os.path.exists(new_query_csv) == True and overwrite == 'True':
            print('It appears that the QRY csv files of chromosome ' + ChrOI + ' have already been written to a single file called "' + new_query_csv + '" ... overwrite was set to "True" and file ' + new_query_csv + ' will be overwritten"')        
            initiating = open(new_query_csv, 'w')
            initiating.close()
            csv_files = [(csv_dir + any_file) for any_file in listdir(csv_dir) if (isfile(join(working_directory  + 'csv_results/', any_file)) and re.match(r'csv_analysis_for_chromosome_P_' + ChrOI + r'(_\d){,2}_polyQs_query\.csv', any_file))]
            print('list of detected QRY csv files to be merged: ' + str(csv_files))
            if csv_files == []:
                print('no macthing query csv files were found at the given location ... ')
            for csv_file in csv_files:
                read_file = open(csv_file, 'r')
                target = open(new_query_csv, 'a')
                for line in read_file:
                    target.write(line)
                read_file.close()
                target.close()
        else:
            raise Exception('HRA_merge_outfiles(ERROR): failed to verify status of merged QRY csv-files for chromosome ' + ChrOI)
        
        
        ### merge reference csv ###
        csv_files = []
        print('working on the merged reference csv for chromosome: ' + ChrOI)
        new_reference_csv = csv_dir + bam_name + '_Chromosome_' + ChrOI + '_merged_reference.csv'
        if os.path.exists(new_reference_csv) == False:
            print('It appears that no file called "' + new_reference_csv + '" exists, therefore csv files of chromosome ' + ChrOI + ' are written to a single file: ' + new_reference_csv)
            initiating = open(new_reference_csv, 'w')
            initiating.close()
            csv_files = [(csv_dir + any_file) for any_file in listdir(csv_dir) if (isfile(join(working_directory  + 'csv_results/', any_file)) and re.match(ref_name + ChrOI + '.fa_P_' + ChrOI + r'(_\d){,2}_polyQs_reference.csv', any_file))]
            print('list of detected REF csv files to be merged: ' + str(csv_files))
            if csv_files == []:
                print('no macthing reference csv files were found at the given location ... ')
            for csv_file in csv_files:
                read_file = open(csv_file, 'r')
                target = open(new_reference_csv, 'a')
                for line in read_file:
                    target.write(line)
                read_file.close()
                target.close()
        elif os.path.exists(new_reference_csv) == True and overwrite == 'True':
            print('It appears that the QRY csv files of chromosome ' + ChrOI + ' have already been written to a single file called "' + new_query_csv + '" ... overwrite was set to "True" and file ' + new_query_csv + ' will be overwritten"')        
            initiating = open(new_reference_csv, 'w')
            initiating.close()
            csv_files = [(csv_dir + any_file) for any_file in listdir(csv_dir) if (isfile(join(working_directory  + 'csv_results/', any_file)) and re.match(ref_name + ChrOI + '.fa_P_' + ChrOI + r'(_\d){,2}_polyQs_reference.csv', any_file))]
            print('list of detected REF csv files to be merged: ' + str(csv_files))
            if csv_files == []:
                print('no macthing reference csv files were found at the given location ... ')
            for csv_file in csv_files:
                read_file = open(csv_file, 'r')
                target = open(new_reference_csv, 'a')
                for line in read_file:
                    target.write(line)
                read_file.close()
                target.close()
        else:
            raise Exception('It appears that the REFERENCE csv files of chromosome merged REF csv-files for chromosome ' + ChrOI)
        
    return True

#@profile
def merge_alignment_pickle_after_HRA_mapping_analyzer(working_directory, bed_name, bam_name):
                                       
    new_pickle_name = working_directory + 'analyze_pickles/' + bed + bam + '_HRA_mapping_analyzer_merged_dict_pickle.txt'
    if os.path.exists(new_pickle_name) == False:
        print('HRA_mapping_analyzer() dict pickles are merged to a single dict pickle: "' + new_pickle_name + '"')
        list_of_dict_pickle_paths = [(working_directory + 'analyze_pickles/' + any_file) for any_file in listdir(working_directory + 'analyze_pickles/') if (isfile(join(working_directory + 'analyze_pickles/', any_file)) and re.match(r'.*align_pickle.txt', any_file))]
        list_of_dicts_loaded_from_pickles = [pickle.load(open(pickled_dict_path)) for pickled_dict_path in list_of_dict_pickle_paths]        
        if len(list_of_dicts_loaded_from_pickles) < 1:
            if len(list_of_dict_pickle_paths) >= 1:
                print('No dicts could be found under the following paths:')
                for path in list_of_dict_pickle_paths:
                    print(path)
                raise Exception('merge_alignment_pickle_after_HRA_mapping_analyzer(ERROR): pickles could not be loaded from given paths')
            else:
                raise Exception('merge_alignment_pickle_after_HRA_mapping_analyzer(ERROR): pickles could not be found under the given criteria')

        else:
            new_alignment_dict = merge_dicts(list_of_dicts_loaded_from_pickles)
            new_pickle = open(new_pickle_name, 'w')
            pickle.dump(new_alignment_dict, new_pickle, protocol=pickle.HIGHEST_PROTOCOL)
            new_pickle.close()
            print('HRA_mapping_analyzer() dict pickles merged and dumped as ' + new_pickle_name)
    
    elif os.path.exists(new_pickle_name) == True:
        print('It appears that seperate HRA_mapping_analyzer() dict pickles have already been written to a single file called "' + new_pickle_name + '" ... if you want to merge pickles anyway remove or rename "' + new_pickle_name + '"')
    return True

#in my thesis referred to as Function 5.2    
#@profile
def HRA_compare_QRY_and_REF(working_directory, bam_name, bed_name, ref_name, ChrOI, GFF3_path, aligner, positions):
        
    t1 = time.time()
    
    if positions != 'relative' and positions != 'absolute':
        raise Exception('HRA_compare_QRY_and_REF(ERROR): use the "position" flag and set it either to "relative" or "absolute" to specify how to handle positions of repeats in each ENST when comparing position changes between the reference and query')
    
    csv_dir = working_directory + 'csv_results/'
    
    comparison_report = csv_dir + bam_name + '_vs_' + bed_name + '_complete_report_on_Chromosome_' + ChrOI + '.csv'
    comparison_report_opened = open(comparison_report, 'w')
    comparison_report_opened.close()

    ENST_dict = {}
    
    QRY = csv_dir + bam_name + '_Chromosome_' + ChrOI + '_merged_query.csv'    
    QRY_exonic_lot_and_sot = csv_to_exonic_lot_and_sot_parser(csv_path = QRY, ChrOI = ChrOI)

    qry_exonic_lot_complete = QRY_exonic_lot_and_sot[0]
    if positions == 'relative':
        print('Positions of repeats will be used to compare their order of occurence between experiments')
        QRY_exonic_lot = [line[0:-2] for line in qry_exonic_lot_complete]
    elif positions == 'absolute':
        QRY_exonic_lot = [line for line in qry_exonic_lot_complete]
        print('Positions of repeats will be used to compare repeats directly between experiments')
    QRY_exonic_ENST_Seq_sot = QRY_exonic_lot_and_sot[1]
    
    REF = csv_dir + bam_name + '_Chromosome_' + ChrOI + '_merged_reference.csv'
    REF_exonic_lot_and_sot = csv_to_exonic_lot_and_sot_parser(csv_path = REF, ChrOI = ChrOI)
    
    ref_exonic_lot_complete = REF_exonic_lot_and_sot[0]
    if positions == 'relative':
        REF_exonic_lot = [line[0:-2] for line in ref_exonic_lot_complete]
    elif positions == 'absolute':
        REF_exonic_lot = [line for line in ref_exonic_lot_complete]
    REF_exonic_ENST_Seq_sot = REF_exonic_lot_and_sot[1]
    
    #as a reminder, here is what a line in ref/qry _exonic_lot_complete looks like:
    #   '4'(chromosome) = line[0], 
    #   'ENSG00000145242.9' (ENSG) = line[1], 
    #   'ENST00000511294.1' (ENST) = line[2], 
    #   '-' (strand) = line[3], 
    #   '5.0'(number of repeats) = line[4], 
    #   '1.0' ('rank') == (placement in all repeats identified in this ENST) = str(position_count),
    #   'TTGCTGCTGCTGCTG' (repeat seq) = line[7],
    #   '66521831' (repeat start) = line[5], 
    #   '66521846' (repeat end) = line[6],
    
    #so REF/QRY _exonic_lot will look like such:
    #   '4'(chromosome) = line[0], 
    #   'ENSG00000145242.9' (ENSG) = line[1], 
    #   'ENST00000511294.1' (ENST) = line[2], 
    #   '-' (strand) = line[3], 
    #   '5.0'(number of repeats) = line[4], 
    #   '1.0' ('rank') == (placement in all repeats identified in this ENST) = str(position_count),
    #   'TTGCTGCTGCTGCTG' (repeat seq) = line[7],
    
    print('QRY_exonic_ENST_Seq_sot: ' + str(QRY_exonic_ENST_Seq_sot))
    for ENST, Seq in QRY_exonic_ENST_Seq_sot:
        ENST_dict['QRY_' + ENST] = Seq
    print('REF_exonic_ENST_Seq_sot: ' + str(REF_exonic_ENST_Seq_sot))
    for ENST, Seq in REF_exonic_ENST_Seq_sot:
        ENST_dict['REF_' + ENST] = Seq
    print('REF_exonic_lot (first 100 characters): ' + str(REF_exonic_lot)[:100])
    print('QRY_exonic_lot (first 100 characters): ' + str(QRY_exonic_lot)[:100])
    
    QRY_exonic_lot_set = set(QRY_exonic_lot)
    REF_exonic_lot_set = set(REF_exonic_lot)
    QRY_unique = QRY_exonic_lot_set - REF_exonic_lot_set
    REF_unique = REF_exonic_lot_set - QRY_exonic_lot_set
    sum_of_uniques = (QRY_exonic_lot_set ^ REF_exonic_lot_set)
    
    unequal_tuples = [(unique_element[1], unique_element[2], unique_element[3]) for unique_element in sum_of_uniques]
    
    unique_tuples = set(unequal_tuples)
    results_log = []
    results_list = []
    
    exon_dict = GFF3_parser(GFF3_path)
    
    print('unique_tuples (complete): ' + str(unique_tuples))
    for unique_tuple in unique_tuples:

        unique_ENST = unique_tuple[1]
        
        if unique_ENST not in results_log:
            print('processing ENST: ' + unique_ENST)
            for exonic_qry_line in QRY_exonic_lot:
                exonic_qry_ENST = exonic_qry_line[2]
                if exonic_qry_ENST == unique_ENST:
                    
                    exonic_qry_chr = exonic_qry_line[0]
                    exonic_qry_ENSG = exonic_qry_line[1]
                    #exonic_qry_ENST = exonic_qry_line[2]
                    exonic_qry_strand = exonic_qry_line[3]
                    exonic_qry_rep_number = exonic_qry_line[4]
                    exonic_qry_rep_position = exonic_qry_line[5]
                    exonic_qry_rep_seq = exonic_qry_line[6]
                    exonic_qry_rep_start = exonic_qry_line[7]
                    exonic_qry_rep_stop = exonic_qry_line[8]
                    exonic_qry_full_seq = ENST_dict['QRY_' + unique_ENST]
                    exonic_AA_Seq = translate_seq_to_AA(exonic_qry_ENST, exonic_qry_strand, exonic_qry_full_seq, exon_dict)
                    if exonic_qry_line in QRY_unique:
                        new_line = ['QRY', exonic_qry_chr, exonic_qry_ENSG, exonic_qry_ENST, exonic_qry_strand, exonic_qry_rep_number, exonic_qry_rep_position, exonic_qry_rep_seq, exonic_qry_rep_start, exonic_qry_rep_stop, exonic_qry_full_seq, '\n']
                    else:
                        new_line = ['qry', exonic_qry_chr, exonic_qry_ENSG, exonic_qry_ENST, exonic_qry_strand, exonic_qry_rep_number, exonic_qry_rep_position, exonic_qry_rep_seq, exonic_qry_rep_start, exonic_qry_rep_stop, exonic_qry_full_seq, '\n']
                    results_list.append(new_line)
                    
                    print('AAS of ' + unique_ENST + ': ' + exonic_AA_Seq[1])
                                        
            for exonic_ref_line in REF_exonic_lot:
                exonic_ref_ENST = exonic_ref_line[2]
                if exonic_ref_ENST == unique_ENST:
                    
                    exonic_ref_chr = exonic_ref_line[0]
                    exonic_ref_ENSG = exonic_ref_line[1]
                    #exonic_ref_ENST = exonic_ref_line[2]
                    exonic_ref_strand = exonic_ref_line[3]
                    exonic_ref_rep_number = exonic_ref_line[4]
                    exonic_ref_rep_position = exonic_ref_line[5]
                    exonic_ref_rep_seq = exonic_ref_line[6]
                    exonic_ref_rep_start = exonic_ref_line[7]
                    exonic_ref_rep_stop = exonic_ref_line[8]
                    exonic_ref_full_seq = ENST_dict['REF_' + unique_ENST]
                    if exonic_ref_line in REF_unique:
                        new_line = ['REF', exonic_ref_chr, exonic_ref_ENSG, exonic_ref_ENST, exonic_ref_strand, exonic_ref_rep_number, exonic_ref_rep_position, exonic_ref_rep_seq, exonic_ref_rep_start, exonic_ref_rep_stop, exonic_ref_full_seq, '\n']
                    else:
                        new_line = ['ref', exonic_ref_chr, exonic_ref_ENSG, exonic_ref_ENST, exonic_ref_strand, exonic_ref_rep_number, exonic_ref_rep_position, exonic_ref_rep_seq, exonic_ref_rep_start, exonic_ref_rep_stop, exonic_ref_full_seq, '\n']
                        
                    results_list.append(new_line)
                    
            results_log.append(unique_ENST)
            
    results_list_sorted_by_rep_position = sorted(results_list, key=return7th)
    results_list_sorted_by_REF_QRY = sorted(results_list_sorted_by_rep_position, key=return1st_case_insensitive)
    results_list_sorted_by_ENST = sorted(results_list_sorted_by_REF_QRY, key=return4th)
    results_list_sorted_by_ENSG = sorted(results_list_sorted_by_ENST, key=return3rd)
    
    result_strings = [[','.join(line)][0] for line in results_list_sorted_by_ENSG]
    final_result = ''.join(result_strings)
    
    comparison_report_opened = open(comparison_report, 'w')
    comparison_report_opened.write(final_result)
    comparison_report_opened.close()
                    
    t2 = time.time()
    return ('alignment took : ' + str((t2 - t1)) + 's')

#in my thesis referred to as Function 5.3
#make the plot from Figure 4
def HRA_trans_intron_polyQ_analyzer(working_directory, gff3, experiment_folder, bam_name):
    #if overwrite_first_call == 'False':
    #    plot_intron_spanning_polyQs_both_sides(working_directory = working_directory)
    #    plot_intron_spanning_polyQs_both_sides_after_side(working_directory = working_directory)
    
    check_if_polyQ_spans_introns(working_directory, gff3, experiment_folder, bam_name)
    plot_intron_spanning_polyQs_both_sides_I(working_directory)
    plot_intron_spanning_polyQs_both_sides_II(working_directory)


def check_if_polyQ_spans_introns(working_directory, gff3, experiment_folder, bam_name):#part of Function 5.3

    random_experiment = experiment_folder
    
    aa_dict = defaultdict(list)

    REF = working_directory + '/' + random_experiment + '/Homo_sapiens.GRCh37.dna_sm.chromosome.'
    GFF3_path = working_directory + gff3
    ChrOI_list = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']
    data_list = []
    hit_list = []

    results = ''
                
    #continue with the first experiment as representative experiment for reference analysis
    complete_results_path = working_directory + '/' + random_experiment + '/csv_results/'
    for ChrOI in ChrOI_list:
        print('working on ChrOI: ' + ChrOI)
        print(REF + ChrOI + r'.fa_P_' + ChrOI + r'_\d{,2}_polyQs_reference.csv')
        ref_csv_files = [(complete_results_path + any_file) for any_file in listdir(complete_results_path) if (isfile(join(complete_results_path, any_file)) and re.match(bam_name + '_Chromosome_' + ChrOI + '_merged_reference.csv', any_file))]
                                                                                                                        
        for ref_csv_file in ref_csv_files:
            ref_complete_results = open(ref_csv_file, 'r')
            for line in ref_complete_results:
                split_line = line.rstrip('\n').split(',')
                print('split_line: ' + str(split_line[:-1]))
                data_list.append(split_line[:7])
    
    GFF3_chromosome_dict = dict()
    for ChrOI in ChrOI_list:
        GFF3_chromosome_dict[ChrOI] = GFF3_parser(GFF3_path + ChrOI + '.gff3')

    for dlist in data_list:
        
        chromosome_in = dlist[0]
        ENSG_in = dlist[1]
        ENST_in = dlist[2]
        strand_in = dlist[3]
        repeat_seq_length = dlist[4]
        repeat_start = dlist[5]
        repeat_stop = dlist[6]
        repeat_start = int(repeat_start)
        repeat_stop = int(repeat_stop)
        current_gff3 = GFF3_chromosome_dict[chromosome_in]
        
        for CDS in current_gff3[ENST_in]:
            cds1_range = 99
            CDS_index = current_gff3[ENST_in].index(CDS)
            CDS_range_start = CDS[0]
            CDS_range_stop = CDS[1]
            if (CDS_range_stop - CDS_range_start) < 99:
                cds1_range = (CDS_range_stop - CDS_range_start) - ((CDS_range_stop - CDS_range_start) % 3)
            sub_feature_phase = int(CDS[2])
            proximity_of_repeat_start = [(repeat_start),(repeat_start - 1),(repeat_start - 2),(repeat_start - 3),(repeat_start - 4),(repeat_start - 5),(repeat_start - 6),(repeat_start - 7),(repeat_start - 8),(repeat_start - 9),(repeat_start - 10)]            
            proximity_of_repeat_stop = [(repeat_stop + 10),(repeat_stop + 9),(repeat_stop + 8),(repeat_stop + 7),(repeat_stop + 6),(repeat_stop + 5),(repeat_stop + 4),(repeat_stop + 3),(repeat_stop + 2),(repeat_stop + 1),(repeat_stop)]            
            
            if CDS_range_start in proximity_of_repeat_start:
                ref_open_for_reading = open(REF + str(chromosome_in) + '.fa', 'r')
                parsed_ref = simple_fasta_parser(fasta_file_object_opened_for_reading = ref_open_for_reading)
                whole_chromosome_seq = '0' + parsed_ref[0][1]
                
                if strand_in == '+':
                    
    ###########################################################################

                    if CDS_index == 0:
                        print('UTR')
                        CDS_range_stop += 1
                        results += 'UTR ... ' + whole_chromosome_seq[(CDS_range_start - 25):CDS_range_start] + '|' + whole_chromosome_seq[CDS_range_start:(CDS_range_start + 50)] + '\n'
                        results += str(chromosome_in) + '|' + str(ENSG_in) + '|' + str(ENST_in) + '|' + str(strand_in) + '|' + str(repeat_seq_length) + '|' + str(repeat_start) + '|' + str(repeat_stop) + '\n\n'
     
                    else:
                        aa_dict['position'].append('<')
                        aa_dict['strand'].append(strand_in)
                        
                        hit_list.append(ENSG_in)
                        prae_CDS = current_gff3[ENST_in][CDS_index - 1]
                        prae_CDS_range_start = prae_CDS[0]
                        prae_CDS_range_stop = prae_CDS[1]
                        prae_sub_feature_phase = int(prae_CDS[2])
                        
                        aa_dict['ENSG'].append(ENSG_in)
                        aa_dict['CDS1_start'].append(CDS_range_start)
                        aa_dict['CDS1_stop'].append(CDS_range_stop)
                        aa_dict['CDS2_start'].append(prae_CDS_range_start)
                        aa_dict['CDS2_stop'].append(prae_CDS_range_stop)  
                        
                        results += '<CDS: ' + str(prae_CDS) + '\n.CDS: ' + str(CDS) + '\n'
                        results += str(prae_sub_feature_phase) + '->' + str(((prae_CDS_range_stop - prae_CDS_range_start) + prae_sub_feature_phase) % 3) + '|' + str(sub_feature_phase) + '->' + str(((CDS_range_stop - CDS_range_start) + sub_feature_phase) % 3) + '\n'
                        cds2_range = 99

                        if (prae_CDS_range_stop - prae_CDS_range_start) < 99:
                            cds2_range = (prae_CDS_range_stop - prae_CDS_range_start) - ((prae_CDS_range_stop - prae_CDS_range_start) % 3)
                                                
                        AA = ''
                        string_triplet = []
                        prae_stop_phase = ((((prae_CDS_range_stop - prae_CDS_range_start) % 3) + prae_sub_feature_phase) % 3)
                        if prae_stop_phase == 0:
                            results += 'intra\n'
                            factor = 1
                            CDS_range_stop += 1
                            prae_CDS_range_stop += 1
                            string_triplet = [whole_chromosome_seq[(prae_CDS_range_stop - (cds2_range + factor)):(prae_CDS_range_stop - factor)], 
                                                                   whole_chromosome_seq[(prae_CDS_range_stop - factor):prae_CDS_range_stop] + whole_chromosome_seq[CDS_range_start:CDS_range_start + (3 - factor)], 
                                                                   whole_chromosome_seq[(CDS_range_start + (3 - factor)):(CDS_range_start + cds1_range  + (3 - factor))]]
                        elif prae_stop_phase == 1:
                            results += 'intra\n'
                            factor = 2
                            CDS_range_stop += 1
                            prae_CDS_range_stop += 1
                            string_triplet = [whole_chromosome_seq[(prae_CDS_range_stop - (cds2_range + factor)):(prae_CDS_range_stop - factor)], 
                                                                   whole_chromosome_seq[(prae_CDS_range_stop - factor):prae_CDS_range_stop] + whole_chromosome_seq[CDS_range_start:CDS_range_start + (3 - factor)], 
                                                                   whole_chromosome_seq[(CDS_range_start + (3 - factor)):(CDS_range_start + cds1_range  + (3 - factor))]]
                        elif prae_stop_phase == 2:
                            results += 'inter\n'
                            print('no split border!')
                            CDS_range_stop += 1
                            prae_CDS_range_stop += 1
                            string_triplet = [whole_chromosome_seq[(prae_CDS_range_stop - cds2_range):(prae_CDS_range_stop)], 
                                                                   '', 
                                                                   whole_chromosome_seq[CDS_range_start:(CDS_range_start + cds1_range)]]
                        AA = simple_translate_seq_to_AA(strand = strand_in, string_triplet = string_triplet)
                        aa_dict['string_triplet'].append(string_triplet)
                        AA = (' ' * int((99 - cds2_range) / 3)) + AA + (' ' * int((99 - cds1_range) / 3))
                        AA = AA.replace('||', ' ').replace('|', '')
                        
                        letter_count = 0
                        for letter in AA:
                            #the 'C' stands for 'column'
                            aa_dict[letter_count].append(letter)
                            letter_count += 1
                                                    
                        results += 'whole_chromosome_seq[' + str(prae_CDS_range_stop - 50) + ':' + str(prae_CDS_range_stop) + "] + '[" + str(prae_stop_phase) + "]|' + whole_chromosome_seq[" + str(prae_CDS_range_stop) + ':' + str(prae_CDS_range_stop + 25) + "] + ' ... ' + whole_chromosome_seq[" + str(CDS_range_start - 25) + ':' + str(CDS_range_start) + "] + '|' + whole_chromosome_seq[" + str(CDS_range_start) + ':' + str(CDS_range_start + 50) + "]\n"
                        results += whole_chromosome_seq[(prae_CDS_range_stop - 50):prae_CDS_range_stop] + '[' + str(prae_stop_phase) + ']|' + whole_chromosome_seq[prae_CDS_range_stop:(prae_CDS_range_stop + 25)] + ' ... ' + whole_chromosome_seq[(CDS_range_start - 25):CDS_range_start] + '|' + whole_chromosome_seq[CDS_range_start:(CDS_range_start + 50)] + '\n'
                        results += str(chromosome_in) + '|' + str(ENSG_in) + '|' + str(ENST_in) + '|' + str(strand_in) + '|' + str(repeat_seq_length) + '|' + str(repeat_start) + '|' + str(repeat_stop) + '\n'
                        results += str(string_triplet) + '\n'
                        results += AA + '\n\n'

                elif strand_in == '-':
                    
    ###########################################################################

                    if CDS_index == 0:
                        print('UTR')
                        CDS_range_stop += 1
                        results += 'UTR ... ' + whole_chromosome_seq[(CDS_range_start - 25):CDS_range_start] + '|' + whole_chromosome_seq[CDS_range_start:(CDS_range_start + 50)] + '\n'
                        results += str(chromosome_in) + '|' + str(ENSG_in) + '|' + str(ENST_in) + '|' + str(strand_in) + '|' + str(repeat_seq_length) + '|' + str(repeat_start) + '|' + str(repeat_stop) + '\n\n'
     
                    else:
                        aa_dict['position'].append('>')
                        aa_dict['strand'].append(strand_in)
                        
                        hit_list.append(ENSG_in)
                        prae_CDS = current_gff3[ENST_in][CDS_index - 1]
                        prae_CDS_range_start = prae_CDS[0]
                        prae_CDS_range_stop = prae_CDS[1]
                        prae_sub_feature_phase = int(prae_CDS[2])
                        
                        aa_dict['ENSG'].append(ENSG_in)
                        aa_dict['CDS1_start'].append(CDS_range_start)
                        aa_dict['CDS1_stop'].append(CDS_range_stop)
                        aa_dict['CDS2_start'].append(prae_CDS_range_start)
                        aa_dict['CDS2_stop'].append(prae_CDS_range_stop)  
                        
                        results += '.CDS: ' + str(CDS) + '\nCDS>: ' + str(prae_CDS) + '\n'
                        results += str(((prae_CDS_range_stop - prae_CDS_range_start) + prae_sub_feature_phase) % 3) + '<-' + str(prae_sub_feature_phase) + '|' + str(((CDS_range_stop - CDS_range_start) + sub_feature_phase) % 3) + '<-' + str(sub_feature_phase) + '\n'

                        
                        
                        cds2_range = 99

                        if (prae_CDS_range_stop - prae_CDS_range_start) < 99:
                            cds2_range = (prae_CDS_range_stop - prae_CDS_range_start) - ((prae_CDS_range_stop - prae_CDS_range_start) % 3)

                        AA = ''
                        string_triplet = []
    
                        prae_stop_phase = prae_sub_feature_phase
                        if prae_stop_phase == 2:
                            results += 'intra\n'
                            factor = 1
                            CDS_range_stop += 1
                            prae_CDS_range_stop += 1
                            string_triplet = [whole_chromosome_seq[(prae_CDS_range_stop - (cds2_range + factor)):(prae_CDS_range_stop - factor)], 
                                                                   whole_chromosome_seq[(prae_CDS_range_stop - factor):prae_CDS_range_stop] + whole_chromosome_seq[CDS_range_start:CDS_range_start + (3 - factor)], 
                                                                   whole_chromosome_seq[(CDS_range_start + (3 - factor)):(CDS_range_start + cds1_range  + (3 - factor))]]
                        elif prae_stop_phase == 1:
                            results += 'intra\n'
                            factor = 2
                            CDS_range_stop += 1
                            prae_CDS_range_stop += 1
                            string_triplet = [whole_chromosome_seq[(prae_CDS_range_stop - (cds2_range + factor)):(prae_CDS_range_stop - factor)], 
                                                                   whole_chromosome_seq[(prae_CDS_range_stop - factor):prae_CDS_range_stop] + whole_chromosome_seq[CDS_range_start:CDS_range_start + (3 - factor)], 
                                                                   whole_chromosome_seq[(CDS_range_start + (3 - factor)):(CDS_range_start + cds1_range  + (3 - factor))]]
                        elif prae_stop_phase == 0:
                            results += 'inter\n'
                            print('no split border!')
                            CDS_range_stop += 1
                            prae_CDS_range_stop += 1
                            string_triplet = [whole_chromosome_seq[(prae_CDS_range_stop - cds2_range):(prae_CDS_range_stop)], 
                                                                   '', 
                                                                   whole_chromosome_seq[CDS_range_start:(CDS_range_start + cds1_range)]]
                        AA = simple_translate_seq_to_AA(strand = strand_in, string_triplet = string_triplet)
                        aa_dict['string_triplet'].append(string_triplet)                       
                        AA = (' ' * int((99 - cds1_range) / 3)) + AA + (' ' * int((99 - cds2_range) / 3))
                        AA = AA.replace('||', ' ').replace('|', '')
                        
                        letter_count = 0
                        for letter in AA:
                            #the 'C' stands for 'column'
                            aa_dict[letter_count].append(letter)
                            letter_count += 1
                        
                        results += 'whole_chromosome_seq[' + str(prae_CDS_range_stop - 50) + ':' + str(prae_CDS_range_stop) + "] + '[" + str(prae_stop_phase) + "]|' + whole_chromosome_seq[" + str(prae_CDS_range_stop) + ':' + str(prae_CDS_range_stop + 25) + "] + ' ... ' + whole_chromosome_seq[" + str(CDS_range_start - 25) + ':' + str(CDS_range_start) + "] + '|' + whole_chromosome_seq[" + str(CDS_range_start) + ':' + str(CDS_range_start + 50) + "]\n"
                        results += whole_chromosome_seq[(prae_CDS_range_stop - 50):prae_CDS_range_stop] + '[' + str(prae_stop_phase) + ']|' + whole_chromosome_seq[prae_CDS_range_stop:(prae_CDS_range_stop + 25)] + ' ... ' + whole_chromosome_seq[(CDS_range_start - 25):CDS_range_start] + '|' + whole_chromosome_seq[CDS_range_start:(CDS_range_start + 50)] + '\n'
                        results += str(chromosome_in) + '|' + str(ENSG_in) + '|' + str(ENST_in) + '|' + str(strand_in) + '|' + str(repeat_seq_length) + '|' + str(repeat_start) + '|' + str(repeat_stop) + '\n'
                        results += str(string_triplet) + '\n'                        
                        results += AA + '\n\n'
                        
                        
            elif CDS_range_stop in proximity_of_repeat_stop:
                ref_open_for_reading = open(REF + str(chromosome_in) + '.fa', 'r')
                parsed_ref = simple_fasta_parser(fasta_file_object_opened_for_reading = ref_open_for_reading)
                whole_chromosome_seq = '0' + parsed_ref[0][1]
                
                if strand_in == '+': 
                
                    if (CDS_index + 1) == len(current_gff3[ENST_in]):
                        print('UTR')
                        CDS_range_stop += 1
                        results += whole_chromosome_seq[(CDS_range_stop - 50):CDS_range_stop] + '|' + whole_chromosome_seq[CDS_range_stop:(CDS_range_stop + 25)] + ' ... UTR' + '\n'
                        results += str(chromosome_in) + '|' + str(ENSG_in) + '|' + str(ENST_in) + '|' + str(strand_in) + '|' + str(repeat_seq_length) + '|' + str(repeat_start) + '|' + str(repeat_stop) + '\n\n'
                    else:
                        aa_dict['position'].append('>')
                        aa_dict['strand'].append(strand_in)
                        
                        hit_list.append(ENSG_in)    
                        post_CDS = current_gff3[ENST_in][CDS_index + 1]
                        post_CDS_range_start = post_CDS[0]
                        post_CDS_range_stop = post_CDS[1]
                        post_sub_feature_phase = int(post_CDS[2])
                        
                        aa_dict['ENSG'].append(ENSG_in)
                        aa_dict['CDS1_start'].append(CDS_range_start)
                        aa_dict['CDS1_stop'].append(CDS_range_stop)
                        aa_dict['CDS2_start'].append(post_CDS_range_start)
                        aa_dict['CDS2_stop'].append(post_CDS_range_stop)                        
                        results += '.CDS: ' + str(CDS) + '\nCDS>: ' + str(post_CDS) + '\n'
                        
                        results += str(sub_feature_phase) + '->' + str(((CDS_range_stop - CDS_range_start) + sub_feature_phase) % 3) + '|' + str(post_sub_feature_phase) + '->' + str(((post_CDS_range_stop - post_CDS_range_start) + post_sub_feature_phase) % 3) + '\n'
                        
                        
                        
                        
                        cds2_range = 99

                        if (post_CDS_range_stop - post_CDS_range_start) < 99:
                            cds2_range = (post_CDS_range_stop - post_CDS_range_start) - ((post_CDS_range_stop - post_CDS_range_start) % 3)
                        AA = ''
                        string_triplet = []
                        post_stop_phase = post_sub_feature_phase
                        if post_stop_phase == 1:
                            results += 'intra\n'
                            factor = 1
                            CDS_range_stop += 1
                            post_CDS_range_stop += 1    
                            string_triplet = [whole_chromosome_seq[(CDS_range_stop - (cds1_range + factor)):(CDS_range_stop - factor)], 
                                              whole_chromosome_seq[(CDS_range_stop - factor):CDS_range_stop] + whole_chromosome_seq[post_CDS_range_start:post_CDS_range_start + (3 - factor)], 
                                              whole_chromosome_seq[(post_CDS_range_start + (3 - factor)):(post_CDS_range_start + cds2_range  + (3 - factor))]
                                              ]
                        elif post_stop_phase == 2:
                            results += 'intra\n'
                            factor = 2
                            CDS_range_stop += 1
                            post_CDS_range_stop += 1    
                            string_triplet = [whole_chromosome_seq[(CDS_range_stop - (cds1_range + factor)):(CDS_range_stop - factor)], 
                                              whole_chromosome_seq[(CDS_range_stop - factor):CDS_range_stop] + whole_chromosome_seq[post_CDS_range_start:post_CDS_range_start + (3 - factor)], 
                                              whole_chromosome_seq[(post_CDS_range_start + (3 - factor)):(post_CDS_range_start + cds2_range  + (3 - factor))]
                                              ]
                        elif post_stop_phase == 0:
                            results += 'inter\n'
                            print('no split border!')
                            CDS_range_stop += 1
                            post_CDS_range_stop += 1    
                            string_triplet = [whole_chromosome_seq[(CDS_range_stop - cds1_range):(CDS_range_stop)], 
                                              '', 
                                              whole_chromosome_seq[post_CDS_range_start:(post_CDS_range_start + cds2_range)]
                                              ]
                        AA = simple_translate_seq_to_AA(strand = strand_in, string_triplet = string_triplet)
                        aa_dict['string_triplet'].append(string_triplet)
                        AA = (' ' * int((99 - cds1_range) / 3)) + AA + (' ' * int((99 - cds2_range) / 3))
                        AA = AA.replace('||', ' ').replace('|', '')

                        letter_count = 0
                        for letter in AA:
                            #the 'C' stands for 'column'
                            aa_dict[letter_count].append(letter)
                            letter_count += 1
                            
                        results += "whole_chromosome_seq[" + str(CDS_range_stop - 50) + ":" + str(CDS_range_stop) + "] + '|' + whole_chromosome_seq[" + str(CDS_range_stop) + ":" + str(CDS_range_stop + 25) + "] + ' ... ' + whole_chromosome_seq[" + str(post_CDS_range_start - 25) + ":" + str(post_CDS_range_start) + "] + '|[" + str(post_stop_phase) + "]' + whole_chromosome_seq[" + str(post_CDS_range_start) + ":" + str(post_CDS_range_start + 50) + "]\n"
                        results += whole_chromosome_seq[(CDS_range_stop - 50):CDS_range_stop] + '|' + whole_chromosome_seq[CDS_range_stop:(CDS_range_stop + 25)] + ' ... ' + whole_chromosome_seq[(post_CDS_range_start - 25):post_CDS_range_start] + '|[' + str(post_stop_phase) + ']' + whole_chromosome_seq[post_CDS_range_start:(post_CDS_range_start + 50)] + '\n'
                        results += str(chromosome_in) + '|' + str(ENSG_in) + '|' + str(ENST_in) + '|' + str(strand_in) + '|' + str(repeat_seq_length) + '|' + str(repeat_start) + '|' + str(repeat_stop) + '\n'
                        results += str(string_triplet) + '\n'
                        results += AA + '\n\n'                    
                    
                elif strand_in == '-':
                    if (CDS_index + 1) == len(current_gff3[ENST_in]):
                        print('UTR')
                        CDS_range_stop += 1
                        results += whole_chromosome_seq[(CDS_range_stop - 50):CDS_range_stop] + '|' + whole_chromosome_seq[CDS_range_stop:(CDS_range_stop + 25)] + ' ... UTR' + '\n'
                        results += str(chromosome_in) + '|' + str(ENSG_in) + '|' + str(ENST_in) + '|' + str(strand_in) + '|' + str(repeat_seq_length) + '|' + str(repeat_start) + '|' + str(repeat_stop) + '\n\n'
                    else:
                        aa_dict['position'].append('<')
                        aa_dict['strand'].append(strand_in)
                        
                        hit_list.append(ENSG_in)    
                        post_CDS = current_gff3[ENST_in][CDS_index + 1]
                        post_CDS_range_start = post_CDS[0]
                        post_CDS_range_stop = post_CDS[1]
                        post_sub_feature_phase = int(post_CDS[2])

                        aa_dict['ENSG'].append(ENSG_in)
                        aa_dict['CDS1_start'].append(CDS_range_start)
                        aa_dict['CDS1_stop'].append(CDS_range_stop)
                        aa_dict['CDS2_start'].append(post_CDS_range_start)
                        aa_dict['CDS2_stop'].append(post_CDS_range_stop)

                        results += '<CDS: ' + str(post_CDS) + '\n.CDS: ' + str(CDS) + '\n'
                        
                        results += str(((CDS_range_stop - CDS_range_start) + sub_feature_phase) % 3) + '<-' + str(sub_feature_phase) + '|' + str(((post_CDS_range_stop - post_CDS_range_start) + post_sub_feature_phase) % 3) + '<-' + str(post_sub_feature_phase) + '\n'
                        
                        cds2_range = 99

                        if (post_CDS_range_stop - post_CDS_range_start) < 99:
                            cds2_range = (post_CDS_range_stop - post_CDS_range_start) - ((post_CDS_range_stop - post_CDS_range_start) % 3)
                        AA = ''
                        string_triplet = []
                            
                        post_stop_phase = ((((post_CDS_range_stop - post_CDS_range_start) % 3) + post_sub_feature_phase) % 3)
                        if post_stop_phase == 1:
                            results += 'intra\n'
                            factor = 1
                            CDS_range_stop += 1
                            post_CDS_range_stop += 1    
                            string_triplet = [whole_chromosome_seq[(CDS_range_stop - (cds1_range + factor)):(CDS_range_stop - factor)], 
                                              whole_chromosome_seq[(CDS_range_stop - factor):CDS_range_stop] + whole_chromosome_seq[post_CDS_range_start:post_CDS_range_start + (3 - factor)], 
                                              whole_chromosome_seq[(post_CDS_range_start + (3 - factor)):(post_CDS_range_start + cds2_range  + (3 - factor))]
                                              ]
                        elif post_stop_phase == 0:
                            results += 'intra\n'
                            factor = 2
                            CDS_range_stop += 1
                            post_CDS_range_stop += 1    
                            string_triplet = [whole_chromosome_seq[(CDS_range_stop - (cds1_range + factor)):(CDS_range_stop - factor)], 
                                              whole_chromosome_seq[(CDS_range_stop - factor):CDS_range_stop] + whole_chromosome_seq[post_CDS_range_start:post_CDS_range_start + (3 - factor)], 
                                              whole_chromosome_seq[(post_CDS_range_start + (3 - factor)):(post_CDS_range_start + cds2_range  + (3 - factor))]
                                              ]
                        elif post_stop_phase == 2:
                            results += 'inter\n'
                            print('no split border!')
                            CDS_range_stop += 1
                            post_CDS_range_stop += 1    
                            string_triplet = [whole_chromosome_seq[(CDS_range_stop - cds1_range):(CDS_range_stop)], 
                                              '', 
                                              whole_chromosome_seq[post_CDS_range_start:(post_CDS_range_start + cds2_range)]
                                              ]
                        AA = simple_translate_seq_to_AA(strand = strand_in, string_triplet = string_triplet)
                        aa_dict['string_triplet'].append(string_triplet)
                        AA = (' ' * int((99 - cds2_range) / 3)) + AA + (' ' * int((99 - cds1_range) / 3))
                        AA = AA.replace('||', ' ').replace('|', '')

                        letter_count = 0
                        for letter in AA:
                            #the 'C' stands for 'column'
                            aa_dict[letter_count].append(letter)
                            letter_count += 1
                            
                        results += "whole_chromosome_seq[" + str(CDS_range_stop - 50) + ":" + str(CDS_range_stop) + "] + '|' + whole_chromosome_seq[" + str(CDS_range_stop) + ":" + str(CDS_range_stop + 25) + "] + ' ... ' + whole_chromosome_seq[" + str(post_CDS_range_start - 25) + ":" + str(post_CDS_range_start) + "] + '|[" + str(post_stop_phase) + "]' + whole_chromosome_seq[" + str(post_CDS_range_start) + ":" + str(post_CDS_range_start + 50) + "]\n"
                        results += whole_chromosome_seq[(CDS_range_stop - 50):CDS_range_stop] + '|' + whole_chromosome_seq[CDS_range_stop:(CDS_range_stop + 25)] + ' ... ' + whole_chromosome_seq[(post_CDS_range_start - 25):post_CDS_range_start] + '|[' + str(post_stop_phase) + ']' + whole_chromosome_seq[post_CDS_range_start:(post_CDS_range_start + 50)] + '\n'
                        results += str(chromosome_in) + '|' + str(ENSG_in) + '|' + str(ENST_in) + '|' + str(strand_in) + '|' + str(repeat_seq_length) + '|' + str(repeat_start) + '|' + str(repeat_stop) + '\n'
                        results += str(string_triplet) + '\n'
                        results += AA + '\n\n'    
                        

    path = working_directory + 'trans_intron_polyQs.txt'
    path_opened = open(path, 'w')
    path_opened.write(results)
    path_opened.close()
    print('hits: ' + str(len(set(hit_list))))
    print('aa_dict: ' + str(aa_dict))
    full_df = pd.DataFrame.from_dict(aa_dict)
    prae_df = full_df[full_df.position == '<']
    post_df = full_df[full_df.position == '>']
    #all_letters_at_position = (set(df['ENSG'].tolist()))
    
    full_df.to_pickle(working_directory + 'trans_intron_polyQ_full_df_pickle.pkl')
    print('whole df: ')
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(full_df)
    
    prae_df.to_pickle(working_directory + 'trans_intron_polyQ_prae_df_pickle.pkl')
    print('prae_df: ')
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(prae_df)
    
    post_df.to_pickle(working_directory + 'trans_intron_polyQ_post_df_pickle.pkl')
    print('post_df: ')
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(post_df)
        
    full_df = pd.read_pickle(working_directory + 'trans_intron_polyQ_full_df_pickle.pkl')
    prae_df = pd.read_pickle(working_directory + 'trans_intron_polyQ_prae_df_pickle.pkl')
    post_df = pd.read_pickle(working_directory + 'trans_intron_polyQ_post_df_pickle.pkl')

    return True

def plot_intron_spanning_polyQs_both_sides_II(working_directory):#part of Function 5.3
    # Intron/Exon boundary plots : characteristics of the exon following a polyQ
    viridis_data = [[0.267004, 0.004874, 0.329415],[0.268510, 0.009605, 0.335427],[0.269944, 0.014625, 0.341379],[0.271305, 0.019942, 0.347269],[0.272594, 0.025563, 0.353093],[0.273809, 0.031497, 0.358853],[0.274952, 0.037752, 0.364543],[0.276022, 0.044167, 0.370164],[0.277018, 0.050344, 0.375715],[0.277941, 0.056324, 0.381191],[0.278791, 0.062145, 0.386592],[0.279566, 0.067836, 0.391917],[0.280267, 0.073417, 0.397163],[0.280894, 0.078907, 0.402329],[0.281446, 0.084320, 0.407414],[0.281924, 0.089666, 0.412415],[0.282327, 0.094955, 0.417331],[0.282656, 0.100196, 0.422160],[0.282910, 0.105393, 0.426902],[0.283091, 0.110553, 0.431554],[0.283197, 0.115680, 0.436115],[0.283229, 0.120777, 0.440584],[0.283187, 0.125848, 0.444960],[0.283072, 0.130895, 0.449241],[0.282884, 0.135920, 0.453427],[0.282623, 0.140926, 0.457517],[0.282290, 0.145912, 0.461510],[0.281887, 0.150881, 0.465405],[0.281412, 0.155834, 0.469201],[0.280868, 0.160771, 0.472899],[0.280255, 0.165693, 0.476498],[0.279574, 0.170599, 0.479997],[0.278826, 0.175490, 0.483397],[0.278012, 0.180367, 0.486697],[0.277134, 0.185228, 0.489898],[0.276194, 0.190074, 0.493001],[0.275191, 0.194905, 0.496005],[0.274128, 0.199721, 0.498911],[0.273006, 0.204520, 0.501721],[0.271828, 0.209303, 0.504434],[0.270595, 0.214069, 0.507052],[0.269308, 0.218818, 0.509577],[0.267968, 0.223549, 0.512008],[0.266580, 0.228262, 0.514349],[0.265145, 0.232956, 0.516599],[0.263663, 0.237631, 0.518762],[0.262138, 0.242286, 0.520837],[0.260571, 0.246922, 0.522828],[0.258965, 0.251537, 0.524736],[0.257322, 0.256130, 0.526563],[0.255645, 0.260703, 0.528312],[0.253935, 0.265254, 0.529983],[0.252194, 0.269783, 0.531579],[0.250425, 0.274290, 0.533103],[0.248629, 0.278775, 0.534556],[0.246811, 0.283237, 0.535941],[0.244972, 0.287675, 0.537260],[0.243113, 0.292092, 0.538516],[0.241237, 0.296485, 0.539709],[0.239346, 0.300855, 0.540844],[0.237441, 0.305202, 0.541921],[0.235526, 0.309527, 0.542944],[0.233603, 0.313828, 0.543914],[0.231674, 0.318106, 0.544834],[0.229739, 0.322361, 0.545706],[0.227802, 0.326594, 0.546532],[0.225863, 0.330805, 0.547314],[0.223925, 0.334994, 0.548053],[0.221989, 0.339161, 0.548752],[0.220057, 0.343307, 0.549413],[0.218130, 0.347432, 0.550038],[0.216210, 0.351535, 0.550627],[0.214298, 0.355619, 0.551184],[0.212395, 0.359683, 0.551710],[0.210503, 0.363727, 0.552206],[0.208623, 0.367752, 0.552675],[0.206756, 0.371758, 0.553117],[0.204903, 0.375746, 0.553533],[0.203063, 0.379716, 0.553925],[0.201239, 0.383670, 0.554294],[0.199430, 0.387607, 0.554642],[0.197636, 0.391528, 0.554969],[0.195860, 0.395433, 0.555276],[0.194100, 0.399323, 0.555565],[0.192357, 0.403199, 0.555836],[0.190631, 0.407061, 0.556089],[0.188923, 0.410910, 0.556326],[0.187231, 0.414746, 0.556547],[0.185556, 0.418570, 0.556753],[0.183898, 0.422383, 0.556944],[0.182256, 0.426184, 0.557120],[0.180629, 0.429975, 0.557282],[0.179019, 0.433756, 0.557430],[0.177423, 0.437527, 0.557565],[0.175841, 0.441290, 0.557685],[0.174274, 0.445044, 0.557792],[0.172719, 0.448791, 0.557885],[0.171176, 0.452530, 0.557965],[0.169646, 0.456262, 0.558030],[0.168126, 0.459988, 0.558082],[0.166617, 0.463708, 0.558119],[0.165117, 0.467423, 0.558141],[0.163625, 0.471133, 0.558148],[0.162142, 0.474838, 0.558140],[0.160665, 0.478540, 0.558115],[0.159194, 0.482237, 0.558073],[0.157729, 0.485932, 0.558013],[0.156270, 0.489624, 0.557936],[0.154815, 0.493313, 0.557840],[0.153364, 0.497000, 0.557724],[0.151918, 0.500685, 0.557587],[0.150476, 0.504369, 0.557430],[0.149039, 0.508051, 0.557250],[0.147607, 0.511733, 0.557049],[0.146180, 0.515413, 0.556823],[0.144759, 0.519093, 0.556572],[0.143343, 0.522773, 0.556295],[0.141935, 0.526453, 0.555991],[0.140536, 0.530132, 0.555659],[0.139147, 0.533812, 0.555298],[0.137770, 0.537492, 0.554906],[0.136408, 0.541173, 0.554483],[0.135066, 0.544853, 0.554029],[0.133743, 0.548535, 0.553541],[0.132444, 0.552216, 0.553018],[0.131172, 0.555899, 0.552459],[0.129933, 0.559582, 0.551864],[0.128729, 0.563265, 0.551229],[0.127568, 0.566949, 0.550556],[0.126453, 0.570633, 0.549841],[0.125394, 0.574318, 0.549086],[0.124395, 0.578002, 0.548287],[0.123463, 0.581687, 0.547445],[0.122606, 0.585371, 0.546557],[0.121831, 0.589055, 0.545623],[0.121148, 0.592739, 0.544641],[0.120565, 0.596422, 0.543611],[0.120092, 0.600104, 0.542530],[0.119738, 0.603785, 0.541400],[0.119512, 0.607464, 0.540218],[0.119423, 0.611141, 0.538982],[0.119483, 0.614817, 0.537692],[0.119699, 0.618490, 0.536347],[0.120081, 0.622161, 0.534946],[0.120638, 0.625828, 0.533488],[0.121380, 0.629492, 0.531973],[0.122312, 0.633153, 0.530398],[0.123444, 0.636809, 0.528763],[0.124780, 0.640461, 0.527068],[0.126326, 0.644107, 0.525311],[0.128087, 0.647749, 0.523491],[0.130067, 0.651384, 0.521608],[0.132268, 0.655014, 0.519661],[0.134692, 0.658636, 0.517649],[0.137339, 0.662252, 0.515571],[0.140210, 0.665859, 0.513427],[0.143303, 0.669459, 0.511215],[0.146616, 0.673050, 0.508936],[0.150148, 0.676631, 0.506589],[0.153894, 0.680203, 0.504172],[0.157851, 0.683765, 0.501686],[0.162016, 0.687316, 0.499129],[0.166383, 0.690856, 0.496502],[0.170948, 0.694384, 0.493803],[0.175707, 0.697900, 0.491033],[0.180653, 0.701402, 0.488189],[0.185783, 0.704891, 0.485273],[0.191090, 0.708366, 0.482284],[0.196571, 0.711827, 0.479221],[0.202219, 0.715272, 0.476084],[0.208030, 0.718701, 0.472873],[0.214000, 0.722114, 0.469588],[0.220124, 0.725509, 0.466226],[0.226397, 0.728888, 0.462789],[0.232815, 0.732247, 0.459277],[0.239374, 0.735588, 0.455688],[0.246070, 0.738910, 0.452024],[0.252899, 0.742211, 0.448284],[0.259857, 0.745492, 0.444467],[0.266941, 0.748751, 0.440573],[0.274149, 0.751988, 0.436601],[0.281477, 0.755203, 0.432552],[0.288921, 0.758394, 0.428426],[0.296479, 0.761561, 0.424223],[0.304148, 0.764704, 0.419943],[0.311925, 0.767822, 0.415586],[0.319809, 0.770914, 0.411152],[0.327796, 0.773980, 0.406640],[0.335885, 0.777018, 0.402049],[0.344074, 0.780029, 0.397381],[0.352360, 0.783011, 0.392636],[0.360741, 0.785964, 0.387814],[0.369214, 0.788888, 0.382914],[0.377779, 0.791781, 0.377939],[0.386433, 0.794644, 0.372886],[0.395174, 0.797475, 0.367757],[0.404001, 0.800275, 0.362552],[0.412913, 0.803041, 0.357269],[0.421908, 0.805774, 0.351910],[0.430983, 0.808473, 0.346476],[0.440137, 0.811138, 0.340967],[0.449368, 0.813768, 0.335384],[0.458674, 0.816363, 0.329727],[0.468053, 0.818921, 0.323998],[0.477504, 0.821444, 0.318195],[0.487026, 0.823929, 0.312321],[0.496615, 0.826376, 0.306377],[0.506271, 0.828786, 0.300362],[0.515992, 0.831158, 0.294279],[0.525776, 0.833491, 0.288127],[0.535621, 0.835785, 0.281908],[0.545524, 0.838039, 0.275626],[0.555484, 0.840254, 0.269281],[0.565498, 0.842430, 0.262877],[0.575563, 0.844566, 0.256415],[0.585678, 0.846661, 0.249897],[0.595839, 0.848717, 0.243329],[0.606045, 0.850733, 0.236712],[0.616293, 0.852709, 0.230052],[0.626579, 0.854645, 0.223353],[0.636902, 0.856542, 0.216620],[0.647257, 0.858400, 0.209861],[0.657642, 0.860219, 0.203082],[0.668054, 0.861999, 0.196293],[0.678489, 0.863742, 0.189503],[0.688944, 0.865448, 0.182725],[0.699415, 0.867117, 0.175971],[0.709898, 0.868751, 0.169257],[0.720391, 0.870350, 0.162603],[0.730889, 0.871916, 0.156029],[0.741388, 0.873449, 0.149561],[0.751884, 0.874951, 0.143228],[0.762373, 0.876424, 0.137064],[0.772852, 0.877868, 0.131109],[0.783315, 0.879285, 0.125405],[0.793760, 0.880678, 0.120005],[0.804182, 0.882046, 0.114965],[0.814576, 0.883393, 0.110347],[0.824940, 0.884720, 0.106217],[0.835270, 0.886029, 0.102646],[0.845561, 0.887322, 0.099702],[0.855810, 0.888601, 0.097452],[0.866013, 0.889868, 0.095953],[0.876168, 0.891125, 0.095250],[0.886271, 0.892374, 0.095374],[0.896320, 0.893616, 0.096335],[0.906311, 0.894855, 0.098125],[0.916242, 0.896091, 0.100717],[0.926106, 0.897330, 0.104071],[0.935904, 0.898570, 0.108131],[0.945636, 0.899815, 0.112838],[0.955300, 0.901065, 0.118128],[0.964894, 0.902323, 0.123941],[0.974417, 0.903590, 0.130215],[0.983868, 0.904867, 0.136897],[0.993248, 0.906157, 0.143936]]

    prae_df = pd.read_pickle(working_directory + 'trans_intron_polyQ_prae_df_pickle.pkl')
    post_df = pd.read_pickle(working_directory + 'trans_intron_polyQ_post_df_pickle.pkl')
    list_of_AAs_single_letters = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z','_',"' '",]
    list_of_dfs_and_names = [(prae_df, 'prae_df'),(post_df, 'post_df'),]
    
    
    
    for df, name_1 in list_of_dfs_and_names:
        print('\n\n\n' + name_1)
        #identify duplicate rows that share the complete subset
        df = df.drop_duplicates(subset=['ENSG', 'CDS1_start', 'CDS1_stop', 'CDS2_start', 'CDS2_stop',], keep='first')#keep options   ->   first : Drop duplicates except for the first occurrence; last : Drop duplicates except for the last occurrence; False : Drop all duplicates.
        #drop unnecessary columns
        df = df.drop('position', 1)# number == axis ... (0 == row; 1 == column)
        df = df.drop('strand', 1)
        df = df.drop('string_triplet', 1)
        df = df.drop('CDS1_start', 1)
        df = df.drop('CDS1_stop', 1)
        df = df.drop('CDS2_start', 1)
        df = df.drop('CDS2_stop', 1)
        df = df.drop('ENSG', 1)
        #show the resulting df
        print('\nsum:')
        print(str(df))
        
        #'prae': polyQ lies after the intron (-> Region of interest is the end of the preceding exon)
        #'post': polyQ lies before the intron (-> Region of interest is the start of the posterior exon)

        #divide into different polyQ distances from the exon boundary        
        if 'prae' in name_1:
            #create 10 sub dfs that cover all shifting frames of where a 4Q-repeat can be found up to 10 amino acids from the exon start
            sub_df_1 = df[(df[ 33 + 1]=='Q') & (df[ 33 + 2]=='Q') & (df[ 33 + 3]=='Q') & (df[ 33 + 4]=='Q')]
            sub_df_2 = df[(df[ 33 + 2]=='Q') & (df[ 33 + 3]=='Q') & (df[ 33 + 4]=='Q') & (df[ 33 + 5]=='Q')]
            sub_df_3 = df[(df[ 33 + 3]=='Q') & (df[ 33 + 4]=='Q') & (df[ 33 + 5]=='Q') & (df[ 33 + 6]=='Q')]
            sub_df_4 = df[(df[ 33 + 4]=='Q') & (df[ 33 + 5]=='Q') & (df[ 33 + 6]=='Q') & (df[ 33 + 7]=='Q')]
            sub_df_5 = df[(df[ 33 + 5]=='Q') & (df[ 33 + 6]=='Q') & (df[ 33 + 7]=='Q') & (df[ 33 + 8]=='Q')]
            sub_df_6 = df[(df[ 33 + 6]=='Q') & (df[ 33 + 7]=='Q') & (df[ 33 + 8]=='Q') & (df[ 33 + 9]=='Q')]
            sub_df_7 = df[(df[ 33 + 7]=='Q') & (df[ 33 + 8]=='Q') & (df[ 33 + 9]=='Q') & (df[ 33 + 10]=='Q')]
            sub_df_8 = df[(df[ 33 + 8]=='Q') & (df[ 33 + 9]=='Q') & (df[ 33 + 10]=='Q') & (df[ 33 + 11]=='Q')]
            sub_df_9 = df[(df[ 33 + 9]=='Q') & (df[ 33 + 10]=='Q') & (df[ 33 + 11]=='Q') & (df[ 33 + 12]=='Q')]
            sub_df_10 = df[(df[ 33 + 10]=='Q') & (df[ 33 + 11]=='Q') & (df[ 33 + 12]=='Q') & (df[ 33 + 13]=='Q')]
        elif 'post' in name_1:
            #create 10 sub dfs that cover all shifting frames of where a 4Q-repeat can be found up to 10 amino acids up to the exon stop
            sub_df_1 = df[(df[ 33 - 1]=='Q') & (df[ 33 - 2]=='Q') & (df[ 33 - 3]=='Q') & (df[ 33 - 4]=='Q')]
            sub_df_2 = df[(df[ 33 - 2]=='Q') & (df[ 33 - 3]=='Q') & (df[ 33 - 4]=='Q') & (df[ 33 - 5]=='Q')]
            sub_df_3 = df[(df[ 33 - 3]=='Q') & (df[ 33 - 4]=='Q') & (df[ 33 - 5]=='Q') & (df[ 33 - 6]=='Q')]
            sub_df_4 = df[(df[ 33 - 4]=='Q') & (df[ 33 - 5]=='Q') & (df[ 33 - 6]=='Q') & (df[ 33 - 7]=='Q')]
            sub_df_5 = df[(df[ 33 - 5]=='Q') & (df[ 33 - 6]=='Q') & (df[ 33 - 7]=='Q') & (df[ 33 - 8]=='Q')]
            sub_df_6 = df[(df[ 33 - 6]=='Q') & (df[ 33 - 7]=='Q') & (df[ 33 - 8]=='Q') & (df[ 33 - 9]=='Q')]
            sub_df_7 = df[(df[ 33 - 7]=='Q') & (df[ 33 - 8]=='Q') & (df[ 33 - 9]=='Q') & (df[ 33 - 10]=='Q')]
            sub_df_8 = df[(df[ 33 - 8]=='Q') & (df[ 33 - 9]=='Q') & (df[ 33 - 10]=='Q') & (df[ 33 - 11]=='Q')]
            sub_df_9 = df[(df[ 33 - 9]=='Q') & (df[ 33 - 10]=='Q') & (df[ 33 - 11]=='Q') & (df[ 33 - 12]=='Q')]
            sub_df_10 = df[(df[ 33 - 10]=='Q') & (df[ 33 - 11]=='Q') & (df[ 33 - 12]=='Q') & (df[ 33 - 13]=='Q')]        
        
        sub_df_1_total_dict = defaultdict(list)
        sub_df_1_relative_dict = defaultdict(list)

        sub_df_2_total_dict = defaultdict(list)
        sub_df_2_relative_dict = defaultdict(list)

        sub_df_3_total_dict = defaultdict(list)
        sub_df_3_relative_dict = defaultdict(list)
        
        sub_df_4_total_dict = defaultdict(list)
        sub_df_4_relative_dict = defaultdict(list)
        
        sub_df_5_total_dict = defaultdict(list)
        sub_df_5_relative_dict = defaultdict(list)
        
        sub_df_6_total_dict = defaultdict(list)
        sub_df_6_relative_dict = defaultdict(list)
        
        sub_df_7_total_dict = defaultdict(list)
        sub_df_7_relative_dict = defaultdict(list)
        
        sub_df_8_total_dict = defaultdict(list)
        sub_df_8_relative_dict = defaultdict(list)
        
        sub_df_9_total_dict = defaultdict(list)
        sub_df_9_relative_dict = defaultdict(list)
        
        sub_df_10_total_dict = defaultdict(list)
        sub_df_10_relative_dict = defaultdict(list)
        
        sub_df_triplet_lot = [(sub_df_1, sub_df_1_total_dict, sub_df_1_relative_dict),
                              (sub_df_2, sub_df_2_total_dict, sub_df_2_relative_dict),
                              (sub_df_3, sub_df_3_total_dict, sub_df_3_relative_dict),
                              (sub_df_4, sub_df_4_total_dict, sub_df_4_relative_dict),
                              (sub_df_5, sub_df_5_total_dict, sub_df_5_relative_dict),
                              (sub_df_6, sub_df_6_total_dict, sub_df_6_relative_dict),
                              (sub_df_7, sub_df_7_total_dict, sub_df_7_relative_dict),
                              (sub_df_8, sub_df_8_total_dict, sub_df_8_relative_dict),
                              (sub_df_9, sub_df_9_total_dict, sub_df_9_relative_dict),
                              (sub_df_10, sub_df_10_total_dict, sub_df_10_relative_dict),]
        count = 0 
        for sub_df, sub_df_total_dict, sub_df_relative_dict in sub_df_triplet_lot:
            count += 1
            #create a total and relative df per frame sub-df
            sub_df_count_row = sub_df.shape[0]  # gives number of row count
            header_list = sub_df.columns.tolist()
            negative_start_value = ((len(header_list) - 1) / (- 2))
            for header in header_list:
                #print('\n\n\nheader: ' + str(header))
                sub_df_total_dict['position'].append(int(negative_start_value + header))
                sub_df_relative_dict['position'].append(int(negative_start_value + header))
                column_list = sub_df[header].tolist()
                for character in list_of_AAs_single_letters:
                    char_counts = int(column_list.count(character))
                    sub_df_total_dict[character].append(char_counts)
                    sub_df_relative_dict[character].append(char_counts / sub_df_count_row)

            
        
        sub_df_1_total_df = pd.DataFrame.from_dict(sub_df_1_total_dict)
        print('\nsub_df_1_total_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_1_total_df)
        sub_df_1_relative_df = pd.DataFrame.from_dict(sub_df_1_relative_dict)
        print('\nsub_df_1_relative_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_1_relative_df)
        
        sub_df_2_total_df = pd.DataFrame.from_dict(sub_df_2_total_dict)
        print('\nsub_df_2_total_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_2_total_df)
        sub_df_2_relative_df = pd.DataFrame.from_dict(sub_df_2_relative_dict)
        print('\nsub_df_2_relative_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_2_relative_df)
        
        sub_df_3_total_df = pd.DataFrame.from_dict(sub_df_3_total_dict)
        print('\nsub_df_3_total_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_3_total_df)
        sub_df_3_relative_df = pd.DataFrame.from_dict(sub_df_3_relative_dict)
        print('\nsub_df_3_relative_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_3_relative_df)
        
        sub_df_4_total_df = pd.DataFrame.from_dict(sub_df_4_total_dict)
        print('\nsub_df_4_total_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_4_total_df)
        sub_df_4_relative_df = pd.DataFrame.from_dict(sub_df_4_relative_dict)
        print('\nsub_df_4_relative_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_4_relative_df)
        
        sub_df_5_total_df = pd.DataFrame.from_dict(sub_df_5_total_dict)
        print('\nsub_df_5_total_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_5_total_df)
        sub_df_5_relative_df = pd.DataFrame.from_dict(sub_df_5_relative_dict)
        print('\nsub_df_5_relative_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_5_relative_df)
        
        sub_df_6_total_df = pd.DataFrame.from_dict(sub_df_6_total_dict)
        print('\nsub_df_6_total_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_6_total_df)
        sub_df_6_relative_df = pd.DataFrame.from_dict(sub_df_6_relative_dict)
        print('\nsub_df_6_relative_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_6_relative_df)
        
        sub_df_7_total_df = pd.DataFrame.from_dict(sub_df_7_total_dict)
        print('\nsub_df_7_total_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_7_total_df)
        sub_df_7_relative_df = pd.DataFrame.from_dict(sub_df_7_relative_dict)
        print('\nsub_df_7_relative_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_7_relative_df)
        
        sub_df_8_total_df = pd.DataFrame.from_dict(sub_df_8_total_dict)
        print('\nsub_df_8_total_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_8_total_df)
        sub_df_8_relative_df = pd.DataFrame.from_dict(sub_df_8_relative_dict)
        print('\nsub_df_8_relative_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_8_relative_df)
        
        sub_df_9_total_df = pd.DataFrame.from_dict(sub_df_9_total_dict)
        print('\nsub_df_9_total_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_9_total_df)
        sub_df_9_relative_df = pd.DataFrame.from_dict(sub_df_9_relative_dict)
        print('\nsub_df_9_relative_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_9_relative_df)
        
        sub_df_10_total_df = pd.DataFrame.from_dict(sub_df_10_total_dict)
        print('\nsub_df_10_total_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_10_total_df)
        sub_df_10_relative_df = pd.DataFrame.from_dict(sub_df_10_relative_dict)
        print('\nsub_df_10_relative_df')
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(sub_df_10_relative_df)
        
        total_and_relative_evaluated_sub_dfs = [(sub_df_1_total_df, sub_df_2_total_df, sub_df_3_total_df, sub_df_4_total_df, sub_df_5_total_df, sub_df_6_total_df, sub_df_7_total_df, sub_df_8_total_df, sub_df_9_total_df, sub_df_10_total_df, '_total'),
                                                (sub_df_1_relative_df, sub_df_2_relative_df, sub_df_3_relative_df, sub_df_4_relative_df, sub_df_5_relative_df, sub_df_6_relative_df, sub_df_7_relative_df, sub_df_8_relative_df, sub_df_9_relative_df, sub_df_10_relative_df, '_relative')]

        num_acids = 33
        number_of_evaluated_Q_positions = 10

        for evaluated_sub_df_1, evaluated_sub_df_2, evaluated_sub_df_3, evaluated_sub_df_4, evaluated_sub_df_5, evaluated_sub_df_6, evaluated_sub_df_7, evaluated_sub_df_8, evaluated_sub_df_9, evaluated_sub_df_10, name_2 in total_and_relative_evaluated_sub_dfs:
            
            print('start bar plot')

            
            print(name_1 + name_2)
            #Intron/Exon boundary bar plots

            if num_acids == 33:
                
                if 'prae' in name_1:
                    #create 10 sub dfs that cover all shifting frames of where a 4Q-repeat can be found up to 10 amino acids up to the exon stop
                    # set height of bar
                    
                    bars_sub_list_1 = [(evaluated_sub_df_1.at[(33 - 0), 'Q']), (evaluated_sub_df_1.at[(33 - 1), 'Q']), (evaluated_sub_df_1.at[(33 - 2), 'Q']), (evaluated_sub_df_1.at[(33 - 3), 'Q']), (evaluated_sub_df_1.at[(33 - 4), 'Q']), (evaluated_sub_df_1.at[(33 - 5), 'Q']), (evaluated_sub_df_1.at[(33 - 6), 'Q']), (evaluated_sub_df_1.at[(33 - 7), 'Q']), (evaluated_sub_df_1.at[(33 - 8), 'Q']), (evaluated_sub_df_1.at[(33 - 9), 'Q']), 
                                       (evaluated_sub_df_1.at[(33 - 10), 'Q']), (evaluated_sub_df_1.at[(33 - 11), 'Q']), (evaluated_sub_df_1.at[(33 - 12), 'Q']), (evaluated_sub_df_1.at[(33 - 13), 'Q']), (evaluated_sub_df_1.at[(33 - 14), 'Q']), (evaluated_sub_df_1.at[(33 - 15), 'Q']), (evaluated_sub_df_1.at[(33 - 16), 'Q']), (evaluated_sub_df_1.at[(33 - 17), 'Q']), (evaluated_sub_df_1.at[(33 - 18), 'Q']), (evaluated_sub_df_1.at[(33 - 19), 'Q']), 
                                       (evaluated_sub_df_1.at[(33 - 20), 'Q']), (evaluated_sub_df_1.at[(33 - 21), 'Q']), (evaluated_sub_df_1.at[(33 - 22), 'Q']), (evaluated_sub_df_1.at[(33 - 23), 'Q']), (evaluated_sub_df_1.at[(33 - 24), 'Q']), (evaluated_sub_df_1.at[(33 - 25), 'Q']), (evaluated_sub_df_1.at[(33 - 26), 'Q']), (evaluated_sub_df_1.at[(33 - 27), 'Q']), (evaluated_sub_df_1.at[(33 - 28), 'Q']), (evaluated_sub_df_1.at[(33 - 29), 'Q']), 
                                       (evaluated_sub_df_1.at[(33 - 30), 'Q']), (evaluated_sub_df_1.at[(33 - 31), 'Q']), (evaluated_sub_df_1.at[(33 - 32), 'Q']), (evaluated_sub_df_1.at[(33 - 33), 'Q']), ][::-1]
                    #bars_sub_list_1 = [(evaluated_sub_df_1.at[(33 - 0), 'Q']), (evaluated_sub_df_1.at[(33 - 1), 'Q']), (evaluated_sub_df_1.at[(33 - 2), 'Q']), (evaluated_sub_df_1.at[(33 - 3), 'Q']), (evaluated_sub_df_1.at[(33 - 4), 'Q']), (evaluated_sub_df_1.at[(33 - 5), 'Q']), (evaluated_sub_df_1.at[(33 - 6), 'Q']), (evaluated_sub_df_1.at[(33 - 7), 'Q']), (evaluated_sub_df_1.at[(33 - 8), 'Q']), (evaluated_sub_df_1.at[(33 - 9), 'Q']), (evaluated_sub_df_1.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_2 = [(evaluated_sub_df_2.at[(33 - 0), 'Q']), (evaluated_sub_df_2.at[(33 - 1), 'Q']), (evaluated_sub_df_2.at[(33 - 2), 'Q']), (evaluated_sub_df_2.at[(33 - 3), 'Q']), (evaluated_sub_df_2.at[(33 - 4), 'Q']), (evaluated_sub_df_2.at[(33 - 5), 'Q']), (evaluated_sub_df_2.at[(33 - 6), 'Q']), (evaluated_sub_df_2.at[(33 - 7), 'Q']), (evaluated_sub_df_2.at[(33 - 8), 'Q']), (evaluated_sub_df_2.at[(33 - 9), 'Q']), 
                                       (evaluated_sub_df_2.at[(33 - 10), 'Q']), (evaluated_sub_df_2.at[(33 - 11), 'Q']), (evaluated_sub_df_2.at[(33 - 12), 'Q']), (evaluated_sub_df_2.at[(33 - 13), 'Q']), (evaluated_sub_df_2.at[(33 - 14), 'Q']), (evaluated_sub_df_2.at[(33 - 15), 'Q']), (evaluated_sub_df_2.at[(33 - 16), 'Q']), (evaluated_sub_df_2.at[(33 - 17), 'Q']), (evaluated_sub_df_2.at[(33 - 18), 'Q']), (evaluated_sub_df_2.at[(33 - 19), 'Q']), 
                                       (evaluated_sub_df_2.at[(33 - 20), 'Q']), (evaluated_sub_df_2.at[(33 - 21), 'Q']), (evaluated_sub_df_2.at[(33 - 22), 'Q']), (evaluated_sub_df_2.at[(33 - 23), 'Q']), (evaluated_sub_df_2.at[(33 - 24), 'Q']), (evaluated_sub_df_2.at[(33 - 25), 'Q']), (evaluated_sub_df_2.at[(33 - 26), 'Q']), (evaluated_sub_df_2.at[(33 - 27), 'Q']), (evaluated_sub_df_2.at[(33 - 28), 'Q']), (evaluated_sub_df_2.at[(33 - 29), 'Q']), 
                                       (evaluated_sub_df_2.at[(33 - 30), 'Q']), (evaluated_sub_df_2.at[(33 - 31), 'Q']), (evaluated_sub_df_2.at[(33 - 32), 'Q']), (evaluated_sub_df_2.at[(33 - 33), 'Q']), ][::-1]
                    #bars_sub_list_2 = [(evaluated_sub_df_2.at[(33 - 0), 'Q']), (evaluated_sub_df_2.at[(33 - 1), 'Q']), (evaluated_sub_df_2.at[(33 - 2), 'Q']), (evaluated_sub_df_2.at[(33 - 3), 'Q']), (evaluated_sub_df_2.at[(33 - 4), 'Q']), (evaluated_sub_df_2.at[(33 - 5), 'Q']), (evaluated_sub_df_2.at[(33 - 6), 'Q']), (evaluated_sub_df_2.at[(33 - 7), 'Q']), (evaluated_sub_df_2.at[(33 - 8), 'Q']), (evaluated_sub_df_2.at[(33 - 9), 'Q']), (evaluated_sub_df_2.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_3 = [(evaluated_sub_df_3.at[(33 - 0), 'Q']), (evaluated_sub_df_3.at[(33 - 1), 'Q']), (evaluated_sub_df_3.at[(33 - 2), 'Q']), (evaluated_sub_df_3.at[(33 - 3), 'Q']), (evaluated_sub_df_3.at[(33 - 4), 'Q']), (evaluated_sub_df_3.at[(33 - 5), 'Q']), (evaluated_sub_df_3.at[(33 - 6), 'Q']), (evaluated_sub_df_3.at[(33 - 7), 'Q']), (evaluated_sub_df_3.at[(33 - 8), 'Q']), (evaluated_sub_df_3.at[(33 - 9), 'Q']), 
                                       (evaluated_sub_df_3.at[(33 - 10), 'Q']), (evaluated_sub_df_3.at[(33 - 11), 'Q']), (evaluated_sub_df_3.at[(33 - 12), 'Q']), (evaluated_sub_df_3.at[(33 - 13), 'Q']), (evaluated_sub_df_3.at[(33 - 14), 'Q']), (evaluated_sub_df_3.at[(33 - 15), 'Q']), (evaluated_sub_df_3.at[(33 - 16), 'Q']), (evaluated_sub_df_3.at[(33 - 17), 'Q']), (evaluated_sub_df_3.at[(33 - 18), 'Q']), (evaluated_sub_df_3.at[(33 - 19), 'Q']), 
                                       (evaluated_sub_df_3.at[(33 - 20), 'Q']), (evaluated_sub_df_3.at[(33 - 21), 'Q']), (evaluated_sub_df_3.at[(33 - 22), 'Q']), (evaluated_sub_df_3.at[(33 - 23), 'Q']), (evaluated_sub_df_3.at[(33 - 24), 'Q']), (evaluated_sub_df_3.at[(33 - 25), 'Q']), (evaluated_sub_df_3.at[(33 - 26), 'Q']), (evaluated_sub_df_3.at[(33 - 27), 'Q']), (evaluated_sub_df_3.at[(33 - 28), 'Q']), (evaluated_sub_df_3.at[(33 - 29), 'Q']), 
                                       (evaluated_sub_df_3.at[(33 - 30), 'Q']), (evaluated_sub_df_3.at[(33 - 31), 'Q']), (evaluated_sub_df_3.at[(33 - 32), 'Q']), (evaluated_sub_df_3.at[(33 - 33), 'Q']), ][::-1]
                    #bars_sub_list_3 = [(evaluated_sub_df_3.at[(33 - 0), 'Q']), (evaluated_sub_df_3.at[(33 - 1), 'Q']), (evaluated_sub_df_3.at[(33 - 2), 'Q']), (evaluated_sub_df_3.at[(33 - 3), 'Q']), (evaluated_sub_df_3.at[(33 - 4), 'Q']), (evaluated_sub_df_3.at[(33 - 5), 'Q']), (evaluated_sub_df_3.at[(33 - 6), 'Q']), (evaluated_sub_df_3.at[(33 - 7), 'Q']), (evaluated_sub_df_3.at[(33 - 8), 'Q']), (evaluated_sub_df_3.at[(33 - 9), 'Q']), (evaluated_sub_df_3.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_4 = [(evaluated_sub_df_4.at[(33 - 0), 'Q']), (evaluated_sub_df_4.at[(33 - 1), 'Q']), (evaluated_sub_df_4.at[(33 - 2), 'Q']), (evaluated_sub_df_4.at[(33 - 3), 'Q']), (evaluated_sub_df_4.at[(33 - 4), 'Q']), (evaluated_sub_df_4.at[(33 - 5), 'Q']), (evaluated_sub_df_4.at[(33 - 6), 'Q']), (evaluated_sub_df_4.at[(33 - 7), 'Q']), (evaluated_sub_df_4.at[(33 - 8), 'Q']), (evaluated_sub_df_4.at[(33 - 9), 'Q']), 
                                       (evaluated_sub_df_4.at[(33 - 10), 'Q']), (evaluated_sub_df_4.at[(33 - 11), 'Q']), (evaluated_sub_df_4.at[(33 - 12), 'Q']), (evaluated_sub_df_4.at[(33 - 13), 'Q']), (evaluated_sub_df_4.at[(33 - 14), 'Q']), (evaluated_sub_df_4.at[(33 - 15), 'Q']), (evaluated_sub_df_4.at[(33 - 16), 'Q']), (evaluated_sub_df_4.at[(33 - 17), 'Q']), (evaluated_sub_df_4.at[(33 - 18), 'Q']), (evaluated_sub_df_4.at[(33 - 19), 'Q']), 
                                       (evaluated_sub_df_4.at[(33 - 20), 'Q']), (evaluated_sub_df_4.at[(33 - 21), 'Q']), (evaluated_sub_df_4.at[(33 - 22), 'Q']), (evaluated_sub_df_4.at[(33 - 23), 'Q']), (evaluated_sub_df_4.at[(33 - 24), 'Q']), (evaluated_sub_df_4.at[(33 - 25), 'Q']), (evaluated_sub_df_4.at[(33 - 26), 'Q']), (evaluated_sub_df_4.at[(33 - 27), 'Q']), (evaluated_sub_df_4.at[(33 - 28), 'Q']), (evaluated_sub_df_4.at[(33 - 29), 'Q']), 
                                       (evaluated_sub_df_4.at[(33 - 30), 'Q']), (evaluated_sub_df_4.at[(33 - 31), 'Q']), (evaluated_sub_df_4.at[(33 - 32), 'Q']), (evaluated_sub_df_4.at[(33 - 33), 'Q']), ][::-1]
                    #bars_sub_list_4 = [(evaluated_sub_df_4.at[(33 - 0), 'Q']), (evaluated_sub_df_4.at[(33 - 1), 'Q']), (evaluated_sub_df_4.at[(33 - 2), 'Q']), (evaluated_sub_df_4.at[(33 - 3), 'Q']), (evaluated_sub_df_4.at[(33 - 4), 'Q']), (evaluated_sub_df_4.at[(33 - 5), 'Q']), (evaluated_sub_df_4.at[(33 - 6), 'Q']), (evaluated_sub_df_4.at[(33 - 7), 'Q']), (evaluated_sub_df_4.at[(33 - 8), 'Q']), (evaluated_sub_df_4.at[(33 - 9), 'Q']), (evaluated_sub_df_4.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_5 = [(evaluated_sub_df_5.at[(33 - 0), 'Q']), (evaluated_sub_df_5.at[(33 - 1), 'Q']), (evaluated_sub_df_5.at[(33 - 2), 'Q']), (evaluated_sub_df_5.at[(33 - 3), 'Q']), (evaluated_sub_df_5.at[(33 - 4), 'Q']), (evaluated_sub_df_5.at[(33 - 5), 'Q']), (evaluated_sub_df_5.at[(33 - 6), 'Q']), (evaluated_sub_df_5.at[(33 - 7), 'Q']), (evaluated_sub_df_5.at[(33 - 8), 'Q']), (evaluated_sub_df_5.at[(33 - 9), 'Q']), 
                                       (evaluated_sub_df_5.at[(33 - 10), 'Q']), (evaluated_sub_df_5.at[(33 - 11), 'Q']), (evaluated_sub_df_5.at[(33 - 12), 'Q']), (evaluated_sub_df_5.at[(33 - 13), 'Q']), (evaluated_sub_df_5.at[(33 - 14), 'Q']), (evaluated_sub_df_5.at[(33 - 15), 'Q']), (evaluated_sub_df_5.at[(33 - 16), 'Q']), (evaluated_sub_df_5.at[(33 - 17), 'Q']), (evaluated_sub_df_5.at[(33 - 18), 'Q']), (evaluated_sub_df_5.at[(33 - 19), 'Q']), 
                                       (evaluated_sub_df_5.at[(33 - 20), 'Q']), (evaluated_sub_df_5.at[(33 - 21), 'Q']), (evaluated_sub_df_5.at[(33 - 22), 'Q']), (evaluated_sub_df_5.at[(33 - 23), 'Q']), (evaluated_sub_df_5.at[(33 - 24), 'Q']), (evaluated_sub_df_5.at[(33 - 25), 'Q']), (evaluated_sub_df_5.at[(33 - 26), 'Q']), (evaluated_sub_df_5.at[(33 - 27), 'Q']), (evaluated_sub_df_5.at[(33 - 28), 'Q']), (evaluated_sub_df_5.at[(33 - 29), 'Q']), 
                                       (evaluated_sub_df_5.at[(33 - 30), 'Q']), (evaluated_sub_df_5.at[(33 - 31), 'Q']), (evaluated_sub_df_5.at[(33 - 32), 'Q']), (evaluated_sub_df_5.at[(33 - 33), 'Q']), ][::-1]
                    #bars_sub_list_5 = [(evaluated_sub_df_5.at[(33 - 0), 'Q']), (evaluated_sub_df_5.at[(33 - 1), 'Q']), (evaluated_sub_df_5.at[(33 - 2), 'Q']), (evaluated_sub_df_5.at[(33 - 3), 'Q']), (evaluated_sub_df_5.at[(33 - 4), 'Q']), (evaluated_sub_df_5.at[(33 - 5), 'Q']), (evaluated_sub_df_5.at[(33 - 6), 'Q']), (evaluated_sub_df_5.at[(33 - 7), 'Q']), (evaluated_sub_df_5.at[(33 - 8), 'Q']), (evaluated_sub_df_5.at[(33 - 9), 'Q']), (evaluated_sub_df_5.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_6 = [(evaluated_sub_df_6.at[(33 - 0), 'Q']), (evaluated_sub_df_6.at[(33 - 1), 'Q']), (evaluated_sub_df_6.at[(33 - 2), 'Q']), (evaluated_sub_df_6.at[(33 - 3), 'Q']), (evaluated_sub_df_6.at[(33 - 4), 'Q']), (evaluated_sub_df_6.at[(33 - 5), 'Q']), (evaluated_sub_df_6.at[(33 - 6), 'Q']), (evaluated_sub_df_6.at[(33 - 7), 'Q']), (evaluated_sub_df_6.at[(33 - 8), 'Q']), (evaluated_sub_df_6.at[(33 - 9), 'Q']), 
                                       (evaluated_sub_df_6.at[(33 - 10), 'Q']), (evaluated_sub_df_6.at[(33 - 11), 'Q']), (evaluated_sub_df_6.at[(33 - 12), 'Q']), (evaluated_sub_df_6.at[(33 - 13), 'Q']), (evaluated_sub_df_6.at[(33 - 14), 'Q']), (evaluated_sub_df_6.at[(33 - 15), 'Q']), (evaluated_sub_df_6.at[(33 - 16), 'Q']), (evaluated_sub_df_6.at[(33 - 17), 'Q']), (evaluated_sub_df_6.at[(33 - 18), 'Q']), (evaluated_sub_df_6.at[(33 - 19), 'Q']), 
                                       (evaluated_sub_df_6.at[(33 - 20), 'Q']), (evaluated_sub_df_6.at[(33 - 21), 'Q']), (evaluated_sub_df_6.at[(33 - 22), 'Q']), (evaluated_sub_df_6.at[(33 - 23), 'Q']), (evaluated_sub_df_6.at[(33 - 24), 'Q']), (evaluated_sub_df_6.at[(33 - 25), 'Q']), (evaluated_sub_df_6.at[(33 - 26), 'Q']), (evaluated_sub_df_6.at[(33 - 27), 'Q']), (evaluated_sub_df_6.at[(33 - 28), 'Q']), (evaluated_sub_df_6.at[(33 - 29), 'Q']), 
                                       (evaluated_sub_df_6.at[(33 - 30), 'Q']), (evaluated_sub_df_6.at[(33 - 31), 'Q']), (evaluated_sub_df_6.at[(33 - 32), 'Q']), (evaluated_sub_df_6.at[(33 - 33), 'Q']), ][::-1]
                    #bars_sub_list_6 = [(evaluated_sub_df_6.at[(33 - 0), 'Q']), (evaluated_sub_df_6.at[(33 - 1), 'Q']), (evaluated_sub_df_6.at[(33 - 2), 'Q']), (evaluated_sub_df_6.at[(33 - 3), 'Q']), (evaluated_sub_df_6.at[(33 - 4), 'Q']), (evaluated_sub_df_6.at[(33 - 5), 'Q']), (evaluated_sub_df_6.at[(33 - 6), 'Q']), (evaluated_sub_df_6.at[(33 - 7), 'Q']), (evaluated_sub_df_6.at[(33 - 8), 'Q']), (evaluated_sub_df_6.at[(33 - 9), 'Q']), (evaluated_sub_df_6.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_7 = [(evaluated_sub_df_7.at[(33 - 0), 'Q']), (evaluated_sub_df_7.at[(33 - 1), 'Q']), (evaluated_sub_df_7.at[(33 - 2), 'Q']), (evaluated_sub_df_7.at[(33 - 3), 'Q']), (evaluated_sub_df_7.at[(33 - 4), 'Q']), (evaluated_sub_df_7.at[(33 - 5), 'Q']), (evaluated_sub_df_7.at[(33 - 6), 'Q']), (evaluated_sub_df_7.at[(33 - 7), 'Q']), (evaluated_sub_df_7.at[(33 - 8), 'Q']), (evaluated_sub_df_7.at[(33 - 9), 'Q']), 
                                       (evaluated_sub_df_7.at[(33 - 10), 'Q']), (evaluated_sub_df_7.at[(33 - 11), 'Q']), (evaluated_sub_df_7.at[(33 - 12), 'Q']), (evaluated_sub_df_7.at[(33 - 13), 'Q']), (evaluated_sub_df_7.at[(33 - 14), 'Q']), (evaluated_sub_df_7.at[(33 - 15), 'Q']), (evaluated_sub_df_7.at[(33 - 16), 'Q']), (evaluated_sub_df_7.at[(33 - 17), 'Q']), (evaluated_sub_df_7.at[(33 - 18), 'Q']), (evaluated_sub_df_7.at[(33 - 19), 'Q']), 
                                       (evaluated_sub_df_7.at[(33 - 20), 'Q']), (evaluated_sub_df_7.at[(33 - 21), 'Q']), (evaluated_sub_df_7.at[(33 - 22), 'Q']), (evaluated_sub_df_7.at[(33 - 23), 'Q']), (evaluated_sub_df_7.at[(33 - 24), 'Q']), (evaluated_sub_df_7.at[(33 - 25), 'Q']), (evaluated_sub_df_7.at[(33 - 26), 'Q']), (evaluated_sub_df_7.at[(33 - 27), 'Q']), (evaluated_sub_df_7.at[(33 - 28), 'Q']), (evaluated_sub_df_7.at[(33 - 29), 'Q']), 
                                       (evaluated_sub_df_7.at[(33 - 30), 'Q']), (evaluated_sub_df_7.at[(33 - 31), 'Q']), (evaluated_sub_df_7.at[(33 - 32), 'Q']), (evaluated_sub_df_7.at[(33 - 33), 'Q']), ][::-1]
                    #bars_sub_list_7 = [(evaluated_sub_df_7.at[(33 - 0), 'Q']), (evaluated_sub_df_7.at[(33 - 1), 'Q']), (evaluated_sub_df_7.at[(33 - 2), 'Q']), (evaluated_sub_df_7.at[(33 - 3), 'Q']), (evaluated_sub_df_7.at[(33 - 4), 'Q']), (evaluated_sub_df_7.at[(33 - 5), 'Q']), (evaluated_sub_df_7.at[(33 - 6), 'Q']), (evaluated_sub_df_7.at[(33 - 7), 'Q']), (evaluated_sub_df_7.at[(33 - 8), 'Q']), (evaluated_sub_df_7.at[(33 - 9), 'Q']), (evaluated_sub_df_7.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_8 = [(evaluated_sub_df_8.at[(33 - 0), 'Q']), (evaluated_sub_df_8.at[(33 - 1), 'Q']), (evaluated_sub_df_8.at[(33 - 2), 'Q']), (evaluated_sub_df_8.at[(33 - 3), 'Q']), (evaluated_sub_df_8.at[(33 - 4), 'Q']), (evaluated_sub_df_8.at[(33 - 5), 'Q']), (evaluated_sub_df_8.at[(33 - 6), 'Q']), (evaluated_sub_df_8.at[(33 - 7), 'Q']), (evaluated_sub_df_8.at[(33 - 8), 'Q']), (evaluated_sub_df_8.at[(33 - 9), 'Q']), 
                                       (evaluated_sub_df_8.at[(33 - 10), 'Q']), (evaluated_sub_df_8.at[(33 - 11), 'Q']), (evaluated_sub_df_8.at[(33 - 12), 'Q']), (evaluated_sub_df_8.at[(33 - 13), 'Q']), (evaluated_sub_df_8.at[(33 - 14), 'Q']), (evaluated_sub_df_8.at[(33 - 15), 'Q']), (evaluated_sub_df_8.at[(33 - 16), 'Q']), (evaluated_sub_df_8.at[(33 - 17), 'Q']), (evaluated_sub_df_8.at[(33 - 18), 'Q']), (evaluated_sub_df_8.at[(33 - 19), 'Q']), 
                                       (evaluated_sub_df_8.at[(33 - 20), 'Q']), (evaluated_sub_df_8.at[(33 - 21), 'Q']), (evaluated_sub_df_8.at[(33 - 22), 'Q']), (evaluated_sub_df_8.at[(33 - 23), 'Q']), (evaluated_sub_df_8.at[(33 - 24), 'Q']), (evaluated_sub_df_8.at[(33 - 25), 'Q']), (evaluated_sub_df_8.at[(33 - 26), 'Q']), (evaluated_sub_df_8.at[(33 - 27), 'Q']), (evaluated_sub_df_8.at[(33 - 28), 'Q']), (evaluated_sub_df_8.at[(33 - 29), 'Q']), 
                                       (evaluated_sub_df_8.at[(33 - 30), 'Q']), (evaluated_sub_df_8.at[(33 - 31), 'Q']), (evaluated_sub_df_8.at[(33 - 32), 'Q']), (evaluated_sub_df_8.at[(33 - 33), 'Q']), ][::-1]
                    #bars_sub_list_8 = [(evaluated_sub_df_8.at[(33 - 0), 'Q']), (evaluated_sub_df_8.at[(33 - 1), 'Q']), (evaluated_sub_df_8.at[(33 - 2), 'Q']), (evaluated_sub_df_8.at[(33 - 3), 'Q']), (evaluated_sub_df_8.at[(33 - 4), 'Q']), (evaluated_sub_df_8.at[(33 - 5), 'Q']), (evaluated_sub_df_8.at[(33 - 6), 'Q']), (evaluated_sub_df_8.at[(33 - 7), 'Q']), (evaluated_sub_df_8.at[(33 - 8), 'Q']), (evaluated_sub_df_8.at[(33 - 9), 'Q']), (evaluated_sub_df_8.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_9 = [(evaluated_sub_df_9.at[(33 - 0), 'Q']), (evaluated_sub_df_9.at[(33 - 1), 'Q']), (evaluated_sub_df_9.at[(33 - 2), 'Q']), (evaluated_sub_df_9.at[(33 - 3), 'Q']), (evaluated_sub_df_9.at[(33 - 4), 'Q']), (evaluated_sub_df_9.at[(33 - 5), 'Q']), (evaluated_sub_df_9.at[(33 - 6), 'Q']), (evaluated_sub_df_9.at[(33 - 7), 'Q']), (evaluated_sub_df_9.at[(33 - 8), 'Q']), (evaluated_sub_df_9.at[(33 - 9), 'Q']), 
                                       (evaluated_sub_df_9.at[(33 - 10), 'Q']), (evaluated_sub_df_9.at[(33 - 11), 'Q']), (evaluated_sub_df_9.at[(33 - 12), 'Q']), (evaluated_sub_df_9.at[(33 - 13), 'Q']), (evaluated_sub_df_9.at[(33 - 14), 'Q']), (evaluated_sub_df_9.at[(33 - 15), 'Q']), (evaluated_sub_df_9.at[(33 - 16), 'Q']), (evaluated_sub_df_9.at[(33 - 17), 'Q']), (evaluated_sub_df_9.at[(33 - 18), 'Q']), (evaluated_sub_df_9.at[(33 - 19), 'Q']), 
                                       (evaluated_sub_df_9.at[(33 - 20), 'Q']), (evaluated_sub_df_9.at[(33 - 21), 'Q']), (evaluated_sub_df_9.at[(33 - 22), 'Q']), (evaluated_sub_df_9.at[(33 - 23), 'Q']), (evaluated_sub_df_9.at[(33 - 24), 'Q']), (evaluated_sub_df_9.at[(33 - 25), 'Q']), (evaluated_sub_df_9.at[(33 - 26), 'Q']), (evaluated_sub_df_9.at[(33 - 27), 'Q']), (evaluated_sub_df_9.at[(33 - 28), 'Q']), (evaluated_sub_df_9.at[(33 - 29), 'Q']), 
                                       (evaluated_sub_df_9.at[(33 - 30), 'Q']), (evaluated_sub_df_9.at[(33 - 31), 'Q']), (evaluated_sub_df_9.at[(33 - 32), 'Q']), (evaluated_sub_df_9.at[(33 - 33), 'Q']), ][::-1]
                    #bars_sub_list_9 = [(evaluated_sub_df_9.at[(33 - 0), 'Q']), (evaluated_sub_df_9.at[(33 - 1), 'Q']), (evaluated_sub_df_9.at[(33 - 2), 'Q']), (evaluated_sub_df_9.at[(33 - 3), 'Q']), (evaluated_sub_df_9.at[(33 - 4), 'Q']), (evaluated_sub_df_9.at[(33 - 5), 'Q']), (evaluated_sub_df_9.at[(33 - 6), 'Q']), (evaluated_sub_df_9.at[(33 - 7), 'Q']), (evaluated_sub_df_9.at[(33 - 8), 'Q']), (evaluated_sub_df_9.at[(33 - 9), 'Q']), (evaluated_sub_df_9.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_10 = [(evaluated_sub_df_10.at[(33 - 0), 'Q']), (evaluated_sub_df_10.at[(33 - 1), 'Q']), (evaluated_sub_df_10.at[(33 - 2), 'Q']), (evaluated_sub_df_10.at[(33 - 3), 'Q']), (evaluated_sub_df_10.at[(33 - 4), 'Q']), (evaluated_sub_df_10.at[(33 - 5), 'Q']), (evaluated_sub_df_10.at[(33 - 6), 'Q']), (evaluated_sub_df_10.at[(33 - 7), 'Q']), (evaluated_sub_df_10.at[(33 - 8), 'Q']), (evaluated_sub_df_10.at[(33 - 9), 'Q']), 
                                       (evaluated_sub_df_10.at[(33 - 10), 'Q']), (evaluated_sub_df_10.at[(33 - 11), 'Q']), (evaluated_sub_df_10.at[(33 - 12), 'Q']), (evaluated_sub_df_10.at[(33 - 13), 'Q']), (evaluated_sub_df_10.at[(33 - 14), 'Q']), (evaluated_sub_df_10.at[(33 - 15), 'Q']), (evaluated_sub_df_10.at[(33 - 16), 'Q']), (evaluated_sub_df_10.at[(33 - 17), 'Q']), (evaluated_sub_df_10.at[(33 - 18), 'Q']), (evaluated_sub_df_10.at[(33 - 19), 'Q']), 
                                       (evaluated_sub_df_10.at[(33 - 20), 'Q']), (evaluated_sub_df_10.at[(33 - 21), 'Q']), (evaluated_sub_df_10.at[(33 - 22), 'Q']), (evaluated_sub_df_10.at[(33 - 23), 'Q']), (evaluated_sub_df_10.at[(33 - 24), 'Q']), (evaluated_sub_df_10.at[(33 - 25), 'Q']), (evaluated_sub_df_10.at[(33 - 26), 'Q']), (evaluated_sub_df_10.at[(33 - 27), 'Q']), (evaluated_sub_df_10.at[(33 - 28), 'Q']), (evaluated_sub_df_10.at[(33 - 29), 'Q']), 
                                       (evaluated_sub_df_10.at[(33 - 30), 'Q']), (evaluated_sub_df_10.at[(33 - 31), 'Q']), (evaluated_sub_df_10.at[(33 - 32), 'Q']), (evaluated_sub_df_10.at[(33 - 33), 'Q']), ][::-1]
                    #bars_sub_list_10 = [(evaluated_sub_df_10.at[(33 - 0), 'Q']), (evaluated_sub_df_10.at[(33 - 1), 'Q']), (evaluated_sub_df_10.at[(33 - 2), 'Q']), (evaluated_sub_df_10.at[(33 - 3), 'Q']), (evaluated_sub_df_10.at[(33 - 4), 'Q']), (evaluated_sub_df_10.at[(33 - 5), 'Q']), (evaluated_sub_df_10.at[(33 - 6), 'Q']), (evaluated_sub_df_10.at[(33 - 7), 'Q']), (evaluated_sub_df_10.at[(33 - 8), 'Q']), (evaluated_sub_df_10.at[(33 - 9), 'Q']), (evaluated_sub_df_10.at[(33 - 10), 'Q']), ][::-1]
                    
                    bar_group_names = ['-33', '-32', '-31', '-30', '-29', '-28', '-27', '-26', '-25', '-24', '-23', '-22', '-21', '-20', '-19', '-18', '-17', '-16', '-15', '-14', '-13', '-12', '-11', '-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', 'split',]
                    #bar_group_names = ['-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', 'split',]
                
                elif 'post' in name_1:
                    #create 10 sub dfs that cover all shifting frames of where a 4Q-repeat can be found up to 10 amino acids from the exon start
                    # set height of bar
                    bars_sub_list_1 = [(evaluated_sub_df_1.at[(33 + 0), 'Q']), (evaluated_sub_df_1.at[(33 + 1), 'Q']), (evaluated_sub_df_1.at[(33 + 2), 'Q']), (evaluated_sub_df_1.at[(33 + 3), 'Q']), (evaluated_sub_df_1.at[(33 + 4), 'Q']), (evaluated_sub_df_1.at[(33 + 5), 'Q']), (evaluated_sub_df_1.at[(33 + 6), 'Q']), (evaluated_sub_df_1.at[(33 + 7), 'Q']), (evaluated_sub_df_1.at[(33 + 8), 'Q']), (evaluated_sub_df_1.at[(33 + 9), 'Q']), 
                                       (evaluated_sub_df_1.at[(33 + 10), 'Q']), (evaluated_sub_df_1.at[(33 + 11), 'Q']), (evaluated_sub_df_1.at[(33 + 12), 'Q']), (evaluated_sub_df_1.at[(33 + 13), 'Q']), (evaluated_sub_df_1.at[(33 + 14), 'Q']), (evaluated_sub_df_1.at[(33 + 15), 'Q']), (evaluated_sub_df_1.at[(33 + 16), 'Q']), (evaluated_sub_df_1.at[(33 + 17), 'Q']), (evaluated_sub_df_1.at[(33 + 18), 'Q']), (evaluated_sub_df_1.at[(33 + 19), 'Q']), 
                                       (evaluated_sub_df_1.at[(33 + 20), 'Q']), (evaluated_sub_df_1.at[(33 + 21), 'Q']), (evaluated_sub_df_1.at[(33 + 22), 'Q']), (evaluated_sub_df_1.at[(33 + 23), 'Q']), (evaluated_sub_df_1.at[(33 + 24), 'Q']), (evaluated_sub_df_1.at[(33 + 25), 'Q']), (evaluated_sub_df_1.at[(33 + 26), 'Q']), (evaluated_sub_df_1.at[(33 + 27), 'Q']), (evaluated_sub_df_1.at[(33 + 28), 'Q']), (evaluated_sub_df_1.at[(33 + 29), 'Q']), 
                                       (evaluated_sub_df_1.at[(33 + 30), 'Q']), (evaluated_sub_df_1.at[(33 + 31), 'Q']), (evaluated_sub_df_1.at[(33 + 32), 'Q']), (evaluated_sub_df_1.at[(33 + 33), 'Q']), ]
                    #bars_sub_list_1 = [(evaluated_sub_df_1.at[(33 + 0), 'Q']), (evaluated_sub_df_1.at[(33 + 1), 'Q']), (evaluated_sub_df_1.at[(33 + 2), 'Q']), (evaluated_sub_df_1.at[(33 + 3), 'Q']), (evaluated_sub_df_1.at[(33 + 4), 'Q']), (evaluated_sub_df_1.at[(33 + 5), 'Q']), (evaluated_sub_df_1.at[(33 + 6), 'Q']), (evaluated_sub_df_1.at[(33 + 7), 'Q']), (evaluated_sub_df_1.at[(33 + 8), 'Q']), (evaluated_sub_df_1.at[(33 + 9), 'Q']), (evaluated_sub_df_1.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_2 = [(evaluated_sub_df_2.at[(33 + 0), 'Q']), (evaluated_sub_df_2.at[(33 + 1), 'Q']), (evaluated_sub_df_2.at[(33 + 2), 'Q']), (evaluated_sub_df_2.at[(33 + 3), 'Q']), (evaluated_sub_df_2.at[(33 + 4), 'Q']), (evaluated_sub_df_2.at[(33 + 5), 'Q']), (evaluated_sub_df_2.at[(33 + 6), 'Q']), (evaluated_sub_df_2.at[(33 + 7), 'Q']), (evaluated_sub_df_2.at[(33 + 8), 'Q']), (evaluated_sub_df_2.at[(33 + 9), 'Q']), 
                                       (evaluated_sub_df_2.at[(33 + 10), 'Q']), (evaluated_sub_df_2.at[(33 + 11), 'Q']), (evaluated_sub_df_2.at[(33 + 12), 'Q']), (evaluated_sub_df_2.at[(33 + 13), 'Q']), (evaluated_sub_df_2.at[(33 + 14), 'Q']), (evaluated_sub_df_2.at[(33 + 15), 'Q']), (evaluated_sub_df_2.at[(33 + 16), 'Q']), (evaluated_sub_df_2.at[(33 + 17), 'Q']), (evaluated_sub_df_2.at[(33 + 18), 'Q']), (evaluated_sub_df_2.at[(33 + 19), 'Q']), 
                                       (evaluated_sub_df_2.at[(33 + 20), 'Q']), (evaluated_sub_df_2.at[(33 + 21), 'Q']), (evaluated_sub_df_2.at[(33 + 22), 'Q']), (evaluated_sub_df_2.at[(33 + 23), 'Q']), (evaluated_sub_df_2.at[(33 + 24), 'Q']), (evaluated_sub_df_2.at[(33 + 25), 'Q']), (evaluated_sub_df_2.at[(33 + 26), 'Q']), (evaluated_sub_df_2.at[(33 + 27), 'Q']), (evaluated_sub_df_2.at[(33 + 28), 'Q']), (evaluated_sub_df_2.at[(33 + 29), 'Q']), 
                                       (evaluated_sub_df_2.at[(33 + 30), 'Q']), (evaluated_sub_df_2.at[(33 + 31), 'Q']), (evaluated_sub_df_2.at[(33 + 32), 'Q']), (evaluated_sub_df_2.at[(33 + 33), 'Q']), ]
                    #bars_sub_list_2 = [(evaluated_sub_df_2.at[(33 + 0), 'Q']), (evaluated_sub_df_2.at[(33 + 1), 'Q']), (evaluated_sub_df_2.at[(33 + 2), 'Q']), (evaluated_sub_df_2.at[(33 + 3), 'Q']), (evaluated_sub_df_2.at[(33 + 4), 'Q']), (evaluated_sub_df_2.at[(33 + 5), 'Q']), (evaluated_sub_df_2.at[(33 + 6), 'Q']), (evaluated_sub_df_2.at[(33 + 7), 'Q']), (evaluated_sub_df_2.at[(33 + 8), 'Q']), (evaluated_sub_df_2.at[(33 + 9), 'Q']), (evaluated_sub_df_2.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_3 = [(evaluated_sub_df_3.at[(33 + 0), 'Q']), (evaluated_sub_df_3.at[(33 + 1), 'Q']), (evaluated_sub_df_3.at[(33 + 2), 'Q']), (evaluated_sub_df_3.at[(33 + 3), 'Q']), (evaluated_sub_df_3.at[(33 + 4), 'Q']), (evaluated_sub_df_3.at[(33 + 5), 'Q']), (evaluated_sub_df_3.at[(33 + 6), 'Q']), (evaluated_sub_df_3.at[(33 + 7), 'Q']), (evaluated_sub_df_3.at[(33 + 8), 'Q']), (evaluated_sub_df_3.at[(33 + 9), 'Q']), 
                                       (evaluated_sub_df_3.at[(33 + 10), 'Q']), (evaluated_sub_df_3.at[(33 + 11), 'Q']), (evaluated_sub_df_3.at[(33 + 12), 'Q']), (evaluated_sub_df_3.at[(33 + 13), 'Q']), (evaluated_sub_df_3.at[(33 + 14), 'Q']), (evaluated_sub_df_3.at[(33 + 15), 'Q']), (evaluated_sub_df_3.at[(33 + 16), 'Q']), (evaluated_sub_df_3.at[(33 + 17), 'Q']), (evaluated_sub_df_3.at[(33 + 18), 'Q']), (evaluated_sub_df_3.at[(33 + 19), 'Q']), 
                                       (evaluated_sub_df_3.at[(33 + 20), 'Q']), (evaluated_sub_df_3.at[(33 + 21), 'Q']), (evaluated_sub_df_3.at[(33 + 22), 'Q']), (evaluated_sub_df_3.at[(33 + 23), 'Q']), (evaluated_sub_df_3.at[(33 + 24), 'Q']), (evaluated_sub_df_3.at[(33 + 25), 'Q']), (evaluated_sub_df_3.at[(33 + 26), 'Q']), (evaluated_sub_df_3.at[(33 + 27), 'Q']), (evaluated_sub_df_3.at[(33 + 28), 'Q']), (evaluated_sub_df_3.at[(33 + 29), 'Q']), 
                                       (evaluated_sub_df_3.at[(33 + 30), 'Q']), (evaluated_sub_df_3.at[(33 + 31), 'Q']), (evaluated_sub_df_3.at[(33 + 32), 'Q']), (evaluated_sub_df_3.at[(33 + 33), 'Q']), ]
                    #bars_sub_list_3 = [(evaluated_sub_df_3.at[(33 + 0), 'Q']), (evaluated_sub_df_3.at[(33 + 1), 'Q']), (evaluated_sub_df_3.at[(33 + 2), 'Q']), (evaluated_sub_df_3.at[(33 + 3), 'Q']), (evaluated_sub_df_3.at[(33 + 4), 'Q']), (evaluated_sub_df_3.at[(33 + 5), 'Q']), (evaluated_sub_df_3.at[(33 + 6), 'Q']), (evaluated_sub_df_3.at[(33 + 7), 'Q']), (evaluated_sub_df_3.at[(33 + 8), 'Q']), (evaluated_sub_df_3.at[(33 + 9), 'Q']), (evaluated_sub_df_3.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_4 = [(evaluated_sub_df_4.at[(33 + 0), 'Q']), (evaluated_sub_df_4.at[(33 + 1), 'Q']), (evaluated_sub_df_4.at[(33 + 2), 'Q']), (evaluated_sub_df_4.at[(33 + 3), 'Q']), (evaluated_sub_df_4.at[(33 + 4), 'Q']), (evaluated_sub_df_4.at[(33 + 5), 'Q']), (evaluated_sub_df_4.at[(33 + 6), 'Q']), (evaluated_sub_df_4.at[(33 + 7), 'Q']), (evaluated_sub_df_4.at[(33 + 8), 'Q']), (evaluated_sub_df_4.at[(33 + 9), 'Q']), 
                                       (evaluated_sub_df_4.at[(33 + 10), 'Q']), (evaluated_sub_df_4.at[(33 + 11), 'Q']), (evaluated_sub_df_4.at[(33 + 12), 'Q']), (evaluated_sub_df_4.at[(33 + 13), 'Q']), (evaluated_sub_df_4.at[(33 + 14), 'Q']), (evaluated_sub_df_4.at[(33 + 15), 'Q']), (evaluated_sub_df_4.at[(33 + 16), 'Q']), (evaluated_sub_df_4.at[(33 + 17), 'Q']), (evaluated_sub_df_4.at[(33 + 18), 'Q']), (evaluated_sub_df_4.at[(33 + 19), 'Q']), 
                                       (evaluated_sub_df_4.at[(33 + 20), 'Q']), (evaluated_sub_df_4.at[(33 + 21), 'Q']), (evaluated_sub_df_4.at[(33 + 22), 'Q']), (evaluated_sub_df_4.at[(33 + 23), 'Q']), (evaluated_sub_df_4.at[(33 + 24), 'Q']), (evaluated_sub_df_4.at[(33 + 25), 'Q']), (evaluated_sub_df_4.at[(33 + 26), 'Q']), (evaluated_sub_df_4.at[(33 + 27), 'Q']), (evaluated_sub_df_4.at[(33 + 28), 'Q']), (evaluated_sub_df_4.at[(33 + 29), 'Q']), 
                                       (evaluated_sub_df_4.at[(33 + 30), 'Q']), (evaluated_sub_df_4.at[(33 + 31), 'Q']), (evaluated_sub_df_4.at[(33 + 32), 'Q']), (evaluated_sub_df_4.at[(33 + 33), 'Q']), ]
                    #bars_sub_list_4 = [(evaluated_sub_df_4.at[(33 + 0), 'Q']), (evaluated_sub_df_4.at[(33 + 1), 'Q']), (evaluated_sub_df_4.at[(33 + 2), 'Q']), (evaluated_sub_df_4.at[(33 + 3), 'Q']), (evaluated_sub_df_4.at[(33 + 4), 'Q']), (evaluated_sub_df_4.at[(33 + 5), 'Q']), (evaluated_sub_df_4.at[(33 + 6), 'Q']), (evaluated_sub_df_4.at[(33 + 7), 'Q']), (evaluated_sub_df_4.at[(33 + 8), 'Q']), (evaluated_sub_df_4.at[(33 + 9), 'Q']), (evaluated_sub_df_4.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_5 = [(evaluated_sub_df_5.at[(33 + 0), 'Q']), (evaluated_sub_df_5.at[(33 + 1), 'Q']), (evaluated_sub_df_5.at[(33 + 2), 'Q']), (evaluated_sub_df_5.at[(33 + 3), 'Q']), (evaluated_sub_df_5.at[(33 + 4), 'Q']), (evaluated_sub_df_5.at[(33 + 5), 'Q']), (evaluated_sub_df_5.at[(33 + 6), 'Q']), (evaluated_sub_df_5.at[(33 + 7), 'Q']), (evaluated_sub_df_5.at[(33 + 8), 'Q']), (evaluated_sub_df_5.at[(33 + 9), 'Q']), 
                                       (evaluated_sub_df_5.at[(33 + 10), 'Q']), (evaluated_sub_df_5.at[(33 + 11), 'Q']), (evaluated_sub_df_5.at[(33 + 12), 'Q']), (evaluated_sub_df_5.at[(33 + 13), 'Q']), (evaluated_sub_df_5.at[(33 + 14), 'Q']), (evaluated_sub_df_5.at[(33 + 15), 'Q']), (evaluated_sub_df_5.at[(33 + 16), 'Q']), (evaluated_sub_df_5.at[(33 + 17), 'Q']), (evaluated_sub_df_5.at[(33 + 18), 'Q']), (evaluated_sub_df_5.at[(33 + 19), 'Q']), 
                                       (evaluated_sub_df_5.at[(33 + 20), 'Q']), (evaluated_sub_df_5.at[(33 + 21), 'Q']), (evaluated_sub_df_5.at[(33 + 22), 'Q']), (evaluated_sub_df_5.at[(33 + 23), 'Q']), (evaluated_sub_df_5.at[(33 + 24), 'Q']), (evaluated_sub_df_5.at[(33 + 25), 'Q']), (evaluated_sub_df_5.at[(33 + 26), 'Q']), (evaluated_sub_df_5.at[(33 + 27), 'Q']), (evaluated_sub_df_5.at[(33 + 28), 'Q']), (evaluated_sub_df_5.at[(33 + 29), 'Q']), 
                                       (evaluated_sub_df_5.at[(33 + 30), 'Q']), (evaluated_sub_df_5.at[(33 + 31), 'Q']), (evaluated_sub_df_5.at[(33 + 32), 'Q']), (evaluated_sub_df_5.at[(33 + 33), 'Q']), ]
                    #bars_sub_list_5 = [(evaluated_sub_df_5.at[(33 + 0), 'Q']), (evaluated_sub_df_5.at[(33 + 1), 'Q']), (evaluated_sub_df_5.at[(33 + 2), 'Q']), (evaluated_sub_df_5.at[(33 + 3), 'Q']), (evaluated_sub_df_5.at[(33 + 4), 'Q']), (evaluated_sub_df_5.at[(33 + 5), 'Q']), (evaluated_sub_df_5.at[(33 + 6), 'Q']), (evaluated_sub_df_5.at[(33 + 7), 'Q']), (evaluated_sub_df_5.at[(33 + 8), 'Q']), (evaluated_sub_df_5.at[(33 + 9), 'Q']), (evaluated_sub_df_5.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_6 = [(evaluated_sub_df_6.at[(33 + 0), 'Q']), (evaluated_sub_df_6.at[(33 + 1), 'Q']), (evaluated_sub_df_6.at[(33 + 2), 'Q']), (evaluated_sub_df_6.at[(33 + 3), 'Q']), (evaluated_sub_df_6.at[(33 + 4), 'Q']), (evaluated_sub_df_6.at[(33 + 5), 'Q']), (evaluated_sub_df_6.at[(33 + 6), 'Q']), (evaluated_sub_df_6.at[(33 + 7), 'Q']), (evaluated_sub_df_6.at[(33 + 8), 'Q']), (evaluated_sub_df_6.at[(33 + 9), 'Q']), 
                                       (evaluated_sub_df_6.at[(33 + 10), 'Q']), (evaluated_sub_df_6.at[(33 + 11), 'Q']), (evaluated_sub_df_6.at[(33 + 12), 'Q']), (evaluated_sub_df_6.at[(33 + 13), 'Q']), (evaluated_sub_df_6.at[(33 + 14), 'Q']), (evaluated_sub_df_6.at[(33 + 15), 'Q']), (evaluated_sub_df_6.at[(33 + 16), 'Q']), (evaluated_sub_df_6.at[(33 + 17), 'Q']), (evaluated_sub_df_6.at[(33 + 18), 'Q']), (evaluated_sub_df_6.at[(33 + 19), 'Q']), 
                                       (evaluated_sub_df_6.at[(33 + 20), 'Q']), (evaluated_sub_df_6.at[(33 + 21), 'Q']), (evaluated_sub_df_6.at[(33 + 22), 'Q']), (evaluated_sub_df_6.at[(33 + 23), 'Q']), (evaluated_sub_df_6.at[(33 + 24), 'Q']), (evaluated_sub_df_6.at[(33 + 25), 'Q']), (evaluated_sub_df_6.at[(33 + 26), 'Q']), (evaluated_sub_df_6.at[(33 + 27), 'Q']), (evaluated_sub_df_6.at[(33 + 28), 'Q']), (evaluated_sub_df_6.at[(33 + 29), 'Q']), 
                                       (evaluated_sub_df_6.at[(33 + 30), 'Q']), (evaluated_sub_df_6.at[(33 + 31), 'Q']), (evaluated_sub_df_6.at[(33 + 32), 'Q']), (evaluated_sub_df_6.at[(33 + 33), 'Q']), ]
                    #bars_sub_list_6 = [(evaluated_sub_df_6.at[(33 + 0), 'Q']), (evaluated_sub_df_6.at[(33 + 1), 'Q']), (evaluated_sub_df_6.at[(33 + 2), 'Q']), (evaluated_sub_df_6.at[(33 + 3), 'Q']), (evaluated_sub_df_6.at[(33 + 4), 'Q']), (evaluated_sub_df_6.at[(33 + 5), 'Q']), (evaluated_sub_df_6.at[(33 + 6), 'Q']), (evaluated_sub_df_6.at[(33 + 7), 'Q']), (evaluated_sub_df_6.at[(33 + 8), 'Q']), (evaluated_sub_df_6.at[(33 + 9), 'Q']), (evaluated_sub_df_6.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_7 = [(evaluated_sub_df_7.at[(33 + 0), 'Q']), (evaluated_sub_df_7.at[(33 + 1), 'Q']), (evaluated_sub_df_7.at[(33 + 2), 'Q']), (evaluated_sub_df_7.at[(33 + 3), 'Q']), (evaluated_sub_df_7.at[(33 + 4), 'Q']), (evaluated_sub_df_7.at[(33 + 5), 'Q']), (evaluated_sub_df_7.at[(33 + 6), 'Q']), (evaluated_sub_df_7.at[(33 + 7), 'Q']), (evaluated_sub_df_7.at[(33 + 8), 'Q']), (evaluated_sub_df_7.at[(33 + 9), 'Q']), 
                                       (evaluated_sub_df_7.at[(33 + 10), 'Q']), (evaluated_sub_df_7.at[(33 + 11), 'Q']), (evaluated_sub_df_7.at[(33 + 12), 'Q']), (evaluated_sub_df_7.at[(33 + 13), 'Q']), (evaluated_sub_df_7.at[(33 + 14), 'Q']), (evaluated_sub_df_7.at[(33 + 15), 'Q']), (evaluated_sub_df_7.at[(33 + 16), 'Q']), (evaluated_sub_df_7.at[(33 + 17), 'Q']), (evaluated_sub_df_7.at[(33 + 18), 'Q']), (evaluated_sub_df_7.at[(33 + 19), 'Q']), 
                                       (evaluated_sub_df_7.at[(33 + 20), 'Q']), (evaluated_sub_df_7.at[(33 + 21), 'Q']), (evaluated_sub_df_7.at[(33 + 22), 'Q']), (evaluated_sub_df_7.at[(33 + 23), 'Q']), (evaluated_sub_df_7.at[(33 + 24), 'Q']), (evaluated_sub_df_7.at[(33 + 25), 'Q']), (evaluated_sub_df_7.at[(33 + 26), 'Q']), (evaluated_sub_df_7.at[(33 + 27), 'Q']), (evaluated_sub_df_7.at[(33 + 28), 'Q']), (evaluated_sub_df_7.at[(33 + 29), 'Q']), 
                                       (evaluated_sub_df_7.at[(33 + 30), 'Q']), (evaluated_sub_df_7.at[(33 + 31), 'Q']), (evaluated_sub_df_7.at[(33 + 32), 'Q']), (evaluated_sub_df_7.at[(33 + 33), 'Q']), ]
                    #bars_sub_list_7 = [(evaluated_sub_df_7.at[(33 + 0), 'Q']), (evaluated_sub_df_7.at[(33 + 1), 'Q']), (evaluated_sub_df_7.at[(33 + 2), 'Q']), (evaluated_sub_df_7.at[(33 + 3), 'Q']), (evaluated_sub_df_7.at[(33 + 4), 'Q']), (evaluated_sub_df_7.at[(33 + 5), 'Q']), (evaluated_sub_df_7.at[(33 + 6), 'Q']), (evaluated_sub_df_7.at[(33 + 7), 'Q']), (evaluated_sub_df_7.at[(33 + 8), 'Q']), (evaluated_sub_df_7.at[(33 + 9), 'Q']), (evaluated_sub_df_7.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_8 = [(evaluated_sub_df_8.at[(33 + 0), 'Q']), (evaluated_sub_df_8.at[(33 + 1), 'Q']), (evaluated_sub_df_8.at[(33 + 2), 'Q']), (evaluated_sub_df_8.at[(33 + 3), 'Q']), (evaluated_sub_df_8.at[(33 + 4), 'Q']), (evaluated_sub_df_8.at[(33 + 5), 'Q']), (evaluated_sub_df_8.at[(33 + 6), 'Q']), (evaluated_sub_df_8.at[(33 + 7), 'Q']), (evaluated_sub_df_8.at[(33 + 8), 'Q']), (evaluated_sub_df_8.at[(33 + 9), 'Q']), 
                                       (evaluated_sub_df_8.at[(33 + 10), 'Q']), (evaluated_sub_df_8.at[(33 + 11), 'Q']), (evaluated_sub_df_8.at[(33 + 12), 'Q']), (evaluated_sub_df_8.at[(33 + 13), 'Q']), (evaluated_sub_df_8.at[(33 + 14), 'Q']), (evaluated_sub_df_8.at[(33 + 15), 'Q']), (evaluated_sub_df_8.at[(33 + 16), 'Q']), (evaluated_sub_df_8.at[(33 + 17), 'Q']), (evaluated_sub_df_8.at[(33 + 18), 'Q']), (evaluated_sub_df_8.at[(33 + 19), 'Q']), 
                                       (evaluated_sub_df_8.at[(33 + 20), 'Q']), (evaluated_sub_df_8.at[(33 + 21), 'Q']), (evaluated_sub_df_8.at[(33 + 22), 'Q']), (evaluated_sub_df_8.at[(33 + 23), 'Q']), (evaluated_sub_df_8.at[(33 + 24), 'Q']), (evaluated_sub_df_8.at[(33 + 25), 'Q']), (evaluated_sub_df_8.at[(33 + 26), 'Q']), (evaluated_sub_df_8.at[(33 + 27), 'Q']), (evaluated_sub_df_8.at[(33 + 28), 'Q']), (evaluated_sub_df_8.at[(33 + 29), 'Q']), 
                                       (evaluated_sub_df_8.at[(33 + 30), 'Q']), (evaluated_sub_df_8.at[(33 + 31), 'Q']), (evaluated_sub_df_8.at[(33 + 32), 'Q']), (evaluated_sub_df_8.at[(33 + 33), 'Q']), ]
                    #bars_sub_list_8 = [(evaluated_sub_df_8.at[(33 + 0), 'Q']), (evaluated_sub_df_8.at[(33 + 1), 'Q']), (evaluated_sub_df_8.at[(33 + 2), 'Q']), (evaluated_sub_df_8.at[(33 + 3), 'Q']), (evaluated_sub_df_8.at[(33 + 4), 'Q']), (evaluated_sub_df_8.at[(33 + 5), 'Q']), (evaluated_sub_df_8.at[(33 + 6), 'Q']), (evaluated_sub_df_8.at[(33 + 7), 'Q']), (evaluated_sub_df_8.at[(33 + 8), 'Q']), (evaluated_sub_df_8.at[(33 + 9), 'Q']), (evaluated_sub_df_8.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_9 = [(evaluated_sub_df_9.at[(33 + 0), 'Q']), (evaluated_sub_df_9.at[(33 + 1), 'Q']), (evaluated_sub_df_9.at[(33 + 2), 'Q']), (evaluated_sub_df_9.at[(33 + 3), 'Q']), (evaluated_sub_df_9.at[(33 + 4), 'Q']), (evaluated_sub_df_9.at[(33 + 5), 'Q']), (evaluated_sub_df_9.at[(33 + 6), 'Q']), (evaluated_sub_df_9.at[(33 + 7), 'Q']), (evaluated_sub_df_9.at[(33 + 8), 'Q']), (evaluated_sub_df_9.at[(33 + 9), 'Q']), 
                                       (evaluated_sub_df_9.at[(33 + 10), 'Q']), (evaluated_sub_df_9.at[(33 + 11), 'Q']), (evaluated_sub_df_9.at[(33 + 12), 'Q']), (evaluated_sub_df_9.at[(33 + 13), 'Q']), (evaluated_sub_df_9.at[(33 + 14), 'Q']), (evaluated_sub_df_9.at[(33 + 15), 'Q']), (evaluated_sub_df_9.at[(33 + 16), 'Q']), (evaluated_sub_df_9.at[(33 + 17), 'Q']), (evaluated_sub_df_9.at[(33 + 18), 'Q']), (evaluated_sub_df_9.at[(33 + 19), 'Q']), 
                                       (evaluated_sub_df_9.at[(33 + 20), 'Q']), (evaluated_sub_df_9.at[(33 + 21), 'Q']), (evaluated_sub_df_9.at[(33 + 22), 'Q']), (evaluated_sub_df_9.at[(33 + 23), 'Q']), (evaluated_sub_df_9.at[(33 + 24), 'Q']), (evaluated_sub_df_9.at[(33 + 25), 'Q']), (evaluated_sub_df_9.at[(33 + 26), 'Q']), (evaluated_sub_df_9.at[(33 + 27), 'Q']), (evaluated_sub_df_9.at[(33 + 28), 'Q']), (evaluated_sub_df_9.at[(33 + 29), 'Q']), 
                                       (evaluated_sub_df_9.at[(33 + 30), 'Q']), (evaluated_sub_df_9.at[(33 + 31), 'Q']), (evaluated_sub_df_9.at[(33 + 32), 'Q']), (evaluated_sub_df_9.at[(33 + 33), 'Q']), ]
                    #bars_sub_list_9 = [(evaluated_sub_df_9.at[(33 + 0), 'Q']), (evaluated_sub_df_9.at[(33 + 1), 'Q']), (evaluated_sub_df_9.at[(33 + 2), 'Q']), (evaluated_sub_df_9.at[(33 + 3), 'Q']), (evaluated_sub_df_9.at[(33 + 4), 'Q']), (evaluated_sub_df_9.at[(33 + 5), 'Q']), (evaluated_sub_df_9.at[(33 + 6), 'Q']), (evaluated_sub_df_9.at[(33 + 7), 'Q']), (evaluated_sub_df_9.at[(33 + 8), 'Q']), (evaluated_sub_df_9.at[(33 + 9), 'Q']), (evaluated_sub_df_9.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_10 = [(evaluated_sub_df_10.at[(33 + 0), 'Q']), (evaluated_sub_df_10.at[(33 + 1), 'Q']), (evaluated_sub_df_10.at[(33 + 2), 'Q']), (evaluated_sub_df_10.at[(33 + 3), 'Q']), (evaluated_sub_df_10.at[(33 + 4), 'Q']), (evaluated_sub_df_10.at[(33 + 5), 'Q']), (evaluated_sub_df_10.at[(33 + 6), 'Q']), (evaluated_sub_df_10.at[(33 + 7), 'Q']), (evaluated_sub_df_10.at[(33 + 8), 'Q']), (evaluated_sub_df_10.at[(33 + 9), 'Q']), 
                                       (evaluated_sub_df_10.at[(33 + 10), 'Q']), (evaluated_sub_df_10.at[(33 + 11), 'Q']), (evaluated_sub_df_10.at[(33 + 12), 'Q']), (evaluated_sub_df_10.at[(33 + 13), 'Q']), (evaluated_sub_df_10.at[(33 + 14), 'Q']), (evaluated_sub_df_10.at[(33 + 15), 'Q']), (evaluated_sub_df_10.at[(33 + 16), 'Q']), (evaluated_sub_df_10.at[(33 + 17), 'Q']), (evaluated_sub_df_10.at[(33 + 18), 'Q']), (evaluated_sub_df_10.at[(33 + 19), 'Q']), 
                                       (evaluated_sub_df_10.at[(33 + 20), 'Q']), (evaluated_sub_df_10.at[(33 + 21), 'Q']), (evaluated_sub_df_10.at[(33 + 22), 'Q']), (evaluated_sub_df_10.at[(33 + 23), 'Q']), (evaluated_sub_df_10.at[(33 + 24), 'Q']), (evaluated_sub_df_10.at[(33 + 25), 'Q']), (evaluated_sub_df_10.at[(33 + 26), 'Q']), (evaluated_sub_df_10.at[(33 + 27), 'Q']), (evaluated_sub_df_10.at[(33 + 28), 'Q']), (evaluated_sub_df_10.at[(33 + 29), 'Q']), 
                                       (evaluated_sub_df_10.at[(33 + 30), 'Q']), (evaluated_sub_df_10.at[(33 + 31), 'Q']), (evaluated_sub_df_10.at[(33 + 32), 'Q']), (evaluated_sub_df_10.at[(33 + 33), 'Q']), ]
                    #bars_sub_list_10 = [(evaluated_sub_df_10.at[(33 + 0), 'Q']), (evaluated_sub_df_10.at[(33 + 1), 'Q']), (evaluated_sub_df_10.at[(33 + 2), 'Q']), (evaluated_sub_df_10.at[(33 + 3), 'Q']), (evaluated_sub_df_10.at[(33 + 4), 'Q']), (evaluated_sub_df_10.at[(33 + 5), 'Q']), (evaluated_sub_df_10.at[(33 + 6), 'Q']), (evaluated_sub_df_10.at[(33 + 7), 'Q']), (evaluated_sub_df_10.at[(33 + 8), 'Q']), (evaluated_sub_df_10.at[(33 + 9), 'Q']), (evaluated_sub_df_10.at[(33 + 10), 'Q']), ]
    
                    bar_group_names = ['split', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33',]
                    #bar_group_names = ['split', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',]
                

            elif num_acids == 10:
                
                if 'prae' in name_1:
                    #create 10 sub dfs that cover all shifting frames of where a 4Q-repeat can be found up to 10 amino acids up to the exon stop
                    # set height of bar
                    
                    bars_sub_list_1 = [(evaluated_sub_df_1.at[(33 - 0), 'Q']), (evaluated_sub_df_1.at[(33 - 1), 'Q']), (evaluated_sub_df_1.at[(33 - 2), 'Q']), (evaluated_sub_df_1.at[(33 - 3), 'Q']), (evaluated_sub_df_1.at[(33 - 4), 'Q']), (evaluated_sub_df_1.at[(33 - 5), 'Q']), (evaluated_sub_df_1.at[(33 - 6), 'Q']), (evaluated_sub_df_1.at[(33 - 7), 'Q']), (evaluated_sub_df_1.at[(33 - 8), 'Q']), (evaluated_sub_df_1.at[(33 - 9), 'Q']), (evaluated_sub_df_1.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_2 = [(evaluated_sub_df_2.at[(33 - 0), 'Q']), (evaluated_sub_df_2.at[(33 - 1), 'Q']), (evaluated_sub_df_2.at[(33 - 2), 'Q']), (evaluated_sub_df_2.at[(33 - 3), 'Q']), (evaluated_sub_df_2.at[(33 - 4), 'Q']), (evaluated_sub_df_2.at[(33 - 5), 'Q']), (evaluated_sub_df_2.at[(33 - 6), 'Q']), (evaluated_sub_df_2.at[(33 - 7), 'Q']), (evaluated_sub_df_2.at[(33 - 8), 'Q']), (evaluated_sub_df_2.at[(33 - 9), 'Q']), (evaluated_sub_df_2.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_3 = [(evaluated_sub_df_3.at[(33 - 0), 'Q']), (evaluated_sub_df_3.at[(33 - 1), 'Q']), (evaluated_sub_df_3.at[(33 - 2), 'Q']), (evaluated_sub_df_3.at[(33 - 3), 'Q']), (evaluated_sub_df_3.at[(33 - 4), 'Q']), (evaluated_sub_df_3.at[(33 - 5), 'Q']), (evaluated_sub_df_3.at[(33 - 6), 'Q']), (evaluated_sub_df_3.at[(33 - 7), 'Q']), (evaluated_sub_df_3.at[(33 - 8), 'Q']), (evaluated_sub_df_3.at[(33 - 9), 'Q']), (evaluated_sub_df_3.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_4 = [(evaluated_sub_df_4.at[(33 - 0), 'Q']), (evaluated_sub_df_4.at[(33 - 1), 'Q']), (evaluated_sub_df_4.at[(33 - 2), 'Q']), (evaluated_sub_df_4.at[(33 - 3), 'Q']), (evaluated_sub_df_4.at[(33 - 4), 'Q']), (evaluated_sub_df_4.at[(33 - 5), 'Q']), (evaluated_sub_df_4.at[(33 - 6), 'Q']), (evaluated_sub_df_4.at[(33 - 7), 'Q']), (evaluated_sub_df_4.at[(33 - 8), 'Q']), (evaluated_sub_df_4.at[(33 - 9), 'Q']), (evaluated_sub_df_4.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_5 = [(evaluated_sub_df_5.at[(33 - 0), 'Q']), (evaluated_sub_df_5.at[(33 - 1), 'Q']), (evaluated_sub_df_5.at[(33 - 2), 'Q']), (evaluated_sub_df_5.at[(33 - 3), 'Q']), (evaluated_sub_df_5.at[(33 - 4), 'Q']), (evaluated_sub_df_5.at[(33 - 5), 'Q']), (evaluated_sub_df_5.at[(33 - 6), 'Q']), (evaluated_sub_df_5.at[(33 - 7), 'Q']), (evaluated_sub_df_5.at[(33 - 8), 'Q']), (evaluated_sub_df_5.at[(33 - 9), 'Q']), (evaluated_sub_df_5.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_6 = [(evaluated_sub_df_6.at[(33 - 0), 'Q']), (evaluated_sub_df_6.at[(33 - 1), 'Q']), (evaluated_sub_df_6.at[(33 - 2), 'Q']), (evaluated_sub_df_6.at[(33 - 3), 'Q']), (evaluated_sub_df_6.at[(33 - 4), 'Q']), (evaluated_sub_df_6.at[(33 - 5), 'Q']), (evaluated_sub_df_6.at[(33 - 6), 'Q']), (evaluated_sub_df_6.at[(33 - 7), 'Q']), (evaluated_sub_df_6.at[(33 - 8), 'Q']), (evaluated_sub_df_6.at[(33 - 9), 'Q']), (evaluated_sub_df_6.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_7 = [(evaluated_sub_df_7.at[(33 - 0), 'Q']), (evaluated_sub_df_7.at[(33 - 1), 'Q']), (evaluated_sub_df_7.at[(33 - 2), 'Q']), (evaluated_sub_df_7.at[(33 - 3), 'Q']), (evaluated_sub_df_7.at[(33 - 4), 'Q']), (evaluated_sub_df_7.at[(33 - 5), 'Q']), (evaluated_sub_df_7.at[(33 - 6), 'Q']), (evaluated_sub_df_7.at[(33 - 7), 'Q']), (evaluated_sub_df_7.at[(33 - 8), 'Q']), (evaluated_sub_df_7.at[(33 - 9), 'Q']), (evaluated_sub_df_7.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_8 = [(evaluated_sub_df_8.at[(33 - 0), 'Q']), (evaluated_sub_df_8.at[(33 - 1), 'Q']), (evaluated_sub_df_8.at[(33 - 2), 'Q']), (evaluated_sub_df_8.at[(33 - 3), 'Q']), (evaluated_sub_df_8.at[(33 - 4), 'Q']), (evaluated_sub_df_8.at[(33 - 5), 'Q']), (evaluated_sub_df_8.at[(33 - 6), 'Q']), (evaluated_sub_df_8.at[(33 - 7), 'Q']), (evaluated_sub_df_8.at[(33 - 8), 'Q']), (evaluated_sub_df_8.at[(33 - 9), 'Q']), (evaluated_sub_df_8.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_9 = [(evaluated_sub_df_9.at[(33 - 0), 'Q']), (evaluated_sub_df_9.at[(33 - 1), 'Q']), (evaluated_sub_df_9.at[(33 - 2), 'Q']), (evaluated_sub_df_9.at[(33 - 3), 'Q']), (evaluated_sub_df_9.at[(33 - 4), 'Q']), (evaluated_sub_df_9.at[(33 - 5), 'Q']), (evaluated_sub_df_9.at[(33 - 6), 'Q']), (evaluated_sub_df_9.at[(33 - 7), 'Q']), (evaluated_sub_df_9.at[(33 - 8), 'Q']), (evaluated_sub_df_9.at[(33 - 9), 'Q']), (evaluated_sub_df_9.at[(33 - 10), 'Q']), ][::-1]
                    
                    bars_sub_list_10 = [(evaluated_sub_df_10.at[(33 - 0), 'Q']), (evaluated_sub_df_10.at[(33 - 1), 'Q']), (evaluated_sub_df_10.at[(33 - 2), 'Q']), (evaluated_sub_df_10.at[(33 - 3), 'Q']), (evaluated_sub_df_10.at[(33 - 4), 'Q']), (evaluated_sub_df_10.at[(33 - 5), 'Q']), (evaluated_sub_df_10.at[(33 - 6), 'Q']), (evaluated_sub_df_10.at[(33 - 7), 'Q']), (evaluated_sub_df_10.at[(33 - 8), 'Q']), (evaluated_sub_df_10.at[(33 - 9), 'Q']), (evaluated_sub_df_10.at[(33 - 10), 'Q']), ][::-1]
                    
                    bar_group_names = ['-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', 'split',]
                
                elif 'post' in name_1:
                    #create 10 sub dfs that cover all shifting frames of where a 4Q-repeat can be found up to 10 amino acids from the exon start
                    # set height of bar
                    bars_sub_list_1 = [(evaluated_sub_df_1.at[(33 + 0), 'Q']), (evaluated_sub_df_1.at[(33 + 1), 'Q']), (evaluated_sub_df_1.at[(33 + 2), 'Q']), (evaluated_sub_df_1.at[(33 + 3), 'Q']), (evaluated_sub_df_1.at[(33 + 4), 'Q']), (evaluated_sub_df_1.at[(33 + 5), 'Q']), (evaluated_sub_df_1.at[(33 + 6), 'Q']), (evaluated_sub_df_1.at[(33 + 7), 'Q']), (evaluated_sub_df_1.at[(33 + 8), 'Q']), (evaluated_sub_df_1.at[(33 + 9), 'Q']), (evaluated_sub_df_1.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_2 = [(evaluated_sub_df_2.at[(33 + 0), 'Q']), (evaluated_sub_df_2.at[(33 + 1), 'Q']), (evaluated_sub_df_2.at[(33 + 2), 'Q']), (evaluated_sub_df_2.at[(33 + 3), 'Q']), (evaluated_sub_df_2.at[(33 + 4), 'Q']), (evaluated_sub_df_2.at[(33 + 5), 'Q']), (evaluated_sub_df_2.at[(33 + 6), 'Q']), (evaluated_sub_df_2.at[(33 + 7), 'Q']), (evaluated_sub_df_2.at[(33 + 8), 'Q']), (evaluated_sub_df_2.at[(33 + 9), 'Q']), (evaluated_sub_df_2.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_3 = [(evaluated_sub_df_3.at[(33 + 0), 'Q']), (evaluated_sub_df_3.at[(33 + 1), 'Q']), (evaluated_sub_df_3.at[(33 + 2), 'Q']), (evaluated_sub_df_3.at[(33 + 3), 'Q']), (evaluated_sub_df_3.at[(33 + 4), 'Q']), (evaluated_sub_df_3.at[(33 + 5), 'Q']), (evaluated_sub_df_3.at[(33 + 6), 'Q']), (evaluated_sub_df_3.at[(33 + 7), 'Q']), (evaluated_sub_df_3.at[(33 + 8), 'Q']), (evaluated_sub_df_3.at[(33 + 9), 'Q']), (evaluated_sub_df_3.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_4 = [(evaluated_sub_df_4.at[(33 + 0), 'Q']), (evaluated_sub_df_4.at[(33 + 1), 'Q']), (evaluated_sub_df_4.at[(33 + 2), 'Q']), (evaluated_sub_df_4.at[(33 + 3), 'Q']), (evaluated_sub_df_4.at[(33 + 4), 'Q']), (evaluated_sub_df_4.at[(33 + 5), 'Q']), (evaluated_sub_df_4.at[(33 + 6), 'Q']), (evaluated_sub_df_4.at[(33 + 7), 'Q']), (evaluated_sub_df_4.at[(33 + 8), 'Q']), (evaluated_sub_df_4.at[(33 + 9), 'Q']), (evaluated_sub_df_4.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_5 = [(evaluated_sub_df_5.at[(33 + 0), 'Q']), (evaluated_sub_df_5.at[(33 + 1), 'Q']), (evaluated_sub_df_5.at[(33 + 2), 'Q']), (evaluated_sub_df_5.at[(33 + 3), 'Q']), (evaluated_sub_df_5.at[(33 + 4), 'Q']), (evaluated_sub_df_5.at[(33 + 5), 'Q']), (evaluated_sub_df_5.at[(33 + 6), 'Q']), (evaluated_sub_df_5.at[(33 + 7), 'Q']), (evaluated_sub_df_5.at[(33 + 8), 'Q']), (evaluated_sub_df_5.at[(33 + 9), 'Q']), (evaluated_sub_df_5.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_6 = [(evaluated_sub_df_6.at[(33 + 0), 'Q']), (evaluated_sub_df_6.at[(33 + 1), 'Q']), (evaluated_sub_df_6.at[(33 + 2), 'Q']), (evaluated_sub_df_6.at[(33 + 3), 'Q']), (evaluated_sub_df_6.at[(33 + 4), 'Q']), (evaluated_sub_df_6.at[(33 + 5), 'Q']), (evaluated_sub_df_6.at[(33 + 6), 'Q']), (evaluated_sub_df_6.at[(33 + 7), 'Q']), (evaluated_sub_df_6.at[(33 + 8), 'Q']), (evaluated_sub_df_6.at[(33 + 9), 'Q']), (evaluated_sub_df_6.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_7 = [(evaluated_sub_df_7.at[(33 + 0), 'Q']), (evaluated_sub_df_7.at[(33 + 1), 'Q']), (evaluated_sub_df_7.at[(33 + 2), 'Q']), (evaluated_sub_df_7.at[(33 + 3), 'Q']), (evaluated_sub_df_7.at[(33 + 4), 'Q']), (evaluated_sub_df_7.at[(33 + 5), 'Q']), (evaluated_sub_df_7.at[(33 + 6), 'Q']), (evaluated_sub_df_7.at[(33 + 7), 'Q']), (evaluated_sub_df_7.at[(33 + 8), 'Q']), (evaluated_sub_df_7.at[(33 + 9), 'Q']), (evaluated_sub_df_7.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_8 = [(evaluated_sub_df_8.at[(33 + 0), 'Q']), (evaluated_sub_df_8.at[(33 + 1), 'Q']), (evaluated_sub_df_8.at[(33 + 2), 'Q']), (evaluated_sub_df_8.at[(33 + 3), 'Q']), (evaluated_sub_df_8.at[(33 + 4), 'Q']), (evaluated_sub_df_8.at[(33 + 5), 'Q']), (evaluated_sub_df_8.at[(33 + 6), 'Q']), (evaluated_sub_df_8.at[(33 + 7), 'Q']), (evaluated_sub_df_8.at[(33 + 8), 'Q']), (evaluated_sub_df_8.at[(33 + 9), 'Q']), (evaluated_sub_df_8.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_9 = [(evaluated_sub_df_9.at[(33 + 0), 'Q']), (evaluated_sub_df_9.at[(33 + 1), 'Q']), (evaluated_sub_df_9.at[(33 + 2), 'Q']), (evaluated_sub_df_9.at[(33 + 3), 'Q']), (evaluated_sub_df_9.at[(33 + 4), 'Q']), (evaluated_sub_df_9.at[(33 + 5), 'Q']), (evaluated_sub_df_9.at[(33 + 6), 'Q']), (evaluated_sub_df_9.at[(33 + 7), 'Q']), (evaluated_sub_df_9.at[(33 + 8), 'Q']), (evaluated_sub_df_9.at[(33 + 9), 'Q']), (evaluated_sub_df_9.at[(33 + 10), 'Q']), ]
                    
                    bars_sub_list_10 = [(evaluated_sub_df_10.at[(33 + 0), 'Q']), (evaluated_sub_df_10.at[(33 + 1), 'Q']), (evaluated_sub_df_10.at[(33 + 2), 'Q']), (evaluated_sub_df_10.at[(33 + 3), 'Q']), (evaluated_sub_df_10.at[(33 + 4), 'Q']), (evaluated_sub_df_10.at[(33 + 5), 'Q']), (evaluated_sub_df_10.at[(33 + 6), 'Q']), (evaluated_sub_df_10.at[(33 + 7), 'Q']), (evaluated_sub_df_10.at[(33 + 8), 'Q']), (evaluated_sub_df_10.at[(33 + 9), 'Q']), (evaluated_sub_df_10.at[(33 + 10), 'Q']), ]
    
                    bar_group_names = ['split', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',]
                



            whole_group_labels = ['dist=1; n=' + str(sub_df_1.shape[0]),'dist=2; n=' + str(sub_df_2.shape[0]),'dist=3; n=' + str(sub_df_3.shape[0]),'dist=4; n=' + str(sub_df_4.shape[0]),'dist=5; n=' + str(sub_df_5.shape[0]),'dist=6; n=' + str(sub_df_6.shape[0]),'dist=7; n=' + str(sub_df_7.shape[0]),'dist=8; n=' + str(sub_df_8.shape[0]),'dist=9; n=' + str(sub_df_9.shape[0]),'dist=10; n=' + str(sub_df_10.shape[0]),]
            
            
            # set width of bar
            barWidth = 0.3
            
            # Set position of bar on X axis
            #r1 = np.arange(len(bars_sub_list_1))
            #r1 = [(x * ((barWidth * number_of_Q_positions) * 1.5)) for x in range(1,(len(bar_group_names) + 1),1)]
            r1 = []
            r = 1
            for p in range(1,(len(bar_group_names) + 1),1):
                r = r + ((barWidth * number_of_evaluated_Q_positions) * 1.2)
                r1.append(r)
            
            r2 = [x + barWidth for x in r1]
            r3 = [x + barWidth for x in r2]
            r4 = [x + barWidth for x in r3]
            r5 = [x + barWidth for x in r4]
            r6 = [x + barWidth for x in r5]
            r7 = [x + barWidth for x in r6]
            r8 = [x + barWidth for x in r7]
            r9 = [x + barWidth for x in r8]
            r10 = [x + barWidth for x in r9]
                     
            # Make the plot
            plt.bar(r1, bars_sub_list_1, color=viridis_data[(25 * 1)], width=barWidth, edgecolor='white', linewidth=(barWidth/10), label=whole_group_labels[0])
            plt.bar(r2, bars_sub_list_2, color=viridis_data[(25 * 2)], width=barWidth, edgecolor='white', linewidth=(barWidth/10), label=whole_group_labels[1])
            plt.bar(r3, bars_sub_list_3, color=viridis_data[(25 * 3)], width=barWidth, edgecolor='white', linewidth=(barWidth/10), label=whole_group_labels[2])
            plt.bar(r4, bars_sub_list_4, color=viridis_data[(25 * 4)], width=barWidth, edgecolor='white', linewidth=(barWidth/10), label=whole_group_labels[3])
            plt.bar(r5, bars_sub_list_5, color=viridis_data[(25 * 5)], width=barWidth, edgecolor='white', linewidth=(barWidth/10), label=whole_group_labels[4])
            plt.bar(r6, bars_sub_list_6, color=viridis_data[(25 * 6)], width=barWidth, edgecolor='white', linewidth=(barWidth/10), label=whole_group_labels[5])
            plt.bar(r7, bars_sub_list_7, color=viridis_data[(25 * 7)], width=barWidth, edgecolor='white', linewidth=(barWidth/10), label=whole_group_labels[6])
            plt.bar(r8, bars_sub_list_8, color=viridis_data[(25 * 8)], width=barWidth, edgecolor='white', linewidth=(barWidth/10), label=whole_group_labels[7])
            plt.bar(r9, bars_sub_list_9, color=viridis_data[(25 * 9)], width=barWidth, edgecolor='white', linewidth=(barWidth/10), label=whole_group_labels[8])
            plt.bar(r10, bars_sub_list_10, color=viridis_data[(25 * 10)], width=barWidth, edgecolor='white', linewidth=(barWidth/10), label=whole_group_labels[9])
                
            # Add xticks on the middle of the group bars
            plt.xlabel((name_1 + name_2), fontweight='bold', verticalalignment='top', horizontalalignment='right')
            plt.xticks([x + ((number_of_evaluated_Q_positions / 2) * barWidth) for x in r1], bar_group_names,)
            #plt.tick_params(labelsize='large')
            # Create legend & Show graphic
            lgd = plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=1)#frameon=False)
            plt.rc('xtick', labelsize=6)     

            #plt.show()
            plt.savefig(working_directory + name_1 + name_2 + '_bars.svg', bbox_extra_artists=(lgd,), bbox_inches='tight')

            plt.close()
            
            #end barplot / start lineplot
            
            #Intron/Exon boundary line plots
            x_start = 0
            x_stop = 0
            
            if num_acids == 10:
                if 'prae' in name_1:
                    bar_group_names = range(-10,1,1)
                    x_start = -10
                    x_stop = 0
                elif 'post' in name_1:
                    bar_group_names = range(0,11,1)
                    x_start = 0
                    x_stop = 10
            elif num_acids == 33:
                if 'prae' in name_1:
                    bar_group_names = range(-33,1,1)
                    x_start = -33
                    x_stop = 0
                elif 'post' in name_1:
                    bar_group_names = range(0,34,1)
                    x_start = 0
                    x_stop = 33
                    
            # Make a data frame
            table = [bars_sub_list_1, bars_sub_list_2, bars_sub_list_3, bars_sub_list_4, bars_sub_list_5, bars_sub_list_6, bars_sub_list_7, bars_sub_list_8, bars_sub_list_9, bars_sub_list_10, bar_group_names]
            whole_group_labels = ['dist=1; n=' + str(sub_df_1.shape[0]),'dist=2; n=' + str(sub_df_2.shape[0]),'dist=3; n=' + str(sub_df_3.shape[0]),'dist=4; n=' + str(sub_df_4.shape[0]),'dist=5; n=' + str(sub_df_5.shape[0]),'dist=6; n=' + str(sub_df_6.shape[0]),'dist=7; n=' + str(sub_df_7.shape[0]),'dist=8; n=' + str(sub_df_8.shape[0]),'dist=9; n=' + str(sub_df_9.shape[0]),'dist=10; n=' + str(sub_df_10.shape[0]), 'position']
            #whole_group_labels = ['polyQ_dist=1','polyQ_dist=2','polyQ_dist=3','polyQ_dist=4','polyQ_dist=5','polyQ_dist=6','polyQ_dist=7','polyQ_dist=8','polyQ_dist=9','polyQ_dist=10', 'position']
            df = pd.DataFrame(table)
            df = df.transpose()
            df.columns = whole_group_labels
            
            ready_to_plot = df#.pivot_table(index='character',columns='position',values='value')
            print('ready_to_plot')
            print(ready_to_plot)
            max_value = sorted((ready_to_plot[whole_group_labels[0:number_of_evaluated_Q_positions]].max()).tolist(), reverse=True)[0]
            #min_value = sorted((ready_to_plot[whole_group_labels[0:number_of_evaluated_Q_positions]].min()).tolist(), reverse=False)[0]
            
            header_list = df.columns.tolist()
            
            # Initialize the figure
            myfig = plt.figure()
            plt.style.use('ggplot')
            # multiple line plot
            num=0
            dropped_positions = ready_to_plot.drop('position', axis=1)
            for column in dropped_positions:
                num+=1
             
                # Find the right spot on the plot
                plt.subplot(5,2, num)
             
                # plot every groups, but discreet
                for v in dropped_positions:
                    plt.plot(ready_to_plot['position'], ready_to_plot[v], marker='', color='grey', linewidth=0.6, alpha=0.3)
             
                # Plot the lineplot
                plt.plot(ready_to_plot['position'], ready_to_plot[column], marker='', color=viridis_data[(num * 20)], linewidth=1.0, alpha=0.9, label=column)
             
                # Same limits for everybody!
                header_list
                plt.xlim(x_start,x_stop)
                if 'relative' in name_2:
                    plt.ylim(0, (1.05))
                elif 'total' in name_2:
                    plt.ylim(-2, (max_value + 1))
                # Not ticks everywhere
                if num in range(9) :
                    plt.tick_params(labelbottom=False)
                if num not in [1,3,5,7,9] :
                    plt.tick_params(labelleft=False)               
                # Add title
                plt.title(column, loc='left', pad=5, fontsize=10, fontweight=0, color='black')
             
            # general title
            plt.suptitle(name, y=5, fontsize=13, fontweight=0, color='black', style='italic',)
            plt.tight_layout()
            #save figure as svg
            myfig.savefig(working_directory + name_1 + name_2 + '_line.svg')
            plt.close()
    return True
            

def plot_intron_spanning_polyQs_both_sides_I(working_directory):#part of Function 5.3
    # Intron/Exon boundary plots : characteristics of the exon with a polyQ (+after)
    viridis_data = [[0.267004, 0.004874, 0.329415],[0.268510, 0.009605, 0.335427],[0.269944, 0.014625, 0.341379],[0.271305, 0.019942, 0.347269],[0.272594, 0.025563, 0.353093],[0.273809, 0.031497, 0.358853],[0.274952, 0.037752, 0.364543],[0.276022, 0.044167, 0.370164],[0.277018, 0.050344, 0.375715],[0.277941, 0.056324, 0.381191],[0.278791, 0.062145, 0.386592],[0.279566, 0.067836, 0.391917],[0.280267, 0.073417, 0.397163],[0.280894, 0.078907, 0.402329],[0.281446, 0.084320, 0.407414],[0.281924, 0.089666, 0.412415],[0.282327, 0.094955, 0.417331],[0.282656, 0.100196, 0.422160],[0.282910, 0.105393, 0.426902],[0.283091, 0.110553, 0.431554],[0.283197, 0.115680, 0.436115],[0.283229, 0.120777, 0.440584],[0.283187, 0.125848, 0.444960],[0.283072, 0.130895, 0.449241],[0.282884, 0.135920, 0.453427],[0.282623, 0.140926, 0.457517],[0.282290, 0.145912, 0.461510],[0.281887, 0.150881, 0.465405],[0.281412, 0.155834, 0.469201],[0.280868, 0.160771, 0.472899],[0.280255, 0.165693, 0.476498],[0.279574, 0.170599, 0.479997],[0.278826, 0.175490, 0.483397],[0.278012, 0.180367, 0.486697],[0.277134, 0.185228, 0.489898],[0.276194, 0.190074, 0.493001],[0.275191, 0.194905, 0.496005],[0.274128, 0.199721, 0.498911],[0.273006, 0.204520, 0.501721],[0.271828, 0.209303, 0.504434],[0.270595, 0.214069, 0.507052],[0.269308, 0.218818, 0.509577],[0.267968, 0.223549, 0.512008],[0.266580, 0.228262, 0.514349],[0.265145, 0.232956, 0.516599],[0.263663, 0.237631, 0.518762],[0.262138, 0.242286, 0.520837],[0.260571, 0.246922, 0.522828],[0.258965, 0.251537, 0.524736],[0.257322, 0.256130, 0.526563],[0.255645, 0.260703, 0.528312],[0.253935, 0.265254, 0.529983],[0.252194, 0.269783, 0.531579],[0.250425, 0.274290, 0.533103],[0.248629, 0.278775, 0.534556],[0.246811, 0.283237, 0.535941],[0.244972, 0.287675, 0.537260],[0.243113, 0.292092, 0.538516],[0.241237, 0.296485, 0.539709],[0.239346, 0.300855, 0.540844],[0.237441, 0.305202, 0.541921],[0.235526, 0.309527, 0.542944],[0.233603, 0.313828, 0.543914],[0.231674, 0.318106, 0.544834],[0.229739, 0.322361, 0.545706],[0.227802, 0.326594, 0.546532],[0.225863, 0.330805, 0.547314],[0.223925, 0.334994, 0.548053],[0.221989, 0.339161, 0.548752],[0.220057, 0.343307, 0.549413],[0.218130, 0.347432, 0.550038],[0.216210, 0.351535, 0.550627],[0.214298, 0.355619, 0.551184],[0.212395, 0.359683, 0.551710],[0.210503, 0.363727, 0.552206],[0.208623, 0.367752, 0.552675],[0.206756, 0.371758, 0.553117],[0.204903, 0.375746, 0.553533],[0.203063, 0.379716, 0.553925],[0.201239, 0.383670, 0.554294],[0.199430, 0.387607, 0.554642],[0.197636, 0.391528, 0.554969],[0.195860, 0.395433, 0.555276],[0.194100, 0.399323, 0.555565],[0.192357, 0.403199, 0.555836],[0.190631, 0.407061, 0.556089],[0.188923, 0.410910, 0.556326],[0.187231, 0.414746, 0.556547],[0.185556, 0.418570, 0.556753],[0.183898, 0.422383, 0.556944],[0.182256, 0.426184, 0.557120],[0.180629, 0.429975, 0.557282],[0.179019, 0.433756, 0.557430],[0.177423, 0.437527, 0.557565],[0.175841, 0.441290, 0.557685],[0.174274, 0.445044, 0.557792],[0.172719, 0.448791, 0.557885],[0.171176, 0.452530, 0.557965],[0.169646, 0.456262, 0.558030],[0.168126, 0.459988, 0.558082],[0.166617, 0.463708, 0.558119],[0.165117, 0.467423, 0.558141],[0.163625, 0.471133, 0.558148],[0.162142, 0.474838, 0.558140],[0.160665, 0.478540, 0.558115],[0.159194, 0.482237, 0.558073],[0.157729, 0.485932, 0.558013],[0.156270, 0.489624, 0.557936],[0.154815, 0.493313, 0.557840],[0.153364, 0.497000, 0.557724],[0.151918, 0.500685, 0.557587],[0.150476, 0.504369, 0.557430],[0.149039, 0.508051, 0.557250],[0.147607, 0.511733, 0.557049],[0.146180, 0.515413, 0.556823],[0.144759, 0.519093, 0.556572],[0.143343, 0.522773, 0.556295],[0.141935, 0.526453, 0.555991],[0.140536, 0.530132, 0.555659],[0.139147, 0.533812, 0.555298],[0.137770, 0.537492, 0.554906],[0.136408, 0.541173, 0.554483],[0.135066, 0.544853, 0.554029],[0.133743, 0.548535, 0.553541],[0.132444, 0.552216, 0.553018],[0.131172, 0.555899, 0.552459],[0.129933, 0.559582, 0.551864],[0.128729, 0.563265, 0.551229],[0.127568, 0.566949, 0.550556],[0.126453, 0.570633, 0.549841],[0.125394, 0.574318, 0.549086],[0.124395, 0.578002, 0.548287],[0.123463, 0.581687, 0.547445],[0.122606, 0.585371, 0.546557],[0.121831, 0.589055, 0.545623],[0.121148, 0.592739, 0.544641],[0.120565, 0.596422, 0.543611],[0.120092, 0.600104, 0.542530],[0.119738, 0.603785, 0.541400],[0.119512, 0.607464, 0.540218],[0.119423, 0.611141, 0.538982],[0.119483, 0.614817, 0.537692],[0.119699, 0.618490, 0.536347],[0.120081, 0.622161, 0.534946],[0.120638, 0.625828, 0.533488],[0.121380, 0.629492, 0.531973],[0.122312, 0.633153, 0.530398],[0.123444, 0.636809, 0.528763],[0.124780, 0.640461, 0.527068],[0.126326, 0.644107, 0.525311],[0.128087, 0.647749, 0.523491],[0.130067, 0.651384, 0.521608],[0.132268, 0.655014, 0.519661],[0.134692, 0.658636, 0.517649],[0.137339, 0.662252, 0.515571],[0.140210, 0.665859, 0.513427],[0.143303, 0.669459, 0.511215],[0.146616, 0.673050, 0.508936],[0.150148, 0.676631, 0.506589],[0.153894, 0.680203, 0.504172],[0.157851, 0.683765, 0.501686],[0.162016, 0.687316, 0.499129],[0.166383, 0.690856, 0.496502],[0.170948, 0.694384, 0.493803],[0.175707, 0.697900, 0.491033],[0.180653, 0.701402, 0.488189],[0.185783, 0.704891, 0.485273],[0.191090, 0.708366, 0.482284],[0.196571, 0.711827, 0.479221],[0.202219, 0.715272, 0.476084],[0.208030, 0.718701, 0.472873],[0.214000, 0.722114, 0.469588],[0.220124, 0.725509, 0.466226],[0.226397, 0.728888, 0.462789],[0.232815, 0.732247, 0.459277],[0.239374, 0.735588, 0.455688],[0.246070, 0.738910, 0.452024],[0.252899, 0.742211, 0.448284],[0.259857, 0.745492, 0.444467],[0.266941, 0.748751, 0.440573],[0.274149, 0.751988, 0.436601],[0.281477, 0.755203, 0.432552],[0.288921, 0.758394, 0.428426],[0.296479, 0.761561, 0.424223],[0.304148, 0.764704, 0.419943],[0.311925, 0.767822, 0.415586],[0.319809, 0.770914, 0.411152],[0.327796, 0.773980, 0.406640],[0.335885, 0.777018, 0.402049],[0.344074, 0.780029, 0.397381],[0.352360, 0.783011, 0.392636],[0.360741, 0.785964, 0.387814],[0.369214, 0.788888, 0.382914],[0.377779, 0.791781, 0.377939],[0.386433, 0.794644, 0.372886],[0.395174, 0.797475, 0.367757],[0.404001, 0.800275, 0.362552],[0.412913, 0.803041, 0.357269],[0.421908, 0.805774, 0.351910],[0.430983, 0.808473, 0.346476],[0.440137, 0.811138, 0.340967],[0.449368, 0.813768, 0.335384],[0.458674, 0.816363, 0.329727],[0.468053, 0.818921, 0.323998],[0.477504, 0.821444, 0.318195],[0.487026, 0.823929, 0.312321],[0.496615, 0.826376, 0.306377],[0.506271, 0.828786, 0.300362],[0.515992, 0.831158, 0.294279],[0.525776, 0.833491, 0.288127],[0.535621, 0.835785, 0.281908],[0.545524, 0.838039, 0.275626],[0.555484, 0.840254, 0.269281],[0.565498, 0.842430, 0.262877],[0.575563, 0.844566, 0.256415],[0.585678, 0.846661, 0.249897],[0.595839, 0.848717, 0.243329],[0.606045, 0.850733, 0.236712],[0.616293, 0.852709, 0.230052],[0.626579, 0.854645, 0.223353],[0.636902, 0.856542, 0.216620],[0.647257, 0.858400, 0.209861],[0.657642, 0.860219, 0.203082],[0.668054, 0.861999, 0.196293],[0.678489, 0.863742, 0.189503],[0.688944, 0.865448, 0.182725],[0.699415, 0.867117, 0.175971],[0.709898, 0.868751, 0.169257],[0.720391, 0.870350, 0.162603],[0.730889, 0.871916, 0.156029],[0.741388, 0.873449, 0.149561],[0.751884, 0.874951, 0.143228],[0.762373, 0.876424, 0.137064],[0.772852, 0.877868, 0.131109],[0.783315, 0.879285, 0.125405],[0.793760, 0.880678, 0.120005],[0.804182, 0.882046, 0.114965],[0.814576, 0.883393, 0.110347],[0.824940, 0.884720, 0.106217],[0.835270, 0.886029, 0.102646],[0.845561, 0.887322, 0.099702],[0.855810, 0.888601, 0.097452],[0.866013, 0.889868, 0.095953],[0.876168, 0.891125, 0.095250],[0.886271, 0.892374, 0.095374],[0.896320, 0.893616, 0.096335],[0.906311, 0.894855, 0.098125],[0.916242, 0.896091, 0.100717],[0.926106, 0.897330, 0.104071],[0.935904, 0.898570, 0.108131],[0.945636, 0.899815, 0.112838],[0.955300, 0.901065, 0.118128],[0.964894, 0.902323, 0.123941],[0.974417, 0.903590, 0.130215],[0.983868, 0.904867, 0.136897],[0.993248, 0.906157, 0.143936]]

    print('plot_intron_spanning_polyQs_both_sides')
    
    full_df = pd.read_pickle(working_directory + 'trans_intron_polyQ_full_df_pickle.pkl')
    prae_df = pd.read_pickle(working_directory + 'trans_intron_polyQ_prae_df_pickle.pkl')
    post_df = pd.read_pickle(working_directory + 'trans_intron_polyQ_post_df_pickle.pkl')
    list_of_AAs_single_letters = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z','_',"' '",]
    list_of_dfs_and_names = []
    ###########################################################################
    full_df_count_row = full_df.shape[0]  # gives number of row count

    full_df = full_df.drop_duplicates(subset=['ENSG', 'CDS1_start', 'CDS1_stop', 'CDS2_start', 'CDS2_stop',], keep='first')#keep options   ->   first : Drop duplicates except for the first occurrence; last : Drop duplicates except for the last occurrence; False : Drop all duplicates.

    full_df = full_df.drop('position', 1)# number == axis ... (0 == row; 1 == column)
    full_df = full_df.drop('strand', 1)
    full_df = full_df.drop('string_triplet', 1)
    full_df = full_df.drop('CDS1_start', 1)
    full_df = full_df.drop('CDS1_stop', 1)
    full_df = full_df.drop('CDS2_start', 1)
    full_df = full_df.drop('CDS2_stop', 1)
    full_df = full_df.drop('ENSG', 1)
        
    full_total_dict = defaultdict(list)
    full_relative_dict = defaultdict(list)

    header_list = full_df.columns.tolist()
    negative_start_value = ((len(header_list) - 1) / (- 2))
    for header in header_list:
        full_total_dict['position'].append(int(negative_start_value + header))
        full_relative_dict['position'].append(int(negative_start_value + header))
        column_list = full_df[header].tolist()
        for character in list_of_AAs_single_letters:
            char_counts = column_list.count(character)
            full_total_dict[character].append(char_counts)
            full_relative_dict[character].append(char_counts / full_df_count_row)
               
    full_total_df = pd.DataFrame.from_dict(full_total_dict)
    full_relative_df = pd.DataFrame.from_dict(full_relative_dict)
    list_of_dfs_and_names.append((full_total_df, 'full_total_df'))
    list_of_dfs_and_names.append((full_relative_df, 'full_relative_df'))
    
    ###########################################################################
    prae_df_count_row = prae_df.shape[0]  # gives number of row count
    
    #identify duplicate rows that share the complete subset
    prae_df = prae_df.drop_duplicates(subset=['ENSG', 'CDS1_start', 'CDS1_stop', 'CDS2_start', 'CDS2_stop',], keep='first')#keep options   ->   first : Drop duplicates except for the first occurrence; last : Drop duplicates except for the last occurrence; False : Drop all duplicates.

    prae_df = prae_df.drop('position', 1)# number == axis ... (0 == row; 1 == column)
    prae_df = prae_df.drop('strand', 1)
    prae_df = prae_df.drop('string_triplet', 1)
    prae_df = prae_df.drop('CDS1_start', 1)
    prae_df = prae_df.drop('CDS1_stop', 1)
    prae_df = prae_df.drop('CDS2_start', 1)
    prae_df = prae_df.drop('CDS2_stop', 1)
    prae_df = prae_df.drop('ENSG', 1)
        
    prae_total_dict = defaultdict(list)
    prae_relative_dict = defaultdict(list)
    
    header_list = prae_df.columns.tolist()
    negative_start_value = ((len(header_list) - 1) / (- 2))
    for header in header_list:
        prae_total_dict['position'].append(int(negative_start_value + header))
        prae_relative_dict['position'].append(int(negative_start_value + header))
        column_list = prae_df[header].tolist()
        for character in list_of_AAs_single_letters:
            char_counts = column_list.count(character)
            prae_total_dict[character].append(char_counts)
            prae_relative_dict[character].append(char_counts / prae_df_count_row)
               
    prae_total_df = pd.DataFrame.from_dict(prae_total_dict)
    prae_relative_df = pd.DataFrame.from_dict(prae_relative_dict)
    list_of_dfs_and_names.append((prae_total_df, 'prae_total_df'))
    list_of_dfs_and_names.append((prae_relative_df, 'prae_relative_df'))
    
    ###########################################################################
    post_df_count_row = post_df.shape[0]  # gives number of row count

    post_df = post_df.drop_duplicates(subset=['ENSG', 'CDS1_start', 'CDS1_stop', 'CDS2_start', 'CDS2_stop',], keep='first')#keep options   ->   first : Drop duplicates except for the first occurrence; last : Drop duplicates except for the last occurrence; False : Drop all duplicates.

    post_df = post_df.drop('position', 1)# number == axis ... (0 == row; 1 == column)
    post_df = post_df.drop('strand', 1)
    post_df = post_df.drop('string_triplet', 1)
    post_df = post_df.drop('CDS1_start', 1)
    post_df = post_df.drop('CDS1_stop', 1)
    post_df = post_df.drop('CDS2_start', 1)
    post_df = post_df.drop('CDS2_stop', 1)
    post_df = post_df.drop('ENSG', 1)
        
    post_total_dict = defaultdict(list)
    post_relative_dict = defaultdict(list)
    
    header_list = post_df.columns.tolist()
    negative_start_value = ((len(header_list) - 1) / (- 2))
    for header in header_list:
        post_total_dict['position'].append(int(negative_start_value + header))
        post_relative_dict['position'].append(int(negative_start_value + header))
        column_list = post_df[header].tolist()
        for character in list_of_AAs_single_letters:
            char_counts = column_list.count(character)
            post_total_dict[character].append(char_counts)
            post_relative_dict[character].append(char_counts / post_df_count_row)
               
    post_total_df = pd.DataFrame.from_dict(post_total_dict)
    post_relative_df = pd.DataFrame.from_dict(post_relative_dict)
    list_of_dfs_and_names.append((post_total_df, 'post_total_df'))
    list_of_dfs_and_names.append((post_relative_df, 'post_relative_df'))
        
    for df, name in list_of_dfs_and_names:
        # Make a data frame
        ready_to_plot = df#.pivot_table(index='character',columns='position',values='value')

        max_value = sorted((ready_to_plot[['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z','_',"' '",]].max()).tolist(), reverse=True)[0]
        #min_value = sorted((ready_to_plot[['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z','_',"' '",]].min()).tolist(), reverse=False)[0]
        
        header_list = df.columns.tolist()
        negative_start_value = -33
        positive_stop_value = 33
        
        # Initialize the figure
        myfig = plt.figure()
        #plt.style.use('seaborn-darkgrid')
        plt.style.use('ggplot')
             
        # multiple line plot
        num=0
        dropped_positions = ready_to_plot.drop('position', axis=1)
        for column in dropped_positions:
            num+=1
         
            # Find the right spot on the plot
            plt.subplot(5,5, num)
            # plot every groups, but discreet
            for v in dropped_positions:
                plt.plot(ready_to_plot['position'], ready_to_plot[v], marker='', color='grey', linewidth=0.6, alpha=0.3)
         
            # Plot the lineplot
            plt.plot(ready_to_plot['position'], ready_to_plot[column], marker='', color=viridis_data[(num * 10)], linewidth=1.0, alpha=0.9, label=column)

            # Same limits for everybody!
            header_list
            #plt.xlim(negative_start_value,positive_stop_value)
            plt.xlim(negative_start_value,positive_stop_value)
            if 'relative' in name:
                plt.ylim(0, (max_value))
            else:
                plt.ylim(-2, (max_value))

            #change grid
            x_ticks = [(-30),(-20),(-10),0,10,20,30]
            x_labels = ["-30","","","","","","30"]
            #plt.xticks(ticks=x_ticks, labels=x_labels)
            # Not ticks everywhere
            if num in range(21) :
                plt.tick_params(labelbottom=False)
            if num not in [1,6,11,16,21] :
                plt.tick_params(labelleft=False)

            
            # Add title
            plt.title(column, loc='left', pad=-10, fontsize=12, fontweight=0, color='black')
         
        # general title
        plt.suptitle(name, y=5, fontsize=13, fontweight=0, color='black', style='italic',)
         
        # Axis title
        plt.text(-260, -30, 'distance from splice-site', ha='center', va='center', fontsize=16,)
        plt.text(-380, 150, 'number of amino acids per position', ha='center', va='center', rotation='vertical', fontsize=16,)
        
        
        #plt.tight_layout()

        #save figure as svg
        myfig.savefig(working_directory + name + '_overview_line.svg')

        plt.close()
        
    return True

#in my thesis referred to as Function 6
def HRA_compare_all_QRYs_and_REF(working_directory, bed_name, REF, gff3, experiment_folders, ChrOI_list, positions, overwrite_old_files, aligner):

    print('Warning: up until now, Ancient_Repeats_Tools.py demanded for a explicitely prepared reference (in form of bed-, sequence-, and gff3-files) for each mapped whole genome (bam file) analyzed. ' + 
          'The upside is, that experiments with different reference genomes can co-exist in the same directory. ' + 
          'This can be justified in analysis of the reference genome not being a computationally heavy task. ' + 
          'The step initiated now will compare all experiment folders (submitted under the "-bam" operator) among each other. ' + 
          'Therefor it is expected, that all their reference genomes, submitted bed- and gff3-files and past analysis-parameters are the same. ' + 
          'As a precaution, it will be checked, whether the reference results stored in the csv_results folder are identical.' +
          'Accordinlgy, out of all experiment folders submitted to this analysis, one (the first) will be chosen for representative extraction of past results ' + 
          'Results extracted will be taken from the HRA_bed_parser_dict and csv_results folder.\nIn addition, if "muscle" was chosen as "aligner" please be aware that difficulties can emerge in using the required biopython wrapper with pypy; if "mafft" was chosen as "aligner" please make sure that mafft is installed and that a subprocess.call() function can spawn a new process safely given the user input')    

    print('before: ' + REF)    
    REF = re.sub(r'\.\w*\.fa$', '.', REF)
    print('after: ' + REF)
    
    print('before: ' + gff3)
    gff3 = re.sub(r'\.\w*\.gff3$', '.', gff3)    
    print('after: ' + gff3)    

    if (os.path.exists(working_directory + 'heatmaps/') == False):
        os.mkdir(working_directory + 'heatmaps/')
    else:
        if overwrite_old_files == 'True':
            shutil.rmtree(working_directory + 'heatmaps/')
            os.mkdir(working_directory + 'heatmaps/')
        else:
            raise Exception('HRA_compare_all_QRYs_and_REF(ERROR): directory ' + working_directory + 'heatmaps/ already exists and "overwrite" is not set to "True" ... ')
    if experiment_folders[0] != 'REF':
        random_experiment = experiment_folders[0]
    elif experiment_folders[0] == 'REF':
        random_experiment = experiment_folders[1]
    
    data_dict = defaultdict(list)
    key_list = set()
    ENSGs = set()
    for experiment in experiment_folders:
        complete_results_path = working_directory + experiment + '/csv_results/'
        for ChrOI in ChrOI_list:
            csv_files = [(complete_results_path + any_file) for any_file in listdir(complete_results_path) if (isfile(join(complete_results_path, any_file)) and re.match(r'csv_analysis_for_chromosome_P_' + ChrOI + r'_.*_polyQs_query.csv', any_file))]
            for csv_file in csv_files:
                complete_results = open(csv_file, 'r')
                if positions == 'absolute':
                    for line in complete_results:
                        
                        split_line = line.rstrip('\n').split(',')
                        ENSG =  split_line[1]
                        data_dict[ENSG + '|' + experiment].append(split_line)
                        key_list.add(ENSG + '|' + experiment)
                        ENSGs.add(ENSG)
                        
                complete_results.close()

    #check that all references (one should have been created per experiment) are the same
    for ChrOI in ChrOI_list:

        if (ChrOI != 'Y') and (ChrOI != 'MT'):
            for ones_pre_dec in ['1','2','3','4','5','6','7','8','9','0']:
                ref_file_name = REF + ChrOI + r'.fa_P_' + ChrOI + '_'+ ones_pre_dec + '_polyQs_reference.csv'
                all_references_read = [working_directory + experiment + '/csv_results/' + ref_file_name for experiment in experiment_folders]
                #print('all_references_read double touple list: ' + str(all_references_read))
                for list_position in range(1,len(experiment_folders),1):
                    if filecmp.cmp(all_references_read[0], all_references_read[list_position]) != True:
                        raise Exception('HRA_compare_all_QRYs_and_REF(ERROR): the reference csv files seem to differ between each experiment which should not be the case: ' + all_references_read[0][1] + ' vs ' + all_references_read[list_position][1])    

        elif (ChrOI == 'Y') or (ChrOI == 'MT'):
            for ones_pre_dec in ['(1|2)','(3|4)','(5|6)','(7|8)','(9|0)']:
                ref_file_name = REF + ChrOI + r'.fa_P_' + ChrOI + '_'+ ones_pre_dec + '_polyQs_reference.csv'
                all_references_read = [working_directory + experiment + '/csv_results/' + ref_file_name for experiment in experiment_folders]
                #print('all_references_read double touple list: ' + str(all_references_read))
                for list_position in range(1,len(experiment_folders),1):
                    if filecmp.cmp(all_references_read[0], all_references_read[list_position]) != True:
                        raise Exception('HRA_compare_all_QRYs_and_REF(ERROR): the reference csv files seem to differ between each experiment which should not be the case: ' + all_references_read[0][1] + ' vs ' + all_references_read[list_position][1])    
                    
    #continue with the first experiment as representative experiment for reference analysis
    complete_results_path = working_directory + random_experiment + '/csv_results/'
    for ChrOI in ChrOI_list:
        ref_csv_files = [(complete_results_path + any_file) for any_file in listdir(complete_results_path) if (isfile(join(complete_results_path, any_file)) and re.match(REF + ChrOI + r'.fa_P_' + ChrOI + r'_.*_polyQs_reference.csv', any_file))]
            
        for ref_csv_file in ref_csv_files:
            ref_complete_results = open(ref_csv_file, 'r')
            if positions == 'absolute':
                for line in ref_complete_results:
                        
                    split_line = line.rstrip('\n').split(',')
                    ENSG =  split_line[1]
                    data_dict[ENSG + '|REF'].append(split_line)
                    key_list.add(ENSG + '|' + experiment)
                    ENSGs.add(ENSG)
    
    print('loading HRA_bed_parser_dict')
    ref_bed_pickle_path = working_directory + random_experiment + '/bed_pickles/'
    print("random experiment's ref_bed_pickle_path: " + ref_bed_pickle_path)
    list_of_HRA_bed_parser_dict_pickle_paths = [(ref_bed_pickle_path + any_file) for any_file in listdir(ref_bed_pickle_path) if (isfile(join(ref_bed_pickle_path, any_file)) and re.match(r'.*_P.*_ENST_split.txt', any_file))]
    list_of_HRA_bed_parser_dicts_loaded_from_pickles = [pickle.load(open(pickled_dict_path)) for pickled_dict_path in list_of_HRA_bed_parser_dict_pickle_paths]        
    if len(list_of_HRA_bed_parser_dicts_loaded_from_pickles) < 1:
        print('No HRA_bed_parser_dict could be found under the following paths:')
        for path in list_of_HRA_bed_parser_dict_pickle_paths:
            print(path)
        raise Exception('HRA_compare_all_QRYs_and_REF(ERROR): HRA_bed_parser_dicts could not be loaded and merged')
    
    experiment_folders = ['REF'] + experiment_folders
    data = defaultdict(list)
    for detected_ENSG in ENSGs:
        for experiment in experiment_folders:
            data_lol = data_dict[detected_ENSG + '|' + experiment]
            for data_list in data_lol:
                if data_list != []:
                    experiment_label = experiment
                    Chromosome = data_list[0]
                    ENSG = data_list[1]
                    ENST = data_list[2]
                    strand = data_list[3]
                    ENST_start = int(data_list[9])#int(merged_HRA_bed_parser_dict[ENST + '|bed_data'][1])
                    ENST_stop = int(data_list[10])#int(merged_HRA_bed_parser_dict[ENST + '|bed_data'][2])
                    repeat_length = int(data_list[4].replace('.0', ''))
                    repeat_start = int(data_list[5])
                    repeat_stop = int(data_list[6])
                    repeat_seq = data_list[7]
                    ENST_seq_plus_introns = data_list[12]
                    by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed = 0
                    by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed = 0
                    repeat_start_in_ENST_seq_plus_introns = (repeat_start - ENST_start)
                    repeat_stop_in_ENST_seq_plus_introns = (repeat_stop - ENST_start) + 1
                    start_minus_50 = repeat_start_in_ENST_seq_plus_introns - 50
                    if start_minus_50 < 0:
                        start_minus_50 = 0
                        by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed = 50 -repeat_start_in_ENST_seq_plus_introns
                    stop_plus_50 = repeat_stop_in_ENST_seq_plus_introns + 50
                    if stop_plus_50 > len(ENST_seq_plus_introns):
                        stop_plus_50 = len(ENST_seq_plus_introns)
                        by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed = stop_plus_50 - len(ENST_seq_plus_introns)
                    repeat_seq_plus_50_bases_left_and_right = ENST_seq_plus_introns[start_minus_50:stop_plus_50]
                    
                    if repeat_seq_plus_50_bases_left_and_right.upper()[0:-1] != data_list[11].upper():
                        print('conflict between these two repeat_seq_plus_50_bases_left_and_right: \n' + repeat_seq_plus_50_bases_left_and_right + '\n' + data_list[11])
                    
                    data['experiment_label'].append(experiment_label)
                    data['Chromosome'].append(Chromosome)
                    data['ENSG'].append(ENSG)
                    data['ENST'].append(ENST)
                    data['strand'].append(strand)
                    #reminder: 
                    #HRA_bed_parser_dict[ENST + '|bed_data'] = (str(Chromosome), str(Start), str(Stop), str(Strand))
                    data['ENST_start'].append(ENST_start)
                    data['ENST_stop'].append(ENST_stop)
                    data['repeat_length'].append(repeat_length)
                    data['repeat_start'].append(repeat_start)
                    data['repeat_stop'].append(repeat_stop)
                    data['repeat_seq'].append(repeat_seq)
                    data['ENST_seq_plus_introns'].append(ENST_seq_plus_introns)
                    data['repeat_seq_plus_50_bases_left_and_right'].append(repeat_seq_plus_50_bases_left_and_right)
                    data['lflank_seq'].append(ENST_seq_plus_introns[start_minus_50:repeat_start_in_ENST_seq_plus_introns])
                    data['start_minus_50'].append(start_minus_50)
                    data['repeat_start_in_ENST_seq_plus_introns'].append(repeat_start_in_ENST_seq_plus_introns) 
                    data['middle_repeat_seq'].append(ENST_seq_plus_introns[repeat_start_in_ENST_seq_plus_introns:repeat_stop_in_ENST_seq_plus_introns])
                    data['rflank_seq'].append(ENST_seq_plus_introns[repeat_stop_in_ENST_seq_plus_introns:stop_plus_50])
                    data['repeat_stop_in_ENST_seq_plus_introns'].append(repeat_stop_in_ENST_seq_plus_introns)
                    data['stop_plus_50'].append(stop_plus_50)
                    data['by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed'].append(by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed)
                    data['by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed'].append(by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed)
                    
                elif data_dict[detected_ENSG + '|' + experiment] == []:
                    data['experiment_label'].append(experiment_label)
                    data['Chromosome'].append('NA')
                    data['ENSG'].append(ENSG)
                    data['ENST'].append('NA')
                    data['strand'].append('NA')
                    #reminder: 
                    #HRA_bed_parser_dict[ENST + '|bed_data'] = (str(Chromosome), str(Start), str(Stop), str(Strand))
                    data['ENST_start'].append('NA')
                    data['ENST_stop'].append('NA')
                    data['repeat_length'].append('0')
                    data['repeat_start'].append('NA')
                    data['repeat_stop'].append('NA')
                    data['repeat_seq'].append('NA')
                    data['ENST_seq_plus_introns'].append('NA')
                    data['repeat_seq_plus_50_bases_left_and_right'].append('NA')               
                    data['lflank_seq'].append('NA')
                    data['start_minus_50'].append('NA')
                    data['repeat_start_in_ENST_seq_plus_introns'].append('NA') 
                    data['middle_repeat_seq'].append('NA')
                    data['rflank_seq'].append('NA')
                    data['repeat_stop_in_ENST_seq_plus_introns'].append('NA')
                    data['stop_plus_50'].append('NA')
                    data['by_how_many_bases_was_the_starting_addition_of_50_bases_curtailed'].append('NA')
                    data['by_how_many_bases_was_the_stop_addition_of_50_bases_curtailed'].append('NA')
                    
                else:
                    message = ('HRA_compare_all_QRYs_and_REF(ERROR): data_dict[' + detected_ENSG + '|' + experiment + '] is neither == [] nor != [] ...')
                    raise Exception(message)
    
    print('creating data_frame_ready_dict_pickle.txt')
    new_pickle = open(working_directory + 'data_frame_ready_dict_pickle.txt', 'w')
    pickle.dump(data, new_pickle, protocol=pickle.HIGHEST_PROTOCOL)
    new_pickle.close()
    
    print("any_seq_changes_matrix = get_statistics(working_directory = working_directory, data_table = data, experiment_list = experiment_folders, desired_information = 'any_seq_changes', list_of_passing_first_lengths = 'any', list_of_passing_second_lengths = 'any', aligner = " + aligner + ', bed_name = bed_name, ref_name = ref_name, gff3 = gff3)')
    any_seq_changes_matrix = get_statistics(working_directory = working_directory, data_table = data, experiment_list = experiment_folders, desired_information = 'any_seq_changes', list_of_passing_first_lengths = 'any', list_of_passing_second_lengths = 'any', aligner = aligner, bed_name = bed_name, ref_name = REF,  gff3 = gff3,)
    return True


#in my thesis referred to as Function 7
def HRA_MSA_analyzer(bed_path, HRA_result_path):
    
    clear_change_alignments_not_redundant = open(HRA_result_path + 'clear_change_alignments_not_redundant.txt', 'w')
    ambiguous_change_alignments_not_redundant = open(HRA_result_path + 'ambiguous_change_alignments_not_redundant.txt', 'w')
    failed_alignments_not_redundant = open(HRA_result_path + 'failed_alignments_not_redundant.txt', 'w')
    less_than_9_genomes_alignments_not_redundant = open(HRA_result_path + 'less_than_9_genomes_alignments_not_redundant.txt', 'w')
    InDel_alignments_redundant = open(HRA_result_path + 'InDel_alignments_redundant.txt', 'w')
    substitution_alignments_redundant = open(HRA_result_path + 'substitution_alignments_redundant.txt', 'w')
    csv = open(HRA_result_path + 'alignments_changes.csv', 'w')
    
    clear_change_alignments_not_redundant.write("'num' => number of different characters of the given alignment position - '1' means that all aligned sequences share the same character at that position, and a number greater than '1' accounts for the number of different characters\n'int' => overall interpretation of the given alignment position - annotated as follows:\n'.' => repeat coding position without any differences\n';' => lower-case bases at same position as same upper-case base were likely mistaken for non-repeat coding position\n'*' => important and distinguished differences at repeat coding position\n'~' => unambigous repeat coding position with differences and gaps\n'-' => unambigous repeat coding position with gaps as only difference\n'@' => ambigious uncertain differences at repeat coding position\n' ' => irrelevant non-repeat coding position\n'?' => unaccounted possibility\n\n\n\n")
    ambiguous_change_alignments_not_redundant.write("'num' => number of different characters of the given alignment position - '1' means that all aligned sequences share the same character at that position, and a number greater than '1' accounts for the number of different characters\n'int' => overall interpretation of the given alignment position - annotated as follows:\n'.' => repeat coding position without any differences\n';' => lower-case bases at same position as same upper-case base were likely mistaken for non-repeat coding position\n'*' => important and distinguished differences at repeat coding position\n'~' => unambigous repeat coding position with differences and gaps\n'-' => unambigous repeat coding position with gaps as only difference\n'@' => ambigious uncertain differences at repeat coding position\n' ' => irrelevant non-repeat coding position\n'?' => unaccounted possibility\n\n\n\n")
    failed_alignments_not_redundant.write("'num' => number of different characters of the given alignment position - '1' means that all aligned sequences share the same character at that position, and a number greater than '1' accounts for the number of different characters\n'int' => overall interpretation of the given alignment position - annotated as follows:\n'.' => repeat coding position without any differences\n';' => lower-case bases at same position as same upper-case base were likely mistaken for non-repeat coding position\n'*' => important and distinguished differences at repeat coding position\n'~' => unambigous repeat coding position with differences and gaps\n'-' => unambigous repeat coding position with gaps as only difference\n'@' => ambigious uncertain differences at repeat coding position\n' ' => irrelevant non-repeat coding position\n'?' => unaccounted possibility\n\n\n\n")
    less_than_9_genomes_alignments_not_redundant.write("'num' => number of different characters of the given alignment position - '1' means that all aligned sequences share the same character at that position, and a number greater than '1' accounts for the number of different characters\n'int' => overall interpretation of the given alignment position - annotated as follows:\n'.' => repeat coding position without any differences\n';' => lower-case bases at same position as same upper-case base were likely mistaken for non-repeat coding position\n'*' => important and distinguished differences at repeat coding position\n'~' => unambigous repeat coding position with differences and gaps\n'-' => unambigous repeat coding position with gaps as only difference\n'@' => ambigious uncertain differences at repeat coding position\n' ' => irrelevant non-repeat coding position\n'?' => unaccounted possibility\n\n\n\n")
    InDel_alignments_redundant.write("'num' => number of different characters of the given alignment position - '1' means that all aligned sequences share the same character at that position, and a number greater than '1' accounts for the number of different characters\n'int' => overall interpretation of the given alignment position - annotated as follows:\n'.' => repeat coding position without any differences\n';' => lower-case bases at same position as same upper-case base were likely mistaken for non-repeat coding position\n'*' => important and distinguished differences at repeat coding position\n'~' => unambigous repeat coding position with differences and gaps\n'-' => unambigous repeat coding position with gaps as only difference\n'@' => ambigious uncertain differences at repeat coding position\n' ' => irrelevant non-repeat coding position\n'?' => unaccounted possibility\n\n\n\n")
    substitution_alignments_redundant.write("'num' => number of different characters of the given alignment position - '1' means that all aligned sequences share the same character at that position, and a number greater than '1' accounts for the number of different characters\n'int' => overall interpretation of the given alignment position - annotated as follows:\n'.' => repeat coding position without any differences\n';' => lower-case bases at same position as same upper-case base were likely mistaken for non-repeat coding position\n'*' => important and distinguished differences at repeat coding position\n'~' => unambigous repeat coding position with differences and gaps\n'-' => unambigous repeat coding position with gaps as only difference\n'@' => ambigious uncertain differences at repeat coding position\n' ' => irrelevant non-repeat coding position\n'?' => unaccounted possibility\n\n\n\n")
    
    clear_change_alignments_not_redundant.close()
    ambiguous_change_alignments_not_redundant.close()
    failed_alignments_not_redundant.close()
    less_than_9_genomes_alignments_not_redundant.close()
    InDel_alignments_redundant.close()
    substitution_alignments_redundant.close()
    csv.close()
    
    opened_path = open(bed_path, 'r')
    
    ENST_ENSG_dict = {}
    for line in opened_path:
        split_line = line.rstrip('\n').split('\t')
        ENST = split_line[4]
        ENSG = split_line[5].split('.')[0]
        ENST_ENSG_dict[ENST] = ENSG
    
    InDel_alignments_redundant_summary_list = []
    substitution_alignments_redundant_summary_list = []
    clear_change_alignments_not_redundant_summary_list = []
    ambiguous_change_alignments_not_redundant_summary_list = []
    failed_alignments_not_redundant_summary_list = []
    less_than_9_genomes_alignments_not_redundant_summary_list = []
    all_alignments_summary_list = []
    
    list_of_repeat_coding_changes = []
    list_of_repeat_coding_changes_minus_gaps = []
    list_of_non_repeat_coding_changes = []
    list_of_non_repeat_coding_changes_minus_gaps = []
    
    alignment_info_dict = defaultdict(list)
    df_dict = defaultdict(list)
    hundred_quantile_dict = defaultdict(list)
    running_repeat_id = 0
    complete_ref_ENSG_GO_list = []
    for each_file in listdir(HRA_result_path):
        strand = '?'
        if re.match(r'ENST\d*\.\d*.cc\d*.mafft.nc.fa.aln.fa.aln.phased_case_change', each_file):# and 'cc2' not in each_file and 'cc3' not in each_file:
            cc_num = each_file.split('.')[2]
            count = 0
            seq_list = []
            rep_string = ''
            summary = ''
            seqs_string = ''
            num_of_changes_string = ''
            
            num_of_repeat_coding_changes = 0
            num_of_non_repeat_coding_changes = 0
            
            total_identity = 0
            non_repeat_identity = 0
            repeat_identity = 0
            #print('\nfound file: ' + each_file)
            summary += 'found file: ' + each_file + '\n'
            
            ENST = each_file.split('.cc')[0]
            ENSG = ENST_ENSG_dict[ENST]
            complete_ref_ENSG_GO_list.append(ENSG)
            summary += 'ENST: ' + ENST + '\n' + 'ENSG: ' + ENSG + '\n'
            
            label_seq_list = []
            seq_line = ''
            label = ''
            orf = ''
            
            x = open(HRA_result_path + each_file, 'r')
            for line in x:
                if 'ORF:' in line:
                    orf = line.rstrip('\n').rstrip('\r').split(':')[1]                
                    label = line.rstrip('\n').rstrip('\r').split(':')[0]
                    orf_line = line.rstrip('\n').rstrip('\r').split(':')[1]
                    if strand == '-':
                        orf_line = seq_line[::-1]
                    
                    seqs_string += label + ':' + orf_line + '\n'
            x.close()
            
            if orf == '':
                print('no orf in ' + str(each_file) + ': ' + str(orf))
                print('note: the reason for this is probably a lower-case "ref:" in the msa')        
                
            else:    
                if '012' in orf:
                    strand = '+'
                elif '210' in orf:
                    strand = '-'
                    orf = orf[::-1]
                else:
                    print('unrecognizable orf in ' + str(each_file) + ': ' + str(orf))
                    print('note: the reason for this is probably a lower-case "ref:" in the msa')
        
                #produce a non-redundant sorted list of all repeat_region associated positions from all genomes in the msa    
                x = open(HRA_result_path + each_file, 'r')
                list_of_positions_of_interest = []
                for line in x:
                    if not 'ORF:' in line:
                        #seqs_string += line
                        label = line.rstrip('\n').rstrip('\r').split(':')[0]
                        seq_line = line.rstrip('\n').rstrip('\r').split(':')[1].upper()
                        if strand == '-':
                            seq_line = seq_line[::-1].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()
                        pattern = re.finditer(r'((CA[AG])+.{3})*(CA[AG]){4,}(.{3}(CA[AG])+)*', seq_line, flags = re.I)
                        for repeat in pattern:
                            repeat_start = repeat.start()
                            repeat_stop = repeat.end()
                            for position_of_interest in range (repeat_start, repeat_stop + 1, 1):
                                list_of_positions_of_interest.append(position_of_interest)
                x.close()
                positions_of_interest = list(set(list_of_positions_of_interest))
                positions_of_interest.sort()
                    
                # iterate over the list of all positions of interest and produce seperate tuples with uninterupted sequences of positions of interest 
                list_of_range_tuples_of_interest = []
                start_pos = positions_of_interest[0]
                last_pos = positions_of_interest[0]
                for current_pos in positions_of_interest:
                    if current_pos == last_pos + 1:
                        last_pos = current_pos
                    elif current_pos > last_pos + 1:
                        list_of_range_tuples_of_interest.append((start_pos, last_pos + 1))
                        start_pos = current_pos
                        last_pos = current_pos
                list_of_range_tuples_of_interest.append((start_pos, last_pos + 1))
                x = open(HRA_result_path + each_file, 'r')                
                
                #iterate over all non-rof lines
                for line in x:
                    if not 'ORF:' in line:
                        label = line.rstrip('\n').rstrip('\r').split(':')[0]
                        seq_line = line.rstrip('\n').rstrip('\r').split(':')[1].upper()
                        
                        #translate minus strands into their reverse compliments
                        if strand == '-':
                            seq_line = seq_line[::-1].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()
                        TROI_count = 0
                        #also iterate over all tuples of interest
                        for tuple_of_interest in list_of_range_tuples_of_interest:
                            tuple_seq_line = seq_line[tuple_of_interest[0]:tuple_of_interest[1]]
                            #warning!
                            #please notice that this regex was modified and tolerates gap-triplets
                            pattern = re.finditer(r'((C(---)*A(---)*[AG])+.{3}(---)*)*((---)*C(---)*A(---)*[AG](---)*){4,}((---)*.{3}(C(---)*A(---)*[AG])+)*', tuple_seq_line, flags = re.I)
                            tuple_orf = orf[tuple_of_interest[0]:tuple_of_interest[1]]
                            
                            #merge all matches within the tuple region of interest into one sequence to read through (long insertions of non Q bases will therefore be included in the analysis if they have an orthologue that features no large insert at the position and which is therefore recognized as one consistent polyQ region)
                            repeat_start = len(tuple_seq_line)
                            repeat_stop = 0
                            for repeat in pattern:
                                if repeat_start > repeat.start():
                                    repeat_start = repeat.start()
                                if repeat_stop < repeat.end():
                                    repeat_stop = repeat.end()
                            
                            codon = ''
                            translation = ''
                            base_count = repeat_start
                            stretch_before_repeat_start = seq_line[:(tuple_of_interest[0] + repeat_start)].replace('-','')[-9:]
                            stretch_after_repeat_stop = seq_line[(tuple_of_interest[0] + repeat_stop):].replace('-','')[:9]
                            if ((len(stretch_before_repeat_start) > 3 and stretch_before_repeat_start.count('N') < 3) or (len(stretch_before_repeat_start) <= 3 and stretch_before_repeat_start.count('N') == 0)) and ((len(stretch_after_repeat_stop) > 3 and stretch_after_repeat_stop.count('N') < 3) or (len(stretch_after_repeat_stop) <= 3 and stretch_after_repeat_stop.count('N') == 0)):
                                assumed_orf_for_gaps = None
                                orf_pos = False
                                codon_pos_1 = repeat_start
                                codon_pos_2 = codon_pos_1 + 1
                                codon_pos_3 = codon_pos_1 + 2
                                EOL_flag = False
                                while EOL_flag == False:
                                    for base_count in range(repeat_start, repeat_stop, 1):
                                        if base_count == codon_pos_1:    
                                            try:
                                                test_letter = tuple_orf[codon_pos_1]
                                            except IndexError:
                                                EOL_flag = True
                                                print('BREAK!')
                                                break
                                            while tuple_orf[codon_pos_1] == ' ' or tuple_orf[codon_pos_1] == '1' or tuple_orf[codon_pos_1] == '2' or tuple_seq_line[codon_pos_1] == '-' or tuple_seq_line[codon_pos_1] == ' ':
                                                codon_pos_1 += 1
                                                try:
                                                    test_letter = tuple_orf[codon_pos_1]
                                                except IndexError:
                                                    EOL_flag = True
                                                    print('BREAK!')
                                                    break
                                            codon_pos_2 = codon_pos_1 + 1
    
                                            try:
                                                test_letter = tuple_orf[codon_pos_2]
                                            except IndexError:
                                                EOL_flag = True
                                                print('BREAK!')
                                                break
                                            while tuple_orf[codon_pos_2] == ' ' or tuple_orf[codon_pos_2] == '0' or tuple_orf[codon_pos_2] == '2' or tuple_seq_line[codon_pos_2] == '-' or tuple_seq_line[codon_pos_2] == ' ':
                                                codon_pos_2 += 1
                                                try:
                                                    test_letter = tuple_orf[codon_pos_2]
                                                except IndexError:
                                                    EOL_flag = True
                                                    print('BREAK!')
                                                    break
                                            codon_pos_3 = codon_pos_2 + 1
    
                                            try:
                                                test_letter = tuple_orf[codon_pos_3]
                                            except IndexError:
                                                EOL_flag = True
                                                print('BREAK!')
                                                break
                                            while tuple_orf[codon_pos_3] == ' ' or tuple_orf[codon_pos_3] == '1' or tuple_orf[codon_pos_3] == '0' or tuple_seq_line[codon_pos_3] == '-' or tuple_seq_line[codon_pos_3] == ' ':
                                                codon_pos_3 += 1
                                                try:
                                                    test_letter = tuple_orf[codon_pos_3]
                                                except IndexError:
                                                    EOL_flag = True
                                                    print('BREAK!')
                                                    break
    
                                            if EOL_flag == True:
                                                break
                                            
                                            
                                            if tuple_orf[codon_pos_1] == '1' or tuple_orf[codon_pos_1] == '2':
                                                print('an issue occurred with this codon: ' + tuple_seq_line[codon_pos_1] + tuple_seq_line[codon_pos_2] + tuple_seq_line[codon_pos_3])
                                                print('an issue occurred with this orf: ' + tuple_orf[codon_pos_1] + tuple_orf[codon_pos_2] + tuple_orf[codon_pos_3])
                                                raise Exception('It seems that orf string positions were skipped for codon_pos_1')
                                            if tuple_orf[codon_pos_2] == '0' or tuple_orf[codon_pos_2] == '2':
                                                print('an issue occurred with this codon: ' + tuple_seq_line[codon_pos_1] + tuple_seq_line[codon_pos_2] + tuple_seq_line[codon_pos_3])
                                                print('an issue occurred with this orf: ' + tuple_orf[codon_pos_1] + tuple_orf[codon_pos_2] + tuple_orf[codon_pos_3])
                                                raise Exception('It seems that orf string positions were skipped for codon_pos_2')
                                            if tuple_orf[codon_pos_3] == '1' or tuple_orf[codon_pos_3] == '0':
                                                print('an issue occurred with this codon: ' + tuple_seq_line[codon_pos_1] + tuple_seq_line[codon_pos_2] + tuple_seq_line[codon_pos_3])
                                                print('an issue occurred with this orf: ' + tuple_orf[codon_pos_1] + tuple_orf[codon_pos_2] + tuple_orf[codon_pos_3])
                                                raise Exception('It seems that orf string positions were skipped for codon_pos_3')
                                            
                                            codon = tuple_seq_line[codon_pos_1] + tuple_seq_line[codon_pos_2] + tuple_seq_line[codon_pos_3]
                                            
                                            if codon[0] in 'ACTG' and codon[1] in 'ACTG' and codon[2] in 'ACTG':
                                                translation = AAS_table[codon]
                                                if codon == 'CAA' or codon == 'CAG':
                                                    codon_type = codon
                                                elif codon in ['AAG', 'GAG', 'AAA', 'CTG', 'CCA', 'CGG', 'CAT', 'CGA', 'TAG', 'TAA', 'CAC', 'CCG', 'GAA', 'CTA']:
                                                    #"sbs" stands for single base substitution away from Q-coding codon
                                                    codon_type = 'sbs'
                                                else:
                                                    #mbs stands for multiple  base substitutions away from Q-coding codon
                                                    codon_type = 'mbs'
                                            else:
                                                translation = 'X'
                                                #aac stands for any ambiguous codon
                                                codon_type = 'aac'
                                            
                                            orf_pos = tuple_orf[codon_pos_1]
                                            if orf_pos == '0':
                                                pass
                                            elif orf_pos == '-':
                                                orf_pos = '0'
                                            else:
                                                print('unexpected orf: ' + orf_pos)
                                                raise Exception('The current orf_pos is neither an annotated in-frame codon start nor is it a gap')
                                            codon_step_count = 0
                                            for base in codon:
                                                    
                                                #file denotes the file of origin
                                                df_dict['file'].append(each_file)
                                                #ENST denotes the ENST
                                                df_dict['ENST'].append(ENST)
                                                #ENSG denotes the ENSG
                                                df_dict['ENSG'].append(ENSG)
                                                #cc_num denotes the running file number which is also found in 'file'
                                                df_dict['cc_num'].append(cc_num)
                                                #label denotes the bam sequence/ref of origin
                                                df_dict['label'].append(label.upper())
                                                #base denotes the base of the read
                                                df_dict['base'].append(base)
                                                #codon denotes the translated base-triplet that the current base contributes to
                                                df_dict['codon'].append(codon)
                                                #codon_type denotes what group the base-triplet is assigned for for further analysis
                                                df_dict['codon_type'].append(codon_type)
                                                #translation denotes the amino acid that the current base contributes to
                                                df_dict['translation'].append(translation)
                                                #orf denotes the bases position in the reading frame
                                                df_dict['orf'].append(orf_pos)
                                                #position denotes the position in the MSA as found under 'file'
                                                if strand == '+':
                                                    df_dict['position'].append(base_count + codon_step_count + tuple_of_interest[0] + repeat_start)
                                                elif strand == '-':
                                                    df_dict['position'].append(len(seq_line) - (base_count + codon_step_count + tuple_of_interest[0] + repeat_start))
                                                #repeat_length denotes the length of the entire repeat
                                                repeat_length = ((repeat_stop - repeat_start) / 3)
                                                df_dict['repeat_length'].append(repeat_length)
                                                #in_repeat_codon_rank denotes the rank of the codon within the detected repeat
                                                in_repeat_codon_rank = int((base_count - (int(orf_pos) + repeat_start)) / 3)
                                                df_dict['in_repeat_codon_rank'].append(in_repeat_codon_rank)
                                                #in_repeat_codon_rank denotes the rank of the codon within the detected repeat
                                                relative_pos = (in_repeat_codon_rank / repeat_length)
                                                df_dict['relative_pos'].append(relative_pos)
                                                #ROI_length denotes the ROI_length of the entire ROI
                                                df_dict['ROI_length'].append(tuple_of_interest[1] - tuple_of_interest[0])
                                                #ROI denotes the uninterrupted stretch of positions in the MSA that contributes to one aligned region of interest
                                                df_dict['ROI'].append(str(tuple_of_interest[0]) + '-' + str(tuple_of_interest[1]))
                                                #TROI_count denotes the MSA's running TROI counter
                                                df_dict['TROI_count'].append(TROI_count)
                                                #running_repeat_id denotes a running read counter of the TROI's covered for each MSA
                                                df_dict['running_repeat_id'].append(running_repeat_id)
                                                df_dict['in_repeat'].append(True)
                                                codon_step_count += 1
                                                
                                            full_percent_range_start = int(relative_pos * 1000) + 1
                                            full_percent_range_stop = int(((in_repeat_codon_rank + 1) / repeat_length) * 1000) + 1                                        
                                            
                                            for full_percent_quantile in range(full_percent_range_start,full_percent_range_stop,1):
                                                hundred_quantile_dict['codon'].append(codon)
                                                hundred_quantile_dict['codon_type'].append(codon_type)
                                                hundred_quantile_dict['label'].append(label)
                                                hundred_quantile_dict['repeat_length'].append(repeat_length)
                                                hundred_quantile_dict['full_percent_quantile'].append(full_percent_quantile)
                                                hundred_quantile_dict['running_repeat_id'].append(running_repeat_id)
                                            
                                            codon_pos_1 = codon_pos_3 + 1                                        
                                            
                                        else:
                                            pass
                                    EOL_flag = True
                                    
                                    
                            else:
                                #filtered out due to unresolved flanking sequence
                                break
                            TROI_count += 1
                            running_repeat_id += 1
                x.close()                                
        
                x = open(HRA_result_path + each_file, 'r')
                orf = ''
                for line in x:        
                    if 'REF:' in line:                    
                        label = line.rstrip('\n').rstrip('\r').split(':')[0]
                        seq_line = line.rstrip('\n').rstrip('\r').split(':')[1]
                        new_line = seq_line
                        seqs_string += label + ':' + seq_line + '\n'
                        seq_list.append(new_line)
    
                    elif 'ORF:' in line:
                        
                        seq_line = line.rstrip('\n').rstrip('\r').split(':')[1]
                        orf = seq_line
                x.close()
                
                x = open(HRA_result_path + each_file, 'r')
                for line in x:
                    if not 'ORF:' in line and not 'REF:' in line:
                        label = line.rstrip('\n').rstrip('\r').split(':')[0]
                        seq_line = line.rstrip('\n').rstrip('\r').split(':')[1]
                        new_line = seq_line
                        seqs_string += label + ':' + seq_line + '\n'
                        seq_list.append(new_line)
                x.close()
         
                ROI_length = len(orf)
                #open dataframe
                df = pd.DataFrame.from_dict(df_dict)
                #only include rows of the current file
                current_df = df.loc[df['cc_num']==cc_num]#['file']['label']['base']['orf']['position']#['ENST']['ENSG']['cc_num']
                running_repeat_id
                if current_df.shape[0] < 4:
                    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
                        print('current_df: ' + str(current_df))
                    print(cc_num)
                    raise Exception()
    
                if len(orf) != len(seq_list[0]):
                    continue
    
    
                pos = -1
                for step in orf:
                    pos += 1
                    #make a lol of all positions per column
                    
                    #the next line doesn't work in alignments with gaps, unfortunately, where gaps lead to a faulty num and int string
                    #pos_list = current_df.loc[current_df['position']== pos]['base'].tolist()
                    pos_list = [seq[pos] for seq in seq_list]
    
                    #iterate through the lol and turn each list into a set                
                    #('C' in pos_list or 'A' in pos_list or 'G' in pos_list or 'T' in pos_list or 'c' in pos_list or 'a' in pos_list or 'g' in pos_list or 't' in pos_list or 'R' in pos_list or 'Y' in pos_list or 'W' in pos_list or 'S' in pos_list or 'M' in pos_list or 'K' in pos_list or 'X' in pos_list or 'N' in pos_list or '-' in pos_list or 'H' in pos_list or 'B' in pos_list or 'D' in pos_list or 'V' in pos_list or 'r' in pos_list or 'y' in pos_list or 'w' in pos_list or 's' in pos_list or 'm' in pos_list or 'k' in pos_list or 'h' in pos_list or 'b' in pos_list or 'd' in pos_list or 'v' in pos_list or 'x')
                    
                    if len(set(pos_list)) > 1 and len(set([element.upper() for element in pos_list])) == 1:
                        rep_string += ';'#==likely mistaken lower case
                        total_identity += 1
                        repeat_identity += 1
                    #>1 & only C A G T
                    elif len(set(pos_list)) > 1 and not ('c' in pos_list or 'a' in pos_list or 'g' in pos_list or 't' in pos_list or 'R' in pos_list or 'Y' in pos_list or 'W' in pos_list or 'S' in pos_list or 'M' in pos_list or 'K' in pos_list or 'X' in pos_list or 'N' in pos_list or '-' in pos_list or 'H' in pos_list or 'B' in pos_list or 'D' in pos_list or 'V' in pos_list or 'r' in pos_list or 'y' in pos_list or 'w' in pos_list or 's' in pos_list or 'm' in pos_list or 'k' in pos_list or 'h' in pos_list or 'b' in pos_list or 'd' in pos_list or 'v' in pos_list or 'x' in pos_list or 'n' in pos_list):
                        rep_string += '*'#==important
                    
                    #>1 & only c a g t r y w s m k h b d v x -
                    elif len(set(pos_list)) > 1 and not ('C' in pos_list or 'A' in pos_list or 'G' in pos_list or 'T' in pos_list or 'R' in pos_list or 'Y' in pos_list or 'W' in pos_list or 'S' in pos_list or 'M' in pos_list or 'K' in pos_list or 'X' in pos_list or 'N' in pos_list or 'H' in pos_list or 'B' in pos_list or 'D' in pos_list or 'V' in pos_list):
                        rep_string += ' '#==irrelevant
                    #>2 & only C A G T -
                    elif len(set(pos_list)) > 2 and not ('c' in pos_list or 'a' in pos_list or 'g' in pos_list or 't' in pos_list or 'R' in pos_list or 'Y' in pos_list or 'W' in pos_list or 'S' in pos_list or 'M' in pos_list or 'K' in pos_list or 'X' in pos_list or 'N' in pos_list or 'H' in pos_list or 'B' in pos_list or 'D' in pos_list or 'V' in pos_list or 'r' in pos_list or 'y' in pos_list or 'w' in pos_list or 's' in pos_list or 'm' in pos_list or 'k' in pos_list or 'h' in pos_list or 'b' in pos_list or 'd' in pos_list or 'v' in pos_list or 'x' in pos_list or 'n' in pos_list):
                        rep_string += '~'#==unambigous SEQ with cange and gap
                    #==2 & only C A G T - OR ==3 and C,c,- or A,a,- or G,g,- or T,t,-
                    elif set(pos_list) == set(['A','-',]) or set(pos_list) == set(['G','-',]) or set(pos_list) == set(['C','-',]) or set(pos_list) == set(['T','-',]) or set(pos_list) == set(['a','A','-',]) or set(pos_list) == set(['c','C','-',]) or set(pos_list) == set(['g','G','-',]) or set(pos_list) == set(['t','T','-',]):
                        rep_string += '-'#==unambigous SEQ with gap
                    #>1
                    elif len(set(pos_list)) > 1:
                        rep_string += '@'#==ambigious
                    #==1
                    elif len(set(pos_list)) == 1 and pos_list[0] in ['C','A','G','T',]:
                        rep_string += '.'#==polyQ coding but nothing changed
                        total_identity += 1
                        repeat_identity += 1
                    elif len(set(pos_list)) == 1:
                        rep_string += ' '#==nothing to report
                        total_identity += 1
                        non_repeat_identity += 1
                    #???
                    else:
                        rep_string += ' '#==unaccounted possibility
                    num_of_changes_string += str(len(set(pos_list)))
                    
                    non_arithmetic_csv_row = ''
                    if 'C' in pos_list or 'A' in pos_list or 'G' in pos_list or 'T' in pos_list: 
                        num_of_repeat_coding_changes += 1
                        list_of_repeat_coding_changes.append(len(set(pos_list)))
                        cg = str(len(set(pos_list)))
                        if len(set(pos_list)) > 1 and '-' in pos_list:
                            list_of_repeat_coding_changes_minus_gaps.append((len(set(pos_list)) -1))
                        else:
                            list_of_repeat_coding_changes_minus_gaps.append((len(set(pos_list))))
                            
                    else:
                        num_of_non_repeat_coding_changes += 1
                        list_of_non_repeat_coding_changes.append(len(set(pos_list)))
                        ncg = str(len(set(pos_list)))
                        if len(set(pos_list)) > 1 and '-' in pos_list:
                            list_of_non_repeat_coding_changes_minus_gaps.append((len(set(pos_list)) -1))
                        else:
                            list_of_non_repeat_coding_changes_minus_gaps.append((len(set(pos_list))))
                
                summary += str(len(seq_list)) + ' sequences total\n'
                seqs_string += '\nnum:' + num_of_changes_string + '\nint:' + rep_string + '\n'
                genomes_counter = 0
                missing_string = '(missing genomes: '
                
                for experiment in ['REF:','BIC:','GB2:','KK1:','LOS:','NE1:','SII:','WC1:','ZVE:']:
                    if experiment in seqs_string:
                        genomes_counter += 1
                    if experiment not in seqs_string:
                        missing_string += experiment.rstrip(':') + ', '
                        
                sum_coding = 0
                for fig in list_of_repeat_coding_changes:
                    sum_coding += fig
                sum_non_coding = 0
                for fig in list_of_non_repeat_coding_changes:
                    sum_non_coding += fig
                sum_coding_gaps = 0
                for fig in list_of_repeat_coding_changes_minus_gaps:
                    sum_coding_gaps += fig
                sum_non_coding_gaps = 0
                for fig in list_of_non_repeat_coding_changes_minus_gaps:
                    sum_non_coding_gaps += fig
        
                average_coding_w_gaps = 'NA'
                if num_of_repeat_coding_changes != 0:
                    average_coding_w_gaps = (sum_coding / num_of_repeat_coding_changes)
                average_coding_w_o_gaps = 'NA'
                if num_of_repeat_coding_changes != 0:
                    average_coding_w_o_gaps = (sum_coding_gaps / num_of_repeat_coding_changes)
                average_non_coding_w_gaps = 'NA'
                if num_of_non_repeat_coding_changes != 0:
                    average_non_coding_w_gaps = (sum_non_coding / num_of_non_repeat_coding_changes)
                average_non_coding_w_o_gaps = 'NA'
                if num_of_non_repeat_coding_changes != 0:
                    average_non_coding_w_o_gaps = (sum_non_coding_gaps / num_of_non_repeat_coding_changes)
                
                list_of_relevant_line_numbers = alignment_info_dict[each_file + '_REF']
                if len(list_of_relevant_line_numbers) > 1:
                    raise Exception('too many REF lines')
                
                summary += ('' +
                            #'Chromosome: ' + chromosome + '\n' + 
                            'Strand: ' + strand + '\n' + 
                            #'ENST start/stop: ' + ENST_start + '/' + ENST_stop + '\n' + 
                            #'repeat start/stop: ' + repeat_start + '/' + repeat_stop + '\n' + 
                            str(genomes_counter) + ' genomes total\n' + 
                            missing_string + ')\n')
                if len(rep_string) == 0:
                    summary += 'total pile-up identity: 0%\n'
                else:
                    summary += 'total pile-up identity: ' + str((total_identity/len(rep_string)) * 100) + '% (' + str(total_identity) + '/' + str(len(rep_string)) + ')\n'
                if rep_string.count(' ') == 0:
                    summary += 'non repeat pile-up identity: 0%\n'
                else:
                    summary += 'non repeat pile-up identity: ' + str((non_repeat_identity/rep_string.count(' ')) * 100) + '% (' + str(non_repeat_identity) + '/' + str(rep_string.count(' ')) + ')\n'
                if (len(rep_string) - rep_string.count(' ')) == 0:
                    summary += 'repeat pile-up identity: 0%\n'
                else:
                    summary += 'repeat pile-up identity: ' + str((repeat_identity/(len(rep_string) - rep_string.count(' '))) * 100) + '% (' + str(repeat_identity) + '/' + str((len(rep_string) - rep_string.count(' '))) + ')\n'
                
                for disease_associated_protein, synonym, disease_associated_ENSG  in [('SBMA', 'AR-001', 'ENSG00000169083'),('HD', 'HTT', 'ENSG00000197386'),('DRPLA', 'DRPLA', 'ENSG00000111676'),('SCA1', 'ATXN1', 'ENSG00000124788'),('SCA2', 'ATXN2', 'ENSG00000204842'),('SCA3', 'ATXN3', 'ENSG00000066427'),('SCA6', 'CACNA1A', 'ENSG00000141837'),('SCA7', 'ATXN7', 'ENSG00000163635'),('SCA17', 'TBP', 'ENSG00000112592'),]:
                    if disease_associated_ENSG == ENSG.split('.')[0]:
                        summary += 'Attention: Disease associated gene!\ndisease_associated_protein: ' + disease_associated_protein + '\nsynonym: '  + synonym + '\ndisease_associated_ENSG: ' + disease_associated_ENSG + '\n'
                summary += ('average number of differences per repeat coding position: ' + str(average_coding_w_gaps) + '\n' + 
                            'average number of differences per repeat coding position (excluding gaps): ' + str(average_coding_w_o_gaps) + '\n' + 
                            'average number of differences per non-repeat coding position: ' + str(average_non_coding_w_gaps) + '\n' + 
                            'average number of differences per non-repeat coding position (excluding gaps): ' + str(average_non_coding_w_o_gaps) + '\n\n' + 
                            seqs_string +
                            '\n\n\n')
        
                if '-' in rep_string or '~' in rep_string:
                    InDel_alignments_redundant_summary_list.append((summary, ENSG))
                    
                if '*' in rep_string or '~' in rep_string or '@' in rep_string:
                    substitution_alignments_redundant_summary_list.append((summary, ENSG))
                    
                if '*' in rep_string or '~' in rep_string:
                    clear_change_alignments_not_redundant_summary_list.append((summary, ENSG))
        
                elif '@' in rep_string or '-' in rep_string:
                    ambiguous_change_alignments_not_redundant_summary_list.append((summary, ENSG))
        
                elif '?' in rep_string and len(seq_list) == 9:
                    failed_alignments_not_redundant_summary_list.append((summary, ENSG))
        
                elif genomes_counter < 9:
                    less_than_9_genomes_alignments_not_redundant_summary_list.append((summary, ENSG))
        
                else:
                    raise Exception('something went wrong')
                
                all_alignments_summary_list .append((summary, ENSG))
    
    InDel_alignments_redundant_summary_list_sorted = sorted(InDel_alignments_redundant_summary_list, key=return2nd)
    substitution_alignments_redundant_summary_list_sorted = sorted(substitution_alignments_redundant_summary_list, key=return2nd) 
    clear_change_alignments_not_redundant_summary_list_sorted = sorted(clear_change_alignments_not_redundant_summary_list, key=return2nd) 
    ambiguous_change_alignments_not_redundant_summary_list_sorted = sorted(ambiguous_change_alignments_not_redundant_summary_list, key=return2nd) 
    failed_alignments_not_redundant_summary_list_sorted = sorted(failed_alignments_not_redundant_summary_list, key=return2nd) 
    less_than_9_genomes_alignments_not_redundant_summary_list_sorted = sorted(less_than_9_genomes_alignments_not_redundant_summary_list, key=return2nd) 
    all_alignments_summary_list = sorted(all_alignments_summary_list, key=return2nd) 
    
    InDel_alignments_ENSG_GO_string = ''
    substitution_alignments_ENSG_GO_string = ''
    clear_change_alignments_ENSG_GO_string = ''
    ambiguous_change_alignments_ENSG_GO_string = ''
    failed_alignments_ENSG_GO_string = ''
    less_than_9_genomes_alignments_ENSG_GO_string = ''
    all_alignments_summary_list_ENSG_GO_string = ''
    complete_ref_ENSG_GO_string = ",".join(complete_ref_ENSG_GO_list)
    
    
    for summary, ENSG in InDel_alignments_redundant_summary_list_sorted:
        InDel_alignments_redundant = open(HRA_result_path + 'InDel_alignments_redundant.txt', 'a')
        InDel_alignments_redundant.write(summary)
        InDel_alignments_redundant.close()
        InDel_alignments_ENSG_GO_string += ENSG + ','
    for summary, ENSG in substitution_alignments_redundant_summary_list_sorted:
        substitution_alignments_redundant = open(HRA_result_path + 'substitution_alignments_redundant.txt', 'a')
        substitution_alignments_redundant.write(summary)
        substitution_alignments_redundant.close()
        substitution_alignments_ENSG_GO_string += ENSG + ','
    for summary, ENSG in clear_change_alignments_not_redundant_summary_list_sorted:
        clear_change_alignments_not_redundant = open(HRA_result_path + 'clear_change_alignments_not_redundant.txt', 'a')
        clear_change_alignments_not_redundant.write(summary)
        clear_change_alignments_not_redundant.close()
        clear_change_alignments_ENSG_GO_string += ENSG + ','
    for summary, ENSG in ambiguous_change_alignments_not_redundant_summary_list_sorted:
        ambiguous_change_alignments_not_redundant = open(HRA_result_path + 'ambiguous_change_alignments_not_redundant.txt', 'a')
        ambiguous_change_alignments_not_redundant.write(summary)
        ambiguous_change_alignments_not_redundant.close()
        ambiguous_change_alignments_ENSG_GO_string += ENSG + ','
    for summary, ENSG in failed_alignments_not_redundant_summary_list_sorted:
        failed_alignments_not_redundant = open(HRA_result_path + 'failed_alignments_not_redundant.txt', 'a')
        failed_alignments_not_redundant.write(summary)
        failed_alignments_not_redundant.close()
        failed_alignments_ENSG_GO_string += ENSG + ','
    for summary, ENSG in less_than_9_genomes_alignments_not_redundant_summary_list_sorted:
        less_than_9_genomes_alignments_not_redundant = open(HRA_result_path + 'less_than_9_genomes_alignments_not_redundant.txt', 'a')
        less_than_9_genomes_alignments_not_redundant.write(summary)
        less_than_9_genomes_alignments_not_redundant.close()
        less_than_9_genomes_alignments_ENSG_GO_string += ENSG + ','
    for summary, ENSG in all_alignments_summary_list:
        all_alignments_summary_list = open(HRA_result_path + 'all_alignments_summary_list.txt', 'a')
        all_alignments_summary_list.write(summary)
        all_alignments_summary_list.close()
        all_alignments_summary_list_ENSG_GO_string += ENSG + ','
    
    for alignment_file, GO_string in [(HRA_result_path + 'InDel_alignments_redundant.txt',InDel_alignments_ENSG_GO_string),
                           (HRA_result_path + 'substitution_alignments_redundant.txt',substitution_alignments_ENSG_GO_string),
                           (HRA_result_path + 'clear_change_alignments_not_redundant.txt',clear_change_alignments_ENSG_GO_string),
                           (HRA_result_path + 'ambiguous_change_alignments_not_redundant.txt',ambiguous_change_alignments_ENSG_GO_string),
                           (HRA_result_path + 'failed_alignments_not_redundant.txt',failed_alignments_ENSG_GO_string),
                           (HRA_result_path + 'less_than_9_genomes_alignments_not_redundant.txt',less_than_9_genomes_alignments_ENSG_GO_string),
                           (HRA_result_path + 'all_alignments_summary_list.txt',all_alignments_summary_list_ENSG_GO_string)]:
    
        open_alignment = open(alignment_file, 'r')
        alignment_count = open_alignment.read().count('found file: ENST')
        open_alignment.close()
        open_alignment = open(alignment_file, 'a')
        open_alignment.write('\n\n\n\n\ntotal number of alignments: ' + str(alignment_count) + '\nComma seperated ENSGs for Gene Ontology Enrichment Analysis: ' + GO_string  + '\nComma seperated reference ENSGs (== all detected ENSGs) for Gene Ontology Enrichment Analysis: ' + complete_ref_ENSG_GO_string + '\nEND OF FILE')
        open_alignment.close()
    
    length_list = [(list_of_repeat_coding_changes, len(list_of_repeat_coding_changes), 'list_of_repeat_coding_changes'), 
                   (list_of_repeat_coding_changes_minus_gaps, len(list_of_repeat_coding_changes_minus_gaps), 'list_of_repeat_coding_changes_minus_gaps'), 
                   (list_of_non_repeat_coding_changes, len(list_of_non_repeat_coding_changes), 'list_of_non_repeat_coding_changes'), 
                   (list_of_non_repeat_coding_changes_minus_gaps, len(list_of_non_repeat_coding_changes_minus_gaps), 'list_of_non_repeat_coding_changes_minus_gaps'),]
    
    sorted_list = sorted(length_list, key=return2nd, reverse=False)
    csv = open(HRA_result_path + 'alignments_changes.csv', 'a')
    for each_list, each_length, each_string in sorted_list:
        csv.write(each_string + ',')
    csv.write('\n')
    longest_length = sorted_list[0][1]
    for step in range(0,longest_length,1):
        for column in sorted_list:
            each_list = column[0]
            each_length = column[1]
            each_string = column[2]
            csv.write('{},'.format(each_list[step]))
        csv.write('\n')
    
    csv.close()
    
    #parse and analyse the df
    
    df = pd.DataFrame.from_dict(df_dict)
    df.to_pickle('main_df.pkl')
    
    
    hundred_quantile_df = pd.DataFrame.from_dict(hundred_quantile_dict)
    hundred_quantile_df.to_pickle('hundred_quantile_df.pkl')
    
    running_repeat_id_list = list(set(df['running_repeat_id'].tolist()))
    cc_num_list = list(set(df['cc_num'].tolist()))
    ENST_list = list(set(df['ENST'].tolist()))
    ENSG_list = list(set(df['ENSG'].tolist()))
    label_list = list(set(df['label'].tolist()))
    label_list.sort()
    
    csv_file = open('codon_types_per_label.csv', 'w')
    csv_file.close()
    csv_file = open('codon_types_per_label.csv', 'a')
    
    for label in label_list:
        csv_file.write(label + '_aac,')
        csv_file.write(label + '_CAA,')
        csv_file.write(label + '_CAG,')
        csv_file.write(label + '_sbs,')
        csv_file.write(label + '_mbs,')
    csv_file.write('\n')
    print('label_list: ' + str(label_list))
    
    violin_plot_dict = defaultdict(list)
    correlation_dict = defaultdict(list)
    for listed_cc_num in cc_num_list:
        for listed_label in label_list:
            print('listed_label: ' + str(listed_label))
            codon_types = df.loc[(df['label'] == listed_label) & (df['cc_num'] == listed_cc_num) & (df['orf'] == '0')]['codon_type'].values.tolist()
            print('codon_types: ' + str(codon_types) + '\nset of codon_types: ' + str(set(codon_types)) + '\n')
            
            if len(codon_types) == 0:
                csv_file.write('0,0,0,')
    
            elif len(codon_types) > 0:
                total_count = len(codon_types)    
                aac_count = codon_types.count('aac')
                CAA_count = codon_types.count('CAA')
                CAG_count = codon_types.count('CAG')
                sbs_count = codon_types.count('sbs')
                mbs_count = codon_types.count('mbs')
    
                aac_percent = np.round((aac_count / total_count), 6)
                CAA_percent = np.round((CAA_count / total_count), 6)
                CAG_percent = np.round((CAG_count / total_count), 6)
                sbs_percent = np.round((sbs_count / total_count), 6)
                mbs_percent = np.round((mbs_count / total_count), 6)
                
                #print('sum of percent: ' + str(np.round((aac_percent + CAA_percent + CAG_percent + sbs_percent + mbs_percent), 5) ) )
                
                if (np.round((aac_percent + CAA_percent + CAG_percent + sbs_percent + mbs_percent), 4)) != 1.0:
                    print(type(np.round((aac_percent + CAA_percent + CAG_percent + sbs_percent + mbs_percent), 5)))
                    print(str(np.round((aac_percent + CAA_percent + CAG_percent + sbs_percent + mbs_percent), 5)))
                    print(type((1.0)))
                    print(str((1.0 )))
                    raise Exception('not 100%')
                
                csv_file.write(str(aac_percent) + ',' + str(CAA_percent) + ',' + str(CAG_percent) + ',' + str(sbs_percent) + ',' + str(mbs_percent) + ',')
                
                violin_plot_dict['codon_type'].append('aac')
                violin_plot_dict['value'].append(aac_percent)
                violin_plot_dict['label'].append(listed_label)
                violin_plot_dict['cc'].append(listed_cc_num)
                
                violin_plot_dict['codon_type'].append('CAA')
                violin_plot_dict['value'].append(CAA_percent)
                violin_plot_dict['label'].append(listed_label)
                violin_plot_dict['cc'].append(listed_cc_num)
                                
                violin_plot_dict['codon_type'].append('CAG')
                violin_plot_dict['value'].append(CAG_percent)
                violin_plot_dict['label'].append(listed_label)
                violin_plot_dict['cc'].append(listed_cc_num)
                
                violin_plot_dict['codon_type'].append('sbs')
                violin_plot_dict['value'].append(sbs_percent)
                violin_plot_dict['label'].append(listed_label)
                violin_plot_dict['cc'].append(listed_cc_num)
                
                violin_plot_dict['codon_type'].append('mbs')
                violin_plot_dict['value'].append(mbs_percent)
                violin_plot_dict['label'].append(listed_label)
                violin_plot_dict['cc'].append(listed_cc_num)
                
                
                if total_count < 9:
                    correlation_dict['quantile'].append(1)
                elif total_count < 14:
                    correlation_dict['quantile'].append(2)
                elif total_count < 19:
                    correlation_dict['quantile'].append(3)
                elif total_count < 24:
                    correlation_dict['quantile'].append(4)
                elif total_count < 29:
                    correlation_dict['quantile'].append(5)
                elif total_count < 34:
                    correlation_dict['quantile'].append(6)
                elif total_count < 39:
                    correlation_dict['quantile'].append(7)
                elif total_count < 44:
                    correlation_dict['quantile'].append(8)
                elif total_count >= 44:
                    correlation_dict['quantile'].append(9)
                else:
                    correlation_dict['quantile'].append(0)
    
                correlation_dict['label'].append(listed_label)
                correlation_dict['aac_count'].append(aac_count)
                correlation_dict['CAA_count'].append(CAA_count)
                correlation_dict['CAG_count'].append(CAG_count)
                correlation_dict['sbs_count'].append(sbs_count)
                correlation_dict['mbs_count'].append(mbs_count)
    
                correlation_dict['aac_percent'].append(aac_percent)
                correlation_dict['CAA_percent'].append(CAA_percent)
                correlation_dict['CAG_percent'].append(CAG_percent)
                correlation_dict['sbs_percent'].append(sbs_percent)
                correlation_dict['mbs_percent'].append(mbs_percent)
                
                correlation_dict['repeat_length'].append(total_count)
                correlation_dict['CAA_by_CAG'].append(CAA_count / CAG_count)
                correlation_dict['impurities_to_CAG'].append((sbs_count + mbs_count) / CAG_count)
                correlation_dict['impurities_and_CAA_to_CAG'].append((sbs_count + mbs_count + CAA_count) / CAG_count)
                correlation_dict['impurities_to_Q'].append((sbs_count + mbs_count) / (CAG_count + CAA_count))
                correlation_dict['sbs_to_CAG'].append((sbs_count) / CAG_count)
                correlation_dict['mbs_to_CAG'].append((mbs_count) / CAG_count)
    
            else:
                print(len(codon_types))
                raise Exception
            
            
        csv_file.write('\n')
        #print('id_df: ')
        #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        #    print(id_df)
        #print('###############################################################################')
    print('loop finished')
    csv_file.close()
    data = pd.read_csv('codon_types_per_label.csv', sep=',', na_values='NA', header=0)
    
    codon_types = (df['codon_type'][df.orf=='0']).values.tolist()
    #print('codon_types: ' + str(codon_types) + '\nset of codon_types: ' + str(set(codon_types)) + '\n')
    
    vp_df = pd.DataFrame.from_dict(violin_plot_dict)
    correlation_df = pd.DataFrame.from_dict(correlation_dict)
    vp_df.to_pickle('vp_df.pkl')
    correlation_df.to_pickle('correlation_df.pkl')
    
    
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #    print('vp_df')
    #    print(vp_df)
    
    nonan_vp_df = vp_df.dropna(axis=0)
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #    print('nonan_vp_df')
    #    print(nonan_vp_df)
    violin_and_box_plot_sample_size_and_median_dict = defaultdict(list)
    
    #make a violinplot
    sns.set(style="whitegrid", palette="pastel", color_codes=True)
    
    f, ax = plt.subplots(figsize=(20, 8))
    #ax.set(ylim=(0.0, 1.0))
    sns.violinplot(x="label", y="value", hue="codon_type", data=nonan_vp_df, linewidth = 0, palette={"aac": "grey", "CAA": "yellow", "CAG": "red", "sbs": "orange", "mbs": "green"}, order=label_list)
    sns.despine(left=True)
    
    f.suptitle('Comparing Codon Usage Bias', fontsize=18, fontweight='bold')
    ax.set_xlabel("genomes",size = 16,alpha=0.7)
    ax.set_ylabel("relative occurence [specific codon / total codons]",size = 16,alpha=0.7)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    plt.legend(bbox_to_anchor=(0.38, 1.05), loc='upper left', ncol=5)
    
    plt.savefig('violinplot.svg')
    plt.close()
    
    
    for listed_codon_type in ['aac','CAA','CAG','sbs','mbs']:
        
        shapiro = open(listed_codon_type + 'shapiro_results.csv', 'w')
        shapiro.close()
        shapiro = open(listed_codon_type + 'shapiro_results.csv', 'a')
        
        two_sided_comp_csv_file_stats = open(listed_codon_type + '_full_unfiltered_comparison_of_codon_types_per_label_two_sided_results_stats.csv', 'w')
        two_sided_comp_csv_file_stats.close()
        two_sided_comp_csv_file_stats = open(listed_codon_type + '_full_unfiltered_comparison_of_codon_types_per_label_two_sided_results_stats.csv', 'a')
        
        two_sided_comp_csv_file_pvalue = open(listed_codon_type + '_full_unfiltered_comparison_of_codon_types_per_label_two_sided_results_pvalue.csv', 'w')
        two_sided_comp_csv_file_pvalue.close()
        two_sided_comp_csv_file_pvalue = open(listed_codon_type + '_full_unfiltered_comparison_of_codon_types_per_label_two_sided_results_pvalue.csv', 'a')
        
        normal_distributions_flag = True
        equal_variances_flag = True
        
        stat_test_dict = dict()
    
        #shapiro test
        shapiro.write('label,shapiro_pvalue,shapiro_W,\n')
        for listed_label in label_list:
            column_array = nonan_vp_df.loc[(nonan_vp_df['codon_type'] == listed_codon_type) & (nonan_vp_df['label'] == listed_label)]['value'].values
            shapiro_W, shapiro_pvalue = stats.shapiro(column_array)
            shapiro.write(listed_label + ',' + str(np.round(shapiro_pvalue, 9)) + ',' + str(np.round(shapiro_W, 9)) + ',\n')
        shapiro.close()
        column_array_list = [nonan_vp_df.loc[(nonan_vp_df['codon_type'] == listed_codon_type) & (nonan_vp_df['label'] == listed_label)]['value'].values for listed_label in label_list]
        median_list = [np.median(nonan_vp_df.loc[(nonan_vp_df['codon_type'] == listed_codon_type) & (nonan_vp_df['label'] == listed_label)]['value'].values) for listed_label in label_list]
        sample_size_list = [nonan_vp_df.loc[(nonan_vp_df['codon_type'] == listed_codon_type) & (nonan_vp_df['label'] == listed_label)]['value'].shape[0] for listed_label in label_list]
        
        violin_and_box_plot_sample_size_and_median_dict[listed_codon_type + '|median_list'] = median_list
        violin_and_box_plot_sample_size_and_median_dict[listed_codon_type + '|sample_size_list'] = sample_size_list
        
        for column_array in column_array_list:
            #test for normal distribution
            try:
                shapiro_W, shapiro_pvalue = stats.shapiro(column_array)
                #print(listed_label + ' ' + listed_codon_type + ' W: ' + str(shapiro_W))
                #print(listed_label + ' ' + listed_codon_type + ' p: ' + str(shapiro_pvalue))
                if shapiro_pvalue < 0.05:
                    #print('maybe\n')
                    pass
                else:
                    #print('nope\n')
                    normal_distributions_flag = False
            except ValueError:
                pass
        
        #test for homogenic variance
        levene_statistic, levene_pvalue = stats.levene(column_array_list[0],column_array_list[1],column_array_list[2],column_array_list[3],column_array_list[4],column_array_list[5],column_array_list[6],column_array_list[7],column_array_list[8],column_array_list[9], center='median')
        print(listed_codon_type)
        print('levene_statistic: ' + str(np.round(levene_statistic, 4)))
        print('levene_pvalue: ' + str(np.round(levene_pvalue, 4)))
        if levene_pvalue < 0.05:
            equal_variances_flag = False
        
        #two sided comparisons
        old_fig1_list = []
        alpha = 0.05
        two_sided_comp_csv_file_stats.write(',')
        two_sided_comp_csv_file_pvalue.write(',')
        for fig1 in range(0,len(label_list),1):
            two_sided_comp_csv_file_stats.write(label_list[fig1] + ',')
            two_sided_comp_csv_file_pvalue.write(label_list[fig1] + ',')    
        two_sided_comp_csv_file_stats.write('\n')
        two_sided_comp_csv_file_pvalue.write('\n')
        for fig1 in range(0,len(label_list),1):
            two_sided_comp_csv_file_stats.write(label_list[fig1] + ',')
            two_sided_comp_csv_file_pvalue.write(label_list[fig1] + ',')
            for fig2 in range(0,len(label_list),1):
                if fig1 != fig2 and fig2 not in old_fig1_list:
                    sample1 = column_array_list[fig1]
                    sample2 = column_array_list[fig2]
                    statistics = np.nan
                    pvalue = np.nan
                    if normal_distributions_flag == True and equal_variances_flag == True:
                        #print('performing independent two sample t-test between ' + label_list[fig1] + ' and ' + label_list[fig2])    
                        try:
                            statistics, pvalue = stats.ttest_ind(sample1, sample2, equal_var=True, nan_policy='omit')
                        #print('statistics: ' + str(statistics))
                        #print('pvalue: ' + str(pvalue))
                        #if pvalue < (alpha / len(label_list)):
                            #print('significant difference: ' + str(pvalue) + ' < ' + str(alpha / len(label_list)))
                            #print(listed_codon_type + ' --- ' + label_list[fig1] + ' and ' + label_list[fig2] + '!!!\n')
                        #else:
                            #print('no significant difference: ' + str(pvalue) + ' < ' + str(alpha / len(label_list)) + '\n')
                        except ValueError:
                            pass
                    elif equal_variances_flag == False and normal_distributions_flag == True:
                        #print("performing two-sided Welch's -test between " + label_list[fig1] + ' and ' + label_list[fig2])
                        try:
                            statistics, pvalue = stats.ttest_ind(sample1, sample2, equal_var=False, nan_policy='omit')
                        #print('statistics: ' + str(statistics))
                        #print('pvalue: ' + str(pvalue))
                        #if pvalue < (alpha / len(label_list)):
                            #print('significant difference: ' + str(pvalue) + ' < ' + str(alpha / len(label_list)))
                            #print(listed_codon_type + '  --- ' + label_list[fig1] + ' and ' + label_list[fig2] + '!!!\n')
                        #else:
                            #print('no significant difference: ' + str(pvalue) + ' < ' + str(alpha / len(label_list)) + '\n')
                        except ValueError:
                            pass
                    elif normal_distributions_flag == False:
                        #print('performing two-sided Mann-Whitney rank test between ' + label_list[fig1] + ' and ' + label_list[fig2])
                        try:
                            statistics, pvalue = stats.mannwhitneyu(sample1, sample2)
                        #print('statistics: ' + str(statistics))
                        #print('pvalue: ' + str(pvalue))
                        #if pvalue < (alpha / len(label_list)):
                            #print('significant difference: ' + str(pvalue) + ' < ' + str(alpha / len(label_list)))
                            #print(listed_codon_type + '  --- ' + label_list[fig1] + ' and ' + label_list[fig2] + '!!!\n')
                        #else:
                            #print('no significant difference: ' + str(pvalue) + ' < ' + str(alpha / len(label_list)) + '\n')
                        except ValueError:
                            pass    
                    two_sided_comp_csv_file_stats.write(str(np.round(statistics, 4)) + ',')
                    two_sided_comp_csv_file_pvalue.write(str(np.round(pvalue, 4)) + ',')
                else:
                    two_sided_comp_csv_file_stats.write(',')
                    two_sided_comp_csv_file_pvalue.write(',')
            two_sided_comp_csv_file_stats.write('\n')
            two_sided_comp_csv_file_pvalue.write('\n')
            old_fig1_list.append(fig1)
        two_sided_comp_csv_file_stats.close()
        two_sided_comp_csv_file_pvalue.close()
    
    #make the boxplot from Figure 2
    sns.set(style="whitegrid", palette="pastel", color_codes=True)
    
    f, ax = plt.subplots(figsize=(20, 8))
    ax.set(ylim=(0.0, 1.0))
    sns.boxplot(x="label", y="value", hue="codon_type", data=nonan_vp_df, palette={"aac": "grey", "CAA": "red", "CAG": "yellow", "sbs": "green", "mbs": "blue"}, order=label_list, notch=True, )
    sns.despine(left=True)
    
    f.suptitle('Comparing Codon Usage Bias', fontsize=18, fontweight='bold')
    ax.set_xlabel("genomes",size = 16,alpha=0.7)
    ax.xaxis.set_label_coords(0.5, -0.075)
    ax.set_ylabel("relative occurrence [specific codon / total codons]",size = 16,alpha=0.7)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    
    median_labels = [str(np.round(s, 2)) for s in median_list]
    pos = range(len(label_list))
    
    for listed_codon_type, offset in [('aac', -0.327),('CAA', -0.1635),('CAG', 0),('sbs', 0.1635),('mbs', 0.327)]:
        
        median_list = violin_and_box_plot_sample_size_and_median_dict[listed_codon_type + '|median_list']    
        median_labels = [str(np.round(s, 3)) for s in median_list]
        sample_size_labels = violin_and_box_plot_sample_size_and_median_dict[listed_codon_type + '|sample_size_list']
    
        for tick,label in zip(pos,ax.get_xticklabels()):
            ax.text(pos[tick], -0.065, 'n=' + sample_size_labels[tick], horizontalalignment='center', size='xsmall', color='b', weight='semibold')
    
    plt.legend(bbox_to_anchor=(0.32, 1.08), loc='upper left', ncol=5)
    
    plt.savefig('boxplot.svg')
    plt.close()
    
    #these heatmaps were of little use in our data analysis, but they might come in handy in the future - so I will leave the code here
    #heatmap_dict = defaultdict(list)
    #
    #overall_max_repeat_codons = df['repeat_length'].max()
    #overall_min_repeat_codons = df['repeat_length'].min()
    #
    #for listed_codon_type in ['aac','CAA','CAG','sbs','mbs']:
    #    
    #    for listed_label in label_list:
    #
    #        global_max = df['repeat_length'][df.label==listed_label].max()
    #        global_min = df['repeat_length'][df.label==listed_label].min()
    #        
    #        for column in range(int(global_min), int(global_max), 1):
    #            COI_rows_at_column = df.loc[(df['in_repeat'] == True) & (df['codon_type'] == listed_codon_type) & (df['orf'] == '0') & (df['label'] == listed_label) & (df['repeat_length'] == column)]
    #            number_of_COI_at_column = COI_rows_at_column.shape[0]
    #            all_seqs_at_column = df.loc[(df['in_repeat'] == True) & (df['codon_type'] != '---') & (df['orf'] == '0') & (df['label'] == listed_label) & (df['repeat_length'] == column)]
    #            all_seqs = df.loc[(df['in_repeat'] == True) & (df['codon_type'] != '---') & (df['orf'] == '0') & (df['label'] == listed_label)]
    #            total_number_of_seqs_at_column = all_seqs_at_column.drop_duplicates(subset=['running_repeat_id'], keep='first').shape[0]
    #            total_number_of_seqs = all_seqs.drop_duplicates(subset=['running_repeat_id'], keep='first').shape[0]
    #            number_of_all_codon_types_in_the_plot = df.loc[(df['in_repeat'] == True) & (df['codon_type'] == listed_codon_type) & (df['orf'] == '0') & (df['label'] == listed_label)].shape[0]
    #            
    #            for rank_step in range(0,(int(global_max)),1):
    #                COI_rows_at_column_and_rank_step = COI_rows_at_column.loc[(COI_rows_at_column['in_repeat_codon_rank'] == rank_step)]
    #                number_of_COI_rows_at_column_and_rank_step = COI_rows_at_column_and_rank_step.shape[0]
    #                if (rank_step + 1) > column or total_number_of_seqs_at_column == 0:
    #                    heatmap_dict['rel_amount_of_codon_types_at_columnand_rank_step'].append(0)
    #                    heatmap_dict['rel_amount_of_codon_types_at_columnand_rank_step_times_log_of_n_columns'].append(0)
    #                    heatmap_dict['rel_amount_of_codon_types_at_columnand_rank_step_divided_by_total_plot_codon_types'].append(0)
    #                    heatmap_dict['total_amount_of_codon_types_at_columnand_rank_step'].append(0)
    #                    heatmap_dict['should_value_be_displayed'].append(True)
    #                    heatmap_dict['codon_type'].append(listed_codon_type)
    #                    heatmap_dict['label'].append(listed_label)
    #                    heatmap_dict['column'].append(column)                
    #                    heatmap_dict['rank_step'].append(rank_step + 1)
    #                    heatmap_dict['column_sample_size'].append(total_number_of_seqs_at_column)
    #                else:
    #                    heatmap_dict['rel_amount_of_codon_types_at_columnand_rank_step'].append((number_of_COI_rows_at_column_and_rank_step / total_number_of_seqs_at_column) * 100)
    #                    heatmap_dict['rel_amount_of_codon_types_at_columnand_rank_step_times_log_of_n_columns'].append((number_of_COI_rows_at_column_and_rank_step / total_number_of_seqs_at_column) * math.log(total_number_of_seqs_at_column))
    #                    heatmap_dict['total_amount_of_codon_types_at_columnand_rank_step'].append((number_of_COI_rows_at_column_and_rank_step))
    #                    if number_of_all_codon_types_in_the_plot != 0:
    #                        heatmap_dict['rel_amount_of_codon_types_at_columnand_rank_step_divided_by_total_plot_codon_types'].append((number_of_COI_rows_at_column_and_rank_step / number_of_all_codon_types_in_the_plot) * 100)
    #                    else:
    #                        heatmap_dict['rel_amount_of_codon_types_at_columnand_rank_step_divided_by_total_plot_codon_types'].append(0)
    #                    heatmap_dict['should_value_be_displayed'].append(False)
    #                    heatmap_dict['codon_type'].append(listed_codon_type)
    #                    heatmap_dict['label'].append(listed_label)
    #                    heatmap_dict['column'].append(column)                
    #                    heatmap_dict['rank_step'].append(rank_step + 1)
    #                    heatmap_dict['column_sample_size'].append(total_number_of_seqs_at_column)
    #
    #pd.set_option('display.float_format', lambda x: '%.3f' % x)
    #
    #heatmap_df = pd.DataFrame.from_dict(heatmap_dict)
    #
    #for listed_codon_type in ['aac','CAA','CAG','sbs','mbs']:
    #    
    #    for listed_label in label_list:
    #        sub_heatmap_df = heatmap_df.loc[(heatmap_df['codon_type'] == listed_codon_type) & (heatmap_df['label'] == listed_label)]
    #        sub_heatmap_df_column_max = sub_heatmap_df['column'].max()
    #        sub_heatmap_df_rank_step_max = sub_heatmap_df['rank_step'].max()
    #        
    #        ready_to_plot = sub_heatmap_df.pivot_table(index='rank_step',columns='column',values='rel_amount_of_codon_types_at_columnand_rank_step')
    #        show_annot_array = sub_heatmap_df.pivot_table(index='rank_step',columns='column',values='should_value_be_displayed')
    #        position_N_array = sub_heatmap_df.pivot_table(index='rank_step',columns='column',values='column_sample_size')
    #        f, ax = plt.subplots(figsize=(10, 8))
    #        hmfig = sns.heatmap(ready_to_plot, xticklabels=1, yticklabels=1, center=0, linewidths=0.5, cmap="Spectral_r", mask=show_annot_array, fmt="d", annot=False)#annot=ready_to_plot, )
    #
    #        plt.title('Heatmap of ' + listed_codon_type + ' codon_type Frequency [%]', fontsize = 20, y=1.04) # title with fontsize 20
    #        plt.xlabel('length of detected repeats [codon_types]', fontsize = 15, verticalalignment="top")
    #        plt.ylabel('ranked position in sequence', fontsize = 15)
    #        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize = 7)
    #        ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize =  7)#, verticalalignment="top")
    #        
    #        sample_size_labels = position_N_array.iloc[0].tolist()
    #        pos = range(len(sample_size_labels))
    #        for tick,label in zip(pos,ax.get_xticklabels()):
    #            ax.text(pos[tick], 0, 'n=' + str(sample_size_labels[tick]), horizontalalignment='center', color='b', weight='semibold', fontsize = 2)#, size='small')
    #        
    #        plt.savefig('heatmap_' + listed_codon_type + '_' + listed_label + '_relative_to_column.svg')
    #        plt.close()
    #        
    #        
    #        #relative to plot
    #        ready_to_plot = sub_heatmap_df.pivot_table(index='rank_step',columns='column',values='rel_amount_of_codon_types_at_columnand_rank_step_divided_by_total_plot_codon_types')
    #        f, ax = plt.subplots(figsize=(10, 8))
    #        hmfig = sns.heatmap(ready_to_plot, xticklabels=1, yticklabels=1, center=0, linewidths=0.5, cmap="Spectral_r", mask=show_annot_array, fmt="d", annot=False)#annot=ready_to_plot, )
    #        plt.title('Heatmap of ' + listed_codon_type + ' codon_type Frequency [%]', fontsize = 20, y=1.04) # title with fontsize 20
    #        plt.xlabel('length of detected repeats [codon_types]', fontsize = 15, verticalalignment="top")
    #        plt.ylabel('ranked position in sequence', fontsize = 15)
    #        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize = 7)
    #        ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize =  7)#, verticalalignment="top")
    #        sample_size_labels = position_N_array.iloc[0].tolist()
    #        pos = range(len(sample_size_labels))
    #        for tick,label in zip(pos,ax.get_xticklabels()):
    #            ax.text(pos[tick], 0, 'n=' + str(sample_size_labels[tick]), horizontalalignment='center', color='b', weight='semibold', fontsize = 2)#, size='small')
    #        plt.savefig('heatmap_' + listed_codon_type + '_' + listed_label + '_relative_to_plot.svg')
    #        plt.close() 
    #            
    #        #relative to plot adjusted for number of sequences in column
    #        ready_to_plot = sub_heatmap_df.pivot_table(index='rank_step',columns='column',values='rel_amount_of_codon_types_at_columnand_rank_step_times_log_of_n_columns')
    #        f, ax = plt.subplots(figsize=(10, 8))
    #        hmfig = sns.heatmap(ready_to_plot, xticklabels=1, yticklabels=1, center=0, linewidths=0.5, cmap="Spectral_r", mask=show_annot_array, fmt="d", annot=False)#annot=ready_to_plot, )
    #        plt.title('Heatmap of ' + listed_codon_type + ' codon_type Frequency times logarithm of sample size [% * log(cn)]', fontsize = 20, y=1.04) # title with fontsize 20
    #        plt.xlabel('length of detected repeats [codon_types]', fontsize = 15, verticalalignment="top")
    #        plt.ylabel('ranked position in sequence', fontsize = 15)
    #        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize = 7)
    #        ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize =  7)#, verticalalignment="top")
    #        sample_size_labels = position_N_array.iloc[0].tolist()
    #        pos = range(len(sample_size_labels))
    #        for tick,label in zip(pos,ax.get_xticklabels()):
    #            ax.text(pos[tick], 0, 'n=' + str(sample_size_labels[tick]), horizontalalignment='center', color='b', weight='semibold', fontsize = 2)#, size='small')
    #        plt.savefig('heatmap_' + listed_codon_type + '_' + listed_label + '_relative_to_plot_adjusted_for_column_n.svg')
    #        plt.close()        
    #        
    #        #total amount per position
    #        ready_to_plot = sub_heatmap_df.pivot_table(index='rank_step',columns='column',values='total_amount_of_codon_types_at_columnand_rank_step')
    #        f, ax = plt.subplots(figsize=(10, 8))
    #        hmfig = sns.heatmap(ready_to_plot, xticklabels=1, yticklabels=1, center=0, linewidths=0.5, cmap="Spectral_r", mask=show_annot_array, fmt="d", annot=False)#annot=ready_to_plot, )
    #        plt.title('Heatmap of ' + listed_codon_type + ' total amount of codon_types per position [codon_types]', fontsize = 20, y=1.04) # title with fontsize 20
    #        plt.xlabel('length of detected repeats [codon_types]', fontsize = 15, verticalalignment="top")
    #        plt.ylabel('ranked position in sequence', fontsize = 15)
    #        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize = 7)
    #        ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize =  7)#, verticalalignment="top")        
    #        sample_size_labels = position_N_array.iloc[0].tolist()
    #        pos = range(len(sample_size_labels))
    #        for tick,label in zip(pos,ax.get_xticklabels()):
    #            ax.text(pos[tick], 0, 'n=' + str(sample_size_labels[tick]), horizontalalignment='center', color='b', weight='semibold', fontsize = 2)#, size='small')        
    #        plt.savefig('heatmap_' + listed_codon_type + '_' + listed_label + '_total.svg')
    #        plt.close()                
    #        
    #for listed_codon_type in ['aac','CAA','CAG','sbs','mbs']:
    #    
    #    sub_heatmap_df = heatmap_df.loc[(heatmap_df['codon_type'] == listed_codon_type)]
    #    sub_heatmap_df_column_max = sub_heatmap_df['column'].max()
    #    sub_heatmap_df_rank_step_max = sub_heatmap_df['rank_step'].max()        
    #    ready_to_plot = sub_heatmap_df.pivot_table(index='rank_step',columns='column',values='rel_amount_of_codon_types_at_columnand_rank_step')    
    #    f, ax = plt.subplots(figsize=(10, 8))
    #    hmfig = sns.heatmap(ready_to_plot, xticklabels=1, yticklabels=1, center=0, linewidths=0.5, cmap="Spectral_r", fmt="d", annot=False)#annot=ready_to_plot,  mask=show_annot_array,)
    #    plt.title('Heatmap of ' + listed_codon_type + ' codon_type Frequency [%]', fontsize = 20, y=1.04) # title with fontsize 20
    #    plt.xlabel('length of detected repeats [codon_types]', fontsize = 15, verticalalignment="top")
    #    plt.ylabel('ranked position in sequence', fontsize = 15)
    #    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize = 7)
    #    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize =  7)#, verticalalignment="top")        
    #    sample_size_labels = position_N_array.iloc[0].tolist()
    #    pos = range(len(sample_size_labels))
    #    for tick,label in zip(pos,ax.get_xticklabels()):
    #        ax.text(pos[tick], 0, 'n=' + str(sample_size_labels[tick]), horizontalalignment='center', color='b', weight='semibold', fontsize = 4)#, size='small')
    #    plt.savefig('heatmap_' + listed_codon_type + '.svg')
    #
    #plt.close()
    
    loaded_quantile_df = pd.read_pickle('hundred_quantile_df.pkl')
    
    #make the plots from Figure 3
    for listed_label in ['BIC','HDE', 'KK1', 'LOS', 'MOT', 'NE1', 'REF', 'SG3', 'WC1', 'ZVE']:
        for listed_codon in ['CAA','sbs','mbs']:
            
            #if listed_label != 'BIC' or listed_codon != 'CAA':
            #    continue
            quantilized_codon_df = loaded_quantile_df.loc[(loaded_quantile_df['codon_type'] == listed_codon) & (loaded_quantile_df['repeat_length'] <= 25) & (loaded_quantile_df['label'] == listed_label)]
            
            if listed_codon == 'CAA':
                codon_color = 'Reds'
            elif listed_codon == 'sbs':
                codon_color = 'Greens'
            elif listed_codon == 'mbs':
                codon_color = 'Blues'
    
            #Bivariate KDE can only use gaussian kernel
            mfig, ax = plt.subplots(figsize=(20, 8))
            ax.set(xlim=(1, 28))
            ax.set(ylim=(0, 1000))
            
            ax = sns.kdeplot(quantilized_codon_df["repeat_length"], quantilized_codon_df["full_percent_quantile"], cmap=codon_color, shade=True, shade_lowest=False, cbar=True, n_levels=20, gridsize=200, bw='silverman')
            mfig.savefig(listed_label + '_' + listed_codon + '_quantile_cbar_gausssilverman_kde_200_grid_20_levels_.svg')
            mfig.clf()
            plt.close()




###############################################################################
#end of major module-functions callable via cmd-agrs
###############################################################################

if __name__ == '__main__':
    ### setting variables ###

    mode = None

    directory = None
    bed = None
    bam = None
    ref = None
    cmd_gff3 = None

    #TODO: fix this chromosome list >_<
    cmd_chr = None
    cmd_tens_pre_dec = None
    cmd_ones_pre_dec = None
    introns = 'NA'
    repeat_positions = 'NA'
    
    for position in xrange(len(sys.argv)):
        if sys.argv[position] == '-m':#mode
            mode = sys.argv[position + 1]
        elif sys.argv[position] == '-d':#directory
            directory = sys.argv[position + 1]
        elif sys.argv[position] == '-b':#bed file
            bed = sys.argv[position + 1]
        elif sys.argv[position] == '-a':#alignments
            bam = sys.argv[position + 1]
        elif sys.argv[position] == '-r':#reference
            ref = sys.argv[position + 1]
        elif sys.argv[position] == '-f':#features
            cmd_gff3 = sys.argv[position + 1]
        elif sys.argv[position] == '-c':#chromosome
            cmd_chr = sys.argv[position + 1]
        elif sys.argv[position] == '-t':#tens predecimal
            cmd_tens_pre_dec = sys.argv[position + 1]
        elif sys.argv[position] == '-o':#ones predecimal
            cmd_ones_pre_dec = sys.argv[position + 1]
        elif sys.argv[position] == '-i':#introns
            exclude_introns = sys.argv[position + 1]
        elif sys.argv[position] == '-p':#positions
            positions = sys.argv[position + 1]
        elif sys.argv[position] == '-w':#write
            overwrite = sys.argv[position + 1]
        elif sys.argv[position] == '-a':#write
            aligner = sys.argv[position + 1]
            
        elif sys.argv[position].startswith('-mode=') or sys.argv[position].startswith('--mode=') or sys.argv[position].startswith('-mode==') or sys.argv[position].startswith('--mode=='):
            mode = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-directory=') or sys.argv[position].startswith('--directory=') or sys.argv[position].startswith('-directory==') or sys.argv[position].startswith('--directory=='):
            directory = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-bed=') or sys.argv[position].startswith('--bed=') or sys.argv[position].startswith('-bed==') or sys.argv[position].startswith('--bed=='):
            bed = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-bam=') or sys.argv[position].startswith('--bam=') or sys.argv[position].startswith('-bam==') or sys.argv[position].startswith('--bam=='):
            bam = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-ref=') or sys.argv[position].startswith('--ref=') or sys.argv[position].startswith('-ref==') or sys.argv[position].startswith('--ref=='):
            ref = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-gff3=') or sys.argv[position].startswith('--gff3=') or sys.argv[position].startswith('-gff3==') or sys.argv[position].startswith('--gff3=='):
            cmd_gff3 = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-chr=') or sys.argv[position].startswith('--chr=') or sys.argv[position].startswith('-chr==') or sys.argv[position].startswith('--chr=='):
            cmd_chr = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-tens_pre_dec=') or sys.argv[position].startswith('--tens_pre_dec=') or sys.argv[position].startswith('-tens_pre_dec==') or sys.argv[position].startswith('--tens_pre_dec=='):
            cmd_tens_pre_dec = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-ones_pre_dec=') or sys.argv[position].startswith('--ones_pre_dec=') or sys.argv[position].startswith('-ones_pre_dec==') or sys.argv[position].startswith('--ones_pre_dec=='):
            cmd_ones_pre_dec = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-introns=') or sys.argv[position].startswith('--introns=') or sys.argv[position].startswith('-introns==') or sys.argv[position].startswith('--introns=='):
            introns = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-positions=') or sys.argv[position].startswith('--positions=') or sys.argv[position].startswith('-positions==') or sys.argv[position].startswith('--positions=='):
            positions = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-overwrite=') or sys.argv[position].startswith('--overwrite=') or sys.argv[position].startswith('-overwrite==') or sys.argv[position].startswith('--overwrite=='):
            overwrite = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
        elif sys.argv[position].startswith('-aligner=') or sys.argv[position].startswith('--aligner=') or sys.argv[position].startswith('-aligner==') or sys.argv[position].startswith('--aligner=='):
            aligner = (sys.argv[position].split('='))[-1].replace('"', '').replace("'", "")
    
    argv_list = [('mode', mode), ('directory', directory), ('bed', bed), ('bam', bam), ('ref', ref), ('cmd_gff3', cmd_gff3), ('cmd_chr', cmd_chr), ('cmd_tens_pre_dec', cmd_tens_pre_dec), ('cmd_ones_pre_dec', cmd_ones_pre_dec), ('introns', introns), ('positions', positions), ('overwrite', overwrite), ('aligner', aligner)]

    
    for name, value in argv_list:
        print(name + ' : ' + str(value))
    
    for name, value in argv_list:
        if value == None:
            print('terminal argument "' +  name + '" = None')
            raise Exception('Error: terminal argument "' +  name + '" = None')
        elif value == '':
            raise Exception('Error: terminal argument "' +  name + '" = "" (empty string)')
        elif type(value) == str and len(value) > 0:
            print('terminal argument "' +  name + '" = None')
        else:
            exception_string = ('Error: uncharacterized terminal argument "' +  name)
            raise Exception(exception_string)
            
    #GFF3_path = directory + cmd_gff3
    bed_file_origin = directory + bed
    bam_file_origin = directory + bam
    
    t1 = time.time()
        
    if mode == 'HRA_bed_parser':
        
        ChrOI_list = []
        
        if ':' in cmd_chr:
            ChrOI_list = cmd_chr.split(':')
        elif '|' in cmd_chr:
            ChrOI_list = cmd_chr.split('|')
        elif ',' in cmd_chr:
            ChrOI_list = cmd_chr.split(',')
        elif ';' in cmd_chr:
            ChrOI_list = cmd_chr.split(';')
        elif '-' in cmd_chr:
            ChrOI_list = cmd_chr.split('-')
        else:
            raise Exception('HRA_bed_parser(ERROR): no list of chromosomes of interest submitted for processing, please submit a list of chromosomes seperated homogenously by one of the following symbols without whitespaces:\n":"    ","    ";"    "|"    "-"')

        HRA_bed_parser(working_directory = directory, bed_name = bed, ChrOI_list = ChrOI_list, overwrite_old_files = overwrite)
        t2 = time.time()
        print('Time spent on HRA_bed_parser: ' + str(t2 - t1))
        
    elif mode == 'HRA_bam_parser':
        HRA_bam_parser(working_directory = directory, bam_name = bam, ChrOI = cmd_chr, overwrite_indexes = False)
        t2 = time.time()
        print('Time spent on HRA_bam_parser: ' + str(t2 - t1))

    elif mode == 'HRA_mapping_analyzer':
        HRA_mapping_analyzer(working_directory = directory, bed_name = bed, bam_name = bam, ChrOI = cmd_chr, tens_pre_dec = cmd_tens_pre_dec, ones_pre_dec = cmd_ones_pre_dec, GFF3_path = cmd_gff3, introns = introns)
        t2 = time.time()
        print('Total time spent on HRA_mapping_analyzer(): ' + str(t2 - t1))
        
    elif mode == 'HRA_ref_analyzer':
        HRA_ref_analyzer(working_directory = directory, bed_name = bed, ChrOI = cmd_chr, ones_pre_dec = cmd_ones_pre_dec, ref_name = ref, GFF3_path = cmd_gff3)
        t2 = time.time()
        print(get_mem('HRA_ref_analyzer'))

        print('Total time spent on HRA_ref_analyzer(): ' + str(t2 - t1))

    elif mode == 'HRA_merge_outfiles':
        
        ChrOI_list = []
        
        if ':' in cmd_chr:
            ChrOI_list = cmd_chr.split(':')
        elif '|' in cmd_chr:
            ChrOI_list = cmd_chr.split('|')
        elif ',' in cmd_chr:
            ChrOI_list = cmd_chr.split(',')
        elif ';' in cmd_chr:
            ChrOI_list = cmd_chr.split(';')
        elif '-' in cmd_chr:
            ChrOI_list = cmd_chr.split('-')
        else:
            raise Exception('HRA_merge_outfiles(ERROR): no list of chromosomes of interest submitted for processing, please submit a list of chromosomes seperated homogenously by one of the following symbols without whitespaces:\n":"    ","    ";"    "|"    "-"')
                
        ### launching clean_up_after_HRA_mapping_analyzer() ###    
        HRA_merge_outfiles(directory, bam, ref, ChrOI_list, overwrite)
        t2 = time.time()
        print('Total time spent on HRA_merge_outfiles(): ' + str(t2 - t1))

    elif mode == 'HRA_compare_QRY_and_REF':
        HRA_compare_QRY_and_REF(working_directory = directory, bam_name = bam, bed_name = bed, ref_name = ref, ChrOI = cmd_chr, GFF3_path = cmd_gff3, aligner = 'NA', positions = positions)
        t2 = time.time()
        print(get_mem('HRA_compare_QRY_and_REF'))
        
        print('Total time spent on HRA_compare_QRY_and_REF(): ' + str(t2 - t1))

    elif mode == 'HRA_trans_intron_polyQ_analyzer':
        
        if ':' in bam:
            experiment_tuple = bam.split(':')
        elif '|' in bam:
            experiment_tuple = bam.split('|')
        elif ',' in bam:
            experiment_tuple = bam.split(',')
        elif ';' in bam:
            experiment_tuple = bam.split(';')
        elif '-' in bam:
            experiment_tuple = bam.split('-')
        else:
            raise Exception('HRA_trans_intron_polyQ_analyzer(ERROR): no list of bam file acronyms of interest plus bam_file name submitted for processing, please submit a the name of an experiment folder of your choice and a contained bam_file (you can chose randomly as only the reference polyQs will be retrieved). Both names ought to be separated by one of the following symbols without whitespaces:\n":"    ","    ";"    "|"    "-", e.g.: "EXP:my_experiment.bam"')

        
        HRA_trans_intron_polyQ_analyzer(working_directory = directory, gff3 = cmd_gff3, experiment_folder = experiment_tuple[0], bam_name = experiment_tuple[1])
        t2 = time.time()
        print('Total time spent on HRA_trans_intron_polyQ_analyzer(): ' + str(t2 - t1))

    elif mode == 'HRA_compare_all_QRYs_and_REF':
        
        experiment_folders = []
        
        if ':' in bam:
            experiment_folders = bam.split(':')
        elif '|' in bam:
            experiment_folders = bam.split('|')
        elif ',' in bam:
            experiment_folders = bam.split(',')
        elif ';' in bam:
            experiment_folders = bam.split(';')
        elif '-' in bam:
            experiment_folders = bam.split('-')
        else:
            raise Exception('HRA_compare_all_QRYs_and_REF(ERROR): no list of bam file acronyms of interest submitted for processing, please submit a list of folders containing each bam files Ancient_Repeats_Tools experiment, seperated homogenously by one of the following symbols without whitespaces:\n":"    ","    ";"    "|"    "-"')

        ChrOI_list = []
        
        if ':' in cmd_chr:
            ChrOI_list = cmd_chr.split(':')
        elif '|' in cmd_chr:
            ChrOI_list = cmd_chr.split('|')
        elif ',' in cmd_chr:
            ChrOI_list = cmd_chr.split(',')
        elif ';' in cmd_chr:
            ChrOI_list = cmd_chr.split(';')
        elif '-' in cmd_chr:
            ChrOI_list = cmd_chr.split('-')
        else:
            raise Exception('HRA_compare_all_QRYs_and_REF(ERROR): no list of chromosomes of interest submitted for processing, please submit a list of chromosomes seperated homogenously by one of the following symbols without whitespaces:\n":"    ","    ";"    "|"    "-"')
        
        HRA_compare_all_QRYs_and_REF(working_directory = directory, bed_name = bed, experiment_folders = experiment_folders, REF = ref, gff3 = cmd_gff3, ChrOI_list = ChrOI_list, positions = positions, overwrite_old_files = overwrite, aligner = aligner)
        t2 = time.time()
        print(get_mem('HRA_compare_all_QRYs_and_REF'))
        
        print('Total time spent on HRA_compare_all_QRYs_and_REF(): ' + str(t2 - t1))

    elif mode == 'HRA_MSA_analyzer':
        HRA_MSA_analyzer(bed_path = directory + '/' + bed, HRA_result_path = directory)
        t2 = time.time()
        print('Total time spent on HRA_MSA_analyzer(): ' + str(t2 - t1))

    else:
        "Choose one of the following main functions you want to execute by submitting them as 5th sys.argv: 'HRA_bed_parser', 'HRA_bam_parser', 'HRA_mapping_analyzer', 'HRA_merge_outfiles', 'HRA_ref_analyzer', 'HRA_compare_QRY_and_REF', 'HRA_trans_intron_polyQ_analyzer', 'HRA_compare_all_QRYs_and_REF', 'HRA_MSA_analyzer'."
