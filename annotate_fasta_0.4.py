#!/usr/bin/env python3
# Given a fasta file and the corresponding blast output this script annotates the fasta with the hit resulting in the longest alignment and covering the beginning of the query
#     Copyright (C) 2018  Jeremias N. Brand

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse
from Bio import SeqIO
from collections import defaultdict
import sys
import os
import subprocess
import fileinput


class BlastRes():
    def __init__(self, qseqid, sseqid, pident, qlen, length, mismatch, 
        gapopen, qstart, qend, sstart, send, evalue, bitscore):
        """
        This class stores the results of a blast search in 
        -outfmt "6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore"
        """

        self.qseqid   = qseqid 
        self.sseqid   = sseqid 
        self.pident   = pident 
        self.qlen     = qlen 
        self.length   = length 
        self.mismatch = mismatch
        self.gapopen  = gapopen 
        self.qstart   = qstart 
        self.qend     = qend  
        self.sstart   = sstart 
        self.send     = send  
        self.evalue   = evalue 
        self.bitscore = bitscore
        
        
def check_hit_position(BlastRes_obj, current_coverage, max_query_start = 1, max_subject_start = 9, max_diff = 0.10, blast_type = 'blastx'):
    """
    Checks if a hit stored in a BlastRes class object meets the hit criteria we give
    returns a boolean, False if the conditions are not met

    """
    assert blast_type in ['blastn', 'blastx'], "Only blastn and blastx is supported."
    res = False
<<<<<<< HEAD
    
    if (int(BlastRes_obj.qstart) <= max_query_start and  int(BlastRes_obj.sstart) <= max_subject_start):
        if blast_type == 'blastn':
            hit_diff = 1-(int(BlastRes_obj.length)/int(BlastRes_obj.qlen))
        else:
            hit_diff = 1-((int(BlastRes_obj.length)*3)/int(BlastRes_obj.qlen))
        if (abs(hit_diff) < max_diff):
            res = True

=======
    if blast_type == 'blastn':
        if (int(BlastRes_obj.qstart) <= max_query_start and  int(BlastRes_obj.sstart) <= max_subject_start):
            hit_diff = 1-(int(BlastRes_obj.length)/int(BlastRes_obj.qlen))
            if (abs(hit_diff) < max_diff):
                res = True
    if blast_type == 'blastx':
        if (int(BlastRes_obj.qstart) <= max_query_start and  int(BlastRes_obj.sstart) <= max_subject_start):
            hit_diff = 1-((int(BlastRes_obj.length)*3)/int(BlastRes_obj.qlen))
            if (abs(hit_diff) < max_diff):
                res = True
>>>>>>> parent of 9e4542c... Cleaned up. Made hit_diff more intuitive.
    return res


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta", type=str, help="Fasta file to annotate")
    parser.add_argument("blast_output", type=str, help="Blast output of input_fasta")
    parser.add_argument("blast_type", type=str, help="Type of blast conducted to produce output: Currently blastn or blastx")
    parser.add_argument("output_file", type=str, help="Name of the outfile")
    parser.add_argument("max_diff", type=float, help="Maximum percent difference between query and subject")
    parser.add_argument("max_query_start", type=int, help="Maximum start position of blast alignment in query")
    parser.add_argument("max_subject_start", type=int, help="Maximum start position of blast alignment in subject")
    parser.add_argument("descriptive_names", type=bool, help="Would you like descriptive gene names? Yes or No.")
    args = parser.parse_args()

    seq_handle = args.input_fasta
    blast_output = args.blast_output
    blast_type = args.blast_type
    output_file = args.output_file
    max_query_start = args.max_query_start
    if blast_type == 'blastn':
        max_subject_start = args.max_subject_start
    if blast_type == 'blastx':
        max_subject_start = (args.max_subject_start)*3
    max_diff = args.max_diff
    name_seqs = args.descriptive_names
    # read all sequences into the file
    record_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))

    # the defaultdict allows to append items
    blast_dict = defaultdict(list)  
    # populate the dict
    with open(blast_output, "r") as blast_f:
        for line in blast_f:
            blast_res = line.strip().split("\t")
            name = blast_res[0]
            blast_dict[name].append(BlastRes(*blast_res))


    with open(output_file, "w") as output_handle:
        for record in record_dict:
            # print(record, record_dict[record],blast_dict[record_dict[record].id] )
            v = blast_dict[record_dict[record].id]
            coverage = 0
            annotation = "No hit"
            if len(v) >= 1:
                for res in v:
                    if blast_type == 'blastn':
                        hit_diff = 1-(int(res.length)/int(res.qlen))
<<<<<<< HEAD
                    else:
                        hit_diff = 1-((int(res.length)*3)/int(res.qlen))
                    if check_hit_position(res, coverage, max_query_start = max_query_start, max_subject_start = max_subject_start, max_diff = max_diff):
                        percent_diff = hit_diff 
                        annotation = res.sseqid + " Diff: " + str(round(percent_diff, 2))
=======
                        if check_hit_position(res, coverage, max_query_start = max_query_start, max_subject_start = max_subject_start, max_diff = max_diff):
                            percent_diff = hit_diff 
                            annotation = res.sseqid + " Diff: " + str(percent_diff)
                    if blast_type == 'blastx':
                        hit_diff = 1-((int(res.length)*3)/int(res.qlen))
                        if check_hit_position(res, coverage, max_query_start = max_query_start, max_subject_start = max_subject_start, max_diff = max_diff):
                            percent_diff = hit_diff
                            annotation = res.sseqid + " Diff: " + str(percent_diff)

>>>>>>> parent of 9e4542c... Cleaned up. Made hit_diff more intuitive.
            else:
                print("No entries for " + record_dict[record].id + "!")


            record_dict[record].description = record_dict[record].id +" Annotation: "+ annotation
            SeqIO.write(record_dict[record], output_handle, "fasta")
    if name_seqs == True:
        make_copy = "cp " + output_file + " Annotations_Named.fa"
        subprocess.call(make_copy, shell=True, executable="/bin/bash")
        get_ids = "grep 'gi' Annotations_Named.fa | cut -f2 -d'|' > gene_ids.tmp"
        subprocess.call(get_ids, shell=True, executable="/bin/bash")
        get_names = "while read ids; do title=$(efetch -id $ids -db protein -format docsum | grep -o -P '(?<=<Title>).*(?=</Title>)'); echo " + "$ids" + "," + "$title " + "; done < gene_ids.tmp > GeneID_Names.out"
        subprocess.call(get_names, shell=True, executable="/bin/bash")
        clean_names = "sed -i 's/ /_/g' GeneID_Names.out"
        subprocess.call(clean_names, shell=True, executable="/bin/bash")
        rename_genes = 'gids=($(cut -d"," -f1 "./GeneID_Names.out")) && names=($(cut -d"," -f2 "./GeneID_Names.out")) && idx=${!gids[*]} && for x in ${idx[@]}; do sed -i "s/gi|${gids[$x]}.*/${names[$x]}/g" Annotations_Named.fa; done'
        subprocess.call(rename_genes, shell=True, executable="/bin/bash")
        cleanup = "rm gene_ids.tmp"
        subprocess.call(cleanup, shell=True, executable="/bin/bash")
        
        

if __name__ == "__main__":
    main()
