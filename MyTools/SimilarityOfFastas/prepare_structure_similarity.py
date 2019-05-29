__author__ = 'Rasto'


from Bio import SeqIO
import os
from shutil import copyfile
from Bio.Emboss.Applications import NeedleCommandline
from Bio import  AlignIO
import collections
from Bio.Emboss.Applications import *
import re

import sys
reload(sys)
sys.setdefaultencoding('ISO-8859-1')

def divide_files_into_folders(dir_name):
    for filename in os.listdir(dir_name):
        if not filename.endswith(".fasta"):
            continue
        handle = open(dir_name + '/' + filename, "rU")
        for record in SeqIO.parse(handle, "fasta") :
            move_to_folder([dir_name, filename, record, len(record.seq)])
            #debug print [dir_name + '/' + filename, record, len(record.seq)]
        handle.close()

def move_to_folder(data): #data is list of [dir_name, filename , fasta, fasta length]
    directories = [{"from":0, "till":50, "folder":"0_50"}, {"from":51, "till":100, "folder":"51_100"}, {"from":101, "till":500, "folder":"101_500"}, {"from":501, "till":1000, "folder":"501_1000"}, {"from":1001, "till":2000, "folder":"1001_2000"}, {"from":2001, "till":10000, "folder":"2001_10000"}]
        #[{"from":0, "till":90000, "folder":"fastasMod"}]#[{"from":0, "till":50, "folder":"0_50"}, {"from":51, "till":100, "folder":"51_100"}, {"from":101, "till":500, "folder":"101_500"}, {"from":501, "till":1000, "folder":"501_1000"}, {"from":1001, "till":2000, "folder":"1001_2000"}, {"from":2001, "till":10000, "folder":"2001_10000"}]
    for dir in directories:
        if data[3] >= dir["from"] and data[3] <= dir["till"]:
            print data[2].id[:4] + "_" + data[2].id[5] + "    " + str(data[3]) + "     " + dir["folder"]
            #copyfile(data[0] + "/" + data[1], data[0] + "/" + dir["folder"] + "/" + data[2].id[:4] + "_" + data[2].id[5] + ".fasta" )
            output_handle = open(data[0] + "/" + dir["folder"] + "/" + data[2].id[:4] + "_" + data[2].id[5] + ".fasta", "w")
            SeqIO.write(data[2], output_handle, "fasta")
            output_handle.close()
            break

def create_similarity_matrix():
    # choose files from which you want to generate pairs
    # very slow if all are chosen
    directories = [
        #{"from":0, "till":50, "folder":"0_50"},
        {"from":51, "till":100, "folder":"51_100"}
        #{"from":101, "till":500, "folder":"101_500"}
        #{"from":501, "till":1000, "folder":"501_1000"},
        #{"from":1001, "till":2000, "folder":"1001_2000"},
        #{"from":2001, "till":10000, "folder":"2001_10000"}
     ]
    for dir in directories:
        heatmap_state =  [[0 for i in range(4)] for j in range(4)]
        buckets_dictionary = {'30-45s0-15g':[],'30-45s15-30g':[],'30-45s30-45g':[],'30-45s45-60g':[],'45-60s0-15g':[],'45-60s15-30g':[],'45-60s30-45g':[],'45-60s45-60g':[],'60-75s0-15g':[],'60-75s15-30g':[],'60-75s30-45g':[],'60-75s45-60g':[],'75-90s0-15g':[],'75-90s15-30g':[],'75-90s30-45g':[],'75-90s45-60g':[] }
        f = open(str(dir['from']) + "_" + str(dir['till']) + ".txt",'w')
        #debug out_file_name = str(dir['from']) + "_" + str(dir['till']) + ".txt"
        for filename_row in os.listdir("fasta/" + dir["folder"]):
            row = {}
            for filename_column in os.listdir("fasta/" + dir["folder"]):
                try:
                    if filename_row == filename_column:
                        row.update({101:filename_column})
                        continue
                    emboss_result = compute_emboss_alignment("fasta/" + dir["folder"] + "/" + filename_row, "fasta/" + dir["folder"] + "/" + filename_column)
                    prepare_data_for_heatmap(heatmap_state, emboss_result[1], emboss_result[0], buckets_dictionary,filename_row, filename_column)
                    if emboss_result[0] in row:
                        row[emboss_result[0]].append(filename_column)
                    else:
                        row.update({emboss_result[0] : [emboss_result[1], filename_column]})
                        #row.update({emboss_result : filename_column})
                except Exception as e:
                    print str(e)
            print >>f, collections.OrderedDict(sorted(row.items(),reverse=True))
            print collections.OrderedDict(sorted(row.items(),reverse=True))
        f.close()
        for x in buckets_dictionary:
            fl = open(name=x,mode='w')
            for y in buckets_dictionary[x]:
                fl.write(y[0] + " " + y[1] + "\n")
            fl.close()
        #debug print buckets_dictionary
        create_heatmap(heatmap_state)

def compute_emboss_alignment(input_file_1, input_file_2):
    #debug print input_file_1 + "       " +  input_file_2
    needle_cli = NeedleCommandline(asequence=input_file_1, \
                               bsequence=input_file_2, \
                               gapopen=10, \
                               gapextend=0.5, \
                               outfile="needle_fname", \
                               )
    needle_cli()
    return get_similarity_from_file("needle_fname")
    #debug aln = AlignIO.read("needle_fname", "emboss")
    #debug print aln

def get_similarity_from_file(file_name):
    file = open(name=file_name,mode='r')
    similarity = 'no_match'
    gap = 'no_match'
    for line in file:
        #debug print line
        matched = re.match('^# Similarity:[^\(]*\((.*)%\)', line)
        matched_gap = re.match('^# Gaps:[^\(]*\((.*)%\)', line)
        if matched:
            similarity = matched.group(1)
            #debug print matched.group(1)
        if matched_gap:
            gap = matched_gap.group(1)
    file.close()
    return [float(similarity), float(gap)]

def prepare_data_for_heatmap(current_state, gap, similarity, buckets_dictionary, fasta1, fasta2):
    # current state is two dimensional array rows means similarity gap and columns similarity buckets
        #row 1
    cell = [fasta1[:4].upper() + fasta1[4:], fasta2[:4].upper() + fasta2[4:]]
    if 30 <= similarity < 45 and 0 <= gap < 15:
        current_state[0][0] = current_state[0][0] + 1
        buckets_dictionary['30-45s0-15g'].append(cell)
    if 30 <= similarity < 45 and 15 <= gap < 30:
        current_state[0][1] = current_state[0][1] + 1
        buckets_dictionary['30-45s15-30g'].append(cell)
    if 30 <= similarity < 45 and 30 <= gap < 45:
        current_state[0][2] = current_state[0][2] + 1
        buckets_dictionary['30-45s30-45g'].append(cell)
    if 30 <= similarity < 45 and 45 <= gap < 60:
        current_state[0][3] = current_state[0][3] + 1
        buckets_dictionary['30-45s45-60g'].append(cell)
        #row 2
    if 45 <= similarity < 60 and 0 <= gap < 15:
        current_state[1][0] = current_state[1][0] + 1
        buckets_dictionary['45-60s0-15g'].append(cell)
    if 45 <= similarity < 60 and 15 <= gap < 30:
        current_state[1][1] = current_state[1][1] + 1
        buckets_dictionary['45-60s15-30g'].append(cell)
    if 45 <= similarity < 60 and 30 <= gap < 45:
        current_state[1][2] = current_state[1][2] + 1
        buckets_dictionary['45-60s30-45g'].append(cell)
    if 45 <= similarity < 60 and 45 <= gap < 60:
        current_state[1][3] = current_state[1][3] + 1
        buckets_dictionary['45-60s45-60g'].append(cell)
        #row 3
    if 60 <= similarity < 75 and 0 <= gap < 15:
        current_state[2][0] = current_state[2][0] + 1
        buckets_dictionary['60-75s0-15g'].append(cell)
    if 60 <= similarity < 75 and 15 <= gap < 30:
        current_state[2][1] = current_state[2][1] + 1
        buckets_dictionary['60-75s15-30g'].append(cell)
    if 60 <= similarity < 75 and 30 <= gap < 45:
        current_state[2][2] = current_state[2][2] + 1
        buckets_dictionary['60-75s30-45g'].append(cell)
    if 60 <= similarity < 75 and 45 <= gap < 60:
        current_state[2][3] = current_state[2][3] + 1
        buckets_dictionary['60-75s45-60g'].append(cell)
        #row 4
    if 75 <= similarity < 90 and 0 <= gap < 15:
        current_state[3][0] = current_state[3][0] + 1
        buckets_dictionary['75-90s0-15g'].append(cell)
    if 75 <= similarity < 90 and 15 <= gap < 30:
        current_state[3][1] = current_state[3][1] + 1
        buckets_dictionary['75-90s15-30g'].append(cell)
    if 75 <= similarity < 90 and 30 <= gap < 45:
        current_state[3][2] = current_state[3][2] + 1
        buckets_dictionary['75-90s30-45g'].append(cell)
    if 75 <= similarity < 90 and 45 <= gap < 60:
        current_state[3][3] = current_state[3][3] + 1
        buckets_dictionary['75-90s45-60g'].append(cell)
    #debug print buckets_dictionary

def create_heatmap(data):
    import matplotlib.pyplot as plt
    import numpy as np

    ylabel= "Similarity"
    xlabel="Gaps"
    column_labels = ['30-45%','45-60%','60-75%','75-90%'] #gaps
    row_labels = ['0-15%','15-30%','30-45%','45-60%'] #similarity
    #data = np.random.rand(4,4)
    data = np.array(data)
    fig, ax = plt.subplots()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    heatmap = ax.pcolor(data, cmap=plt.cm.RdYlGn)

    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            plt.text(x + 0.5, y + 0.5, '%.0f' % data[y, x],
                     horizontalalignment='center',
                     verticalalignment='center',
                    )

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

    ax.set_xticklabels(row_labels, minor=False)
    ax.set_yticklabels(column_labels, minor=False)
    print data
    plt.show()




#divide_files_into_folders("fasta")
#compute_emboss_alignment("fasta/51_100/5BTM_B.fasta","fasta/51_100/4ZNP_B.fasta")
create_similarity_matrix()

