import time
import EmbossNeedle
import Helper
import superposition
from Bio.PDB import *
from Bio import SeqIO
import os
import config


class Gap(object):
    gap_from = -1  # including
    gap_till = -1  # including +1
    gap_length = -1
    similarity = -1
    template_name = ''  # XXXX_A
    template_from = ''  # including
    template_till = ''  # including
    template_sequence = ''
    target_sequence = ''
    extended_gap_from = -1 # including
    extended_gap_till = -1 # including +1
    extended_gap_length = extended_gap_till - extended_gap_from
    target_connector_ids = []
    template_connector_ids = []
    template_fill_gap_ids = []
    template_chain_id = ''
    template_exists = False

    def fill_gap_object(self, gap_from, gap_till, gap_length):
        self.gap_from = gap_from
        self.gap_till = gap_till
        self.gap_length = gap_length


class MultipleTemplateModul:
    """Fill long gaps in predicted structure with residues chains from another structure."""
    longGaps = []
    enrichedGaps = []
    CONST_PDB_DIR = "../pdbs/"
    modelID = 0
    CONST_FASTAS_DIR = "../fastas/"
    CONST_FASTA_FILE_TYPE = ".fasta"
    #LONG_GAP_LOW_LIMIT = 9
    NUMBER_OF_RES_TO_EXTEND_GAP_FROM_EACH_SIDE = 5
    MIN_SIMILARITY_TO_ACCEPT_ALN = 70
    MAX_GLOBAL_SIMILARITY_TO_ACCEPT_ALN = 93
    # if zero then there is gap on that index in target structure if 1 there is residue
    target_structure_gap_map = []

    def __init__(self, target_structure, target_name, target_structure_chain_id):
        """pdb file with path to fill long gaps, XXXX_X fasta name of target structure, chain_id of target_structure (usually from primary template)"""
        self.target_name = target_name
        self.target_structure_chain_id = target_structure_chain_id
        self.target_structure = target_structure # structure with gaps to fill
        self.target_fasta_length = self.__get_target_fasta_length__()
        self.target_fasta_seq = self.__get_target_fasta_seq__()
        self.gaps = []
        self.indexed_fastas = list()

    def identify_gaps_in_structure(self):
        """returns list of [from(including), till(including + 1), gapLength]"""
        parser_pdb = PDBParser()
        predicted_part_structure = parser_pdb.get_structure("origReadyForFARFARStr", self.target_structure)
        gaps = []
        last_res_id = 0
        self.target_structure_gap_map = [0 for i in xrange(0, self.target_fasta_length+1)]
        predicted_part_structure[self.modelID][self.target_structure_chain_id].child_list.sort(key=lambda x: x.id[1])
        for res in predicted_part_structure[self.modelID][self.target_structure_chain_id]:
            self.target_structure_gap_map[res.id[1]] = 1
            if res.id[1] > last_res_id + 1:
                gap = [last_res_id + 1, res.id[1],
                       res.id[1] - (last_res_id + 1)]  # from(including), till(including + 1), gapLength
                if gap[2] >= config.configMinGabLengthToOwnPrediction: #self.LONG_GAP_LOW_LIMIT:
                    gap_object = Gap()
                    gap_object.fill_gap_object(gap[0], gap[1], gap[2])
                    self.gaps.append(gap_object)
                    gaps.append(gap)
            last_res_id = res.id[1]
        if last_res_id < self.target_fasta_length:
            gap = [last_res_id + 1, self.target_fasta_length,
                   self.target_fasta_length - (last_res_id + 1)]  # from(including), till(including + 1), gapLength
            if gap[2] >= config.configMinGabLengthToOwnPrediction: #self.LONG_GAP_LOW_LIMIT:
                gap_object = Gap()
                gap_object.fill_gap_object(gap[0], gap[1], gap[2])
                self.gaps.append(gap_object)
                gaps.append(gap)
        return gaps

    def index_fastas(self):
        """when finding suitable templates align only those long enough and making any sense"""
        for fasta in os.listdir(self.CONST_FASTAS_DIR):
            records = list(SeqIO.parse(self.CONST_FASTAS_DIR + fasta, "fasta"))
            fasta_seq = ""
            invalid_residues_count = 0
            for r in records[0]:
                fasta_seq = fasta_seq + r
                if r != 'A' and r != 'G' and r != 'C' and r != 'U':
                    invalid_residues_count += 1
            if len(fasta_seq) == 0:
                suitable_for_rna_template = 1
            else:
                suitable_for_rna_template = float(invalid_residues_count)/len(fasta_seq)
            self.indexed_fastas.append((fasta, len(fasta_seq), suitable_for_rna_template))

    def find_suitable_alignments_for_gaps(self):
        for gap in self.gaps:
            print "Working on gap from " + str(gap.gap_from) + ", gap till " + str(gap.gap_till)
            self.__get_extended_gap__(gap)
            gap_sequence = self.target_fasta_seq[gap.extended_gap_from-1:gap.extended_gap_till]
            temp_file_name = "temp_multiple_template_module.fasta"
            temp_seq_name = "TEMP_T"
            self.__save_temp_fasta_file__(gap_sequence, temp_file_name, temp_seq_name[:-2])
            for possible_template in self.indexed_fastas:
                if possible_template[1] < gap.gap_till - gap.gap_from:
                    # print "Skip" + possible_template[0]
                    continue
                # 50% of different residues than A,G,C,T
                if possible_template[2] > 0.5:
                    # print "Skip" + possible_template[0]
                    continue
                if possible_template[0][:4] == self.target_name.upper()[:4]:
                    # this is the predicted structure - does not make sense to use it
                    continue

                with EmbossNeedle.MyEmboss(temp_seq_name, possible_template[0][:-6],
                                           typeOfAlignemnt='global', target_fasta_file=temp_file_name) as emboss_aln:
                    a_seq, b_seq, b_start, b_end, similarity = emboss_aln.ParseSemiGlobalAlignment()

                if similarity is not None and float(similarity)*100 > self.MIN_SIMILARITY_TO_ACCEPT_ALN:
                    if not Helper.is_fast_and_pdb_ordered_correctly(self.CONST_PDB_DIR + Helper.TrimPDBName(possible_template[0]) + '.pdb', Helper.GetChainID(possible_template[0]), self.CONST_FASTAS_DIR + possible_template[0]):
                        print "Possible template fasta and pdp does not match sequence or file does not exists = " + possible_template[0] \
                              + ". Continue with another template."
                        continue
                    with EmbossNeedle.MyEmboss(self.target_name, possible_template[0][:-6]) as emboss_aln:
                        if float(emboss_aln.GetSimilarity()) >= self.MAX_GLOBAL_SIMILARITY_TO_ACCEPT_ALN:
                            print "Structures " + possible_template[
                                0] + " and " + self.target_name + " are probably duplicates (global similarity = " + str(emboss_aln.GetSimilarity()) + "). There is no reason to predict one using another. Skip..."
                            continue

                    gap.similarity = similarity
                    gap.template_from = b_start
                    gap.template_till = b_end
                    gap.template_name = possible_template[0][:-8]
                    gap.template_chain_id = possible_template[0][-7]
                    gap.template_sequence = b_seq
                    gap.target_sequence = a_seq
                    target_alignment_indexes = []
                    counter = 0
                    # number target alignment
                    for i in gap.target_sequence:
                        if i != '-':
                            target_alignment_indexes.append(gap.extended_gap_from+counter)
                            counter += 1
                        else:
                            target_alignment_indexes.append(i)
                    counter = 0
                    template_alignment_indexes = []
                    # number template alignment
                    for i in gap.template_sequence:
                        if i != '-':
                            template_alignment_indexes.append(gap.template_from+counter)
                            counter += 1
                        else:
                            template_alignment_indexes.append(i)
                    gap.target_connector_ids = []
                    gap.template_connector_ids = []
                    gap.template_fill_gap_ids = []
                    gap.template_exists = True
                    # prepare input for superposition function
                    for i in xrange(0, len(target_alignment_indexes)):
                        if gap.target_sequence[i] == gap.template_sequence[i]:
                            if target_alignment_indexes[i] < gap.gap_from:
                                # tento if a jeden po d9m sa stara o to, aby sedel pocet spajajucich resides
                                if self.target_structure_gap_map[target_alignment_indexes[i]] == 1 and gap.template_sequence[i] != 'N':
                                    gap.target_connector_ids.append(target_alignment_indexes[i])
                                    gap.template_connector_ids.append(template_alignment_indexes[i])
                            else:
                                if target_alignment_indexes[i] < gap.gap_till:
                                    # inserting tuples (index of residue in template structure, target_alignment_residue -> for renumbering of passed residues)
                                    gap.template_fill_gap_ids.append((template_alignment_indexes[i], target_alignment_indexes[i]))
                                else: # target_alignment_indexes[i] < gap.extended_gap_till
                                    if self.target_structure_gap_map[target_alignment_indexes[i]] == 1 and gap.template_sequence[i] != 'N':
                                        gap.target_connector_ids.append(target_alignment_indexes[i])
                                        gap.template_connector_ids.append(template_alignment_indexes[i])
                    print "Match for gap found with similarity " + str(similarity) + " and template " + possible_template[0]
                    # TODO: najdi najlepssiu zhodu a nie prvu
                    break


    def callSuperpositionMethod(self):
            for gap in self.gaps:
                try:
                    if gap.template_exists:
                        superposition.Merge(target_name=self.target_structure,
                                            template_name="../pdbs/" + str.lower(gap.template_name) + ".pdb",
                                            fill_gap_template_ids=gap.template_fill_gap_ids,
                                            target_connector_residue_ids=gap.target_connector_ids,
                                            template_connector_residue_ids=gap.template_connector_ids,
                                            template_chain_id=gap.template_chain_id,
                                            target_chain_id=self.target_structure_chain_id)
                except Exception as e:
                    print("ERROR in superposition.py: " + str(e))


    def __get_extended_gap__(self, gap):
        if gap.gap_from <= self.NUMBER_OF_RES_TO_EXTEND_GAP_FROM_EACH_SIDE:
            start_aln_residue = 1
        else:
            start_aln_residue = gap.gap_from - self.NUMBER_OF_RES_TO_EXTEND_GAP_FROM_EACH_SIDE
        if gap.gap_till + self.NUMBER_OF_RES_TO_EXTEND_GAP_FROM_EACH_SIDE >= self.target_fasta_length:
            end_aln_residue = self.target_fasta_length - 1
        else:
            end_aln_residue = gap.gap_till + self.NUMBER_OF_RES_TO_EXTEND_GAP_FROM_EACH_SIDE
        gap.extended_gap_from = start_aln_residue
        gap.extended_gap_till = end_aln_residue

    def __get_target_fasta_length__(self):
        fasta_file = open(self.CONST_FASTAS_DIR + self.target_name.upper() + self.CONST_FASTA_FILE_TYPE, 'rU')
        for rec in SeqIO.parse(fasta_file, 'fasta'):
            seq_len = len(rec)
        fasta_file.close()
        return seq_len

    def __get_target_fasta_seq__(self):
        fasta_file = open(self.CONST_FASTAS_DIR + self.target_name.upper() + self.CONST_FASTA_FILE_TYPE, 'rU')
        for rec in SeqIO.parse(fasta_file, 'fasta'):
            seq = rec
        fasta_file.close()
        return seq

    @staticmethod
    def __save_temp_fasta_file__(sequence, temp_file_name, temp_seq_name):
        with open(temp_file_name, 'w') as output_file_handler:
            output_file_handler.writelines(">"+temp_seq_name+":T|PDBID|CHAIN|SEQUENCE\n")
            output_file_handler.writelines(sequence)


def test_identify_gaps_in_structure(given_gaps):
    mtm = MultipleTemplateModul("test_data/1p9x.pdb", "1p9x_0", "0")
    identified_gaps = mtm.identify_gaps_in_structure()
    print "Starting... test_identify_gaps_in_structure"
    # return differences
    diff = [x for x in given_gaps if x not in identified_gaps] + [x for x in identified_gaps if x not in given_gaps]
    if len(diff) == 0:
        print "Test passed successfully."
    else:
        print diff
    print "End... test_identify_gaps_in_structure"


def run_module(target_structure, target_name, target_structure_chain_id):
    print "MULTIPLE_TEMPLATES_MODULE: start"
    #mtm = MultipleTemplateModul("test_data/4znp.pdb", "4ZNP_B", "B")
    mtm = MultipleTemplateModul(target_structure, target_name, target_structure_chain_id)
    mtm.index_fastas()
    identified_gaps = mtm.identify_gaps_in_structure()
    mtm.find_suitable_alignments_for_gaps()
    mtm.callSuperpositionMethod()
    print "MULTIPLE_TEMPLATES_MODULE: end"


def run_all_tests():
    test_identify_gaps_in_structure([[1, 26, 25], [600, 701, 101], [404, 430, 26], [249, 292, 43], [2850, 2880, 30]])


#run_module("test_data/1p9x.pdb", "1p9x_0", "0")

