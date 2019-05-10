from Bio.PDB import *
import os
import EmbossNeedle

CONST_SIMILARITY_LIMIT = 99
CONST_FILTER_SIMILAR_STRUCTURES = False


def save_fasta_file(sequence, fasta_file_name, sequence_name, chain_id):
    with open(fasta_file_name, 'w') as output_file_handler:
        output_file_handler.writelines(">" + sequence_name + ":" + chain_id + "|PDBID|CHAIN|SEQUENCE\n")
        output_file_handler.writelines(sequence)


def get_fasta_seq(fasta):
    from Bio import SeqIO
    records = list(SeqIO.parse(fasta, "fasta"))
    fasta_seq = ""
    for r in records[0]:
        fasta_seq = fasta_seq + r
    return fasta_seq


def create_fasta_from_pdb(pdb):
    parser_pdb = PDBParser()
    structure = parser_pdb.get_structure('self', './pdbs/' + pdb)
    model = structure[0]
    for chain in model:
        fasta_file_name = './fastas/' + str.upper(pdb[:4]) + '_' + str.upper(chain.id) + '.fasta'
        downloaded_fasta_name = './fastas_downloaded/' + str.upper(pdb[:4]) + '_' + str.upper(chain.id) + '.fasta'
        continuity_counter = 0
        fasta_seq = ""
        invalid_residue_counter = 0

        # if PDB is valid for downloaded fasta sequence - use this fasta sequence
        try:
            downloaded_fasta_seq = get_fasta_seq(downloaded_fasta_name)
            pdb_fasta_mismatch = False
            for res in chain:
                if res.resname[2] != 'A' and res.resname[2] != 'G' and res.resname[2] != 'C' and res.resname[2] != 'U':
                    continue
                if len(downloaded_fasta_seq) > res.id[1] and downloaded_fasta_seq[res.id[1]-1] != res.resname[2]:
                    pdb_fasta_mismatch = True
                    break
        except BaseException as e:
            print e
            pdb_fasta_mismatch = True

        if pdb_fasta_mismatch:
            # else create fasta sequence from PDB
            for res in chain:
                # count invalid residues in pdb
                if res.resname[2] != 'A' and res.resname[2] != 'G' and res.resname[2] != 'C' and res.resname[2] != 'U':
                    invalid_residue_counter += 1
                    continue
                continuity_counter += 1
                # set X into gaps in numbering
                while res.id[1] > continuity_counter:
                    fasta_seq += 'X'
                    invalid_residue_counter += 1
                    continuity_counter += 1
                # residue index should only increase
                if res.id[1] < continuity_counter:
                    continue
                fasta_seq += res.resname[2]
            # filter out too short sequences
            if continuity_counter < 10:
                continue
            if invalid_residue_counter/float(continuity_counter) > 0.1:
                continue
        else:
            fasta_seq = downloaded_fasta_seq

        # filter out too similar sequences (we do not want to have different named sequences, which are the same)
        if CONST_FILTER_SIMILAR_STRUCTURES:
            similarity = 0
            temp_file = './possible_candidate.fasta'
            save_fasta_file(fasta_seq, temp_file, pdb[:4], chain.id)
            parsed_similarity = 0
            for fasta_file in os.listdir('./fastas/'):
                with EmbossNeedle.MyEmboss(pdb[:4] + "_" + chain.id, fasta_file[:-6],
                                       typeOfAlignemnt='global', target_fasta_file=temp_file) as emboss_aln:
                    similarity = emboss_aln.GetSimilarity()
                    try:
                        parsed_similarity = float(similarity)
                        if parsed_similarity > CONST_SIMILARITY_LIMIT:
                            break
                    except:
                        print "Could not get similarity of " + pdb[:4] + "_" + chain.id + ' and ' + fasta_file[:-6]

            if parsed_similarity > CONST_SIMILARITY_LIMIT:
                print 'Similarity too high (' + similarity + '). Skip + ' + fasta_file_name
                continue
        save_fasta_file(fasta_seq, fasta_file_name, pdb[:4], chain.id)


for pdb_file in os.listdir('./pdbs/'):
    create_fasta_from_pdb(pdb_file)
