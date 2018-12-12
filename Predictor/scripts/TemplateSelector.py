import os
import Helper
import EmbossNeedle

# returns the best suitable templates in array
# targetFastaName like "46XY_A"
def SelectTemplate(targetFastaName, similarity_min = 50, similarity_max = 90, can_return_itself = False, chose_only_not_similar_templates = True, maximal_similarity_between_templates = 90, number_of_templates_to_return = 3, debug = False):
    suitable_templates = []
    for fastaFile in os.listdir('../fastas/'):
        potential_template = Helper.GetFastaNAmeFromFileName(fastaFile)
        similarity = __get_similarity__(targetFastaName, potential_template, debug)

        # error getting similarity - continue
        if similarity == -1:
            continue
        if similarity < similarity_min:
            continue
        if similarity > similarity_max:
            continue
        if targetFastaName == potential_template and not can_return_itself:
            continue
        if chose_only_not_similar_templates:
            is_suitable = __compare_with_already_chosen_templates__(suitable_templates, potential_template, maximal_similarity_between_templates, debug)
            if not is_suitable:
                continue
        Helper.modify_pdb('../pdbs/' + potential_template + ".pdb", Helper.GetChainID(potential_template))
        if not Helper.check_order_of_fasta_and_pdb('./template.pdb', Helper.GetChainID(potential_template), '../fastas/' + fastaFile):
            continue

        suitable_templates.append(potential_template)
        if debug:
            print Helper.bcolors.OKGREEN + "Chosen template " + potential_template + " for target " + targetFastaName + "with similarity " + str(similarity) + "%"

        if len(suitable_templates) >= number_of_templates_to_return:
            break
    return suitable_templates

# returns similarity of two sequences(expects seq. names like "46XY_A")
# in case of error returns -1
def __get_similarity__(sequence_A, sequence_B, debug):
    seq_A = Helper.GetFastaNAmeFromFileName(sequence_A)
    seq_B = Helper.GetFastaNAmeFromFileName(sequence_B)
    with EmbossNeedle.MyEmboss(seq_A, seq_B) as instance:
        try:
            similarity = instance.GetSimilarity()
            if debug:
                print Helper.bcolors.OKBLUE + "Similarity of sequence " + sequence_A + " and sequence " + sequence_B + " is: " + str(similarity)
            return float(similarity)
        except:
            if debug:
                print Helper.bcolors.WARNING + "Error when getting similarity of sequence " + sequence_A + " and sequence " + sequence_B
            return -1


# check if the potential template is different enough to be useful
def __compare_with_already_chosen_templates__(suitable_templates, potential_template, max_similarity, debug):
    for template in suitable_templates:
        try:
            similarity = __get_similarity__(potential_template, template, False)
            # ToDo: what to do in case of error (-1)
            if similarity > max_similarity:
                if debug:
                    print Helper.bcolors.WARNING + "Template  " + potential_template + " wont be used because of high similarity ("+ str(similarity) + "%) with already chosen template " + template
                return False
        except:
            print Helper.bcolors.WARNING + "Error when getting similarity of sequence " + potential_template + " and sequence " + template + " when comparing suitable templates."
    return True


#SelectTemplate("1P9X_0", debug=True)
