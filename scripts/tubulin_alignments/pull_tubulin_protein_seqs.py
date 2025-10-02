import os
import time

# get the date in yyyymmmdd format

date = time.strftime("%Y%m%d", time.localtime())

out_dir = f"data/proc/tubulin_alignments/{date}_tubulin_protein_seqs_{date}"

# check if the directory exists, if not create it
if not os.path.exists(out_dir):
    print(f"Creating output directory: {out_dir}")
    os.makedirs(out_dir)


ce_prot_file = "data/blast_data/20230324_blast/c_elegans/ce_all_prot_flat.fa"
cb_prot_file = "data/blast_data/20230324_blast/c_briggsae/cb_all_prot_flat.fa"
ct_prot_file = "data/blast_data/20230324_blast/c_tropicalis/ct_all_prot_flat.fa"


# out_file = 'test_pull.fa'
# out_file_cb = 'test_pull_cb.fa'
# out_file_ce = 'test_pull_ce.fa'


# function that will pull the line with the sequence id and the next line with the sequence
def pull_protein_sequence(protein_id, prot_file, out_file):
    print(f"Pulling protein sequence for ID: {protein_id} from file: {prot_file}")
    with open(prot_file, "r") as f:
        lines = f.readlines()
        print(f"Total lines read from {prot_file}: {len(lines)}")
        for i, line in enumerate(lines):
            if protein_id in line:
                print(f"Found protein ID {protein_id} at line {i}")
                with open(out_file, "a") as out_f:
                    out_f.write(line)
                    out_f.write(lines[i + 1])
                print(f"Protein sequence written to {out_file}")
                break
        else:
            print(f"Protein ID {protein_id} not found in {prot_file}")


# pull_protein_sequence('QX1410.13336.1', cb_prot_file, out_file_cb)
# pull_protein_sequence('C54C6.2.1', ce_prot_file, out_file_ce)

#### ben_1 ####
# ben_1 = {
#     'cb': 'QX1410.13336.1',
#     'ce': 'C54C6.2.1',
#     'ct': 'NIC58.15504.1'
# }


def pull_and_rename_protein_sequence(protein_id, formal_name, prot_file, out_file):
    print(f"Pulling protein sequence for ID: {protein_id} from file: {prot_file}")
    with open(prot_file, "r") as f:
        lines = f.readlines()
        print(f"Total lines read from {prot_file}: {len(lines)}")
        for i, line in enumerate(lines):
            if protein_id in line:
                print(f"Found protein ID {protein_id} at line {i}")
                with open(out_file, "a") as out_f:
                    out_f.write(
                        f">{formal_name}\n"
                    )  # Write the formal name as the new header
                    out_f.write(lines[i + 1])  # Write the sequence


def process_protein_sequences(
    protein_dict, ce_prot_file, cb_prot_file, ct_prot_file, out_file
):
    for key, protein_info in protein_dict.items():
        protein_id = protein_info[0]  # Extract the transcript ID
        formal_name = protein_info[1]  # Extract the formal name
        if key == "ce":
            pull_and_rename_protein_sequence(
                protein_id, formal_name, ce_prot_file, out_file
            )
        elif key == "cb":
            pull_and_rename_protein_sequence(
                protein_id, formal_name, cb_prot_file, out_file
            )
        elif key == "ct":
            pull_and_rename_protein_sequence(
                protein_id, formal_name, ct_prot_file, out_file
            )
        else:
            print(f"Invalid species key {key}")


# use a list for each entry to include transcript ID and the formal name
ben_1 = {
    "ce": ["C54C6.2.1", "CE_ben-1"],
    "cb": ["QX1410.13336.1", "CB_ben-1"],
    "ct": ["NIC58.15504.1", "CT_ben-1"],
}

out_file = f"{out_dir}/ben_1.fa"

# check if the file exists and contains data, if so, delete it
if os.path.exists(out_file):
    print(f"Deleting existing file: {out_file}")
    os.remove(out_file)

process_protein_sequences(ben_1, ce_prot_file, cb_prot_file, ct_prot_file, out_file)

tbb_1 = {
    "ce": ["K01G5.7.1", "CE_tbb-1"],
    "cb": ["QX1410.12556.1", "CB_tbb-1"],
    "ct": ["NIC58.14724.4", "CT_tbb-1"],
}

out_file = f"{out_dir}/tbb_1.fa"

# check if the file exists and contains data, if so, delete it
if os.path.exists(out_file):
    print(f"Deleting existing file: {out_file}")
    os.remove(out_file)

process_protein_sequences(tbb_1, ce_prot_file, cb_prot_file, ct_prot_file, out_file)


tbb_2 = {
    "ce": ["C36E8.5.1", "CE_tbb-2"],
    "cb": ["QX1410.14085.1", "CB_tbb-2"],
    "ct": ["NIC58.15977.2", "CT_tbb-2"],
}

out_file = f"{out_dir}/tbb_2.fa"

# check if the file exists and contains data, if so, delete it
if os.path.exists(out_file):
    print(f"Deleting existing file: {out_file}")
    os.remove(out_file)

process_protein_sequences(tbb_2, ce_prot_file, cb_prot_file, ct_prot_file, out_file)
