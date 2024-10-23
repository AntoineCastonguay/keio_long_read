import subprocess
import os
import pathlib
import gzip
from glob import glob
from Bio import SeqIO
import csv
import pandas as pd
from concurrent.futures import ThreadPoolExecutor


class Methods(object):
    accepted_extensions = ['.fq', '.fq.gz',
                           '.fastq', '.fastq.gz',
                           '.fasta', '.fasta.gz',
                           '.fa', '.fa.gz',
                           '.fna', '.fna.gz']

    @staticmethod
    def check_input(my_input):
        if not os.path.exists(my_input):
            raise Exception('Please select an existing file or folder as input.')

        # Check if folder
        if os.path.isdir(my_input):
            file_list = os.listdir(my_input)  # List content of folder
        else:  # elif os.path.isfile(my_input):
            file_list = [my_input]

        # if folder is not empty and all files have the accepted extensions
        if not all([f.endswith(tuple(Methods.accepted_extensions)) for f in file_list]):
            raise Exception('Make sure files in input folder end with {}'.format(Methods.accepted_extensions))

    @staticmethod
    def check_ref(ref):
        if not os.path.isfile(ref):
            raise Exception('The reference file provided does not exist')

        with gzip.open(ref, 'rt') if ref.endswith('.gz') else open(ref, 'r') as f:
            first_header = f.readline()
            first_character = first_header[0]
            if first_character != '>':
                raise Exception('The reference file provided does not appear to be a valid fasta file.')

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
    
    @staticmethod
    def find_sam_files(folder):
        sam_files = []
        for filename in os.listdir(folder):
            if filename.endswith('.sam'):
                sam_files.append(os.path.join(folder, filename))
        return sam_files

    @staticmethod
    def list_to_file(my_list, output_file):
        with open(output_file, 'wt') as f:
            for l in my_list:
                f.write('{}\n'.format(l))

    @staticmethod
    def list_files_in_folder(folder, extension):
        return glob(folder + '/*' + extension)

    @staticmethod
    def flag_done(flag_file):
        with open(flag_file, 'w') as f:
            pass
    
    @staticmethod
    def process_convert_file(f, output_folder):
        if "_barcode" in f:
            barcode = f.split('_')[-1].split('.')[0]
            print(f'\t{barcode}')
            output_file = os.path.join(output_folder, f"{barcode}.fasta")

            with gzip.open(f, "rt") as fastq_file:
                with open(output_file, "w") as fasta_file:
                    for record in SeqIO.parse(fastq_file, "fastq"):
                        fasta_file.write(f">{record.id}\n")
                        fasta_file.write(f"{str(record.seq)}\n")

            return {barcode: output_file}

    @staticmethod
    def fastq_to_fasta(fastq_folder, output_folder):
        Methods.make_folder(output_folder)
        
        files = Methods.list_files_in_folder(fastq_folder, 'fastq.gz')
        my_dict = {}

        # Utilisation du ThreadPoolExecutor pour paralléliser le traitement
        with ThreadPoolExecutor(max_workers=10) as executor:
            results = list(executor.map(lambda f: Methods.process_convert_file(f, output_folder), files))

        # Fusion des résultats dans my_dict
        for result in results:
            if result:
                my_dict.update(result)

        return my_dict

    @staticmethod
    def blast(read, ref, output_folder):
        Methods.make_folder(output_folder)
        for key,f in read.items():
            print(f'\t{key}')
            out = f'{output_folder}{key}/'
            Methods.make_folder(out)
            makeblastdb_cmd = ['makeblastdb', '-in', f, '-dbtype', 'nucl', '-out', f'{out}{key}_reads_db']
            subprocess.run(makeblastdb_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)

            out_res = f'{output_folder}all_res/'
            Methods.make_folder(out_res)
            blast_cmd = ['blastn', '-query', ref, '-db', f'{out}{key}_reads_db', '-out', f'{out_res}{key}_results.txt', '-outfmt', '6', '-max_hsps', '1']
            subprocess.run(blast_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)

    
    @staticmethod
    def extract(res, read, output_folder):
        Methods.make_folder(output_folder)

        for key,f in res.items():
            print(f'\t{key}')

            dict_match = {}

            with open(f, "r") as file:
                reader = csv.reader(file, delimiter="\t")
                for row in reader:
                    if float(row[2]) > float(90) and int(row[3]) > 100:
                        if key not in dict_match:
                            dict_match[key] = {}

                        if row[1] not in dict_match[key]:
                            dict_match[key][row[1]] = [row[3], row[8], row[9]]
                        elif row[3] > dict_match[key][row[1]][0]:
                            dict_match[key][row[1]] = [row[3], row[8], row[9]]
        
            output_file_l = os.path.join(output_folder, f"{key}_l_selected_sequences.fasta")
            output_file_r = os.path.join(output_folder, f"{key}_r_selected_sequences.fasta")

            with open(output_file_l, "w") as fasta_output, open(output_file_r, "w") as fasta_output2:
                # Parcours des séquences dans le fichier FASTA
                found = False  # Variable pour vérifier si au moins une séquence a été trouvée
                for record in SeqIO.parse(read[key], "fasta"):
                    if record.id in dict_match[key]:
                        # Écrire la séquence au format FASTA
                        a,b = 1,2
                        if dict_match[key][record.id][a] > dict_match[key][record.id][b]:
                            b,a = 1,2
                        fasta_output.write(f">{record.id}\n")
                        fasta_output.write(f"{record.seq[1:int(dict_match[key][record.id][a])]}\n")
                        fasta_output2.write(f">{record.id}\n")
                        fasta_output2.write(f"{record.seq[int(dict_match[key][record.id][b]):len(record.seq)]}\n")
                        found = True

                if not found:
                    print(f"\tError : {key}_selected_sequences.fasta")
                else:
                    print(f"\tFichier FASTA created : {key}_selected_sequences.fasta")

    @staticmethod
    def run_blast2(key, pos, h, ref_db, output_folder):
        print(f'\t{key} {pos}')

        out_res = f'{output_folder}all_res/'
        Methods.make_folder(out_res)

        # Commande BLAST
        blast_cmd = ['blastn', '-query', h, '-db', ref_db, '-out', f'{out_res}{key}_{pos}_results.txt', '-outfmt', '6', '-max_hsps', '1']
        subprocess.run(blast_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)

        # Lecture du fichier de résultats avec pandas
        output_file_txt = {}
        df = pd.read_csv(f"{out_res}{key}_{pos}_results.txt", sep="\t", header=None)

        # Vérification si le fichier de résultats est vide
        if df.empty:
            return key, output_file_txt
        
        # Lecture des colonnes pertinentes (colonne 0 pour query_id, colonnes 8 et 9 pour les positions sujet)
        for _, row in df.iterrows():
            query_id = row[0]
            pos_1_subject = int(row[8])
            pos_2_subject = int(row[9])
            list_var = [pos_1_subject, pos_2_subject]

            if query_id not in output_file_txt:
                output_file_txt[query_id] = {}
            output_file_txt[query_id][pos] = list_var

        return key, output_file_txt

    @staticmethod
    def blast2(align, ref, output_folder):
        Methods.make_folder(output_folder)

        # Création de la base de données BLAST
        makeblastdb_cmd = ['makeblastdb', '-in', ref, '-dbtype', 'nucl', '-out', f'{output_folder}ref_db']
        subprocess.run(makeblastdb_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)

        ref_db = f'{output_folder}ref_db'
        all_results = {}

        # Parallélisation des commandes BLAST avec ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = []
            for key, f in align.items():
                for pos, h in f.items():
                    futures.append(executor.submit(Methods.run_blast2, key, pos, h, ref_db, output_folder))

            # Collecte des résultats une fois les tâches terminées
            for future in futures:
                key, output_file_txt = future.result()
                if key not in all_results:
                    all_results[key] = {}
                all_results[key].update(output_file_txt)

        # Sauvegarder les résultats dans un fichier final avec pandas
        out_output = f'{output_folder}all_output/'
        Methods.make_folder(out_output)
        
        for key, output_file_txt in all_results.items():
            rows = []  # Liste pour stocker les lignes de résultats
            for query_id, sub_dict in output_file_txt.items():
                test = []
                for pos, list_var in sub_dict.items():
                    test.append(pos)
                    test.extend(list_var)

                if len(test) == 6:
                    rows.append([query_id, int(test[1]), int(test[2]), int(test[4]), int(test[5]), int(test[4])-int(test[2])])
                elif len(test) == 3 and test[0] == 'l':
                    rows.append([query_id, int(test[1]), int(test[2]), 'nd', 'nd', 'nd'])
                elif len(test) == 3 and test[0] == 'r':
                    rows.append([query_id, 'nd', 'nd', int(test[1]), int(test[2]), 'nd'])
                else:
                    rows.append([query_id, 'error'])

            # Convertir la liste des résultats en DataFrame et écrire dans un fichier
            df_result = pd.DataFrame(rows, columns=["query_id", "pos_1_l_subject", "pos_2_l_subject", "pos_1_r_subject", "pos_2_r_subject", "length_gene"])
            df_result.to_csv(f'{out_output}{key}_aln_output.txt', sep='\t', index=False)

    @staticmethod
    def process_resultat_file(key, f, ecoli_positif, output_folder):
        rows = []
        read = pd.read_csv(f, sep='\t')

        for i in range(len(read)):
            if read['pos_2_l_subject'][i] != 'nd' and read['pos_1_r_subject'][i] != 'nd':
                try:
                    pos_2_l = float(read['pos_2_l_subject'][i])
                    pos_1_r = float(read['pos_1_r_subject'][i])
                except ValueError:
                    continue

                # Définir 'ens1'
                if pos_2_l > pos_1_r:
                    ens1 = range(int(pos_1_r), int(pos_2_l) + 1)
                    var = 'l'
                else:
                    ens1 = range(int(pos_2_l), int(pos_1_r) + 1)
                    var = 'r'

                # Boucle à travers 'ecoli_positif'
                for j in range(len(ecoli_positif)):
                    ens2 = range(int(ecoli_positif['first_pos'][j]), int(ecoli_positif['second_pos'][j]) + 1)
                    if len(ens1) > len(ens2):
                        continue
                    
                    tab = pd.Series([item in ens1 for item in ens2]).value_counts()

                    if len(tab) == 2:
                        #print(tab)
                        #print(ecoli_positif['gene'][j])
                        if var == 'l':
                            pos_r_l = read['pos_2_l_subject'][i]
                            pos_g_l = ecoli_positif['second_pos'][j]
                            pos_r_r = read['pos_1_r_subject'][i]
                            pos_g_r = ecoli_positif['first_pos'][j]
                            res_l = float(read['pos_2_l_subject'][i]) - ecoli_positif['second_pos'][j]
                            res_r = float(read['pos_1_r_subject'][i]) - ecoli_positif['first_pos'][j]
                        else:
                            pos_r_l = read['pos_2_l_subject'][i]
                            pos_g_l = ecoli_positif['first_pos'][j]
                            pos_r_r = read['pos_1_r_subject'][i]
                            pos_g_r = ecoli_positif['second_pos'][j]
                            res_l = float(read['pos_2_l_subject'][i]) - ecoli_positif['first_pos'][j]
                            res_r = float(read['pos_1_r_subject'][i]) - ecoli_positif['second_pos'][j]

                        new_row = {
                            "query_id": read['query_id'][i],
                            "match": tab.get(True, 0),
                            "miss_match": tab.get(False, 0),
                            "gene": ecoli_positif['gene'][j],
                            "pos_l_read": pos_r_l,
                            "pos_l_gene": pos_g_l,
                            "diff_l": res_l,
                            "pos_r_read": pos_r_r,
                            "pos_r_gene": pos_g_r,
                            "diff_r": res_r
                        }
                        rows.append(new_row)

        # Convertir en DataFrame et sauvegarder le résultat
        df = pd.DataFrame(rows)
        output_file = f"{output_folder}{key}_resultats.csv"
        df.to_csv(output_file, index=False)

    @staticmethod
    def resultat(out, pos, output_folder):
        Methods.make_folder(output_folder)
        ecoli_positif = pd.read_csv(pos)

        # Utilisation de ThreadPoolExecutor pour le parallélisme
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = [executor.submit(Methods.process_resultat_file, key, f, ecoli_positif, output_folder) for key, f in out.items()]

        # Attendre que toutes les tâches soient terminées
        for future in futures:
            future.result()
