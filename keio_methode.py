import subprocess
import os
import pathlib
import gzip
from glob import glob
from Bio import SeqIO
import csv
import pandas as pd
import concurrent.futures


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
    def list_files_in_folder(folder, extensions):
        if isinstance(extensions, str):
            extensions = [extensions]
        files = []
        for ext in extensions:
            files.extend(glob(f"{folder}/*.{ext}"))
        return files

    @staticmethod
    def flag_done(flag_file):
        with open(flag_file, 'w') as f:
            pass
    
    @staticmethod
    def fastq_to_fasta(files, output_folder):
        Methods.make_folder(output_folder)

        for f in files:
            if "_barcode" in f:
                barcode = f.split('_')[-1].split('.')[0]
                print(f'\t{barcode}')
                output_file = os.path.join(output_folder, f"{barcode}.fasta")

                #with open(output_file, 'w') as out_f:
                #    seqtk_cmd = ['seqtk', 'seq', '-A', f]
                #    subprocess.run(seqtk_cmd, stdout=out_f, stderr=subprocess.STDOUT, check=True)

                with gzip.open(f, "rt") as fastq_file:
                    with open(output_file, "w") as fasta_file:
                       for record in SeqIO.parse(fastq_file, "fastq"):
                            fasta_file.write(f">{record.id}\n")
                            fasta_file.write(f"{str(record.seq)}\n")

    @staticmethod
    def blast(sequence, ref, output_folder):
        Methods.make_folder(output_folder)
        for key,f in sequence.items():
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
    def extract(blast_kan, sequence, output_folder):
        Methods.make_folder(output_folder)

        for key,f in blast_kan.items():
            print(f'\t{key}')

            dict_match = {}

            with open(f, "r") as file:
                reader = csv.reader(file, delimiter="\t")
                for row in reader:
                    if float(row[2]) > float(90) and int(row[3]) > 300:
                        if key not in dict_match:
                            dict_match[key] = {}

                        if row[1] not in dict_match[key]:
                            dict_match[key][row[1]] = [row[3], row[8], row[9]]
                        elif row[3] > dict_match[key][row[1]][0]:
                            dict_match[key][row[1]] = [row[3], row[8], row[9]]
        
            output_file_l = os.path.join(output_folder, f"{key}_l_select_seq.fasta")
            output_file_r = os.path.join(output_folder, f"{key}_r_select_seq.fasta")

            with open(output_file_l, "w") as fasta_output, open(output_file_r, "w") as fasta_output2:
                # Parcours des séquences dans le fichier FASTA
                found = False  # Variable pour vérifier si au moins une séquence a été trouvée
                for record in SeqIO.parse(sequence[key], "fasta"):
                    if record.id in dict_match[key]:
                        # Écrire la séquence au format FASTA
                        a,b = 1,2
                        if dict_match[key][record.id][a] > dict_match[key][record.id][b]:
                            b,a = 1,2
                        pos_first = int(dict_match[key][record.id][a])
                        pos_second = int(dict_match[key][record.id][b])
                                        
                        fasta_output.write(f">{record.id}\n")
                        if pos_first > 300:
                            fasta_output.write(f"{record.seq[pos_first-300:pos_first]}\n")
                        else:
                            fasta_output.write(f"{record.seq[1:pos_first]}\n")
                        fasta_output2.write(f">{record.id}\n")
                        if len(record.seq)-pos_second > 300:
                            fasta_output2.write(f"{record.seq[pos_second:pos_second+300]}\n")
                        else:
                            fasta_output2.write(f"{record.seq[pos_second:len(record.seq)]}\n")
                        found = True

                if not found:
                    print(f"\tError : {key}_selected_sequences.fasta")

    @staticmethod
    def blast2(align, ref, output_folder):
        Methods.make_folder(output_folder)

        makeblastdb_cmd = ['makeblastdb', '-in', ref, '-dbtype', 'nucl', '-out', f'{output_folder}ref_db']
        subprocess.run(makeblastdb_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)
        
        for key,f in align.items():
            align_dict = {}
            for pos,h in f.items():
                print(f'\t{key} {pos}')

                out_res = f'{output_folder}all_res/'
                Methods.make_folder(out_res)
                blast_cmd = ['blastn', '-query', h, '-db', f'{output_folder}ref_db', '-out', f'{out_res}{key}_{pos}_results.txt', '-outfmt', '6', '-max_hsps', '1']
                subprocess.run(blast_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)
                
                with open(f"{out_res}{key}_{pos}_results.txt", 'r') as file:
                    reader = csv.reader(file, delimiter="\t")
                    for row in reader:
                        query_id = row[0]
                        pos_1_subject = int(row[8])
                        pos_2_subject = int(row[9])
                        
                        list_var = [pos_1_subject,pos_2_subject]

                        if query_id not in align_dict:
                            align_dict[query_id] = {}
                        align_dict[query_id][pos] = list_var
            
            out_output = f'{output_folder}all_output/'
            Methods.make_folder(out_output)
            with open(f'{out_output}{key}_aln_output.txt', 'w') as f:
                f.write(f"query_id\tpos_1_l_subject\tpos_2_l_subject\tpos_1_r_subject\tpos_2_r_subject\tlenth_gene\n")
                for query_id, sub_dict in align_dict.items():
                    row_align = []
                    for pos, list_var in sub_dict.items():
                        row_align.append(pos)
                        row_align.extend(list_var)
                    if len(row_align) == 6 and row_align[0] == 'l':
                        f.write(f"{query_id}\t{int(row_align[1])}\t{int(row_align[2])}\t{int(row_align[4])}\t{int(row_align[5])}\t{int(row_align[4])-int(row_align[2])}\n")
                    elif len(row_align) == 6 and row_align[0] == 'r':
                        f.write(f"{query_id}\t{int(row_align[4])}\t{int(row_align[5])}\t{int(row_align[1])}\t{int(row_align[2])}\t{int(row_align[5])-int(row_align[1])}\n")
                    elif len(row_align) == 3 and row_align[0] == 'l':
                        f.write(f"{query_id}\t{row_align[1]}\t{row_align[2]}\tnd\tnd\tnd\n")
                    elif len(row_align) == 3 and row_align[0] == 'r':
                        f.write(f"{query_id}\tnd\tnd\t{row_align[1]}\t{row_align[2]}\tnd\n")
                    else:
                        f.write(f"{query_id} error \n")
         
    @staticmethod
    def run_res(key, align, ecoli_positif, output_folder):
        rows = []  # Utilisation d'une liste temporaire pour stocker les lignes
        df_align = pd.read_csv(align, sep='\t')
        for i in range(len(df_align)):
            w = 0
            if df_align['pos_2_l_subject'][i] != 'nd' and df_align['pos_1_r_subject'][i] != 'nd':
                try:
                    pos_2_l = float(df_align['pos_2_l_subject'][i])
                    pos_1_r = float(df_align['pos_1_r_subject'][i])
                    w = 0
                except ValueError:
                    continue
            elif df_align['pos_2_l_subject'][i] == 'nd':
                try:
                    pos_1_r = float(df_align['pos_1_r_subject'][i])
                    if float(df_align['pos_1_r_subject'][i]) > float(df_align['pos_2_r_subject'][i]):
                        pos_2_l = float(df_align['pos_1_r_subject'][i]) + 40.
                    else:
                        pos_2_l = float(df_align['pos_1_r_subject'][i]) - 40
                    w = 1
                except ValueError:
                    continue
            elif df_align['pos_1_r_subject'][i] == 'nd':
                try:
                    pos_2_l = float(df_align['pos_2_l_subject'][i])
                    if float(df_align['pos_2_l_subject'][i]) > float(df_align['pos_1_l_subject'][i]):
                        pos_1_r = float(df_align['pos_2_l_subject'][i]) + 40
                    else:
                        pos_1_r = float(df_align['pos_2_l_subject'][i]) - 40
                    w = 2
                except ValueError:
                    continue
            else:
                continue

            #print(key)
            #print(df_align['query_id'][i])            

            # Vérifier et créer 'ens1'
            if pos_2_l > pos_1_r:
                ens1 = range(int(pos_1_r), int(pos_2_l) + 1)
                var = 'l'
            else:
                ens1 = range(int(pos_2_l), int(pos_1_r) + 1)
                var = 'r'

            # Boucle à travers 'ecoli_positif'
            for j in range(len(ecoli_positif)):
                ens2 = range(int(ecoli_positif['first_pos'][j]), int(ecoli_positif['second_pos'][j]) + 1)

                # Vérifier la longueur des séquences
                if len(ens1) > len(ens2):
                    continue

                # Créer une table avec un comptage des correspondances et des erreurs
                tab = pd.Series([item in ens1 for item in ens2]).value_counts()

                if len(tab) == 2:
                    #print(tab)
                    #print(ecoli_positif['gene'][j])
                    # #if tab[True] > abs(int(df_align['lenth_gene'][i]))*0.8:
                    if var == 'l':
                        if w == 1:
                            pos_r_l, res_l = 0, 0
                        else: 
                            pos_r_l = df_align['pos_2_l_subject'][i]
                            res_l = float(df_align['pos_2_l_subject'][i]) - ecoli_positif['second_pos'][j]

                        if w == 2:
                            pos_r_r, res_r = 0, 0
                        else:
                            pos_r_r = df_align['pos_1_r_subject'][i]
                            res_r = float(df_align['pos_1_r_subject'][i]) - ecoli_positif['first_pos'][j]                        

                        pos_g_l = ecoli_positif['second_pos'][j]
                        pos_g_r = ecoli_positif['first_pos'][j]
                    else:
                        if w == 1:
                            pos_r_l, res_l = 0, 0
                        else: 
                            pos_r_l = df_align['pos_2_l_subject'][i]
                            res_l = float(df_align['pos_2_l_subject'][i]) - ecoli_positif['first_pos'][j]

                        if w == 2:
                            pos_r_r, res_r = 0, 0
                        else:
                            pos_r_r = df_align['pos_1_r_subject'][i]
                            res_r = float(df_align['pos_1_r_subject'][i]) - ecoli_positif['second_pos'][j]                        

                        pos_g_l = ecoli_positif['first_pos'][j]
                        pos_g_r = ecoli_positif['second_pos'][j]

                    # Ajouter une ligne au DataFrame
                    if var == 'l':
                        new_row = {
                            "query_id": df_align['query_id'][i],
                            "match": tab.get(True, 0),
                            "offset": tab.get(False, 0),
                            "gene": ecoli_positif['gene'][j],
                            "pos_l_read": pos_r_r,
                            "pos_l_gene": pos_g_r,
                            "diff_l": res_r,
                            "pos_r_read": pos_r_l,
                            "pos_r_gene": pos_g_l,
                            "diff_r": res_l
                        }    
                    else:
                        new_row = {
                            "query_id": df_align['query_id'][i],
                            "match": tab.get(True, 0),
                            "offset": tab.get(False, 0),
                            "gene": ecoli_positif['gene'][j],
                            "pos_l_read": pos_r_l,
                            "pos_l_gene": pos_g_l,
                            "diff_l": res_l,
                            "pos_r_read": pos_r_r,
                            "pos_r_gene": pos_g_r,
                            "diff_r": res_r
                        }                                            

                    rows.append(new_row)

        # Convertir la liste de lignes en DataFrame
        df = pd.DataFrame(rows)

        # Sauvegarder dans un fichier CSV
        output_file = os.path.join(output_folder, f"{key}_resultats.csv")
        df.to_csv(output_file, index=False)

    def resultat(align, position, output_folder):
        Methods.make_folder(output_folder)
        ecoli_positif = pd.read_csv(position)

        # Utilisation de ThreadPoolExecutor pour paralléliser la boucle
        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as executor:
            futures = []
            for key, f in align.items():
                print(f'\t{key}')
                # Exécution en parallèle de la fonction process_file pour chaque clé
                futures.append(executor.submit(Methods.run_res, key, f, ecoli_positif, output_folder))

            # Attendre que toutes les tâches soient terminées
            concurrent.futures.wait(futures)

    def final(file, pos_ref, output_folder): 
        # Créer une liste pour stocker les lignes du fichier final
        rows = []
        
        for key, f in file.items():
            # Lire le fichier CSV
            df = pd.read_csv(f)
            
            # Compter les occurrences et calculer les pourcentages pour la colonne 'gene'
            compte = df["gene"].value_counts()
            pourcentage = (compte / compte.sum()) * 100
            
            # Créer un dictionnaire avec 'key' et les colonnes gene_1, gene_1_pourcentage, ...
            row = {'key': key}
            
            # Ajouter les gènes et leurs pourcentages dans les colonnes
            for i, (gene, pct) in enumerate(pourcentage.items(), start=1):
                if i > 10:  # Limiter à un maximum de 10 gènes
                    break
                row[f'gene_{i}'] = gene
                row[f'gene_{i}_pourcentage'] = pct
            
            # Ajouter la ligne au tableau final
            rows.append(row)
        
        # Convertir les lignes en DataFrame
        result_df = pd.DataFrame(rows)
        
        # Écrire le DataFrame dans le fichier CSV de sortie
        output_file = os.path.join(output_folder, f"resultats_all.csv")
        result_df.to_csv(output_file, index=False, sep=';')
