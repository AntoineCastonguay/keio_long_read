import subprocess
import os
import pathlib
import gzip
from glob import glob
from Bio import SeqIO
import csv


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
    def fastq_to_fasta(fastq_folder, output_folder):
        Methods.make_folder(output_folder)

        my_dict = {}

        files = Methods.list_files_in_folder(fastq_folder, 'fastq.gz')
        for f in files:
            if "_barcode" in f:
                barcode = f.split('_')[-1].split('.')[0]
                print(f'\t{barcode}')
                output_file = os.path.join(output_folder, f"{barcode}.fasta")
                my_dict[barcode] = output_file

                with gzip.open(f, "rt") as fastq_file:
                    with open(output_file, "w") as fasta_file:
                        for record in SeqIO.parse(fastq_file, "fastq"):
                            fasta_file.write(f">{record.id}\n")
                            fasta_file.write(f"{str(record.seq)}\n")
        
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
            blast_cmd = ['blastn', '-query', ref, '-db', f'{out}{key}_reads_db', '-out', f'{out_res}{key}_results.txt', '-outfmt', '6']
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
    def blast2(align, ref, output_folder):
        Methods.make_folder(output_folder)

        makeblastdb_cmd = ['makeblastdb', '-in', ref, '-dbtype', 'nucl', '-out', f'{output_folder}ref_db']
        subprocess.run(makeblastdb_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)
        
        for key,f in align.items():
            for pos,h in f.items():
                print(f'\t{key} {pos}')

                out_res = f'{output_folder}all_res/'
                Methods.make_folder(out_res)
                blast_cmd = ['blastn', '-query', h, '-db', f'{output_folder}ref_db', '-out', f'{out_res}{key}_{pos}_results.txt', '-outfmt', '6']
                subprocess.run(blast_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)


    @staticmethod
    def alignment(align,ref, output_folder):
        Methods.make_folder(output_folder)

        # Alignment BWA
        BWA_index_cmd = ['bwa', 'index', ref]
        subprocess.run(BWA_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        for key,f in align.items():
            for pos,h in f.items():
                if pos == 'l':
                    left = h
                elif pos == 'r':
                    right = h
                else:
                    print('Error')
            print(f'\t{key}')
            BWA_cmd = ['bwa', 'mem', ref, left, right]
            with open(f'{output_folder}/{key}BWA.sam', 'w') as outfile, open(os.devnull, 'w') as errfile:
                    subprocess.run(BWA_cmd, stdout=outfile, stderr=errfile)