import subprocess
import os
import pathlib
import gzip
from glob import glob


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
    def format_sequence(sequence, max_length=80):
        return [sequence[i:i + max_length] for i in range(0, len(sequence), max_length)]
    
    @staticmethod
    def fastq_to_fasta(fastq, output_folder):
        Methods.make_folder(output_folder)

        my_dict = {}

        file = Methods.list_files_in_folder(fastq, 'fastq.gz')
        for f in file:
            if "_barcode" in f:
                barcode = f.split('_')[-1].split('.')[0]
                output_file = output_folder + barcode + '.fasta'
                my_dict[barcode] = output_file
                with gzip.open(f, "rt") as fastq_file, open(output_file, "w") as fasta_file:
                    while True:
                        # Lire 4 lignes à la fois (une entrée FASTQ)
                        header = fastq_file.readline().strip()
                        if not header:  # Vérifier si c'est la fin du fichier
                            break
                        sequence = fastq_file.readline().strip()
                        plus_line = fastq_file.readline().strip()  # Ligne '+'
                        quality = fastq_file.readline().strip()  # Ligne de qualité

                        # Écrire en format FASTA
                        fasta_file.write(f">{header[1:]}\n")
                        formatted_sequences = Methods.format_sequence(sequence, max_length=80)
                        for seq_line in formatted_sequences:
                            fasta_file.write(f"{seq_line}\n")
        
        return my_dict

    @staticmethod
    def blast(read, ref, output_folder):
        Methods.make_folder(output_folder)
        for key,f in read.items():
            print('/t'+key)
            out = f'{output_folder}{key}/'
            Methods.make_folder(out)
            makeblastdb_cmd = ['makeblastdb', '-in', f, '-dbtype', 'nucl', '-out', f'{out}{key}_reads_db']
            subprocess.run(makeblastdb_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)

            out_res = f'{output_folder}all_res/'
            Methods.make_folder(out_res)
            blast_cmd = ['blastn', '-query', ref, '-db', f'{out}{key}_reads_db', '-out', f'{out_res}{key}_results.txt', '-outfmt', '6']
            subprocess.run(blast_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)

