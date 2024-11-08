import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from keio_methode import Methods


__author__ = 'castonguaya'
__version__ = '0.1'


class Keio(object):
    def __init__(self, args):
        # I/O
        self.input = os.path.abspath(args.input)
        self.output_folder = os.path.abspath(args.output)
        self.ref = os.path.abspath(args.reference)
        self.genome = os.path.abspath(args.genome)
        self.position = os.path.abspath(args.position)

        self.run()

    def run(self):
        print('Checking a few things...')
        Methods.check_input(self.input)
        Methods.check_ref(self.ref)

        done_convert = self.output_folder + '/done_convert'
        done_blast = self.output_folder + '/done_blast'
        done_extract = self.output_folder + '/done_extract'
        done_alignment = self.output_folder + '/done_alignment'
        done_resultat = self.output_folder + '/done_resultat'
        done_final = self.output_folder + '/done_final'

        convert_folder = os.path.join(self.output_folder, '1_convert/')
        blast_folder = os.path.join(self.output_folder, '2_blast/')
        extract_folder = os.path.join(self.output_folder, '3_extract/')
        alignment_folder = os.path.join(self.output_folder, '4_alignment/')
        resultat_folder = os.path.join(self.output_folder, '5_resultat/')

        Methods.make_folder(self.output_folder)
        print('\tAll good!')

        # Aller chercher les fichier fastq ou fasta du input
        extention_fasta = ['fasta', 'fasta.gz','fa', 'fa.gz']
        extension_fastq = ['fq', 'fq.gz','fastq', 'fastq.gz']
        files_fasta = Methods.list_files_in_folder(self.input, extention_fasta)
        files_fastq = Methods.list_files_in_folder(self.input, extension_fastq)
        sequence = {}

        # S'il y a des fichiers fasta ajout au dictionnaire sequence
        if files_fasta != []:
            for f in files_fasta:
                barcode = f.split('/')[-1].split('.')[0].split('_')[-1]
                sequence[barcode] = f

        # S'il y a des fichier fastq convertir en fasta
        if files_fastq != []:
            # convert
            if not os.path.exists(done_convert):
                print('Convert fastq to fasta...')
                Methods.fastq_to_fasta(files_fastq, convert_folder)
                Methods.flag_done(done_convert)
            else:
                print('Skipping convert. Already done.')
        else:
            Methods.make_folder(convert_folder)
        
        # Ajout les fichier fastq convertie dans le dictionnaire sequence
        files_convert = Methods.list_files_in_folder(convert_folder, 'fasta')
        for f in files_convert:
            barcode = f.split('/')[-1].split('.')[0]
            sequence[barcode] = f

        # blast kanamicine
        if not os.path.exists(done_blast):
            print('Blast...')
            Methods.blast(sequence, self.ref, blast_folder)
            Methods.flag_done(done_blast)
        else:
            print('Skipping blast. Already done.')

        # Ajout les fichier des résultat du blast de la kanamicine sur les reads dans le dictionnaire blast_kan
        file_blast_kan = Methods.list_files_in_folder(blast_folder + 'all_res/', 'txt')
        blast_kan = {}
        for f in file_blast_kan:
            barcode = f.split('/')[-1].split('_')[0]
            blast_kan[barcode] = f

        # extract kan
        if not os.path.exists(done_extract):
            print('Extract...')
            Methods.extract(blast_kan, sequence, extract_folder)
            Methods.flag_done(done_extract)
        else:
            print('Skipping extract. Already done.')

        # Ajout les fichier left et right apres l'extract de la kanamicine dans le dictionnaire extract_kan
        file_extract_kan = Methods.list_files_in_folder(extract_folder, 'fasta')
        extract_kan = {}
        for f in file_extract_kan:
            barcode,pos = f.split('/')[-1].split('_')[0:2]
            if barcode not in extract_kan:
                extract_kan[barcode] = {}
            extract_kan[barcode][pos] = f
        
        # alignment
        if not os.path.exists(done_alignment):
            print('Alignment...')
            Methods.blast2(extract_kan,self.genome, alignment_folder)
            Methods.flag_done(done_alignment)
        else:
            print('Skipping alignment. Already done.')

        # Ajout les fichier du resultat de l'alignement des portion left et right sur le génome dnas le dictionnaire align
        file_align = Methods.list_files_in_folder(f"{alignment_folder}all_output/", 'txt')
        align = {}
        for f in file_align:
            barcode = f.split('/')[-1].split('_')[0]
            align[barcode] = f

        # resultat
        if not os.path.exists(done_resultat):
            print('Resultat...')
            Methods.resultat(align,self.position,resultat_folder)
            Methods.flag_done(done_resultat)
        else:
            print('Skipping resultat. Already done.')

        # Ajout les résultats de l'analyse pour faire une synthese sur un document
        file_res = Methods.list_files_in_folder(f"{resultat_folder}",'csv')
        res = {}
        for f in file_res:
            barcode = f.split('/')[-1].split('_')[0]
            res[barcode] = f        
        
        # synthese
        if not os.path.exists(done_final):
            print('Output...')
            Methods.final(res,self.position,self.output_folder)

        print('DONE!')

if __name__ == "__main__":

    parser = ArgumentParser(description='blast, assemble and compare Nanopore reads matching a reference sequence.')
    parser.add_argument('-r', '--reference', 
                        required=True, 
                        help='Reference kanamycine')
    parser.add_argument('-g', '--genome', 
                        required=True, 
                        help='Reference genome')
    parser.add_argument('-p', '--position', 
                        required=True, 
                        help='Position des gènes dans un file.csv')    
    parser.add_argument('-i', '--input', 
                        required=True, 
                        help='Folder that contains the fasta files or individual fasta file.') 
    parser.add_argument('-o', '--output', 
                        required=True, 
                        help='Folder to hold the extract files.')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    
    arguments = parser.parse_args()
    Keio(arguments)