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

        convert_folder = os.path.join(self.output_folder, '1_convert/')
        blast_folder = os.path.join(self.output_folder, '2_blast/')
        extract_folder = os.path.join(self.output_folder, '3_extract/')
        alignment_folder = os.path.join(self.output_folder, '4_alignment/')
        resultat_folder = os.path.join(self.output_folder, '5_resultat/')

        Methods.make_folder(self.output_folder)
        print('\tAll good!')

        # convert
        if not os.path.exists(done_convert):
            print('Convert fastq to fasta...')
            read = Methods.fastq_to_fasta(self.input, convert_folder)
            Methods.flag_done(done_convert)
        else:
            print('Skipping convert. Already done.')
            file = Methods.list_files_in_folder(convert_folder, 'fasta')
            read = {}
            for f in file:
                barcode = f.split('/')[-1].split('.')[0]
                read[barcode] = f

        # blast
        if not os.path.exists(done_blast):
            print('Blast...')
            Methods.blast(read, self.ref, blast_folder)
            Methods.flag_done(done_blast)
        else:
            print('Skipping blast. Already done.')

        file = Methods.list_files_in_folder(blast_folder + 'all_res/', 'txt')
        res = {}
        for f in file:
            barcode = f.split('/')[-1].split('_')[0]
            res[barcode] = f

        # extract kan
        if not os.path.exists(done_extract):
            print('Extract...')
            Methods.extract(res, read, extract_folder)
            Methods.flag_done(done_extract)
        else:
            print('Skipping extract. Already done.')

        file = Methods.list_files_in_folder(extract_folder, 'fasta')
        align = {}
        for f in file:
            barcode,pos = f.split('/')[-1].split('_')[0:2]
            if barcode not in align:
                align[barcode] = {}
            align[barcode][pos] = f
        
        # alignment
        if not os.path.exists(done_alignment):
            print('Alignment...')
            Methods.blast2(align,self.genome, alignment_folder)
            Methods.flag_done(done_alignment)
        else:
            print('Skipping alignment. Already done.')

        file = Methods.list_files_in_folder(f"{alignment_folder}all_output/", 'txt')
        out = {}
        for f in file:
            barcode = f.split('/')[-1].split('_')[0]
            out[barcode] = f

        # resultat
        if not os.path.exists(done_resultat):
            print('Resultat...')
            Methods.resultat(out,self.position,resultat_folder)
            Methods.flag_done(done_resultat)
        else:
            print('Skipping resultat. Already done.')

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
                        help='Folder that contains the fasta files or individual fasta file. Mandatory.') 
    parser.add_argument('-o', '--output', 
                        required=True, 
                        help='Folder to hold the extract files. Mandatory.')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    
    arguments = parser.parse_args()
    Keio(arguments)