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

        self.run()

    def run(self):
        print('Checking a few things...')
        Methods.check_input(self.input)
        Methods.check_ref(self.ref)

        done_convert = self.output_folder + '/done_convert'
        done_blast = self.output_folder + '/done_blast'
        done_result = self.output_folder + '/done_result'

        convert_folder = os.path.join(self.output_folder, '1_convert/')
        blast_folder = os.path.join(self.output_folder, '2_blast/')
        result_folder = os.path.join(self.output_folder, '3_result/')

        Methods.make_folder(self.output_folder)
        print('\tAll good!')

        # convert
        if not os.path.exists(done_convert):
            print('Convert fastq to fasta...')
            read = Methods.fastq_to_fasta(self.input, convert_folder)
            #Methods.flag_done(done_convert)
        else:
            print('Skipping convert. Already done.')

        # blast
        if not os.path.exists(done_blast):
            print('Blast...')
            Methods.blast(read, self.ref, blast_folder)
            #Methods.flag_done(done_blast)
        else:
            print('Skipping blast. Already done.')

        print('DONE!')

if __name__ == "__main__":

    parser = ArgumentParser(description='blast, assemble and compare Nanopore reads matching a reference sequence.')
    parser.add_argument('-r', '--reference', 
                        required=True, 
                        help='Reference reference for primer mapping. Mandatory.')
    parser.add_argument('-i', '--input', 
                        required=True, 
                        help='Folder that contains the fasta files or individual fasta file. Mandatory.') 
    parser.add_argument('-o', '--output', 
                        required=True, 
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    
    arguments = parser.parse_args()
    Keio(arguments)