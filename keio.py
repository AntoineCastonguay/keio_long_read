import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
from keio_methode import Methods


__author__ = 'castonguaya'
__version__ = '0.1'


class Keio(object):
    def __init__(self, args):
        # I/O
        self.input = os.path.abspath(args.input)
        self.output_folder = os.path.abspath(args.output)
        self.ref_genome = os.path.abspath(args.genome)

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel
        self.mem = args.memory

        self.run()

    def run(self):
        print('Checking a few things...')
        Methods.check_input(self.input)
        Methods.check_ref(self.ref_genome)

        # Check if number of CPU and memory requested are valid
        self.cpu, self.parallel = Methods.check_cpus(self.cpu, self.parallel)
        self.mem = Methods.check_mem(self.mem)

        done_extract = self.output_folder + '/done_extract'
        done_result = self.output_folder + '/done_result'

        extract_folder = os.path.join(self.output_folder, '1_extract/')
        result_folder = os.path.join(self.output_folder, '2_result/')

        Methods.make_folder(self.output_folder)
        print('\tAll good!')

        # Extraction
        if not os.path.exists(done_extract):
            Methods.alignment(self.ref_genome, self.input, extract_folder)
            #Methods.flag_done(done_extract)
        else:
            print('Skipping extract. Already done.')

        sam_file = Methods.find_sam_files(extract_folder)

        # Manipulation
        if not os.path.exists(done_result):
            position = Methods.extract_primer_positions(sam_file[0], self.essentiel)
            Methods.write_result(position, result_folder)
            #Methods.flag_done(done_result)
        else:
            print('Skipping result. Already done.')

        print('DONE!')

if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Extract, assemble and compare Nanopore reads matching a reference sequence.')
    parser.add_argument('-g', '--genome', 
                        required=True, 
                        help='Reference genome for primer mapping. Mandatory.')
    parser.add_argument('-i', '--input', 
                        required=True, 
                        help='Folder that contains the fasta files or individual fasta file. Mandatory.') 
    parser.add_argument('-o', '--output', 
                        required=True, 
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-p', '--parallel', metavar='2',
                    required=False,
                    type=int, default=2,
                    help='Number of samples to process in parallel. Default is 2. Optional.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                    required=False,
                    type=int, default=max_cpu,
                    help='Number of threads. Default is maximum available({}). Optional.'.format(max_cpu))
    parser.add_argument('-m', '--memory', metavar=str(max_mem),
                        required=False,
                        type=int, default=max_mem,
                        help='Memory in GB. Default is 85%% of total memory ({})'.format(max_mem))
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    
    arguments = parser.parse_args()
    Ecoli(arguments)