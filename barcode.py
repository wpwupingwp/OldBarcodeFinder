#!/usr/bin/python3

import argparse
from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from collections import defaultdict
from glob import glob
from multiprocessing import cpu_count
from subprocess import run


def find_longest(fasta_files):
    avg_length = list()
    for fasta in fasta_files:
        length = list()
        raw = SeqIO.parse(fasta, 'fasta')
        for sequence in raw:
            length.append(len(sequence))
        avg_length.append([fasta, sum(length)/len(length)])
    avg_length.sort(key=lambda i: i[1])
    return [i[0] for i in avg_length]


def get_sample(fasta_files, target):
    new_fasta_files = list()
    for fasta in fasta_files:
        output = fasta.replace('.fasta', '_{0}.fasta'.format(target))
        raw = SeqIO.parse(fasta, 'fasta')
        with open(output, 'w') as output_file:
            for n in range(target):
                SeqIO.write(next(raw), output_file, 'fasta')
        new_fasta_files.append(output)
    return new_fasta_files


def makeblastdb(db_file):
    db_name = db_file.replace('.fasta', '')
    run('makeblastdb -in {0} -title {1} -out {1} -dbtype nucl'.format(
        db_file, db_name), shell=True)
    return db_name


def blast(query_file, db_file, output_file='BLASTResult.xml'):
    """Here we use "max_hsps" to restrict only first hsp.
    """
    cmd = nb(num_threads=cpu_count(),
             query=query_file,
             db=db_file,
             task='blastn',
             outfmt=5,
             out=output_file)
    stdout, stderr = cmd()
    return output_file


def parse(blast_results, limit):
    handle = open('stats.tmp', 'w')
    for blast_result in blast_results:
        result = SearchIO.parse(blast_result, 'blast-xml')
        handle.write(
            'query_id,hit_id,query_start,query_end,hit_start,hit_end,score\n')
        for query in result:
            for hit in query:
                for hsp in hit:
                    if hsp.bitscore >= limit:
                        yield hsp.query_start, hsp.query_end, hsp.bitscore
                        line = [hsp.query_id, hsp.hit_id, hsp.query_start,
                                hsp.query_end, hsp.hit_start, hsp.hit_end,
                                hsp.bitscore]
                        line = [str(_) for _ in line]
                        handle.write(','.join(line)+'\n')
        handle.write(
            '#############################################################\n')


def main():
    """This program will try to find out single-copy barcode to devide
    different species while ignore distinction among subspecies level.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('-p', '--path', default='.',
                        help='target path, default is present directory')
    parser.add_argument('-d', '--db', default=None, help='''fasta file to make blast
    database, which contains longest sequence''')
    parser.add_argument('-s', '--sample', default=5, type=int,
                        help='sample numbers')
    parser.add_argument('-m', '--min_length', default=200, type=int,
                        help='minium barcode length')
    parser.print_help()
    arg = parser.parse_args()
    fasta_files = glob(arg.path+'/*.fasta')
    if arg.sample is not None:
        fasta_files = get_sample(fasta_files, arg.sample)
    if arg.db is None:
        *query, db = find_longest(fasta_files)
    else:
        db = arg.db
        query = set(fasta_files) - db
    db_name = makeblastdb(db)
    blast_result = list()
    for fasta in query:
        result_file = fasta.replace('.fasta', '.xml')
        blast_result.append(blast(fasta, db_name, result_file))
    count = defaultdict(lambda: 0)
    count_2 = defaultdict(lambda: 0)
    for line in parse(blast_result, arg.min_length):
        for n in range(line[0], line[1]+1):
            count[n] += line[2]
            count[n] += 1
    print(len(count_2))
    for i in count_2.keys():
        print(i, count[i])


if __name__ == '__main__':
    main()
