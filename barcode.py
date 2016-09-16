#!/usr/bin/python3

import argparse
from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from collections import Counter
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
    cmd = nb(num_threads=cpu_count(),
             query=query_file,
             db=db_file,
             task='blastn',
             outfmt=5,
             out=output_file)
    stdout, stderr = cmd()
    return output_file


def parse(blast_results, length, samples, evalue):
    for blast_result in blast_results:
        raw = list()
        result = SearchIO.parse(blast_result, 'blast-xml')
        for query in result:
            for hit in query:
                for hsp in hit:
                    hsp_query_length = hsp.query_end - hsp.query_start
                    if hsp_query_length < length or hsp.evalue > evalue:
                        # if hsp.bitscore < length:
                        continue
                    line = [hsp.query, hsp.hit, hsp.query_start,
                            hsp.query_end, hsp.hit_start, hsp.hit_end,
                            hsp.bitscore, hsp.evalue]
                    raw.append(line)
    return raw


def remove_multicopy(raw, length, samples):
    """raw:
    hsp.query_id, hsp.hit_id, hsp.query_start, hsp.query_end,
    hsp.hit_start, hsp.hit_end, hsp.bitscore, hsp.evalue
    """
    query_start = [i[2] for i in raw]
    query_start = Counter(query_start)
    to_remove = set()
    tmp = list(query_start.keys())
    tmp.sort()
    for n, key in enumerate(tmp):
        if query_start[key] != samples:
            to_remove.add(key)
        if n == 0:
            continue
        last_key = tmp[n-1]
        if key - last_key < length:
            to_remove.add(key)
            to_remove.add(last_key)
    raw = [i for i in raw if i[2] not in to_remove]
    hit_start = [i[4] for i in raw]
    hit_start = Counter(hit_start)
    to_remove = set()
    tmp = list(hit_start.keys())
    tmp.sort()
    for n, key in enumerate(tmp):
        if hit_start[key] != samples:
            to_remove.add(key)
        if n == 0:
            continue
        last_key = tmp[n-1]
        if key - last_key < length:
            to_remove.add(key)
            to_remove.add(last_key)
    raw = [i for i in raw if i[4] not in to_remove]
    raw.sort(key=lambda i: i[-1])
    return raw


def extract(singlecopy):
    pass


def main():
    """This program will try to find out single-copy barcode to devide
    different species while ignore distinction among subspecies level.
    Notice that this program assuming that the sequence length of every record
    in each input fasta file has little difference.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('-p', '--path', default='.',
                        help='target path, default is present directory')
    parser.add_argument('-d', '--db', default=None, help='''fasta file to make blast
    database, which contains longest sequence''')
    parser.add_argument('sample', type=int, help='sample numbers')
    parser.add_argument('-m', '--min_length', default=200, type=int,
                        help='minium barcode length')
    parser.add_argument('-e', '--evalue', default=1e-20, type=float,
                        help='evalue for BLAST')
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
    # db_name = db
    blast_result = list()
    for fasta in query:
        result_file = fasta.replace('.fasta', '.xml')
        blast_result.append(blast(fasta, db_name, result_file))
    raw_result = parse(blast_result, arg.min_length, arg.sample, arg.evalue)
    singlecopy = remove_multicopy(raw_result, arg.min_length, arg.sample)
    for i in singlecopy:
        print(i[2], i[4], i[0].id, i[1].id, i[-2], i[-1])


if __name__ == '__main__':
    main()
