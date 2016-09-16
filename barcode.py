#!/usr/bin/python3

import argparse
from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from collections import Counter, defaultdict
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


def parse(blast_result, **arg):
    raw = list()
    result = SearchIO.parse(blast_result, 'blast-xml')
    for query in result:
        for hit in query:
            for hsp in hit:
                hsp_query_length = hsp.query_end - hsp.query_start
                if (hsp_query_length < arg['length'] or
                        hsp.evalue > arg['evalue']):
                    continue
                line = [hsp.query, hsp.hit, hsp.query_start, hsp.hit_start]
                raw.append(line)
    return raw


def remove_multicopy(raw, **arg):
    """raw:
    hsp.query, hsp.hit, hsp.query_start, hsp.hit_start
    """
    query_info = ['{0}SPLIT{1}SPLIT{2}'.format(
        i[0].id, i[1].id, i[2]) for i in raw]
    query_info = Counter(query_info)
    to_remove = set()
    tmp = list(query_info.keys())
    tmp.sort()
    for n, key in enumerate(tmp):
        if query_info[key] != 1:
            to_remove.add(key)
    hit_info = ['{0}SPLIT{1}SPLIT{2}'.format(
        i[0].id, i[1].id, i[3]) for i in raw]
    hit_info = Counter(hit_info)
    tmp = list(hit_info.keys())
    tmp.sort()
    for n, key in enumerate(tmp):
        if hit_info[key] != 1:
            to_remove.add(key)
    if arg['strict']:
        limit = arg['sample'] ** 2
        count_info = [i[3] for i in raw]
        count_info = Counter(count_info)
        for hit in count_info.keys():
            if count_info[hit] < limit:
                to_remove.add(hit)
    singlecopy = list()
    for i in raw:
        if (i[3] not in to_remove and
                ('{0}SPLIT{1}SPLIT{2}'.format(
                    i[0].id, i[1].id, i[3]) not in to_remove and
                 '{0}SPLIT{1}SPLIT{2}'.format(
                     i[0].id, i[1].id, i[2]) not in to_remove)):
            singlecopy.append(i)
    singlecopy.sort(key=lambda i: i[3])
    return singlecopy


def extract(query_file, singlecopy):
    query_db = makeblastdb(query_file)
    hits = defaultdict(lambda: list())
    for record in singlecopy:
        hits[record[3]].append(record[0:2])
    for n, hit in enumerate(hits.keys()):
        # only use first record
        query_to_blast = hits[hit][0][0]
        query_output = 'query-{0}.fasta'.format(n)
        SeqIO.write(query_to_blast, query_output, 'fasta')
        query_blast_result = blast(query_output, query_db,
                                   query_output.replace('.fasta', 'xml'))
        query_seq = parse(query_blast_result)
        # hit only extract once
        if n != 0:
            continue
        hit_to_blast = hits[hit][0][1].seq
        hit_output = 'hit-{0}.fasta'.format(n)
        SeqIO.write(hit_to_blast, hit_output, 'fasta')
        pass


def main():
    """This program will try to find out single-copy barcode to devide
    different species in given two fasta files while ignore distinction
    among subspecies level.
    Notice that this program assuming that the sequence length of every record
    in each input fasta file has slight difference.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('sample', default=15, type=int,
                        help='sample numbers')
    parser.add_argument('-p', '--path', default='.',
                        help='target path, default is present directory')
    parser.add_argument('-d', '--db', default=None, help='''fasta file to make blast
    database, which contains longest sequence''')
    parser.add_argument('-l', '--length', default=200, type=int,
                        help='minium barcode length')
    parser.add_argument('-e', '--evalue', default=1e-20, type=float,
                        help='evalue for BLAST')
    parser.add_argument('-s', '--strict', action='store_false',
                        help='barcode location among subspecies must be same')
    arg = parser.parse_args()
    arg = vars(arg)
    fasta_files = glob(arg['path']+'/*.fasta')
    fasta_files = get_sample(fasta_files, arg['sample'])
    if arg['db'] is None:
        *query, db = find_longest(fasta_files)
    else:
        db = arg[db]
        query = set(fasta_files) - db
    db_name = makeblastdb(db)
    blast_result = list()
    for fasta in query:
        result_file = fasta.replace('.fasta', '.xml')
        blast_result = blast(fasta, db_name, result_file)
    # to be continue
    raw_result = parse(blast_result, **arg)
    singlecopy = remove_multicopy(raw_result, **arg)
    for i in singlecopy:
        print(i)


if __name__ == '__main__':
    main()
