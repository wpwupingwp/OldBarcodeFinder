#!/usr/bin/python3

import argparse
import sys
from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from collections import Counter
from functools import wraps
from glob import glob
from multiprocessing import cpu_count
from os import path, mkdir
from random import shuffle
from subprocess import run
from timeit import default_timer as timer


def print_time(function):
    @wraps(function)
    def wrapper(*args, **kargs):
        start = timer()
        result = function(*args, **kargs)
        end = timer()
        print('='*80)
        print('The function {0} costed {1:.3f}s.'.format(
            function.__name__, end-start))
        return result
    return wrapper


@print_time
def check_dependence():
    """Check dependent programs by printing version info of them. Since
    BLAST suite does not use traditional "--version", here are two loops
    to handle the difference.
    """
    for program in ('makeblastdb', 'blastn'):
        check = run('{0} -version'.format(program), shell=True)
        if check.returncode != 0:
            raise Exception('{0} missing, see README for install'.format(
                program))
    for program in ('mafft', ):
        check = run('{0} --version'.format(program), shell=True)
        if check.returncode != 0:
            print(check)
            raise Exception('{0} missing, see README for install'.format(
                program))


@print_time
def merge_fasta(fasta_files):
    merge_file = path.join(arg.path, 'merge')
    with open(merge_file, 'w') as merge:
        for fasta in fasta_files:
            with open(fasta, 'r') as f:
                merge.write(f.read())
    merge_db = makeblastdb(merge_file)
    return merge_db


@print_time
def get_sample(fasta, target):
    output = path.join(arg.tempdir,
                       '{0}-{1}'.format(target, path.basename(fasta)))
    raw = SeqIO.index(fasta, 'fasta')
    target_list = list(raw.keys())
    shuffle(target_list)
    target_list = target_list[:arg.sample]
    with open(output, 'w') as output_file:
        for index in target_list:
            SeqIO.write(raw[index], output_file, 'fasta')
    return output


@print_time
def makeblastdb(db_file):
    db_name = db_file.replace('.fasta', '')
    run('makeblastdb -in {0} -out {1} -logfile {2} -dbtype nucl'.format(
        db_file, db_name, db_name+'.log'), shell=True)
    return db_name


@print_time
def blast(query_file, db_file, output_file='BLASTResult.xml'):
    cmd = nb(num_threads=cpu_count(),
             query=query_file,
             db=db_file,
             task='blastn',
             outfmt=5,
             out=output_file)
    stdout, stderr = cmd()
    return output_file


@print_time
def parse(blast_result):
    raw = list()
    result = SearchIO.parse(blast_result, 'blast-xml')
    for query in result:
        for hit in query:
            for hsp in hit:
                hsp_query_length = hsp.query_end - hsp.query_start
                if (hsp_query_length < arg.length or
                        hsp.evalue > arg.evalue):
                    continue
                line = [hsp.query, hsp.hit, hsp.query_start, hsp.hit_start,
                        hsp.query_end, hsp.hit_end]
                raw.append(line)
    return raw


@print_time
def remove_multicopy(raw):
    """raw:
    hsp.query, hsp.hit, hsp.query_start, hsp.hit_start, hsp.query_end,
    hsp.hit_end
    """
    query_info = ['{0}SPLIT{1}SPLIT{2}'.format(
        i[0].id, i[1].id, i[2]) for i in raw]
    query_info = Counter(query_info)
    to_remove = set()
    for key in query_info.keys():
        if query_info[key] != 1:
            to_remove.add(key)
    hit_info = ['{0}SPLIT{1}SPLIT{2}'.format(
        i[0].id, i[1].id, i[3]) for i in raw]
    hit_info = Counter(hit_info)
    for key in hit_info.keys():
        if hit_info[key] != 1:
            to_remove.add(key)
    if arg.strict:
        limit = arg.sample ** 2
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


@print_time
def extract(db, singlecopy, n_query):
    """Extract barcode sequence with BLAST against merged query files to
    ensure the validation of the barcode.
    """
    barcode = list()
    hits = dict()
    for record in singlecopy:
        hits[record[3]] = record[0]
    n = 1
    for hit in hits.keys():
        # only use first record
        hit_to_blast = hits[hit]
        hit_sample = path.join(arg.tempdir, 'hit-{0}.fasta'.format(n))
        SeqIO.write(hit_to_blast, hit_sample, 'fasta')
        hit_blast_result = blast(hit_sample, db,
                                 hit_sample.replace('.fasta', '.xml'))
        hit_seq = parse(hit_blast_result)
        hit_seq = [i[1] for i in hit_seq]
        for record in hit_seq:
            record.id = record.id.replace(' ', '_')
            record.description = record.description.replace(' ', '_')
            record.id = '-'.join([record.id, record.description])
            record.description = ''
        barcode_output = path.join(arg.tempdir,
                                   'barcode-{0}-{1}.fasta'.format(n_query, n))
        SeqIO.write(hit_seq, barcode_output, 'fasta')
        barcode.append(barcode_output)
        n += 1
    return barcode


@print_time
def mafft(barcode_file):
    barcode_aln = list()
    for barcode in barcode_file:
        aln_file = barcode.replace('.fasta', '.aln')
        aln_file = path.join(arg.output, path.basename(aln_file))
        run('mafft --quiet --thread {0} {1} > {2}'.format(
            cpu_count(), barcode, aln_file), shell=True)
        barcode_aln.append(aln_file)
    return barcode_aln


def main():
    """This program will try to find out single-copy barcode to devide
    different species in given two fasta files while ignore distinction
    among subspecies level.
    Notice that this program assuming that the sequence length of every record
    in each input fasta file has slight difference.
    """
    times = dict()
    times['start'] = timer()
    check_dependence()
    sys.stderr = open('logfile.txt', 'w')
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('-n', '--sample', default=3, type=int,
                        help='sample numbers')
    parser.add_argument('-p', '--path', default='.',
                        help='target path, default is present directory')
    parser.add_argument('-l', '--length', default=200, type=int,
                        help='minium barcode length')
    parser.add_argument('-e', '--evalue', default=1e-20, type=float,
                        help='evalue for BLAST')
    parser.add_argument('-s', '--strict', action='store_false',
                        help='barcode location among subspecies must be same')
    parser.add_argument('-o', '--output', default='out', help='output path')
    parser.add_argument('-t', '--tempdir', default='tmp',
                        help='temp file directory')
    global arg
    arg = parser.parse_args()
    if not path.exists(arg.tempdir):
        mkdir(arg.tempdir)
    if not path.exists(arg.output):
        mkdir(arg.output)
    fasta_files = glob(path.join(arg.path, '*.fasta'))
    fasta_files = [get_sample(i, arg.sample) for i in fasta_files]
    merge_db = merge_fasta(fasta_files)
    *query, db = fasta_files
    db_name = makeblastdb(db)
    for n_query, fasta in enumerate(query):
        result_file = fasta.replace('.fasta', '.xml')
        blast_result = blast(fasta, db_name, result_file)
    # to be continue
        raw_result = parse(blast_result)
        singlecopy = remove_multicopy(raw_result)
        barcode = extract(merge_db, singlecopy, n_query)
        mafft(barcode)
    times['end'] = timer()
    print('''\n\nFinished with {0:.3f}s. You can find barcodes as aligned fasta
          format in the output folder "{1}".\n'''.format(
              times['end']-times['start'], arg.output))


if __name__ == '__main__':
    main()
