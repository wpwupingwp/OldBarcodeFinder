#!/usr/bin/python3

import argparse
from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from collections import Counter
from functools import wraps
from glob import glob
from multiprocessing import cpu_count
from os import path, mkdir
from random import sample
from subprocess import run
from tempfile import mkdtemp
from timeit import default_timer as timer


def print_time(function):
    @wraps(function)
    def wrapper(*args, **kargs):
        start = timer()
        result = function(*args, **kargs)
        end = timer()
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
def makeblastdb(db_file):
    db_name = db_file.replace('.fasta', '')
    db_name = path.join(tmp, db_name)
    run('makeblastdb -in {0} -out {1} -logfile {2} -dbtype nucl'.format(
        db_file, db_name, db_name+'.log'), shell=True)
    return db_name


@print_time
def merge_and_split(fasta_files, target):
    merge_file = path.join(tmp, 'merge')
    count = 0
    sample_list = list()
    with open(merge_file, 'w') as merge:
        for fasta in fasta_files:
            output = path.join(tmp, '{0}-{1}'.format(
                target, path.basename(fasta)))
            raw = list(SeqIO.parse(fasta, 'fasta'))
            count += len(raw)
            target_list = sample(raw, target)
            SeqIO.write(target_list, output, 'fasta')
            sample_list.append(output)
            SeqIO.write(raw, merge, 'fasta')
    return count, merge_file, sample_list


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
def remove_multicopy(raw, is_merge):
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
    if not is_merge:
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
def extract(db, singlecopy, count, n_query):
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
        hit_sample = path.join(tmp, 'hit-{0}.fasta'.format(n))
        SeqIO.write(hit_to_blast, hit_sample, 'fasta')
        hit_blast_result = blast(hit_sample, db,
                                 hit_sample.replace('.fasta', '.xml'))
        hit_seq = parse(hit_blast_result)
        hit_seq = remove_multicopy(hit_seq, True)
        hit_seq = [i[1] for i in hit_seq]
        cover = len(hit_seq) / count
        if cover < arg.cover:
            print('''The coverage of this barcode candidate is too small
({0:.3f}), drop it.'''.format(cover))
            continue
        for record in hit_seq:
            record.id = record.id.replace(' ', '_')
            record.description = record.description.replace(' ', '_')
            record.id = '-'.join([record.id, record.description])
            record.description = ''
        barcode_output = path.join(tmp,
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


def find_barcode():
    fasta_files = glob(path.join(arg.path, '*.fasta'))
    count, merge_file, sample_list = merge_and_split(fasta_files, arg.sample)
    merge_db = makeblastdb(merge_file)
    *query, db = sample_list
    db_name = makeblastdb(db)
    for n_query, fasta in enumerate(query):
        result_file = fasta.replace('.fasta', '.xml')
        blast_result = blast(fasta, db_name, result_file)
    # to be continue
        raw_result = parse(blast_result)
        singlecopy = remove_multicopy(raw_result, False)
        barcode = extract(merge_db, singlecopy, count, n_query)
        mafft(barcode)
    return len(barcode)


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
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('-n', '--sample', default=3, type=int,
                        help='sample numbers')
    parser.add_argument('-p', '--path', default='.',
                        help='target path, default is present directory')
    parser.add_argument('-l', '--length', default=200, type=int,
                        help='minium barcode length')
    parser.add_argument('-e', '--evalue', default=1e-20, type=float,
                        help='evalue for BLAST')
    parser.add_argument('-c', '--cover', default=0.8, type=float,
                        help='coverage of barcode among query data')
    parser.add_argument('-o', '--output', default='out', help='output path')
    global arg
    arg = parser.parse_args()
    global tmp
    tmp = mkdtemp()
    if not path.exists(arg.output):
        mkdir(arg.output)
    round = 1
    while True:
        print('='*80)
        print('ROUND {0}:'.format(round))
        n_barcode = find_barcode()
        if n_barcode != 0:
            break
        print('\n\nNo good barcode found this turn, restarting now ...\n\n''')
        round += 1
    times['end'] = timer()
    print('''\n\nFinished with {0:.3f}s. You can find {1} barcodes as aligned fasta
          format in the output folder "{2}".\n'''.format(
              times['end']-times['start'], n_barcode, arg.output))


if __name__ == '__main__':
    main()
