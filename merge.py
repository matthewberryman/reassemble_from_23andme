#!/usr/local/bin/python3
import sys

if (len(sys.argv)!=3):
    print('Usage: python merge.py genome_file /path/to/HCR_fasta_files/')
    print('Requires files from ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRC.p13_chr*.fa.gz to be downloaded and gunzipped to the location specified in the second command line parameter.')
    quit(0)

from Bio import SeqIO, SeqRecord

# Store the two copies as 'a' and 'b'

a_genotype = {'X': {}, 'Y': {}, 'MT': {}}
b_genotype = {'X': {}}

path = str(sys.argv[2])

for i in range(1,23):
    a_genotype[str(i)]={}
    b_genotype[str(i)]={}

# read 23andme data
with open(str(sys.argv[1])) as f:
    for line in f:
        if (not line.startswith('#')):
            (rsid, chromosome, position, genotype) = tuple(line.split(None,4))
            
            a_genotype[str(chromosome)][position] = genotype[0]
            if (not (chromosome == 'MT' or chromosome == 'Y')):
                try:
                    b_genotype[str(chromosome)][position] = genotype[1]
                except IndexError:
                    if chromosome != 'X':
                        raise

male = 'Y' in a_genotype
# we then use this to throw away the PAR regions listed against the X chromosome, since it's too hard to
# reconstruct the actual position on the Y from the listed position on the X.

# Also discard the Y chromosome key from females so we don't iterate over it later.
if (not male):
    a_genotype.pop('Y')

for chromosome in a_genotype.keys():
    x = path + 'hs_ref_GRCh37.p13_chr' + chromosome + '.fa'
    a = path + 'hs_ref_GRCh37.p13_chr' + chromosome + '_a.fa'
    a_hdl = open(a,'w')
    new_record_a = ''
    if (not (chromosome == 'MT' or chromosome == 'Y')):
        b = path + 'hs_ref_GRCh37.p13_chr' + chromosome + '_b.fa'
        b_hdl = open(b,'w')
        new_record_b = ''

    
    for seq_record in SeqIO.parse(x,'fasta'):
        # Convert them to MutableSeq objects so that we can edit them:
        seq_a=seq_record.seq.tomutable()
        seq_b=seq_record.seq.tomutable()
        for key in a_genotype[chromosome].keys():
            # Seq objects use indexing starting from 0, but biology convention in 23andme file starts from 1, so use int(key)-1 to convert

            seq_a[int(key)-1] = a_genotype[chromosome][key]

            #quit()
            if (not (chromosome == 'MT' or chromosome == 'Y')):
                try:
                    seq_b[int(key)-1] = b_genotype[chromosome][key]
                except:
                    # Some out of bound errors still occur for X chromosome, I guess related to the PAR region
                    print ('OOB in Chr ' +chromosome + ': ' + key + ' > ' + str(seq_b.__len__()))
        new_record_a = SeqRecord.SeqRecord(seq_a,id=seq_record.id,name=seq_record.name,description=seq_record.description)
        if (not (chromosome == 'MT' or chromosome == 'Y')):
            new_record_b = SeqRecord.SeqRecord(seq_b,id=seq_record.id,name=seq_record.name,description=seq_record.description)

    SeqIO.write([new_record_a],a_hdl,'fasta')
    # last boolean here throws away PAR data as discussed above:
    if (not (chromosome == 'MT' or chromosome == 'Y' or (male and chromosome=='X'))):
        SeqIO.write([new_record_b],b_hdl,'fasta')
    a_hdl.close()

