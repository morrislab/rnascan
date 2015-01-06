import sys,csv

from math import log


def read_pfm(pfmfile):
	pfm = {}
	with open(pfmfile) as f:
		reader = csv.reader(f, delimiter='\t')
		headerline = reader.next()
		alphabet = headerline[1:]
		for base in alphabet:
			pfm[base] = []
		for row in reader:
			for i in range(len(row)):
				if i==0:
					continue
				pfm[alphabet[i-1]].append(float(row[i]))
	
	f.close()
	return(pfm)

def write_pfm(pfm,pfmoutfile):
	of = open(pfmoutfile, 'w')
		
	alphabet = sorted(pfm.keys())
	pfm_length = len(pfm[alphabet[0]])
	
	of.write('PO')
	for base in alphabet:
		of.write('\t' + base)
	of.write('\n')

	for pos in range(0,pfm_length):
		if pfm[alphabet[0]][pos] is None:
			break
		of.write(str(pos))
		for base in alphabet:
			of.write('\t' + str(pfm[base][pos]))
		of.write('\n')

	of.close()

def norm_pfm(pfm):
	alphabet = sorted(pfm.keys())
	pfm_length = len(pfm[alphabet[0]])

	pfm_normalized = {}
	for base in alphabet:
		pfm_normalized[base] = [None] * (pfm_length)

	for pos in range(0,pfm_length):
		sum = 0
		for base in alphabet:
			sum = sum + pfm[base][pos]
		for base in alphabet:
			normalized = pfm[base][pos] / float(sum)
			pfm_normalized[base][pos] = normalized
	return pfm_normalized


def pfm_to_pwm(pfm, num_sites):
	# assume equal background model for all bases (ie .25/.25/.25/.25 for nucleotides)
	alphabet = sorted(pfm.keys())
	num_bases = len(alphabet)
	
	pfm_length = len(pfm[alphabet[0]])
	
	pwm = {}
	for base in alphabet:
		pwm[base] = [None] * (pfm_length)
	
	for pos in range(0,pfm_length):
		for base in alphabet:
			corrected_p = (float(pfm[base][pos])*num_sites + 1./num_bases) / (num_sites+1)
			pwm[base][pos] = log(corrected_p / (1./num_bases), 2)
	
	return pwm

def pwm_scan_fwd(pwm,seq):
	#score forward strand only
	alphabet = sorted(pwm.keys())
	pwm_length = len(pwm[alphabet[0]])
	
	pos = 0
	FwdScores = []
	
	while pos < (len(seq) - pwm_length + 1):
			score = 0
			subseq = seq[pos:(pos + pwm_length)]
			for subpos in range(0, len(subseq)):
					base = subseq[subpos]
					score = score + pwm[base][subpos]
			FwdScores.append(score)
			pos = pos + 1
	return FwdScores
