### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.31 ###

__usage__ = """
					python contigs2pseudochromosomes.py
					--in <CONTIG_FILE>
					--ref <REFERENCE_GENOME_SEQ_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					--gap <INT, size of gaps to place between contigs>[100]
					"""
import os, sys
from operator import itemgetter

# --- end of imports --- #

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_best_blast_hit( blast_result_file ):
	"""! @brief load best blast hit per query """
	
	best_hits = {}
	
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				data = best_hits[ parts[0] ]
				if float( parts[-1] ) > data['score']:
					del best_hits[ parts[0] ]
					best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'line': line } } )
			except:
				best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'line': line } } )
			line = f.readline()
	return best_hits


def get_best_hit_per_chr( tmp ):
	"""! @brief get best hit per chr """
	
	# --- sort data by chromosome --- #
	hits = {}
	best_line_per_chr = {}
	for each in tmp:
		parts = each.strip().split('\t')
		try:
			hits[ parts[1] ] += float( parts[-1] )
		except KeyError:
			hits.update( { parts[1]: float( parts[-1] ) } )
			best_line_per_chr.update( { parts[1]: each } )
	
	# --- find best chromosome --- #
	best_chr = False
	best_score = 0
	for hit in hits.keys():
		if hits[ hit ] > best_score:
			best_chr = hit
			best_score = hits[ hit ]
	return { 'score': float( best_line_per_chr[ best_chr ].strip().split('\t')[-1] ), 'line': best_line_per_chr[ best_chr ] }


def load_best_blast_hit2( result_file ):
	"""! @brief load best chromosome hit """
	
	sim_cutoff = 95
	len_cutoff = 20000
	
	best_hits = {}
	
	with open( result_file, "r" ) as f:
		
		line = f.readline()
		prev_query = line.split('\t')[0]
		tmp = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_query:
				if len( tmp ) > 0:
					best_hits.update( { prev_query: get_best_hit_per_chr( tmp ) } )
				tmp = [ line ]
				prev_query = parts[0]
			else:
				if float( parts[2] ) > sim_cutoff:
					if int( parts[3] ) > len_cutoff:
						tmp.append( line )
			line = f.readline()
		if len( tmp ) > 0:
			best_hits.update( { prev_query: get_best_hit_per_chr( tmp ) } )
	return best_hits


def load_seq_lengths( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	seq_lens = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split( " " )[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					seq_lens.update( { header: len( "".join( seq ) ) } )
					header = line.strip()[1:].split( " " )[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		seq_lens.update( { header: len( "".join( seq ) ) } )	
	return seq_lens


def get_pos_per_contig( best_hit_lines, contig_lens, chr_lens ):
	"""! @brief get position of each contig on reference """
	
	contig_positions = []
	best_hits = []
	for line in best_hit_lines.values():
		best_hits.append( line['line'] )
	for hit in best_hits:
		parts = hit.strip().split('\t')
		length = contig_lens[ parts[0] ]
		qstart, qend = int( parts[6] ), int( parts[7] )
		sstart, send = int( parts[8] ), int( parts[9] )
		if sstart < send:
			sstart = sstart - qstart
			send = send + ( length - qend )
			orientation = True
		else:
			sstart = sstart + qstart
			send = send - ( length - qstart )
			orientation = False
		contig_positions.append( { 'ID': parts[0], 'chr': parts[1],
													'start': max( [ 1, min( [ sstart, send ] ) ] ),
													'end': min( [ max( [ sstart, send ] ), chr_lens[ parts[1] ] ] ),
													'contig_len': contig_lens[ parts[0] ],
													'orientation': orientation
													} )
	return contig_positions


def calculate_overlap( sequences, contig ):
	"""! @brief calculate overlap length with all contigs """
	
	counter = 0
	for seq in sequences:
		if seq['start'] < contig['end']:
			if seq['end'] > contig['start']:
				pos = sorted( [ seq['start'], seq['end'], contig['start'], contig['end'] ] )
				counter += pos[2]-pos[1]
	return counter


def get_contig_cumlen( contig_set ):
	"""! @brief calculate cumulative length of all contigs in list """
	
	counter = 0
	for contig in contig_set:
		counter += contig['contig_len']
	return counter


def group_contigs_to_haplophases( pos_per_contig, chr_lens ):
	"""! @brief group all contigs by haplophases """
	
	# --- sort contigs by chromosome --- #
	contigs_by_chromosomes = {}
	for contig in pos_per_contig:
		try:
			contigs_by_chromosomes[ contig['chr'] ].append( contig )
		except KeyError:
			contigs_by_chromosomes.update( { contig['chr']: [ contig ] } )
	
	# --- run analysis per chromosome --- #
	haplo1 = {}
	haplo2 = {}
	for chromosome in contigs_by_chromosomes.keys():
		contigs = sorted( contigs_by_chromosomes[ chromosome ], key=itemgetter('contig_len') )[::-1]
		x1 = []	#haplophase 1
		x2 = []	#haplophase 2
		for contig in contigs:
			x1_overlap = calculate_overlap( x1, contig )
			x2_overlap = calculate_overlap( x2, contig )
			if x2_overlap < x1_overlap:
				x2.append( contig )
			else:
				x1.append( contig )
		x1_len = get_contig_cumlen( x1 )
		x2_len = get_contig_cumlen( x2 )
		haplo1.update( { chromosome: x1 } )
		haplo2.update( { chromosome: x2 } )
		print chromosome + " len: " + str( chr_lens[ chromosome ] ) + "\tx1: " + str( x1_len ) + "\tx2: " + str( x2_len )
	return haplo1, haplo2


def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] ).upper()


def construct_agp_and_fasta( haplo, agp, fasta, prefix, contigs, gap_len ):
	"""! @brief construct AGP and FASTA file """
	
	with open( fasta, "w" ) as seq_out:
		with open( agp, "w" ) as agp_out:
			#{'start': 9639543, 'chr': 'chr6', 'end': 13598982, 'ID': 'tig2391', 'contig_len': 3813575}
			for key in haplo.keys():
				sorted_contigs = sorted( haplo[ key ], key=itemgetter( 'start' ) )
				seq = []
				start_pos = 1
				for idx, contig in enumerate( sorted_contigs ):
					if idx > 0:	#add spacer/gap line
						agp_out.write( "\t".join( map( str, [ prefix+key,
																				start_pos,
																				start_pos+contig['contig_len'],
																				(2*idx),
																				"N",
																				gap_len,
																				"scaffold",
																				"yes",
																				"reference-based"
																				] ) ) + '\n' )
						start_pos += gap_len+1
					if not contig['orientation']:
						seq.append( revcomp( contigs[ contig['ID'] ] ) )
						agp_out.write( "\t".join( map( str, [ prefix+key,
																				start_pos,
																				start_pos+contig['contig_len'],
																				(2*idx)+1,
																				"W",
																				contig['ID'],
																				"1",
																				contig['contig_len'],
																				"+"
																				] ) ) + '\n' )
					else:
						seq.append( contigs[ contig['ID'] ] )
						agp_out.write( "\t".join( map( str, [ prefix+key,
																				start_pos,
																				start_pos+contig['contig_len'],
																				(2*idx)+1,
																				"W",
																				contig['ID'],
																				"1",
																				contig['contig_len'],
																				"-"
																				] ) ) + '\n' )
					start_pos += contig['contig_len']+1
				seq_out.write( '>' + key + '\n' + ("N"*gap_len).join( seq ) + '\n' )	#N gap bug fixed in v0.31


def main( arguments ):
	
	contig_file = arguments[ arguments.index('--in')+1 ]
	ref_file = arguments[ arguments.index('--ref')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	if '--gap' in arguments:
		gap_len = int( arguments[ arguments.index('--gap')+1 ] )
	else:
		gap_len = 100
	
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- generate BLAST results --- #
	result_file = output_folder + "blast_results.txt"
	blast_db = output_folder + "blastdb"
	
	if not os.path.isfile( result_file ):
		os.popen( "makeblastdb -in " +ref_file + " -out " + blast_db + " -dbtype nucl" )
		os.popen( "blastn -query " + contig_file + " -db " + blast_db + " -out " + result_file + " -outfmt 6 -evalue 0.0000000001 -word_size 50 -num_threads 20" )
	
	# --- get position per contig --- #
	#best_hit_lines = load_best_blast_hit( result_file )	#get best BLAST hit per contig
	best_hit_lines = load_best_blast_hit2( result_file )	#get best BLAST hit per contig
	contig_lens = load_seq_lengths( contig_file )
	chr_lens = load_seq_lengths( ref_file )
	pos_per_contig = get_pos_per_contig( best_hit_lines, contig_lens, chr_lens )	#extend best BLAST hit based on contig length
	print len( pos_per_contig )
	
	# --- sort contigs into haplophases --- #
	haplo1, haplo2 = group_contigs_to_haplophases( pos_per_contig, chr_lens )
	
	# --- generate AGP and FASTA files --- #
	contigs = load_sequences( contig_file )
	agp1 = output_folder + "haplo1.agp"
	fasta1 = output_folder + "haplo1.fasta"
	prefix1 = "hap1."
	construct_agp_and_fasta( haplo1, agp1, fasta1, prefix1, contigs, gap_len )
	
	agp2 = output_folder + "haplo2.agp"
	fasta2 = output_folder + "haplo2.fasta"
	prefix2 = "hap2."
	construct_agp_and_fasta( haplo2, agp2, fasta2, prefix2, contigs, gap_len )


if '--in' in sys.argv and '--ref' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
