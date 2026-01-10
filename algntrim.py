### Boas Pucker ###
### b.pucker@tu-bs.de ###

### some functions derived from KIPEs (Pucker et al., 2020: https://doi.org/10.3390/plants9091103) ###

__description__ = "Removes all columns from an alignment that exceed a certain percentage of gaps. This shortens the alignment and reduces the computational costs for inferring a phylogenetic tree"

__version__ = "v0.15"

__reference__ = """Pucker, B. & Iorizzo M. (2022). Apiaceae FNS I originated from F3H through tandem gene duplication. PLoS ONE 18(1): e0280155. https://doi.org/10.1371/journal.pone.0280155"""

__usage__ = """
					python3 algntrim.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					
					optional:
					--occ <MINIMAL_OCCUPANCY>
					--sort <ACTIVATES_ALPHANUMERICAL_SORTING>
					
					bug reports and feature requests: b.pucker@tu-bs.de
					"""

import os, sys

# --- end of imports --- #


def load_alignment( aln_file ):
	"""! @brief load alignment from input file """
	
	sequences, names = {}, []
	with open( aln_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				sequences.update( { header: "".join( seq ) } )
				names.append( header )
				header = line.strip()[1:]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )
		names.append( header )
	return sequences, names


def aln_len_check( alignment ):
	"""! @brief check if all sequences in alignment have the same length """
	
	lengths = []
	for key in list( alignment.keys() ):
		lengths.append( len( alignment[ key ] ) )
	if len( list( set( lengths ) ) ) != 1:	#not all sequences have same length
		for key in list( alignment.keys() ):
			sys.stdout.write( key + " - len: " + str( len( alignment[ key ] ) ) + "\n" )
		sys.stdout.flush()
		sys.exit( "ERROR: sequences do not have the same length." )		


def alignment_trimming( aln_file, cln_aln_file, occupancy, sorting ):
	"""! @brief remove all alignment columns with insufficient occupancy """
	
	alignment, names = load_alignment( aln_file )
	
	if sorting:	#alphanumerical sorting of names if activated
		names = sorted( names )
	
	# --- if there is an alignment (expected case) --- #
	if len( names ) > 0:
		aln_len_check( alignment )	#perform alignment check (same sequence lengths?)
		
		# --- identify valid residues in aligned sequences (columns with sufficient occupancy) --- #
		valid_index = []
		for idx, aa in enumerate( list(alignment.values())[0] ):
			counter = 0
			for key in names:
				if alignment[ key ][ idx ] != "-":
					counter += 1
			if counter / float( len( list(alignment.keys()) ) ) >= occupancy:
				valid_index.append( idx )
		
		# --- generate new sequences --- #
		with open( cln_aln_file, "w" ) as out:
			for key in names:
				seq = alignment[ key ]
				new_seq = []
				for idx in valid_index:
					new_seq.append( seq[ idx ] )
				new_seq  = "".join( new_seq )
				if new_seq.count('-') == len( new_seq ):	#exclude sequences that contain only gaps after trimming
					sys.stdout.write( "WARNING: only gaps remaining in sequence - " + key + " (sequence not included in output)\n" )
					sys.stdout.flush()
				else:
					out.write( ">" + key + '\n' + new_seq + '\n' )
	
	# --- just in case the alignment file is empty ---#
	else:
		sys.stdout.write( "WARNING: input file was empty (" + aln_file + ")\n" )
		sys.stdout.flush()
		with open( cln_aln_file, "w" ) as out:
			out.write( "" )


def main( arguments ):
	"""! @brief run everything """
	
	aln_file = arguments[ arguments.index('--in')+1 ]
	cln_aln_file = arguments[ arguments.index('--out')+1 ]
	if '--occ' in arguments:
		occupancy = float( arguments[ arguments.index('--occ')+1 ] )
	else:
		occupancy = 0.1
	
	if '--sort' in arguments:
		sorting = True
	else:
		sorting = False
	
	alignment_trimming( aln_file, cln_aln_file, occupancy, sorting )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
elif '--version' in sys.argv:
	sys.exit( __version__ )
elif '--reference' in sys.argv or '--cite' in sys.argv:
	sys.exit( __reference__ )
else:
	sys.exit( __description__ + "\n" + __version__ + "\n" + __usage__ + "\n" + __reference__ )
