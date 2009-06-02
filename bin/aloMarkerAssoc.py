#!/usr/local/bin/python

# Name: aloMarkerAssoc.py
# Purpose: handle associations between ALOs (Allele-Like Objects) and markers,
#	and related data.  See USAGE statement below for more information.
# Author: jsb
# Assumes: exclusive write access to the database
# Alters Data in these Database Tables:
#	SEQ_GeneTrap, SEQ_Allele_Assoc, ALL_Marker_Assoc, ALL_Allele,
#	MGI_Note, MGI_NoteChunk
# Environment variables:
#	required: MGD_DBSERVER, MGD_DBNAME, MGD_DBUSER, and (MGD_DBPASSWORD or MGD_DBPASSWORDFILE)
#	optional: LOGDIR, LOG_DEBUG
# Logs: writes log files to the current working directory, or to an alternate
#	logging directory specified by LOGDIR

import os
import sys
import time
import getopt
import tempfile
import traceback
import db
import runCommand
import aloMarkerLogger

###-----------------------###
###--- usage statement ---###
###-----------------------###

USAGE='''Usage: %s [-f]
    Purpose:
    	To update associations between ALOs (Allele-Like Objects) and markers,
	and related data.  This process includes:
	    1. setting the point coordinate for gene trap sequences
	    2. flag alleles which are from mixed cell lines
	    3. flag the representative sequence for each ALO
	    4. establish ALO/marker associations for alleles without curated
	    	marker associations
	    5. update ALO symbols to reflect marker nomenclature
	    6. create/update molecular notes for ALOs
    Flags:
	-f : write a log file of SQL commands which were used to update the db
    Notes:
	1. Requires several environment variables to be set.  See header
	   comments at the top of the script if you need the list.
	3. Using -f will write to the 'sql.log' file the statements that were
	   executed to make the changes to the database.  (Any read-only
	   statements are not written.)  This is just an extra logging
	   mechanism; the script still executes the commands inline.
	4. Does not do a drop-and-reload of any data; uses inline SQL to
	   apply changes to the database after computing a diff with the
	   current data state.
	5. Given #4, all changes are computed using other data in the
	   database.  So, if a failure occurs, the script can be re-run as
	   needed once any source data problems have been resolved.
	6. Assumes exclusive write-access to the database.
	7. Writes logs to the directory given in the LOGDIR environment
	   variable.  If that is not defined, then the current working
	   directory is used instead.
	8. If we cannot write curator and diagnostic logs to the log dir (see
	   item #7), we will instead write them to /tmp with an error message
	   to stdout and to stderr.
''' % sys.argv[0]

###------------------------###
###--- global variables ---###
###------------------------###

START_TIME = time.time()		# time (in seconds) at start of run
VERBOSE = False				# run in verbose mode?
SQL_FILE = None				# file pointer to output SQL cmds
TEMP_FILES = []				# temporary files to delete when done
UPDATED_TABLES = []			# list of tables with updated records
LOGDIR = os.getcwd()			# logging directory (start in cwd)

# for collecting log items
LOGGER = aloMarkerLogger.AloMarkerLogManager()

# We will need to adjust the molecular notes for alleles, based on changes
# made by this script.  There are five standard molecular notes:
#	A. alleles associated with one marker by load
#	B. unassociated alleles with a single good hit to the genome
#	C. alleles with zero or multiple good hits to the genome
#	D. alleles from mixed lines
#	E. alleles where the marker association was invalidated by a curator
#		(E is never set by the load)

NOTE_A = None	# initialized later by setupMolecularNotes() function
NOTE_B = None
NOTE_C = None
NOTE_D = None

###-----------------###
###--- functions ---###
###-----------------###

###------------------------------------------------------------------------###

def bailout (
	msg,			# string; error message to write to user
	showUsage = True	# boolean; give usage statement or not?
	):
	# Purpose: to give an error message to the user and exit with an error
	#	code
	# Returns: nothing
	# Assumes: nothing
	# Effects: writes to stderr, exits to OS
	# Throws: SystemExit, to force exit to OS

	if showUsage:
		sys.stderr.write (USAGE)
	sys.stderr.write ('Error: %s\n' % msg)
	sys.exit(1)

###------------------------------------------------------------------------###

def optimizeTables ():
	# Purpose: run 'update statistics' on any database tables that we
	#	altered, to ensure optimal query performance
	# Returns: nothing
	# Assumes: nothing
	# Effects: gives SQL commands to the database
	# Throws: propagates any exceptions from db.sql()

	debug ('in optimizeTables()')
	for table in UPDATED_TABLES:
		update ('UPDATE STATISTICS %s' % table)
		LOGGER.log ('diag', 'Updated table statistics on %s' % table)
	return

###------------------------------------------------------------------------###

def flagTable (table):
	# Purpose: note that we have altered data in the given 'table'
	# Returns: nothing
	# Assumes: nothing
	# Effects: alters global 'UPDATED_TABLES'
	# Throws: nothing

	global UPDATED_TABLES
	if table not in UPDATED_TABLES:
		UPDATED_TABLES.append (table)
	return

###------------------------------------------------------------------------###

def debug (
	s			# string; message to be given to user
	):
	# Purpose: if running in verbose mode, writes the given message
	#	to stdout with a timestamp
	# Returns: nothing
	# Assumes: nothing
	# Effects: writes to stdout
	# Throws: nothing

	if VERBOSE:
		LOGGER.log ('debug', s)
		print 'debug: ' + s
	return

###------------------------------------------------------------------------###

def getTempFile (
	prefix, 		# string; prefix for desired filename
	suffix			# string; suffix for desired filename
	):
	# Purpose: generate a new temp file named according to the given
	#	'prefix' and 'suffix'
	# Returns: string -- path to uniquely-named temporary file
	# Assumes: nothing
	# Effects: creates and closes file in file system, so its name will
	#	be reserved and not scooped up by another process
	# Throws: propagates any exceptions from tempfile.mkstemp()

	global TEMP_FILES

	(fp, path) = tempfile.mkstemp (suffix, prefix, text = True)
	TEMP_FILES.append (path)
	return path

###------------------------------------------------------------------------###

def update (
	cmd		# string; SQL command to update data in the database
	):
	# Purpose: excecute the given 'cmd' command against the database to
	#	update some data
	# Returns: results returned from db.sql(cmd) call, or None if no
	#	results (which is most common for updates)
	# Assumes: db module has been properly initialized
	# Effects: executes 'cmd' as a SQL command against the database
	# Throws: propagates any exceptions from db.sql()
	# Notes: If excecuting in SQL-file mode (-f), will print 'cmd' to
	#	SQL_FILE, in addition to executing it against the database.
	#	This provides us a record of changes made by the load.

	# if we need to log the update commands, do so before executing it
	if SQL_FILE:
		SQL_FILE.write ('%s\n' % cmd)
	results = db.sql (cmd, 'auto')
	return results

###------------------------------------------------------------------------###

def processCommandLine():
	# Purpose: process the command-line parameters and update various
	#	global variables.  also handles the db module setup.
	# Returns: nothing
	# Assumes: nothing
	# Effects: updates global variables, initializes db module, does
	#	test query against the database
	# Throws: propagates SystemExit from bailout() in case of errors

	global VERBOSE, SQL_FILE, LOGDIR

	# extract options from command-line

	try:
		(options, args) = getopt.getopt (sys.argv[1:], 'f')
	except getopt.GetoptError:
		bailout ('Invalid command-line')

	# pull out options and update global variables

	sql_filename = None
	for (option, value) in options:
		if option == '-f':
			sql_filename = 'sql.log'
		else:
			bailout ('Unsupported option: %s' % option)

	# pull necessary values out of environment variables

	environ = os.environ
	for var in [ 'MGD_DBSERVER', 'MGD_DBNAME', 'MGD_DBUSER' ]:
		if not environ.has_key(var):
			bailout ('Missing environment variable: %s' % var)

	server = environ['MGD_DBSERVER']
	database = environ['MGD_DBNAME']
	user = environ['MGD_DBUSER']

	if environ.has_key('LOGDIR'):
		LOGDIR = environ['LOGDIR']

	if environ.has_key('LOG_DEBUG'):
		if environ['LOG_DEBUG'].lower() != 'false':
			VERBOSE = True

	if environ.has_key('MGD_DBPASSWORD'):
		password = environ['MGD_DBPASSWORD']
	elif environ.has_key('MGD_DBPASSWORDFILE'):
		try:
			fp = open (environ['MGD_DBPASSWORDFILE'], 'r')
			password = fp.readline().strip()
			fp.close()
		except:
			bailout ('Cannot read password from: %s' % value)
	else:
		bailout('Must specify either MGD_DBPASSWORD or MGD_DBPASSWORDFILE' + \
			' environment variable')

	# set up SQL logging, if necessary

	if not os.path.exists (LOGDIR):
		bailout ('Log directory does not exist: %s' % LOGDIR)
	if not os.path.isdir (LOGDIR):
		bailout ('Log directory is not a directory: %s' % LOGDIR)

	if sql_filename:
		sql_filename = os.path.join (LOGDIR, sql_filename)
		try:
			SQL_FILE = open (sql_filename, 'w')
		except:
			bailout ('Cannot write SQL log to: %s' % sql_filename)

	# verify database-access parameters

	if (server and database and user and password):
		db.set_sqlLogin (user, password, server, database)
		db.useOneConnection(1)
		try:
			db.sql ('SELECT COUNT(1) FROM MGI_dbInfo', 'auto')
		except:
			bailout ('Cannot query database')
	else:
		bailout ('Missing database access parameter(s)')

	# debugging output of parameters

	debug ('Log directory: %s' % LOGDIR)
	debug ('Database: %s..%s  User: %s' % (server, database, user))
	if sql_filename:
		debug ('Write to SQL File: %s' % sql_filename)
	if VERBOSE:
		debug ('Verbose output: True')
	
	LOGGER.log ('diag', 'Processed command-line')
	return

###------------------------------------------------------------------------###

def conditional (comparison, trueValue, falseValue):
	# Purpose: choose either 'trueValue' or 'falseValue' based on whether
	#	'comparison' evaluates to true or false
	# Returns: 'trueValue' if 'comparison' evaluates to true, or
	#	'falseValue' otherwise
	# Assumes: nothing
	# Effects: nothing
	# Throws: nothing

	if comparison:
		return trueValue
	return falseValue

###------------------------------------------------------------------------###

def setPointCoordinates ():
	# Purpose: assign the point coordinate for each gene trap sequence,
	#	where possible (some will not have point coordinates)
	# Returns: nothing
	# Assumes: we can write to the database
	# Effects: updates SEQ_GeneTrap.pointCoordinate field in the database
	# Throws: propagates any exceptions raised by db.sql() function

	debug ('in setPointCoordinates()')

	race5 = lookupTerm ('Sequence Tag Method', "5' Race")
	race3 = lookupTerm ('Sequence Tag Method', "3' Race")

	upstream = lookupTerm ('Gene Trap Vector End', 'upstream')
	downstream = lookupTerm ('Gene Trap Vector End', 'downstream')

	revComp = lookupTerm ('Reverse Complement', 'yes')

	cmd = '''SELECT gt._Sequence_key,
			gt.pointCoordinate,
			cc.startCoordinate,
			cc.endCoordinate,
			cc.strand,
			gt._TagMethod_key,
			gt._VectorEnd_key,
			gt._ReverseComp_key
		FROM SEQ_GeneTrap gt,
			SEQ_Coord_Cache cc
		WHERE gt._Sequence_key *= cc._Sequence_key'''

	results = db.sql (cmd, 'auto')
	LOGGER.log ('diag', 'Retrieved %d sequence coordinates' % len(results))

	toSet = []
	toNull = []
	for row in results:
		pointCoord = None
		lower = row['startCoordinate']

		if lower != None:
			upper = row['endCoordinate']
			strand = row['strand']
			tagMethod = row['_TagMethod_key']

			# if DNA-based sequence tag method...
			if 3983002 <= tagMethod <= 3983007:
				if (row['_VectorEnd_key'] == upstream) and \
				    (row['_ReverseComp_key'] == revComp):
					pointCoord = conditional (
						strand == '+', upper, lower)
				else:
					pointCoord = conditional (
						strand == '+', lower, upper)

			# two RNA-based sequence tag methods...
			elif tagMethod == race5:
				pointCoord = conditional (
					strand == '+', upper, lower)
			elif tagMethod == race3:
				pointCoord = conditional (
					strand == '+', lower, upper)

		# if this is a change for the existing point coordinate, then
		# we need to update the field

		if pointCoord != row['pointCoordinate']:
			if pointCoord == None:
				toNull.append (row['_Sequence_key'])
			else:
				toSet.append (
					(row['_Sequence_key'], pointCoord) )

	LOGGER.log ('diag', 'Identified new point coords')

	cmd1 = '''UPDATE SEQ_GeneTrap SET pointCoordinate = %d 
		WHERE _Sequence_key = %d'''
	for (seqKey, pointCoord) in toSet:
		update (cmd1 % (pointCoord, seqKey))
	LOGGER.log ('diag', 'Updated point coordinate for %d sequences' % \
		len(toSet))

	cmd2 = '''UPDATE SEQ_GeneTrap SET pointCoordinate = null 
		WHERE _Sequence_key = %d'''
	for seqKey in toNull:
		update (cmd2 % seqKey)
	LOGGER.log ('diag',
		'Changed to null point coordinate for %d sequences' % \
		len(toNull))

	LOGGER.log ('diag',
		'Left point coordinates as-is for %d sequences' % \
		(len(results) - len(toSet) - len(toNull)) )

	if toSet or toNull:
		flagTable ('SEQ_GeneTrap')
	return

###------------------------------------------------------------------------###

def initLogger ():
	# Purpose: initialize the logger module by retrieving some key/value
	#	mappings from the database and passing them on to the logger
	#	module
	# Returns: nothing
	# Assumes: we can query the database
	# Effects: queries the database, initializes aloMarkerLogger.py
	# Throws: propagates all exceptions from db.sql()

	# allele keys and symbols

	debug ('in initLogger()')

	cmd1 = 'SELECT _Allele_key, symbol FROM ALL_Allele'
	results1 = db.sql (cmd1, 'auto')

	alleles = {}
	for row in results1:
		alleles[row['_Allele_key']] = row['symbol']
	aloMarkerLogger.setAlleleSymbols(alleles)
	LOGGER.log ('diag', 'Retrieved %d allele symbols (for reporting)' % \
		len(alleles))

	# marker keys and symbols (current and interim)

	cmd2 = '''SELECT _Marker_key, symbol
			FROM MRK_Marker
			WHERE _Organism_key = 1
				AND _Marker_Status_key IN (1,3)'''
	results2 = db.sql (cmd2, 'auto')

	markers = {}
	for row in results2:
		markers[row['_Marker_key']] = row['symbol']
	aloMarkerLogger.setMarkerSymbols(markers)
	LOGGER.log ('diag', 'Retrieved %d marker symbols (for reporting)' % \
		len(markers))

	# reference keys and J: numbers for mixed references

	cmd3a = '''SELECT _RefAssocType_key
		FROM MGI_RefAssocType
		WHERE _MGIType_key = 11		-- Allele
			AND assocType = "Mixed" '''
	results3a = db.sql (cmd3a, 'auto')

	if len(results3a) > 0:
		refAssocType = results3a[0]['_RefAssocType_key']

		cmd3b = '''SELECT r._Refs_key, b.jnumID
			FROM MGI_Reference_Assoc r,
				BIB_Citation_Cache b
			WHERE r._Refs_key = b._Refs_key
				AND r._RefAssocType_key = %d''' % refAssocType
		results3 = db.sql (cmd3b, 'auto')

		jnums = {}
		for row in results3:
			jnums[row['_Refs_key']] = row['jnumID']
		aloMarkerLogger.setReferences (jnums)
		LOGGER.log ('diag', 'Retrieved %d J: numbers (for reporting)'\
			% len(jnums))

	# sequence keys and IDs for sequences associated with alleles

	cmd4 = '''SELECT s._Sequence_key, a.accID
		FROM SEQ_Allele_Assoc s,
			ACC_Accession a
		WHERE s._Sequence_key = a._Object_key
			AND a._MGIType_key = 19
			AND a.private = 0
			AND a.preferred = 1'''
	results4 = db.sql (cmd4, 'auto')

	seqs = {}
	for row in results4:
		seqs[row['_Sequence_key']] = row['accID']
	aloMarkerLogger.setSeqIDs (seqs)
	LOGGER.log ('diag', 'Retrieved %d sequence IDs (for reporting)' % \
		len(seqs))
	return

###------------------------------------------------------------------------###

def writeGffFile (fp,		# file pointer, opened for writing
	rows, 			# list of dictionaries, one for each object
	source, 		# string; value for sequence source
	objectType, 		# string; type of object contained in 'rows'
	idField			# string; name of unique ID field in
				# ...dictionary in 'rows'
	):
	# Purpose: write contents of 'rows' into a roughly GFF-formatted file
	# Returns: nothing
	# Assumes: each dictionary in rows contains fields: 'startCoordinate',
	#	'chromosome', 'endCoordinate', 'strand', and idField
	# Effects: writes to 'fp', then closes 'fp' when done
	# Throws: any exceptions raised when writing to (or closing) 'fp'

	for row in rows:
		if row['startCoordinate'] != None:
			fp.write ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
				row['chromosome'],
				source,
				objectType,
				int(row['startCoordinate']),
				int(row['endCoordinate']),
				'-',				# unused
				row['strand'],
				'-',				# unused
				row[idField]) )
	fp.close()
	return

###------------------------------------------------------------------------###

def getSequences():
	# Purpose: retrieve relevant sequence data from the database
	# Returns: string -- path to GFF-formatted file of sequences
	# Assumes: db module has been initialized, can write to file system
	# Effects: queries database, writes to file system
	# Throws: propagates all exceptions from querying database or from
	#	writing to file system

	debug ('in getSequences()')
	upstream = lookupTerm ('Gene Trap Vector End', 'upstream')
	notRevComp = lookupTerm ('Reverse Complement', 'no')

	# we are seeking gene trap sequences which have coordinates

	cmd = '''SELECT gt._Sequence_key,
			cc.chromosome,
			gt.pointCoordinate,
			cc.startCoordinate,
			cc.endCoordinate,
			cc.strand,
			gt._TagMethod_key,
			gt._VectorEnd_key,
			gt._ReverseComp_key
		FROM SEQ_GeneTrap gt,
			SEQ_Coord_Cache cc
		WHERE gt._Sequence_key = cc._Sequence_key
			AND cc.startCoordinate != null'''

	results = db.sql (cmd, 'auto')
	LOGGER.log ('diag', 'Got %d gene trap sequences with coordinates' % \
		len(results))

	flipCt = 0	# number of sequences where we flipped the strand
	pcCt = 0	# number of sequences with point coordinates

	for row in results:
		# if we have a DNA-based sequence tag method, and if the
		# vector end is 'upstream', then we need to reverse the strand
		# for the sake of comparison with markers

		if 3983002 <= row['_TagMethod_key'] <= 3983007:
			if (row['_VectorEnd_key'] == upstream) and \
			    (row['_ReverseComp_key'] == notRevComp):
				if row['strand'] == '+':
					row['strand'] = '-'
				else:
					row['strand'] = '+'
				flipCt = flipCt + 1

		# if we have a point coordinate, we should use that for the
		# comparisons

		if row['pointCoordinate'] != None:
			row['startCoordinate'] = row['pointCoordinate']
			row['endCoordinate'] = row['pointCoordinate']
			pcCt = pcCt + 1

	LOGGER.log ('diag', 'Needed to swap strand on ' + \
		'%d sequences for comparison with markers' % flipCt)
	LOGGER.log ('diag', 'Using point coordinate for comparison ' + \
		'on %d sequences' % pcCt)
	LOGGER.log ('diag', 'Using start/end coordinates for comparison ' + \
		'on %d sequences' % (len(results) - pcCt) )

	# write sequences to temporary file

	path = getTempFile ('sequences.', '.gff')
	seqFile = open (path, 'w')
	writeGffFile (seqFile, results, '-', 'Sequence', '_Sequence_key')
	LOGGER.log ('diag', 'Wrote file of sequence data to: %s' % path)
	return path

###------------------------------------------------------------------------###

def getMarkers():
	# Purpose: retrieve relevant marker data from the database
	# Returns: string -- path to GFF-formatted file of markers
	# Assumes: db module has been initialized, can write to file system
	# Effects: queries database, writes to file system
	# Throws: propagates all exceptions from querying database or from
	#	writing to file system

	debug ('in getMarkers()')

	# seeking all markers with coordinate information

	cmd = '''SELECT _Marker_key,
			chromosome,
			startCoordinate,
			endCoordinate,
			strand
		FROM MRK_Location_Cache
		WHERE startCoordinate != null'''
	results = db.sql (cmd, 'auto')
	LOGGER.log ('diag', 'Retrieved %d markers with coordinates' % \
		len(results))

	# write markers to temporary GFF-formatted file

	path = getTempFile ('markers.', '.gff')
	mrkFile = open (path, 'w')
	writeGffFile (mrkFile, results, '-', 'Marker', '_Marker_key')
	LOGGER.log ('diag', 'Wrote file of marker data to: %s' % path)
	return path

###------------------------------------------------------------------------###

def fjoin (
	f1,		# string; path to GFF-formatted file of sequences
	f2		# string; path to GFF-formatted file of markers
	):
	# Purpose: execute fjoin script against files 'f1' and 'f2', to
	#	determine where sequences and markers overlap
	# Returns: string -- path to temp file with output from fjoin
	# Assumes: can read 'f1' and 'f2', and can write a new temp file
	# Effects: reads and writes from file system
	# Throws: propagates exceptions from reading/writing files and from
	#	running fjoin via runCommand

	debug ('in fjoin()')
	path = getTempFile ('fjoin.', '.out')

	# if a 'sort' implementation exists, it is a bit faster than just
	# letting fjoin sort the files in python

	if os.path.exists ('/usr/local/bin/sort'):
		gnuSort = ' --gnuSort /usr/local/bin/sort'
	else:
		gnuSort = ''

	# command necessary to run fjoin for overlap of 1bp, which does
	# sorting of input files, and which considers strand when looking for
	# overlaps

	cmd = '/usr/local/mgi/live/lib/python/fjoin.py ' + \
			'-1 %s ' % f1 + \
			'-2 %s ' % f2 + \
			'-s both ' + \
			'-k 1 ' + \
			'-o %s%s' % (path, gnuSort)
	debug ('Running: %s' % cmd)

	(stdout, stderr, exitcode) = runCommand.runCommand (cmd)
	if (exitcode):
		bailout ('fJoin failed with exit code: %d -- stderr: %s' % (
			exitcode, stderr), False)

	# if running in verbose mode, echo fjoin's reporting output

	debug ('fJoin output: <PRE>%s</PRE>' % stderr)
	LOGGER.log ('diag', 'Ran fjoin on markers and sequences')
	return path

###------------------------------------------------------------------------###

def readFjoinFile (
	f			# string; path to file of output from fjoin
	):
	# Purpose: read and parse fjoin output file 'f'
	# Returns: list of two-item tuples, each of which is an integer
	#	database key
	# Assumes: can read from 'f'; unique IDs for items in 'f' are integer
	#	database keys
	# Effects: reads from file system
	# Throws: propagates any exceptions from reading the file

	debug ('in readFjoinFile()')
	try:
		fp = open (f, 'r')
		lines = fp.readlines()
		fp.close()
	except:
		bailout ('Cannot read fJoined file: %s' % f)

	pairs = []
	for line in lines:
		items = line[:-1].split('\t')
		if len(items) >= 18:
			# remember that column 0 contains the amount of
			# overlap, then the columns from the two files are
			# to the right.  last column for each file is the
			# unique ID (which is all we care about).

			pairs.append ( (int(items[9]), int(items[18])) )
		else:
			bailout ('fjoin row with too few columns: %s' % line)

	LOGGER.log ('diag',
		'Fjoin found %d pairs of overlapping sequences and markers' \
		% len(pairs) )
	return pairs

###------------------------------------------------------------------------###

def groupby (
	pairs,		# list of 2-item tuples from readFjoinFile()
	i = 0,		# which item in tuple to group by
	j = 1		# which item in tuple to collect for each item
	):
	# Purpose: collate the items in 'pairs' such that for each tuple (i,j)
	#	we collect all the 'j' items which share a common 'i'
	# Returns: dictionary; keys are 'i' values, each referring to a list
	#	of 'j' values
	# Assumes: nothing
	# Effects: nothing
	# Throws: nothing

	debug ('in groupby()')
	dict = {}
	for pair in pairs:
		key = pair[i]
		if dict.has_key(key):
			dict[key].append (pair[j])
		else:
			dict[key] = [ pair[j] ]
	return dict

###------------------------------------------------------------------------###

def splitSingles (
	dict		# dictionary, as produced by groupby()
	):
	# Purpose: extract those items from 'dict' that only have one item
	#	in their list of items for a given key
	# Returns: two dictionaries; one with keys that refer to 1-item lists,
	#	and one with keys that refer to multi-item lists.
	# Assumes: nothing
	# Effects: nothing
	# Throws: nothing

	debug ('in splitSingles()')
	singles = {}
	multiples = {}

	for (key, values) in dict.items():
		if len(values) == 1:
			singles[key] = values
		else:
			multiples[key] = values

	LOGGER.log ('diag',
		'Identified %d sequences which overlap a single marker' % \
		len(singles))
	LOGGER.log ('diag',
		'Identified %d sequences which overlap multiple markers' % \
		len(multiples))
	return singles, multiples

###------------------------------------------------------------------------###

def findSequencesNoMarkers (
	dict,			# dictionary, as produced by groupby()
	lengths			# dictionary, mapping seq key to seq length
	):
	# Purpose: find gene trap sequences which overlap no markers
	# Returns: dictionary, each key of which is a sequence key which does
	#	not overlap a marker
	# Assumes: keys of 'dict' are sequence keys; findSequenceLengths() has
	#	already been run
	# Effects: nothing
	# Throws: nothing

	debug ('in findSequenceNoMarkers()')
	noMarkers = {}				# dictionary to return

	# if a given gene trap sequence has a length, but does not appear in
	# the dictionary of seq-to-marker overlaps, then it has no overlap

	for k in lengths.keys():
		if not dict.has_key(k):
			noMarkers[k] = 1

	LOGGER.log ('diag',
		'Identified %d sequences which overlap no markers' % \
		len(noMarkers))
	return noMarkers

###------------------------------------------------------------------------###

def findSequenceLengths():
	# Purpose: retrieve the lengths of the various gene trap sequences
	#	from the database
	# Returns: dictionary, mapping seq keys to their seq lengths
	# Assumes: db module has been initialized
	# Effects: queries the database
	# Throws: propagates any exceptions from db.sql()

	debug ('in findSequenceLengths()')

	# sequence-to-allele associations are maintained by the ALO load and
	# by curators.  We'll pull lengths for all of them.

	cmd = '''SELECT aa._Sequence_key, ss.length
		FROM SEQ_Allele_Assoc aa,
			SEQ_Sequence ss
		WHERE aa._Sequence_key = ss._Sequence_key'''
	results = db.sql (cmd, 'auto')

	# update them in a dictionary that stores the length of each sequence

	lengths = {}
	for row in results:
		lengths[row['_Sequence_key']] = row['length']
	LOGGER.log ('diag', 'Retrieved %d sequence lengths' % len(lengths))
	return lengths

###------------------------------------------------------------------------###

def lookupTerm (
	vocab, 			# string; name of vocabulary
	term			# string; name of vocabulary term
	):
	# Porpoise: find the term key for the given 'term' in the given
	#	'vocab'  :)
	# Returns: integer -- term key
	# Assumes: db module has been initialized
	# Effects: queries the database
	# Throws: propagates any exceptions raised in db.sql()

	cmd = '''SELECT t._Term_key
		FROM VOC_Term t, VOC_Vocab v
		WHERE t._Vocab_key = v._Vocab_key
			AND t.term = "%s"
			AND v.name = "%s"''' % (term, vocab)
	results = db.sql (cmd, 'auto')

	if len(results) == 0:
		bailout ('Cannot find "%s" in vocab "%s"' % (term, vocab),
			False)
	return results[0]['_Term_key']

###------------------------------------------------------------------------###

def getRepSeq (associations,	# list of dictionaries, each an association
	singleMarkers, 		# dictionary; seq keys map to lists w/1 mrk
	multiMarkers, 		# dictionary; seq keys map to lists w/2+ mrk
	noMarkers,		# dictionary; seq keys which map to no mrk
	mixedAlleles,		# dictionary; allele keys which are mixed
	lengths			# dictionary; seq keys map to seq lengths
	):
	# Purpose: identify the representative sequence for this allele
	# Returns: dictionary; item from 'associations' for the representative
	#	sequence
	# Assumes: each dictionary in 'associations' has the following keys:
	#	'_Allele_key', '_Sequence_key'; also, 'associations' has at
	#	least one dictionary
	# Effects: nothing
	# Throws: nothing

	alleleKey = associations[0]['_Allele_key']

	# the easy cases first...
	# 1. if is a mixed allele, designate no representative sequence
	# 2. if allele has only one sequence, designate that one as rep

	if mixedAlleles.has_key(alleleKey):
		return None
	if len(associations) == 1:
		return associations[0]

	# At this point, we have an allele with multiple sequences and which
	# is not mixed.  We must examine overlaps with markers and lengths for
	# the sequences.  The next cases...
	# 2b. prefer DNA sequence tag methods over RNA ones
	# 3. longest sequence which overlaps only one marker
	# 4. longest sequence which overlaps more than one marker
	# 5. longest sequence which overlaps no markers
	# 6. longest sequence which has no coordinates

	seqs = []
	for assoc in associations:
		seqKey = assoc['_Sequence_key']

		if lengths.has_key(seqKey):
			length = lengths[seqKey]
		else:
			length = 0			# should not happen

		# if DNA-based sequence tag method, prefer over RNA-based
		if 3983002 <= row['_TagMethod_key'] <= 3983007:
			isDNA = 1
		else:
			isDNA = 0

		if singleMarkers.has_key (seqKey):
			tpl = (3, isDNA, length, seqKey, assoc)
		elif multiMarkers.has_key (seqKey):
			tpl = (2, isDNA, length, seqKey, assoc)
		elif noMarkers.has_key (seqKey):
			tpl = (1, isDNA, length, seqKey, assoc)
		else:						# no coords
			tpl = (0, isDNA, length, seqKey, assoc)

		seqs.append (tpl)

	seqs.sort()
	return seqs[-1][-1]		# last one is the best one

###------------------------------------------------------------------------###

def splitList (
	original, 		# list; the list to be divided into sublists
	length = 100		# integer; max length of each sublist
	):
	# Purpose: break 'original' into smaller lists, each with a maximum
	#	given 'length'
	# Returns: list of lists, containing the same items as 'original'
	#	in the same order, but broken into smaller chunks
	# Assumes: nothing
	# Effects: nothing
	# Throws: nothing

	sublists = []
	lenOriginal = len(original)
	i = 0
	j = 0
	while j < lenOriginal:
		i = j
		j = j + length
		sublists.append (original[i:j])
	return sublists

###------------------------------------------------------------------------###

def updateRepSeqs (
	singleMarkers, 		# dictionary; seq keys map to lists w/1 mrk
	multiMarkers, 		# dictionary; seq keys map to lists w/2+ mrk
	noMarkers,		# dictionary; seq keys which map to no mrk
	mixedAlleles,		# dictionary; allele keys which are mixed
	lengths			# dictionary; seq keys map to seq lengths
	):
	# Purpose: update which sequence is flagged as 'representative' for
	#	each allele
	# Returns: dictionary which maps from each sequence key to a list of
	#	its allele keys
	# Assumes: db module has been initialized
	# Effects: queries the database, updates data in the database
	# Throws: propagates any exceptions from db.sql()

	debug ('in updateRepSeqs()')

	# term key for 'representative' qualifier and N/A qualifier
	repQualifier = lookupTerm ("Sequence Allele Association Qualifier",
		"Representative")
	notAppQualifier = lookupTerm ("Sequence Allele Association Qualifier",
		"Not Applicable")

	# get the set of sequence/allele associations, and the current
	# representative flags

	cmd = '''SELECT a._Assoc_key,
			a._Allele_key,
			a._Sequence_key,
			a._Qualifier_key,
			g._TagMethod_key
		FROM SEQ_Allele_Assoc a,
			SEQ_GeneTrap g
		WHERE a._Sequence_key *= g._Sequence_key'''
	results = db.sql (cmd, 'auto')

	LOGGER.log ('diag', 'Retrieved %d sequence/allele associations' % \
		len(results))

	alleles = {}		# allele key -> list of rows from results
	seqToAllele = {}	# seq key -> list of allele keys

	for row in results:
		alleleKey = row['_Allele_key']
		if alleles.has_key(alleleKey):
			alleles[alleleKey].append (row)
		else:
			alleles[alleleKey] = [ row ]

		seqKey = row['_Sequence_key']
		if seqToAllele.has_key(seqKey):
			seqToAllele[seqKey].append (alleleKey)
		else:
			seqToAllele[seqKey] = [ alleleKey ]

	LOGGER.log ('diag', 'Collated seq/allele assoc')

	asIs = 0	# integer; count of representative seqs left as-is
	demoted = []	# list of integer assoc keys to remove rep flag
	promoted = []	# list of integer assoc keys to add rep flag

	# walk through each allele and its associated sequences

	for (alleleKey, associations) in alleles.items():
		oldRep = None	# previous representative sequence association
		newRep = None	# new representative sequence association

		# find the previously-flagged representative seq association

		for assoc in associations:
			if assoc['_Qualifier_key'] == repQualifier:
				oldRep = assoc
				break

		# find the new representative sequence association

		newRep = getRepSeq (associations, singleMarkers,
			multiMarkers, noMarkers, mixedAlleles, lengths)

		# if they're different, we need to promote one and demote the
		# other (assuming they're non-null)

		if oldRep != newRep:
			if oldRep != None:
				demoted.append (oldRep['_Assoc_key'])
			if newRep != None:
				promoted.append (newRep['_Assoc_key'])
		else:
			asIs = asIs + 1

	LOGGER.log ('diag', 'Determined rep seqs')

	# sort the lists, so that changes are grouped for associations stored
	# close together in the database.  This should result in fewer pages
	# swapped from disk and be slightly more efficient.

	# un-flag all the newly-demoted associations

	demoted.sort()
	demotedSublists = splitList(demoted)
	cmd = '''UPDATE SEQ_Allele_Assoc
		SET _Qualifier_key = %d
		WHERE _Assoc_key IN (%s)'''
	for sublist in demotedSublists:
		update (cmd % (notAppQualifier, ','.join(map (str, sublist))))
	LOGGER.log ('diag',
	  'Demoted %d previously representative sequence/allele associations'\
		% len(demoted))

	# flag all the newly-promoted associations

	promoted.sort()
	promotedSublists = splitList(promoted)
	cmd = '''UPDATE SEQ_Allele_Assoc
		SET _Qualifier_key = %d
		WHERE _Assoc_key IN (%s)''' % (repQualifier, '%s')
	for sublist in promotedSublists:
		update (cmd % ','.join(map (str, sublist)) )
	LOGGER.log ('diag',
		'Promoted %d new representative sequence/allele associations'\
		% len(promoted))

	LOGGER.log ('diag',
		'Left %d sequence/allele associations as-is' % asIs)

	if demoted or promoted:
		flagTable ('SEQ_Allele_Assoc')
	return seqToAllele

###------------------------------------------------------------------------###

GAPS = []		# list of unused _Assoc_key values in ALL_Marker_Assoc
NEXT_GAP = None		# index into GAPS of the next _Assoc_key value to use
MAX_ASSOC_KEY = None	# integer; max value of ALL_Marker_Assoc._Assoc_key

def nextAssocKey():
	# Purpose: get the next _Assoc_key to use when inserting a record
	#	into ALL_Marker_Assoc
	# Returns: integer
	# Assumes: db module has been initialized
	# Effects: queries the database
	# Throws: propagates any exceptions from db.sql()

	global GAPS, NEXT_GAP, MAX_ASSOC_KEY

	# Look up the set of association keys, so we can see which ones have
	# been deleted previously.  Let's fill in the gaps, to avoid just
	# having an ever-increasing association key value...  (not required;
	# just programmer-preference)

	if NEXT_GAP == None:
		cmd = '''SELECT _Assoc_key 
			FROM ALL_Marker_Assoc 
			ORDER BY _Assoc_key'''
		results = db.sql (cmd, 'auto')

		lastKey = 0
		for row in results:
			assocKey = row['_Assoc_key']
			if (assocKey - 1) > lastKey:
				GAPS = GAPS + range(lastKey + 1, assocKey)
			lastKey = assocKey

		LOGGER.log ('diag',
			'Found %d unused _Assoc_key values in ALL_Marker_Assoc that we can re-use' % len(GAPS))
		NEXT_GAP = 0

	if NEXT_GAP < len(GAPS):
		key = GAPS[NEXT_GAP]
		NEXT_GAP = NEXT_GAP + 1
		return key

	# no more gaps to use, so we just need to keep incrementing the max
	# assoc key from the database and add records to the end

	if MAX_ASSOC_KEY == None:
		cmd = 'SELECT MAX(_Assoc_key) FROM ALL_Marker_Assoc'
		results = db.sql (cmd, 'auto')

		if results:
			MAX_ASSOC_KEY = results[0]['']
		else:
			MAX_ASSOC_KEY = 0

	MAX_ASSOC_KEY = MAX_ASSOC_KEY + 1
	return MAX_ASSOC_KEY

###------------------------------------------------------------------------###

def updateMarkerAssoc (
	singleMarkers, 		# dictionary; seq keys map to lists w/1 mrk
	seqToAllele,		# dictionary; seq keys to list of allele keys
	mixedAlleles		# dictionary; allele keys which are mixed
	):
	# Purpose: update the allele/marker associations
	# Returns: list of integer allele keys which had changed marker
	#	associations
	# Assumes: db module has been initialized
	# Effects: updates data in the database
	# Throws: propagates any exceptions from db.sql()

	debug ('in updateMarkerAssoc()')

	# get a couple term keys that we will need for the associations

	notApp = lookupTerm ('Marker-Allele Association Qualifier',
		'Not Specified')
#		'Not Applicable')
	loadedTerm = lookupTerm ('Marker-Allele Association Status', 'Loaded')

	# get set of alleles associated with nomen markers (the load should 
	# not alter marker associations for these alleles)

	nomenAlleles = {}
	symbols = {}
	cmd = '''SELECT _Allele_key, symbol, nomenSymbol
		FROM ALL_Allele 
		WHERE nomenSymbol != null'''
	results = db.sql(cmd, 'auto')
	for row in results:
		nomenAlleles[row['_Allele_key']] = row['nomenSymbol']
		symbols[row['_Allele_key']] = row['symbol']

	# build mappings from alleles to markers, and vice versa...

	# dictionaries for allele/marker mappings as determined by sequence
	# overlaps:

	all2mrk = {}	# maps allele key to marker key
	mrk2all = {}	# maps marker key to list of allele keys

	for (seqKey, markers) in singleMarkers.items():
		allKeys = seqToAllele[seqKey]
		mrkKey = markers[0]

		for alleleKey in allKeys:
			# 1. if this allele is associated with a nomen marker
			# then skip it; 2. if we already have this allele key
			# and it's for another marker, then report the
			# conflict and just leave the original marker

			if nomenAlleles.has_key(alleleKey):
				LOGGER.log ('Has Nomen Marker',
					(alleleKey, symbols[alleleKey],
					nomenAlleles[alleleKey],
					mrkKey) )

			elif not all2mrk.has_key(alleleKey):
				all2mrk[alleleKey] = mrkKey
			
				# collect all allele keys for each marker

				if not mrk2all.has_key(mrkKey):
					mrk2all[mrkKey] = [ alleleKey ]
				elif alleleKey not in mrk2all[mrkKey]:
					mrk2all[mrkKey].append (alleleKey)

			elif all2mrk[alleleKey] != mrkKey:
				LOGGER.log ('Multiple Markers',
					(alleleKey, symbols[alleleKey],
					all2mrk[alleleKey],
					mrkKey) )

	# get the set of all alleles currently associated with each marker

	cmd = '''SELECT a._Assoc_key,
			a._Allele_key,
			a._Marker_key,
			a._Refs_key,
			t.term
		FROM ALL_Marker_Assoc a,
			VOC_Term t
		where a._Status_key = t._Term_key'''
	results = db.sql (cmd, 'auto')

	allAssoc = {}			# all associated alleles
	curAssoc = {}			# associations tagged by curator
	loadAssoc = {}			# load-associated alleles

	for row in results:
		# we assume that each allele is associated with -- at most --
		# one marker, so we only worry about collecting one row for
		# each allele

		allAssoc[row['_Allele_key']] = row

		# if the status of the association is either 'Curated' or
		# 'Curator Invalidated', the load should not alter those

		if row['term'][:5] == 'Curat':
			curAssoc[row['_Allele_key']] = row
		else:
			loadAssoc[row['_Allele_key']] = row

	LOGGER.log ('diag', 'Retrieved %d marker/allele associations' % \
		len(allAssoc))

	# Now, reconcile the currently-stored associations with those we need
	# to have.

	revisedAlleles = {}	# dictionary of allele keys w/revisions

	# 1. If a load association does not appear in singleMarkers, or if the
	#	allele appears in mixedAlleles, remove it from the database.

	toDelete = []		# list of assoc keys to delete

	msg = 'mixed allele %d assoc w/marker %d with refs %d'
	for (alleleKey, assoc) in loadAssoc.items():
		if mixedAlleles.has_key(alleleKey):
			toDelete.append (assoc['_Assoc_key'])
			revisedAlleles[alleleKey] = 1
		else:
			markerKey = assoc['_Marker_key']
			if (not mrk2all.has_key(markerKey)) or \
				(alleleKey not in mrk2all[markerKey]):
					toDelete.append (assoc['_Assoc_key'])
					revisedAlleles[alleleKey] = 1

	# 2. If a curated association appears in mixedAlleles, report it.
	
	for (alleleKey, assoc) in curAssoc.items():
		if mixedAlleles.has_key(alleleKey):
			LOGGER.log ('Curated Mixed',
				(alleleKey, assoc['_Marker_key']) )

	# 3. If an association appears in singleMarkers that is not in the
	#	set of all associations, and if the allele is not in 
	#	mixedAlleles, add it as a load association.

	toAdd = []	# list of (allele key, marker key) tuples to add

	for (alleleKey, markerKey) in all2mrk.items():
		if mixedAlleles.has_key(alleleKey):
			continue

		# if this is a novel association, we must add it

		if not allAssoc.has_key(alleleKey):
			toAdd.append ( (alleleKey, markerKey) )
			revisedAlleles[alleleKey] = 1
		else:
			# if the allele is already associated in the database
			# with a different marker, we should report it

			if allAssoc[alleleKey]['_Marker_key'] != markerKey:
				LOGGER.log ('Marker Mismatch',
					(alleleKey,
					allAssoc[alleleKey]['_Marker_key'],
					markerKey,
					allAssoc[alleleKey]['term']) )

	# apply the deletions

	sublists = splitList (toDelete)
	cmd = 'DELETE FROM ALL_Marker_Assoc WHERE _Assoc_key IN (%s)'
	for sublist in sublists:
		update (cmd % ','.join (map (str, sublist)) )
	LOGGER.log ('diag',
		'Removed %d previously loaded allele/marker associations' % \
		len(toDelete))

	LOGGER.log ('diag',
		'Left %d previously loaded allele/marker associations as-is' \
		% (len(loadAssoc) - len(toDelete)) )

	# apply the additions, getting a new assoc key at each step

	cmd1 = '''INSERT ALL_Marker_Assoc
			(_Assoc_key, _Allele_key, _Marker_key, _Qualifier_key,
			_Status_key)
		VALUES (%d, %d, %d, %d, %d)'''

	for (alleleKey, markerKey) in toAdd:
		assocKey = nextAssocKey()
		update (cmd1 % (assocKey, alleleKey, markerKey, notApp,
			loadedTerm) )

	LOGGER.log ('diag', 'Added %d new allele/marker associations' % \
		len(toAdd))

	changed = 0
	for (alleleKey, markerKey) in toAdd:
		if setMolecularNote (alleleKey, NOTE_A):
			changed = changed + 1

	LOGGER.log ('diag',
	'Set %d type A molecular notes (loaded association based on overlap)'\
		% changed)

	# if we applied any updates, we need to refresh the cached marker keys
	# in the allele table

	if toAdd or toDelete:
		update ('EXEC ALL_cacheMarker')
		LOGGER.log ('diag', 'Ran ALL_cacheMarker stored procedure')
		flagTable ('ALL_Marker_Assoc')
		flagTable ('ALL_Allele')

	return revisedAlleles.keys()

###------------------------------------------------------------------------###

def updateSymbols (
	alleles		# list of allele keys which need updated symbols
	):
	# Purpose: refresh the stored allele symbols to reflect current
	#	marker associations
	# Returns: nothing
	# Assumes: db module has been initialized
	# Effects: updates symbol field in ALL_Allele
	# Throws: propagates any exceptions from db.sql()

	debug ('in updateSymbols()')

	current = []	# list of dictionaries for current symbols (results
			# from 'cmd' below)

	cmd = '''SELECT a._Allele_key,
			a.symbol AS alleleSymbol,
			m.symbol AS markerSymbol,
			t.term AS status
		FROM ALL_Allele a,
			ALL_Marker_Assoc ma,
			MRK_Marker m,
			VOC_Term t
		WHERE a._Allele_key IN (%s)
			AND ma._Status_key *= t._Term_key
			AND a._Allele_key *= ma._Allele_key
			AND ma._Marker_key *= m._Marker_key'''

	# break 'alleles' into smaller lists so we don't go beyond Sybase
	# limits; run the queries, and join the results together in 'current'

	sublists = splitList(alleles)
	for sublist in sublists:
		results = db.sql (cmd % ','.join (map (str, sublist)), 'auto')
		current = current + results

	LOGGER.log ('diag', 'Retrieved %d allele symbols' % len(current))

	cmd = 'UPDATE ALL_Allele SET symbol = "%s" WHERE _Allele_key = %d'

	# now, loop through all the alleles and update each individually

	updated = 0	# count of alleles with revised symbols
	asIs = 0	# count of alleles with same symbols after revisions

	for row in current:
		alleleKey = row['_Allele_key']
		alleleSymbol = row['alleleSymbol']
		markerSymbol = row['markerSymbol']
		status = row['status']
		originalSymbol = alleleSymbol

		# We want to get down to the base allele symbol, so trim off
		# any trailing bracket and any marker symbol leading up to the
		# leading bracket.

		if alleleSymbol[-1] == '>':
			alleleSymbol = alleleSymbol[:-1]
		startPos = alleleSymbol.find('<')
		if startPos >= 0:
			alleleSymbol = alleleSymbol[startPos+1:]

		# If we have a valid marker association, prepend the marker
		# symbol and superscript the allele symbol portion.

		if (status != 'Curator Invalidated') and markerSymbol:
			alleleSymbol = '%s<%s>' % (markerSymbol, alleleSymbol)

		if (originalSymbol != alleleSymbol):
			update (cmd % (alleleSymbol, alleleKey) )
			updated = updated + 1
		else:
			asIs = asIs + 1

	if updated:
		flagTable ('ALL_Allele')
	LOGGER.log ('diag', 'Updated %d allele symbols' % updated)
	LOGGER.log ('diag', 'Left %d allele symbols as-is (for alleles with altered marker associations)' % asIs)
	return

###------------------------------------------------------------------------###

def anySeparate (
	seqs		# list of dictionaries; associated seqs for an allele
	):
	# Purpose: determine if any of the sequences given in 'seqs' do not
	#	overlap any of the other sequences within 50bp
	# Returns: True if there are any non-overlapping sequences, False if
	#	all are overlapping
	# Assumes: each dictionary in 'seqs' has these keys: 'chromosome',
	#	'startCoordinate', 'endCoordinate'
	# Effects: nothing
	# Throws: nothing

	# We need to do pairwise comparisons of all sequences in 'seqs'.  To
	# do this, we will walk two counters through 'seqs' in a nested loop.

	i = 0			# outer counter; first seq to compare
	threshold = 25		# amount to pad the ends of each sequence and
				# ...have them considered to be overlapping
	lenSeqs = len(seqs)	# number of sequences

	while i < lenSeqs:
		# information about first (leftmost) sequence
		seq1 = seqs[i]
		start1 = seq1['startCoordinate'] - threshold
		end1 = seq1['endCoordinate'] + threshold
		chrom1 = seq1['chromosome']

		# now walk through sequences to its right

		j = i + 1	# inner counter; steps through seqs after i
		while j < lenSeqs:
			seq2 = seqs[j]

			# easy case; if different chromosomes, no overlap
			if chrom1 != seq2['chromosome']:
				return True

			# check for non-overlapping coordinates
			start2 = seq2['startCoordinate'] - threshold
			end2 = seq2['endCoordinate'] + threshold
			if (start1 > end2) or (start2 > end1):
				return True

			j = j + 1	# try next seq to right
		i = i + 1		# move to next seq on left side

	return False

###------------------------------------------------------------------------###

def flagMixedAlleles():
	# Purpose: update the isMixed flag in the database
	# Returns: dictionary; keys are allele keys which have isMixed = 1
	#	in the database
	# Assumes: db module has been initialized
	# Effects: updates isMixed field in ALL_Allele table in the database
	# Throws: propagates any exceptions from db.sql()

	debug ('in flagMixedAlleles()')

	# Get alleles which have a curated isMixed flag; we identify these by
	# a "Mixed" reference which is only set by curators (not by the load).

	cmd = '''SELECT a._Allele_key,
			a.isMixed,
			ra._Refs_key
		FROM MGI_Reference_Assoc ra,
			MGI_RefAssocType rat,
			ALL_Allele a
		WHERE rat._MGIType_key = 11
			AND rat.assocType = "Mixed"
			AND rat._RefAssocType_key = ra._RefAssocType_key
			AND ra._Object_key = a._Allele_key'''

	curated = {}		# allele keys with curated isMixed flags
	mixed = {}		# allele keys for alleles with isMixed = 1

	results = db.sql (cmd, 'auto')
	for row in results:
		alleleKey = row['_Allele_key']
		curated[alleleKey] = 1

		# we will remember the curator-flagged mixed alleles first,
		# then pick up the load-flagged ones later on

		if row['isMixed'] == 1:
			mixed[alleleKey] = 1
			LOGGER.log ('Mixed Refs', 
				(row['_Allele_key'], row['_Refs_key']) )

	LOGGER.log ('diag', 'Found %d alleles with curated isMixed flags' % \
		len(curated))

	# Get all alleles, their associated sequences, the coordinates for
	# those sequences, and their gene trap tag methods.

	cmd = '''SELECT a._Allele_key,
			a.isMixed,
			cc._Sequence_key,
			cc.chromosome,
			cc.startCoordinate,
			cc.endCoordinate,
			gt._TagMethod_key
		FROM ALL_Allele a,
			SEQ_Allele_Assoc aa,
			SEQ_Coord_Cache cc,
			SEQ_GeneTrap gt
		WHERE a._Allele_key = aa._Allele_key
			AND aa._Sequence_key = gt._Sequence_key
			AND gt._Sequence_key = cc._Sequence_key'''
	results = db.sql (cmd, 'auto')

	alleles = {}
	for row in results:
		alleleKey = row['_Allele_key']

		# we want to ignore alleles which have curated isMixed flags;
		# these are ones with a Mixed reference

		if curated.has_key(alleleKey):
			continue

		# we only want to consider sequences with DNA-based seq tags

		if 3983002 <= row['_TagMethod_key'] <= 3983007:
			if alleles.has_key (alleleKey):
				alleles[alleleKey].append (row)
			else:
				alleles[alleleKey] = [ row ]

	LOGGER.log ('diag', 'Retrieved coordinates for %d alleles and their %d associated sequences' % (len(alleles), len(results)) )

	# An allele is mixed if (using only DNA-based seq tags) it...
	# 1. has >1 sequence, and
	# 2. at least one pair of sequences does not overlap within 50bp

	clearMixed = []		# list of allele keys to clear the isMixed bit
	setMixed = []		# list of allele keys to set the isMixed bit
	asIs = 0		# count of allele keys w/ correct isMixed bits

	for (allele, seqs) in alleles.items():
		mixedInDb = seqs[0]['isMixed']

		# alleles with only one sequence cannot be mixed; clear the
		# bit if it's already set in the database

		if len(seqs) == 1:
			if mixedInDb:
				clearMixed.append (allele)
			else:
				asIs = asIs + 1
		else:
			# allele has multiple sequences, so see if any of them
			# are non-overlapping.  Update the isMixed bit as
			# needed.

			mixedNow = anySeparate (seqs)
			if mixedNow and not mixedInDb:
				setMixed.append (allele)
			elif mixedInDb and not mixedNow:
				clearMixed.append (allele)
			else:
				asIs = asIs + 1

			# if this allele is now mixed, we need to remember it

			if mixedNow:
				mixed[allele] = 1
				LOGGER.log ("Mixed",
					(allele, seqs) )

	cmd = 'UPDATE ALL_Allele SET isMixed = %d WHERE _Allele_key IN (%s)'

	# update alleles which now need to be flagged as mixed

	for sublist in splitList (setMixed):
		update (cmd % (1, ','.join (map (str, sublist)) ) )
	LOGGER.log ('diag', 'Flagged %d alleles as now being mixed' % \
		len(setMixed) )

	changed = 0
	for alleleKey in setMixed:
		if setMolecularNote(alleleKey, NOTE_D):
			changed = changed + 1
	LOGGER.log ('diag', 'Set %d type D molecular notes (mixed allele)' % \
		changed)

	# clear alleles which are no longer mixed

	for sublist in splitList (clearMixed):
		update (cmd % (0, ','.join (map (str, sublist)) ) )
	LOGGER.log ('diag', 'Removed isMixed flag for %d alleles' % \
		len(clearMixed))

	LOGGER.log ('diag', 'Left %d isMixed flags as-is' % asIs)
	LOGGER.log ('diag', 'Found total of %d mixed alleles' % len(mixed))
	if setMixed or clearMixed:
		flagTable ('ALL_Allele')
	return mixed

###------------------------------------------------------------------------###

def setupMolecularNotes():
	# Purpose: initialize the various standard molecular notes that reside
	#	in global variables
	# Returns: nothing
	# Assumes: nothing
	# Effects: initializes globals NOTE_A, NOTE_B, NOTE_C, and NOTE_D
	# Throws: nothing

	global NOTE_A, NOTE_B, NOTE_C, NOTE_D

	debug ('in setupMolecularNotes()')

	# NOTE_A is set whenever we make a new ALO/marker association in the
	# load, and where there is not already a curated molecular note
	NOTE_A = splitString('''Trapped gene association determined by genome coordinate overlap and strand.''')

	# NOTE_B is set whenever an allele has a single good hit to the genome
	# but where there was not an overlap with a single marker (and where
	# there is not a curated molecular note)
	NOTE_B = splitString('''Sequence tags for this gene trap insertion do not overlap a single known gene in the trapped gene orientation. It is possible that sequence tags for this gene trap overlap multiple genes in the trapped gene orientation, one or more genes in reverse orientation (not disrupted by the gene trap mechanism), or no known genes (using current gene model coordinates).''')

	# NOTE_C is set whenever an allele has the number of good hits != 1
	# and where there is not already a curated molecular note
	NOTE_C = splitString('''Sequence tags for this gene trap insertion either do not align with the reference mouse genome, or align in multiple locations.''')

	# NOTE_D is set whenever an allele has been flagged as a mixed line,
	# and where there is not already a curated molecular note
	NOTE_D = splitString('''Mixed cell lines contain multiple mutant cell line clones mixed together in a single culture well. Mixed lines are characterized by multiple sequence tags for the same cell line that align to different locations in the genome.''')

	LOGGER.log ('diag', 'Initialized 4 in-memory molecular notes')
	return

###------------------------------------------------------------------------###

CURATED_NOTES = {}	# has allele key for each allele with curated mol.note
LOADED_NOTES = {}	# allele key -> (note key, note chunks) for alleles w/
			# ...loaded mol.notes
NEXT_NOTE_KEY = None	# next _Note_key to use when inserting new notes

def setupNoteCache ():
	# Purpose: cache all existing molecular notes for alleles in-memory
	# Returns: nothing
	# Assumes: we can query the database
	# Effects: queries the database, initializes globals CURATED_NOTES,
	#	LOADED_NOTES, and NEXT_NOTE_KEY
	# Throws: propagates all exceptions from db.sql()

	global CURATED_NOTES, LOADED_NOTES, NEXT_NOTE_KEY

	debug ('in setupNoteCache()')

	cmds = [
		# curated notes
		'''SELECT _Object_key
		FROM MGI_Note
		WHERE _NoteType_key = 1021
			AND _MGIType_key = 11
			AND _ModifiedBy_key NOT IN (SELECT _User_key
				FROM MGI_User
				WHERE login LIKE "%load")''',

		# loaded notes
		'''SELECT n._Object_key, n._Note_key, c.note
		FROM MGI_Note n,
			MGI_NoteChunk c
		WHERE n._NoteType_key = 1021
			AND n._MGIType_key = 11
			AND n._ModifiedBy_key IN (SELECT _User_key
				FROM MGI_User
				WHERE login LIKE "%load")
			AND n._Note_key = c._Note_key
		ORDER BY c._Note_key, c.sequenceNum''',

		# maximum note key stored so far
		'SELECT MAX(_Note_key) FROM MGI_Note',
		]

	[ results1, results2, results3 ] = db.sql (cmds, 'auto')

	# just set a flag if an allele has a curated molecular note; we must
	# not change these ones
	for row in results1:
		CURATED_NOTES[row['_Object_key']] = 1

	# if an allele has a loaded molecular note, then we should collect
	# the contents of the note to see if we need to update it or if it
	# is still correct
	for row in results2:
		objectKey = row['_Object_key']
		if LOADED_NOTES.has_key(objectKey):
			(key, notes) = LOADED_NOTES[objectKey]
			notes.append (row['note'])
		else:
			LOADED_NOTES[objectKey] = (row['_Note_key'],
				[ row['note'] ])

	# remember the next note key value that we should use as we add more
	# molecular notes
	NEXT_NOTE_KEY = results3[0][''] + 1
	LOGGER.log ('diag', 'Retrieved %d existing molecular notes' % \
		(len(CURATED_NOTES) + len(LOADED_NOTES)) )
	return

###------------------------------------------------------------------------###

def setMolecularNote (alleleKey, noteChunks):
	# Purpose: set the given molecular note (in 'noteChunks') for the
	#	specified 'alleleKey' in the database
	# Returns: True if we set the note for 'alleleKey', or False if we
	#	did not need to (because it matched the already-stored note)
	# Assumes: setupNoteCache() has been called
	# Effects: updates a molecular note in the database
	# Throws: propagates any exceptions from db.sql()

	global NEXT_NOTE_KEY

	# if this allele has a curated molecular note, skip the update
	if CURATED_NOTES.has_key(alleleKey):
		return False

	# if this allele already has a loaded molecular note, we can either
	# keep it (if the same) or delete its chunks (if different)
	if LOADED_NOTES.has_key(alleleKey):
		if LOADED_NOTES[alleleKey][1] == noteChunks:
			return False

		noteKey = LOADED_NOTES[alleleKey][0]
		update ('DELETE FROM MGI_NoteChunk WHERE _Note_key = %d' % \
			noteKey)
		createNote = False	# can re-use the existing note
	else:
		noteKey = NEXT_NOTE_KEY
		NEXT_NOTE_KEY = NEXT_NOTE_KEY + 1
		createNote = True	# need a new note

	if createNote:
		cmd = '''INSERT MGI_Note (_Note_key, _Object_key,
			_MGIType_key, _NoteType_key)
			VALUES (%d, %d, 11, 1021)''' % (noteKey, alleleKey)
		update(cmd)
		flagTable ('MGI_Note')

	i = 1
	for chunk in noteChunks:
		cmd = '''INSERT MGI_NoteChunk (_Note_key, sequenceNum, note)
			VALUES (%d, %d, "%s")''' % (noteKey, i, chunk)
		update(cmd)
		i = i + 1
	flagTable ('MGI_NoteChunk')
	return True

###------------------------------------------------------------------------###

def splitString (s):
	# Purpose: split string 's' into multiple chunks of length <= 255
	# Returns: list of strings, each of 255 characters or less; only the
	#	last string will be less than 255 characters (any before the
	#	last will be the full 255)
	# Assumes: nothing
	# Effects: nothing
	# Throws: nothing

	chunks = []
	while len(s) > 255:
		chunks.append (s[:255])
		s = s[255:]
	if s:
		chunks.append(s)
	return chunks

###------------------------------------------------------------------------###

def fixMolecularNotes():
	# Purpose: assign type B and C molecular notes for alleles, which are
	#	for alleles without marker associations and are based on the
	#	goodHitCount for the allele's associated sequence
	# Returns: nothing
	# Assumes: updates have been applied to ALL_Marker_Assoc, not just
	#	collected in an output file of SQL
	# Effects: updates molecular notes in the database
	# Throws: propagates any exceptions from db.sql()

	debug ('in fixMolecularNotes()')

	cmd = '''SELECT a._Allele_key, t.goodHitCount
		FROM SEQ_GeneTrap t,
			SEQ_Allele_Assoc a
		WHERE t._Sequence_key = a._Sequence_key
			AND NOT EXISTS (SELECT 1 FROM ALL_Marker_Assoc m
				WHERE a._Allele_key = m._Allele_key)'''
	results = db.sql (cmd, 'auto')

	countB = 0
	countC = 0

	for row in results:
		if row['goodHitCount'] == 1:
			if setMolecularNote (row['_Allele_key'], NOTE_B):
				countB = countB + 1
		elif setMolecularNote (row['_Allele_key'], NOTE_C):
			countC = countC + 1
			
	LOGGER.log ('diag', 'Set %d type B molecular notes (single good hit to genome, no marker association)' % countB)
	LOGGER.log ('diag', 'Set %d type C molecular notes (zero or multiple good hits to genome)' % countC)
	return

###------------------------------------------------------------------------###

def reportMixed(seqMrkDict):
	# Purpose: do reporting of mixed alleles, their associated sequences,
	#	and the markers which overlap those sequences
	# Returns: nothing
	# Assumes: updates have been applied to ALL_Marker_Assoc, not just
	#	collected in an output file of SQL
	# Effects: queries the database
	# Throws: propagates any exceptions from db.sql()

	debug ('in reportMixed()')

	cmd = '''SELECT a._Allele_key,
			s._Sequence_key
		FROM ALL_Allele a,
			SEQ_Allele_Assoc s
		WHERE a._Allele_key = s._Allele_key
			AND a.isMixed = 1'''
	results = db.sql (cmd, 'auto')

	alleles = {}
	for row in results:
		alleleKey = row['_Allele_key']
		if alleles.has_key(alleleKey):
			alleles[alleleKey].append (row['_Sequence_key'])
		else:
			alleles[alleleKey] = [ row['_Sequence_key'] ]

	for (alleleKey, seqKeys) in alleles.items():
		seqs = []
		for seqKey in seqKeys:
			if seqMrkDict.has_key(seqKey):
				seqs.append ( (seqKey, seqMrkDict[seqKey]) )
			else:
				seqs.append ( (seqKey, []) )

		LOGGER.log ('Mixed With Markers', (alleleKey, seqs) )

	LOGGER.log ('diag',
		'Logged sequence and marker overlaps for mixed alleles')
	return

###------------------------------------------------------------------------###

def writeException (
	excInfo,	# three-item tuple, as from sys.exc_info()
	fp		# file pointer to which to write
	):
	# Purpose: record information about the given exception to 'fp'
	# Returns: nothing
	# Assumes: nothing
	# Effects: writes to 'fp'
	# Throws: any exceptions raised by writing to 'fp'

	fp.write ('<HR>\n')
	fp.write ('<H3>Error: %s : %s</H3>' % (excInfo[0], excInfo[1]))
	fp.write ('<B>Traceback</B>:<PRE>\n')
	traceback.print_tb (excInfo[2], None, fp)
	fp.write ('</PRE>\n')
	fp.write ('<HR>\n')
	return

###------------------------------------------------------------------------###

def writeLogFile (
	excInfo,	 	# 3-item tuple, as from sys.exc_info(), to be
				# ...recorded in the log file
	filename, 		# string; filename to which to write
	logs, 			# list of tuples: (log name, title, writer)
	title			# string; overall title for page to write
	):
	# Purpose: to write a log file, containing various log tables, to
	#	a given filename
	# Returns: nothing
	# Assumes: we can write to 'filename'
	# Effects: writes to 'filename' in the file system
	# Throws: any exceptions raised as we try to write to 'filename'

	file = os.path.join (LOGDIR, filename)
	try:
		fp = open (file, 'w')
	except:
		msg = 'Cannot write log file to: %s' % file
		print msg
		sys.stderr.write (msg + '\n')

		file = os.path.join ('/tmp/', filename)
		fp = open (file, 'w')

		msg2 = 'Wrote instead to: %s' % file
		print msg2
		sys.stderr.write (msg2 + '\n')

	fp.write ('<HTML><HEAD><TITLE>%s</TITLE></HEAD>\n' % title)
	fp.write ('<BODY><H2>%s</H2>\n' % title)

	if excInfo:
		writeException (excInfo, fp)

	fp.write ('<OL>\n')
	for (logName, logTitle, logWriter) in logs:
		fp.write ('<LI><A HREF="#%s">%s</A></LI>\n' % (
			logTitle, logTitle) )
	fp.write ('</OL>\n')
	fp.write ('<HR>\n')
	for (logName, logTitle, logWriter) in logs:
		LOGGER.writeLog (logName, logTitle, logWriter, fp)
		fp.write ('<HR>\n')

	fp.write ('</BODY></HTML>\n')
	fp.close()
	return

###------------------------------------------------------------------------###

def writeCuratorLogs(excInfo = None):
	# Purpose: write the curatorial logs for this process, including any
	#	exception (in 'excInfo') that caused its termination
	# Returns: nothing
	# Assumes: we can write to the file system
	# Effects: writest to the file system
	# Throws: any exceptions raised by writing to the file system

	logs = [
		('Has Nomen Marker',
		'Alleles with Nomen Markers which overlap a current marker',
		aloMarkerLogger.NomenMarkerWriter() ),

		('Multiple Markers',
		'Alleles with multiple Marker associations by overlap',
		aloMarkerLogger.MultiMarkerWriter() ),

		('Curated Mixed',
		'Mixed Alleles with curated Marker associations',
		aloMarkerLogger.CuratedAssocMixedWriter() ),

		('Marker Mismatch',
		'Alleles associated with a Marker other than one they overlap',
		aloMarkerLogger.MultiMarkerStatusWriter() ),

		('Mixed Refs',
		'Alleles with curated isMixed flag',
		aloMarkerLogger.MixedRefsWriter() ),

		('Mixed',
		'Alleles flagged by load as being mixed',
		aloMarkerLogger.MixedSeqWriter() ),

		('Mixed With Markers',
		'Mixed Alleles with their sequences and overlapping markers',
		aloMarkerLogger.MixedMarkerWriter() ),
		]

	writeLogFile (excInfo, 'alomrkload.cur.html', logs,
		'ALO/Marker Assoc Curator Logs')
	return

###------------------------------------------------------------------------###

def writeDiagLog(excInfo = None):
	# Purpose: write the diagnostic logs for this process, including any
	#	exception (in 'excInfo') that caused its termination
	# Returns: nothing
	# Assumes: we can write to the file system
	# Effects: writest to the file system
	# Throws: any exceptions raised by writing to the file system

	logs = [
		('diag', 'Diagnostic Log', 
			aloMarkerLogger.ProfilingLogWriter() ),
		('debug', 'Debug Log', 
			aloMarkerLogger.AloMarkerLogWriterHTML() ),
		]

	writeLogFile (excInfo, 'alomrkload.diag.html', logs, 
		'ALO/Marker Assoc Diagnostic Logs')
	return

###------------------------------------------------------------------------###

def writeLogs(errorString = None):
	# Purpose: write the various logs for this process
	# Returns: nothing
	# Assumes: we can write to the file system
	# Effects: writes files to the file system
	# Throws: propagates any exceptions raised while writing to the file
	#	system

	debug ('in writeLogs()')
	writeCuratorLogs(errorString)
	writeDiagLog(errorString)
	return

###------------------------------------------------------------------------###

def main():
	# Purpose: main program
	# Returns: nothing
	# Assumes: nothing
	# Effects: nothing
	# Throws: propagates all exceptions

	processCommandLine()			# process command-line args
	initLogger()				# initialize logger module
	setPointCoordinates()			# set seqs' point coordinates
	seqLengths = findSequenceLengths()	# get all GT seq lengths
	setupMolecularNotes()			# intialize global notes
	setupNoteCache()			# cache existing mol.notes

	# find overlaps between sequences and markers (using fjoin for
	# performance)

	seqFile = getSequences()
	mrkFile = getMarkers()
	fjoinFile = fjoin (seqFile, mrkFile)

	# read fjoin results as a list of (seq, marker) pairs, then convert
	# to a dictionary mapping sequence key to list of marker keys

	seqMrkPairs = readFjoinFile (fjoinFile)
	seqMrkDict = groupby (seqMrkPairs)

	# get dictionaries of sequence keys which map to single markers,
	# multiple markers, and no markers

	singleMarkers, multiMarkers = splitSingles (seqMrkDict)
	noMarkers = findSequencesNoMarkers (seqMrkDict, seqLengths)

	# update the isMixed flag for all alleles and get a dictionary
	# of allele keys for mixed alleles; also handles setting molecular
	# note type D

	mixedAlleles = flagMixedAlleles() 
	reportMixed (seqMrkDict)

	# update the choice of representative sequence for each allele

	seqToAllele = updateRepSeqs (singleMarkers, multiMarkers, noMarkers,
		mixedAlleles, seqLengths)

	# update the marker/allele associations, based on sequence overlaps;
	# also handle molecular note type A

	revisedAlleles = updateMarkerAssoc (singleMarkers, seqToAllele,
		mixedAlleles)

	# update allele symbols to incorporate the new marker associations

	updateSymbols (revisedAlleles)

	# set the B and C molecular notes for alleles which are not associated
	# with a marker

	fixMolecularNotes()

	# update statistics on any tables we altered, to boost performance

	optimizeTables()

	# close of the output file of SQL statements, if we have one

	if SQL_FILE:
		SQL_FILE.close()

	# delete all temp files, to clean up after ourselves

	for file in TEMP_FILES:
		os.remove (file)
	LOGGER.log ('diag', 'Removed temp files')

	# write logs...

	writeLogs()
	return

###--------------------###
###--- main program ---###
###--------------------###

if __name__ == '__main__':
	try:
		main()
	except:
		writeLogs (sys.exc_info())
