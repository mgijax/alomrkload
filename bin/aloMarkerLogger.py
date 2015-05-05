# Name: aloMarkerLogger.py
# Purpose: Python library to handle logging for the ALO/Marker Association
#	load.  Deals with collecting log items, timestamping them, and
#	providing a flexible mechanism for writing them.
# Author: jsb
# Assumes: nothing

import mgi_html		# used for its doSubSupTags() function
import time		# used for its time.time() function (for timestamps)
import sys

###------------------------------------------------------------------------###

###------------------------###
###--- Global Variables ---###
###------------------------###

DEBUG = False			# set to true for debugging output to stderr
START_TIME = time.time()	# time at which the module is imported
ALLELE_SYMBOLS = {}		# maps allele key to allele symbol
MARKER_SYMBOLS = {}		# maps marker key to marker symbol
SEQ_IDS = {}			# maps sequence key to seq ID
JNUMS = {}			# maps reference key to J: number

# relevant URLs for linking to MGI detail pages; append a key to each to make
# a full, valid URL

ALLELE_DETAIL = 'http://www.informatics.jax.org/allele/key/'
MARKER_DETAIL = 'http://www.informatics.jax.org/marker/key/'
REFERENCE_DETAIL = 'http://www.informatics.jax.org/reference/key/'
SEQUENCE_DETAIL = 'http://www.informatics.jax.org/sequence/key/'

###------------------------------------------------------------------------###

###-----------------###
###--- Functions ---###
###-----------------###

def setAlleleSymbols (dict):
	# Purpose: to remember 'dict' as the mapping from allele key to symbol
	# Returns: nothing
	# Assumes: keys of 'dict' are integer allele keys, and the values are
	#	string allele symbols
	# Effects: updated global 'ALLELE_SYMBOLS'
	# Throws: nothing

	global ALLELE_SYMBOLS
	ALLELE_SYMBOLS = dict
	return

def setMarkerSymbols (dict):
	# Purpose: to remember 'dict' as the mapping from marker key to symbol
	# Returns: nothing
	# Assumes: keys of 'dict' are integer marker keys, and the values are
	#	string marker symbols
	# Effects: updated global 'MARKER_SYMBOLS'
	# Throws: nothing

	global MARKER_SYMBOLS
	MARKER_SYMBOLS = dict
	return

def setSeqIDs (dict):
	# Purpose: to remember 'dict' as the mapping from sequence key to ID
	# Returns: nothing
	# Assumes: keys of 'dict' are integer sequence keys, and the values
	#	are string seq IDs
	# Effects: updated global 'SEQ_IDS'
	# Throws: nothing

	global SEQ_IDS
	SEQ_IDS = dict
	return

def setReferences (dict):
	# Purpose: to remember 'dict' as the mapping from reference key to its
	#	J: number ID
	# Returns: nothing
	# Assumes: keys of 'dict' are integer reference keys, and the values
	#	are string J: numbers
	# Effects: updated global 'JNUMS'
	# Throws: nothing

	global JNUMS
	JNUMS = dict
	return

def alleleSymbol (key):
	# Purpose: get the allele symbol for the given integer allele 'key'
	# Returns: string (allele symbol) or integer (allele key, if symbol is
	#	not known)
	# Assumes: setAlleleSymbols() has been called to set up the cache of
	#	allele key/symbol mappings
	# Effects: nothing
	# Throws: nothing

	if ALLELE_SYMBOLS.has_key(key):
		return mgi_html.doSubSupTags(ALLELE_SYMBOLS[key])
	return key

def markerSymbol (key):
	# Purpose: get the marker symbol for the given integer marker 'key'
	# Returns: string (marker symbol) or integer (marker key, if symbol is
	#	not known)
	# Assumes: setMarkerSymbols() has been called to set up the cache of
	#	marker key/symbol mappings
	# Effects: nothing
	# Throws: nothing

	if MARKER_SYMBOLS.has_key(key):
		return MARKER_SYMBOLS[key]
	return key

def seqID (key):
	# Purpose: get the seq ID for the given integer sequence 'key'
	# Returns: string (seq ID) or integer (sequence key, if seq ID is
	#	not known)
	# Assumes: setSeqIDs() has been called to set up the cache of
	#	sequence key/ID mappings
	# Effects: nothing
	# Throws: nothing

	if SEQ_IDS.has_key(key):
		return SEQ_IDS[key]
	return key

def jnum (key):
	# Purpose: get the J: number for the given integer reference 'key'
	# Returns: string (J: number) or integer (reference key, if the J:
	#	number is not known)
	# Assumes: setReferences() has been called to set up the cache of
	#	reference key/J: number mappings
	# Effects: nothing
	# Throws: nothing

	if JNUMS.has_key(key):
		return JNUMS[key]
	return key

###------------------------------------------------------------------------###

###---------------###
###--- Classes ---###
###---------------###

class AloMarkerLog:
	# Is: an ordered collection of timestamped log entries
	# Has: an ordered collection of timestamped log entries; these entries
	#	are not required to be of a particular type
	# Does: accepts entries and timestamps them; returns the current set
	#	of entries with their timestamps

	def __init__ (self):
		# Purpose: constructor
		# Returns: nothing
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing

		# list of tuples: (elapsed time in seconds, log entry)
		# ordered by ascending time
		self.entries = []
		return

	def addEntry (self, entry):
		# Purpose: add an entry to the log
		# Returns: nothing
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing

		# timestamp the entry and add it to the end of the list
		self.entries.append ( (time.time() - START_TIME, entry) )
		return

	def getEntries (self):
		# Purpose: retrieve the entries collected so far
		# Returns: list of (elapsed time, log entry) tuples
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing
		# Notes: A reference to the actual list is returned; if you
		#	modify the list, then you may break the AloMarkerLog.
		#	Better to treat it as being read-only.

		return self.entries

###------------------------------------------------------------------------###

class AloMarkerLogWriter:
	# Is: a writer for an AloMarkerLog
	# Has: nothing
	# Does: writes entries from an AloMarkerLog to a given file pointer
	#	with a given format and a given title (see write() method)
	# Notes: This base class simply writes the log entries in a plain
	#	text format with no timestamping, one entry per line.

	def __init__ (self):
		# Purpose: constructor
		# Returns: nothing
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing

		return

	def write (self,
		fp,	# file pointer to which to write the log
		title,	# string; title for the log
		log	# AloMarkerLog; which has the entries to write
		):
		# Purpose: write the given 'log' to 'fp' with a given 'title'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes to 'fp'
		# Throws: propagates all exceptions from writing to 'fp'

		fp.write ('%s\n' %  title)
		fp.write ('%s\n' % ('-' * len(title)) )
		if log:
			for (timestamp, entry) in log.getEntries():
				self.writeEntry (fp, entry)
		else:
			self.write ('No entries\n')
		return

	def writeEntry (self,
		fp,		# file pointer to which to write the log
		entry		# one entry to write
		):
		# Purpose: (private) write the given log 'entry' to 'fp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes to 'fp'
		# Throws: propagates all exceptions from writing to 'fp'

		fp.write ('%s\n' % str(entry))
		return

###------------------------------------------------------------------------###

class AloMarkerLogWriterHTML (AloMarkerLogWriter):
	# Is: an AloMarkerLogWriter which specializes in writing the log
	#	entries in an HTML table (for inclusion in an HTML document)
	# Has: see AloMarkerLogWriter
	# Does: see AloMarkerLogWriter

	def write (self,
		fp,	# file pointer to which to write the log
		title,	# string; title for the log
		log	# AloMarkerLog; which has the entries to write
		):
		# Purpose: write the given 'log' to 'fp' with a given 'title'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes to 'fp'
		# Throws: propagates all exceptions from writing to 'fp'

		self.preProcess(log)
		fp.write ('<A NAME="%s"><H3>%s</H3></A><P>\n' % (title, title))
		if log:
			fp.write ('<TABLE BORDER="1">\n')
			self.writeHeader(fp)
			for (timestamp, entry) in log.getEntries():
				self.writeTimedEntry (fp, timestamp, entry)
			fp.write ('</TABLE>\n')
		else:
			fp.write ('No entries<P>\n')
		return

	def preProcess(self, log):
		# Purpose: (private) do any pre-processing necessary for the
		#	given 'log'
		# Returns: nothing
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing
		# Notes: This is just a stub, to be filled in as-needed by
		#	subclasses.

		return

	def writeHeader (self, fp):
		# Purpose: (private) write an HTML table header row to 'fp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'

		fp.write ('<TR><TH>Entry</TH></TR>\n')
		return

	def writeTimedEntry (self, fp, timestamp, entry):
		# Purpose: (private) write a log 'entry' as an HTML table row
		#	to 'fp' with any necessary info from 'timestamp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'
		# Notes: In this class, this is just a wrapper over the
		#	writeEntry() method.  It will be filled out as
		#	needed in subclasses.

		self.writeEntry (fp, entry)
		return

	def writeEntry (self, fp, entry):
		# Purpose: (private) write a log 'entry' to 'fp' as an HTML
		#	table row
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'

		fp.write ('<TR><TD>%s</TD></TR>\n' % str(entry))
		return

###------------------------------------------------------------------------###

class ProfilingLogWriter (AloMarkerLogWriterHTML):
	# Is: a subclass of AloMarkerLogWriterHTML which analyzes the
	#	timestamps for the log entries and provides simple profiling
	#	output
	# Has: see AloMarkerLogWriter
	# Does: see AloMarkerLogWriter

	def __init__ (self):
		# Purpose: constructor
		# Returns: nothing
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing

		AloMarkerLogWriterHTML.__init__ (self)
		self.lastTime = 0.0
		self.perMark = 5.0
		return

	def preProcess (self, log):
		# Purpose: (private) do any pre-processing necessary for the
		#	given 'log'
		# Returns: nothing
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing
		# Notes: We must look for the maximum elapsed time between any
		#	two log entries, so we can use that time difference to
		#	compute the scale for our histogram.  (to know how
		#	many seconds per asterisk, if we use 20 asterisks for
		#	the maximum difference)

		if log == None:
			return

		lastTime = None
		diff = 0.0
		maxTime = 0.0
		for (timestamp, entry) in log.getEntries():
			if lastTime == None:
				diff = timestamp
			else:
				diff = timestamp - lastTime
			lastTime = timestamp
			if diff > maxTime:
				maxTime = diff

		# maximum of twenty asterisks per line
		self.perMark = maxTime / 20.0

		# prevent divide-by-0 errors
		if self.perMark == 0:
			self.perMark = 1
		return

	def writeHeader (self, fp):
		# Purpose: (private) write an HTML table header row to 'fp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'

		fp.write ('<TR><TH>Total (sec)</TH><TH>Elapsed (sec)</TH>')
		fp.write ('<TH>Histogram</TH><TH>Entry</TH></TR>\n')
		return

	def writeTimedEntry (self, fp, timestamp, entry):
		# Purpose: (private) write a log 'entry' as an HTML table row
		#	to 'fp' with any necessary info from 'timestamp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'
		# Notes: includes four columns (elapsed time from start of
		#	script, elapsed time since last log entry, histogram
		#	for value from column 2, the log entry itself)

		diff = timestamp - self.lastTime
		self.lastTime = timestamp
		histogram = int(round(diff / self.perMark)) * '*'
		if histogram == '':
			histogram = '&nbsp;'
		fp.write ('<TR><TD ALIGN="right">%8.3f</TD>' % timestamp)
		fp.write ('<TD ALIGN="right">%8.3f</TD>' % diff)
		fp.write ('<TD ALIGN="left">%s</TD>' % histogram)
		fp.write ('<TD>%s</TD></TR>\n' % entry) 
		return

###------------------------------------------------------------------------###

class NomenMarkerWriter (AloMarkerLogWriterHTML):
	# Is: an AloMarkerLogWriterHTML subclass that specializes in writing
	#	an HTML table showing alleles which have an association with
	#	a nomen marker (a pre-publication marker) and which overlaps
	#	an already-public marker according to its coordinates
	# Has: see AloMarkerLogWriter
	# Does: see AloMarkerLogWriter
	# Assumes: log entries are tuples with four items: allele key,
	#	allele symbol, nomen marker symbol, marker key

	def writeHeader (self, fp):
		# Purpose: (private) write an HTML table header row to 'fp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'

		fp.write ('<TR><TH>Allele</TH><TH>Nomen Marker</TH>')
		fp.write ('<TH>Overlap Marker</TR>\n')
		return

	def writeEntry (self, fp, entry):
		# Purpose: (private) write a log 'entry' to 'fp' as an HTML
		#	table row
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'
		# Notes: writes a three-column row (allele symbol linked to
		#	allele detail page, nomen symbol, marker symbol linked
		#	to marker detail page)

		fp.write ('<TR><TD><A HREF="%s%s">%s</A></TD>' % (
			ALLELE_DETAIL, entry[0],
			mgi_html.doSubSupTags(entry[1]) ) )
		fp.write ('<TD>%s</TD><TD><A HREF="%s%s">%s</A></TD></TR>' % (
			entry[2],
			MARKER_DETAIL, entry[3], markerSymbol(entry[3])) )
		return

###------------------------------------------------------------------------###

class MultiMarkerWriter (AloMarkerLogWriterHTML):
	# Is: an AloMarkerLogWriterHTML subclass that specializes in writing
	#	an HTML table showing alleles which would map to multiple
	#	markers, according to its coordinates
	# Has: see AloMarkerLogWriter
	# Does: see AloMarkerLogWriter
	# Assumes: log entries are tuples with four items:  allele key,
	#	sequence key, list of marker keys

	def writeHeader (self, fp):
		# Purpose: (private) write an HTML table header row to 'fp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'

		fp.write ('<TR><TH>Allele</TH><TH>Representative Sequence</TH>')
		fp.write ('<TH>Overlapped Markers</TH></TR>\n')
		return

	def writeEntry (self, fp, entry):
		# Purpose: (private) write a log 'entry' to 'fp' as an HTML
		#	table row
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'
		# Notes: writes a three-column row (allele symbol linked to
		#	allele detail page, sequence ID linked to sequence
		#	detail page, overlapped marker symbols linked to
		#	marker detail pages)

		fp.write ('<TR><TD><A HREF="%s%s">%s</TD>' % (
			ALLELE_DETAIL, entry[0], alleleSymbol(entry[0])) )
		fp.write ('<TD><A HREF="%s%s">%s</A></TD>' % (
			SEQUENCE_DETAIL, entry[1], seqID(entry[1])) )
		fp.write ('<TD>')
		for mrkKey in entry[2]:
			fp.write ('<A HREF="%s%s">%s</A>' % (
				MARKER_DETAIL, mrkKey, markerSymbol(mrkKey)) )
			if mrkKey != entry[2][-1]:
				fp.write ('<BR>')
		fp.write ('</TD></TR>\n')
		return

###------------------------------------------------------------------------###

class MultiMarkerStatusWriter (AloMarkerLogWriterHTML):
	# Is: an AloMarkerLogWriterHTML subclass that specializes in writing
	#	an HTML table showing alleles which have an association with
	#	an existing marker (with the status of the association) and
	#	which overlaps a different marker according to its coordinates
	# Has: see AloMarkerLogWriter
	# Does: see AloMarkerLogWriter
	# Assumes: log entries are tuples with four items:  allele key,
	#	existing marker key, overlapped marker key, association status

	def writeHeader (self, fp):
		# Purpose: (private) write an HTML table header row to 'fp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'

		fp.write ('<TR><TH>Allele</TH>')
		fp.write ('<TH>Status : Current Marker</TH>')
		fp.write ('<TH>Overlap Marker</TH></TR>\n')
		return

	def writeEntry (self, fp, entry):
		# Purpose: (private) write a log 'entry' to 'fp' as an HTML
		#	table row
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'
		# Notes: writes a three-column row (allele symbol linked to
		#	allele detail page, association status with existing
		#	associated marker symbol linked to marker detail page,
		#	overlapped marker symbol linked to marker detail page)

		fp.write ('<TR><TD><A HREF="%s%s">%s</TD>' % (
			ALLELE_DETAIL, entry[0], alleleSymbol(entry[0])) )
		fp.write ('<TD>%s : <A HREF="%s%s">%s</A></TD>' % (
			entry[3],
			MARKER_DETAIL, entry[1], markerSymbol(entry[1])) )
		fp.write ('<TD><A HREF="%s%s">%s</A></TD></TR>\n' % (
			MARKER_DETAIL, entry[2], markerSymbol(entry[2])) )
		return

###------------------------------------------------------------------------###

class CuratedAssocMixedWriter (AloMarkerLogWriterHTML):
	# Is: an AloMarkerLogWriterHTML subclass that specializes in writing
	#	an HTML table showing alleles which have a curated association
	#	with an existing marker (often used for mixed alleles)
	# Has: see AloMarkerLogWriter
	# Does: see AloMarkerLogWriter
	# Assumes: log entries are tuples with two items:  allele key,
	#	marker key

	def writeHeader (self, fp):
		# Purpose: (private) write an HTML table header row to 'fp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'

		fp.write ('<TR><TH>Allele</TH><TH>Curated Marker</TH></TR>')
		return

	def writeEntry (self, fp, entry):
		# Purpose: (private) write a log 'entry' to 'fp' as an HTML
		#	table row
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'
		# Notes: writes a two-column row (allele symbol linked to
		#	allele detail page, marker symbol linked to marker
		#	detail page)

		fp.write ('<TR><TD><A HREF="%s%s">%s</A></TD>' % (
			ALLELE_DETAIL, entry[0], alleleSymbol(entry[0])) )
		fp.write ('<TD><A HREF="%s%s">%s</A></TD></TR>\n' % (
			MARKER_DETAIL, entry[1], markerSymbol(entry[1])) )
		return

###------------------------------------------------------------------------###

class MixedRefsWriter (AloMarkerLogWriterHTML):
	# Is: an AloMarkerLogWriterHTML subclass that specializes in writing
	#	an HTML table showing alleles which have a curated mixed flag
	# Has: see AloMarkerLogWriter
	# Does: see AloMarkerLogWriter
	# Assumes: log entries are tuples with two items:  allele key,
	#	reference key

	def writeHeader (self, fp):
		# Purpose: (private) write an HTML table header row to 'fp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'

		fp.write ('<TR><TH>Allele</TH><TH>Mixed Reference</TH></TR>')
		return

	def writeEntry (self, fp, entry):
		# Purpose: (private) write a log 'entry' to 'fp' as an HTML
		#	table row
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'
		# Notes: writes a two-column row (allele symbol linked to
		#	allele detail page, J: number linked to reference
		#	detail page)

		fp.write ('<TR><TD><A HREF="%s%s">%s</TD>' % (
			ALLELE_DETAIL, entry[0], alleleSymbol(entry[0])) )
		fp.write ('<TD><A HREF="%s%s">%s</A></TD></TR>\n' % (
			REFERENCE_DETAIL, entry[1], jnum(entry[1])) )
		return

###------------------------------------------------------------------------###

class MixedSeqWriter (AloMarkerLogWriterHTML):
	# Is: an AloMarkerLogWriterHTML subclass that specializes in writing
	#	an HTML table showing alleles which have been flagged by the
	#	load as being mixed, along with their associated sequences
	# Has: see AloMarkerLogWriter
	# Does: see AloMarkerLogWriter
	# Assumes: log entries are tuples with two items:  allele key,
	#	list of dictionaries with sequence key, chromosome, coords

	def writeHeader (self, fp):
		# Purpose: (private) write an HTML table header row to 'fp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'

		fp.write ('<TR><TH>Allele</TH><TH>Sequences</TH></TR>')
		return

	def writeEntry (self, fp, entry):
		# Purpose: (private) write a log 'entry' to 'fp' as an HTML
		#	table row
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'
		# Notes: writes a two-column row (allele symbol linked to
		#	allele detail page, seq ID linked to its sequence
		#	detail page, plus seq coordinates in the same column)

		fp.write ('<TR><TD><A HREF="%s%s">%s</TD><TD>' % (
			ALLELE_DETAIL, entry[0], alleleSymbol(entry[0])) )
		for row in entry[1]:
			fp.write ('<A HREF="%s%s">%s</A> : chr%s : %s-%s<BR>'\
				% (SEQUENCE_DETAIL, row['_Sequence_key'],
				seqID(row['_Sequence_key']),
				row['chromosome'],
				row['startCoordinate'], row['endCoordinate']))
		fp.write ('</TD></TR>\n')
		return

###------------------------------------------------------------------------###

class MixedMarkerWriter (AloMarkerLogWriterHTML):
	# Is: an AloMarkerLogWriterHTML subclass that specializes in writing
	#	an HTML table showing alleles which have been flagged as being
	#	mixed, along with their associated sequences and overlapping
	#	markers
	# Has: see AloMarkerLogWriter
	# Does: see AloMarkerLogWriter
	# Assumes: log entries are tuples with two items:  allele key,
	#	list of tuples (seq key, list of overlapping marker keys)

	def writeHeader (self, fp):
		# Purpose: (private) write an HTML table header row to 'fp'
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'

		fp.write ('<TR><TH>Allele</TH><TH>Seq : Markers</TH></TR>')
		return

	def writeEntry (self, fp, entry):
		# Purpose: (private) write a log 'entry' to 'fp' as an HTML
		#	table row
		# Returns: nothing
		# Assumes: we can write to 'fp'
		# Effects: writes a line to 'fp'
		# Throws: propagates any exceptions from writing to 'fp'
		# Notes: writes a two-column row (allele symbol linked to
		#	allele detail page, seq ID linked to its sequence
		#	detail page, plus marker symbols linked to their
		#	marker detail pages in the same column)

		fp.write ('<TR><TD><A HREF="%s%s">%s</TD><TD>' % (
			ALLELE_DETAIL, entry[0], alleleSymbol(entry[0])) )
		for (seqKey, markers) in entry[1]:
			fp.write ('<A HREF="%s%s">%s</A> : ' % (
				SEQUENCE_DETAIL, seqKey, seqID(seqKey)) )
			mrkList = []
			for markerKey in markers:
				mrkList.append ('<A HREF="%s%s">%s</A>' % \
					(MARKER_DETAIL, markerKey,
						markerSymbol(markerKey)) )
			fp.write (', '.join (mrkList))
			fp.write ('<BR>\n')
		fp.write ('</TD></TR>\n')
		return

###------------------------------------------------------------------------###

class AloMarkerLogManager:
	# Is: a managing class for AloMarkerLog objects, providing easy access
	#	to multiple logs
	# Has: an unordered collection of AloMarkerLog objects
	# Does: provides methods to write to a named log and to write the
	#	contents of a named log

	def __init__ (self):
		# Purpose: constructor
		# Returns: nothing
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing

		self.logs = {}
		return

	def log (self,
		logName, 	# string; name of log to which to send 'entry'
		entry		# any type; the entry to be logged
		):
		# Purpose: send a given log 'entry' to a named AloMarkerLog
		#	object
		# Returns: nothing
		# Assumes: nothing
		# Effects: nothing
		# Throws: nothing
		# Notes: The type of 'entry' is dictated by what sort of
		#	AloMarkerLogWriter you will want to use to write the
		#	log later on.

		# if we don't already have a log by that name, create new one

		if not self.logs.has_key(logName):
			self.logs[logName] = AloMarkerLog()
		self.logs[logName].addEntry (entry)
		if DEBUG:
			sys.stderr.write ('%9.3f : %s : %s\n' % (
				time.time() - START_TIME, logName, entry) )
		return

	def writeLog (self,
		logName, 	# string; name of log to write
		title, 		# string; title for the output table
		logWriter, 	# AloMarkerLogWriter; to use to build output
		fp		# file pointer; where to write the output
		):
		# Purpose: to build a table of output from the log entries in
		#	the log specified by 'logName', formatted by
		#	'logWriter', with a given 'title', and to write that
		#	table out to 'fp'
		# Returns: nothing
		# Assumes: the contents of the log entries for 'logName' are
		#	appropriate for the writer specified as 'logWriter';
		#	if not, exceptions will likely result
		# Effects: writes several lines to 'fp'
		# Throws: propagates all exceptions if we cannot write to 'fp'
		#	of if the assumptions are violated

		# look up the log corresponding to 'logName' and pass it along
		# to the specified writer

		myLog = None
		if self.logs.has_key(logName):
			myLog = self.logs[logName]
		logWriter.write (fp, title, myLog)
		return
