# Purpose: Python library to handle logging for the ALO/Marker Association
#	load.  Deals with collecting log items, timestamping them, and
#	providing a flexible mechanism for writing them.
# Author: jsb
# Assumes: nothing

import time		# used for its time.time() function (for timestamps)
import sys

###------------------------------------------------------------------------###

###------------------------###
###--- Global Variables ---###
###------------------------###

DEBUG = False			# set to true for debugging output to stderr
START_TIME = time.time()	# time at which the module is imported

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
                title,	# str. title for the log
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
                title,	# str. title for the log
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
                logName, 	# str. name of log to which to send 'entry'
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

                if logName not in self.logs:
                        self.logs[logName] = AloMarkerLog()
                self.logs[logName].addEntry (entry)
                if DEBUG:
                        sys.stderr.write ('%9.3f : %s : %s\n' % (
                                time.time() - START_TIME, logName, entry) )
                return

        def writeLog (self,
                logName, 	# str. name of log to write
                title, 		# str. title for the output table
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
                if logName in self.logs:
                        myLog = self.logs[logName]
                logWriter.write (fp, title, myLog)
                return
