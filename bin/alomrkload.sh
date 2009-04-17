#!/bin/sh 
#
#  alomrkload.sh
###########################################################################
#
#  Purpose:  This script determines associations between ALOs 
#	(Allele-Like Objects) and markers, and related data.
#
#  Usage:
#
   Usage=alomrkload.sh
#
#  Env Vars:
#
#      See the configuration file
#
#  Inputs:
#
#      - Common configuration file (common.config.sh)
#      - Load configuration file (alomrkload.config)
#      - MGD database
#
#  Outputs:
#
#      - An archive file
#      - Log files defined by the environment variables ${LOG_PROC},
#        ${LOG_DIAG}, ${LOG_CUR} and ${LOG_VAL}
#      - BCP files for each database table to be loaded
#      - Records written to the database tables
#      - inline updates to the database
#      - Exceptions written to standard error
#      - Configuration and initialization errors are written to a log file
#        for the shell script
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Fatal error occurred
#      2:  Non-fatal error occurred
#
#  Assumes:  Nothing
#
#  Implementation:  Description
#
#  Notes:  None
#
###########################################################################

#
#  Set up a log file for the shell script in case there is an error
#  during configuration and initialization.
#
cd `dirname $0`/..
LOG=`pwd`/alomrkload.log
rm -f ${LOG}

#
#  Verify the argument(s) to the shell script.
#
if [ $# -ne 0 ]
then
    echo "Usage: $0" | tee -a ${LOG}
    exit 1
fi

#
#  Establish the configuration file names.
#
LOAD_CONFIG=`pwd`/alomrkload.config

#
#  Make sure the configuration files are readable.
#
if [ ! -r ${LOAD_CONFIG} ]
then
    echo "Cannot read configuration file: ${LOAD_CONFIG}" | tee -a ${LOG}
    exit 1
fi

#
#  Source the load configuration file.
#

. ${LOAD_CONFIG}
echo "curator log ${LOG_CUR}"

#
#  Source the common DLA functions script.
#
if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}"
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined."
    exit 1
fi

##################################################################
##################################################################
#
# main
#
##################################################################
##################################################################
#
#  Perform pre-load tasks (archive, start logs, print config environ
#                          get and set jobstream key)
#
preload

#
# optionally rm all files and subdirs dirs of directories on the command line
# Note: archiving does not remove them
#
cleanDir ${OUTPUTDIR}

#
#  Run the alo marker association load
#

echo "\n`date`" >> ${LOG_PROC}
echo "Running alo marker association  load" >> ${LOG_PROC}
${ALOMRKLOAD}/bin/aloMarkerAssoc.py ${DBDEBUG} >> ${LOG_DIAG} 2>&1
STAT=$?
if [ ${STAT} -ne 0 ]
then
    echo "ALO Marker Association Load failed.  Return status: ${STAT}" >> ${LOG_PROC}
    shutDown
    exit 1
fi
echo "ALO Marker Association Load completed successfully" >> ${LOG_PROC}

shutDown

exit 0

