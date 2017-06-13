# This file is part of nit-picker.
#
# nit-picker is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# nit-picker is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with nit-picker. If not, see <http://www.gnu.org/licenses/>.

# Version 1.0


import sys
import os
import getopt
import nit_picker

SoftwareVersion = "nit-picker Version 1.0"

# TODO Fix usage
def usage():
    print("usage:\n" + 
    "\tThis script is written for python 2.7.11\n" + 
    "\tDo this instruction because it is wrong right now.\n" + 
    "\t\tOR\n" + 
    "\t\tSubdirectories containing .fast5 reads (Barcoded Reads, ex. /BCO1)\n\n" + 
    "\tThe output directory will be filled with .fasta and .fastq reads, sorted by subfolders.\n\n" + 

    "\tOptions:\n" +  
    "\t-i\t--idir   \tInput Directory (required)\n" +  
    "\t-o\t--odir   \tOutput Directory (required)\n" +  
    "\t-m\t--minlen \tMinimum Read Length Filter\n" +  
    "\t-M\t--maxlen \tMaximum Read Length Filter\n" +  
    "\t-h\t--help   \tPrint this message\n" +   
    "\t-r\t--rundate\tSequencing Date or other information which will be included in the extract filename.\n" +   
    
    "\n\tSee README.MD for instructions on how to set up an anaconda environment for this script\n"
    )   
    # Usage: you must provide a read input and output directory. 
    # reads are fastq or a directory with fastq.
    # The rest are optional.

# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    # Default to None.  So I can easily check if they were not passed in.
   
    global readInput
    global outputResultDirectory
    global barcodeFileLocation    
    global minimumReadLength
    global maximumReadLength
    global minimumQuality
    global maximumQuality    
    global barcodeSampleMapFilename
    global sampleID
            
    readInput                = None
    outputResultDirectory    = None
    barcodeFileLocation      = None    
    minimumReadLength        = None
    maximumReadLength        = None
    minimumQuality           = None
    maximumQuality           = None    
    barcodeSampleMapFilename = None
    sampleID                 = None

 

    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt.getopt(sys.argv[1:]
            ,"m:M:q:Q:hvbo:r:s:"
            ,["minlen=", "maxlen=", "minqual=", "maxqual=", "help", "version","barcode=","outputdir=","reads=", "sampleid="])

        for opt, arg in opts:

            if opt in ('-h', '--help'):
                print (SoftwareVersion)
                usage()
                return False

            elif opt in ('-v', '--version'):
                print (SoftwareVersion)
                return False

            elif opt in ("-o", "--outputdir"):
                outputResultDirectory = arg
            elif opt in ("-r", "--reads"):
                readInput = arg
                
            elif opt in ("-b", "--barcode"):
                barcodeFileLocation = arg
                
                
            elif opt in ("-m", "--minlen"):
                minimumReadLength = int(arg)
            elif opt in ("-M", "--maxlen"):
                maximumReadLength = int(arg)
            
            elif opt in ("-q", "--minqual"):
                minimumQuality = int(arg)   
                
            elif opt in ("-q", "--maxqual"):
                maximumQuality = int(arg)    
            
            elif opt in ("-s", "--sampleid"):
                sampleID = arg
                
            else:
                print('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
                raise Exception('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
            
        if(len(sys.argv) < 3):
            print ('I don\'t think you have enough arguments.\n')
            usage()
            return False     

    except getopt.GetoptError, errorMessage:
        print ('Something seems wrong with your commandline parameters.')
        print (errorMessage)
        usage()
        return False

    
    # Sanity Checks    
    if (os.path.isfile(readInput)):
        print ('Read input is a file that exists.')
    elif (os.path.isdir(readInput)):
        print ('Read input is a directory that exists.')
    else :
        print ('I don\'t understand the read input specified, it is not a file or directory:' + readInput)
        return False
    
    # This output directory should exist
    if not os.path.exists(outputResultDirectory):
        os.makedirs(outputResultDirectory)
        
    # TODO: trim out spaces and special characters from the sample id.
    # Check on all the values.


    return True




if __name__=='__main__':
    
    try:    
        if(readArgs()):
            print('Commandline arguments look fine.\n Now I will prepare the reads and calculate quality statistics.')
            
            nit_picker.prepareReads(readInput, outputResultDirectory, sampleID, barcodeFileLocation, minimumReadLength, maximumReadLength, minimumQuality, maximumQuality )
            
            print ('Done with nit-picker for now. Have a nice day.')    
        else:
            print('\nI\'m giving up because I was not satisfied with your commandline arguments.')  
            
    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print 'Fatal problem during read extraction:'
        print sys.exc_info()
        raise
   