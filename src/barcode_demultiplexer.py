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
from StringIO import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from os import listdir
from os.path import isfile, join, dirname

inputDirectory=''
outputDirectory=''

# Read Barcodes from an input file
def readBarcodes():
    try:
        #runningDirectory = os.path.dirname(sys.argv[0])
        #barcodeFileNameWithPath = join(runningDirectory, 'barcodes.txt')
        barcodeFileNameWithPath = join(inputDirectory, 'barcodes.txt')
        # print ('running directory:' + runningDirectory)
        print ('barcodeFileName=' + barcodeFileNameWithPath)

        barcodes = {}
        barcodeFile = open(barcodeFileNameWithPath, 'r')
        for i, line in enumerate(barcodeFile):
            if (len(line.strip()) > 0):
                # Filter comments.
                if not (line.strip()[0:1]=='#'):
                    tokens = line.split()
                    barcodes[tokens[0]] = tokens[1]

        barcodeFile.close()
        return barcodes

    except Exception, e:
        print ('Problem reading barcode file: ' + str(e))
        raise
        
def splitByBarcode():
    print('Searching for barcodes in ' + inputDirectory)
    barcodeList = readBarcodes()

    #List of the files in this directory
    fileNames = [f for f in listdir(inputDirectory) if isfile(join(inputDirectory, f))]
    fileNames.sort()

    print('I found ' + str(len(fileNames)) + ' files in the input directory.')
    
    print('THe first file looks like this:' + fileNames[0])

    print('Preparing output files in ' + outputDirectory)
    #twoDFileOutputs = {}
    allReadFileOutputs = {}

    for key in barcodeList.keys():
    #    twoDFileName = join(outputDirectory, key + '_2D_reads.fastq')
    #   twoDFileOutputs[key] = open(twoDFileName, 'w')
        allReadFileName = join(outputDirectory, key + '_all_reads.fastq')
        allReadFileOutputs[key] = open(allReadFileName, 'w')
        
    unbarcodedOutput = open(join(outputDirectory, 'unbarcoded.fastq'), 'w')

    print('Splitting by barcodes, Just a second...')

    # get fastq records from the input file.
    inputReadFileNameWithPath = join(inputDirectory, 'reads.fastq')
    parsedInputReads = SeqIO.parse(inputReadFileNameWithPath, 'fastq') 
    
    #print('I opened the input file. It has this many fastq records:' + str(len(parsedInputReads)))

    print('This is the barcode sequence for BC01:' + barcodeList['BC01'])
    print('2 nested loops. First loop through each fastq record. Then loop thru each barcode.')
    
    
    # TODO: Maybe I can parse the reads in one step, and write in another
    # Writing a sequence to the hard drive for each iteration is quite slow.
    # But if i keep all the data in memory, it also might be slow.
    for rec in parsedInputReads:
        
        print('Searching Barcodes,id=' + rec.id)
        barcodeFound = False
    
        # Loop through each barcode.
        for barcodeKey in barcodeList.keys():
            barcodeText = barcodeList[barcodeKey]
            
            revcomBarcode = Seq(barcodeText).reverse_complement()
            #print ('barcode:' + barcodeText)
            #print ('reverse complement barcode:' + revcomBarcode)
            # Here is where I check if the barcode exists in the reads.  
            if barcodeText in rec:
                print('I found ' + barcodeKey + ' in the record with id:' + rec.id)
                SeqIO.write([rec], allReadFileOutputs[barcodeKey], 'fastq')
                barcodeFound = True
                #if('TwoDir' in fast5Key):
                #    SeqIO.write([rec], twoDFileOutputs[barcodeKey], 'fastq')
                 
            elif revcomBarcode in rec:
                #beforeQualities = rec.letter_annotations["phred_quality"] 
                #beforeSequence = rec.seq
                print('I found ' + barcodeKey + ' (REVCOM) in the record with id:' + rec.id)
                
                barcodeFound = True
                
                #calculatedRevcom = rec.reverse_complement(id=True, name=True, description=True)
                #print('revcom type:' + str(type(calculatedRevcom)))
                
                
                rec = rec.reverse_complement(id=rec.id+"_reverse_complement", name=True, description=True)
                SeqIO.write([rec], allReadFileOutputs[barcodeKey], 'fastq')
                #afterQualities = rec.letter_annotations["phred_quality"] 
                #print ('I reverse complemented it. the new id is:' + str(rec.id))
                #print ('before qualities:' + str(beforeQualities))
                #print ('after qualities:' + str(afterQualities))
            
                    
            else:
                #print('Barcode not found in the record with id' + rec.id)
                #Don't write the file here. You're too deep in the loop.
                #SeqIO.write([rec], unbarcodedOutput, 'fastq')
                pass
   
        if(barcodeFound):
            # Great.  Do nothing.
            pass
        else:
            SeqIO.write([rec], unbarcodedOutput, 'fastq')
            
    #close output files.
    for key in barcodeList.keys():
        #twoDFileOutputs[key].close()
        allReadFileOutputs[key].close()
        
    unbarcodedOutput

    print('Done.')

if __name__=='__main__':

    if(len(sys.argv)==2) and (
        sys.argv[1].lower() == '-v' or 
        sys.argv[1].lower() == '--version' or 
        sys.argv[1].lower() == '-version'
    ):
        print (SoftwareVersion)

    elif(len(sys.argv)==3): 
        
        inputDirectory = sys.argv[1]
        outputDirectory = sys.argv[2]

        # This output directory should exist
        if not os.path.exists(outputDirectory):
            os.makedirs(outputDirectory)

        splitByBarcode()

        print ('Done!')

    else:
        print("usage:\n" + 
            "\tRun this program using a standard python call, with two parameters specifying input and output:\n" + 
            "\t$python search_barcode.py inputDirectory outputDirectory\n" + 
            "\tSee README.MD for instructions on how to set up an anaconda environment for this script\n"
        )        

