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
from Bio import SeqIO
from Bio.Seq import Seq
from os import listdir
from os.path import isfile, join
import minion_read_collection
from minion_read_collection import *#minionReadCollection


# Read Barcodes from an input file
def readBarcodes(barcodeFileNameWithPath):
    try:
        print ('Reading Barcodes:' + barcodeFileNameWithPath)

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
        
def splitByBarcode(currentReadCollection, outputDirectory, barcodeFileNameWithPath, sampleID):
    #print('Searching for barcodes in ' + inputDirectory)
    barcodeList = readBarcodes(barcodeFileNameWithPath)


    print('Preparing output files in ' + outputDirectory)

    allReadFileOutputs = {}    
    allBarcodeCollections = {}

    for key in barcodeList.keys():
        #allReadFileName = join(outputDirectory, key + '_pass_reads.fastq')
        #allReadFileOutputs[key] = open(allReadFileName, 'w')
        allBarcodeCollections[key] = minionReadCollection([])
        
    #unbarcodedOutput = open(join(outputDirectory, 'unbarcoded_pass_reads.fastq'), 'w')
    unbarcodedReadCollection = minionReadCollection([])

    print('Splitting by barcodes, Just a second...')


    #print('This is the barcode sequence for BC01:' + barcodeList['BC01'])
    #print('2 nested loops. First loop through each fastq record. Then loop thru each barcode.')
    
    
    # TODO: Maybe I can parse the reads in one step, and write in another
    # Writing a sequence to the hard drive for each iteration is quite slow.
    # But if i keep all the data in memory, it also might be slow.
    for rec in currentReadCollection.readCollection:
        
        #print('Searching Barcodes,id=' + rec.id)
        barcodeFound = False
    
        # Loop through each barcode.
        for barcodeKey in barcodeList.keys():
            barcodeText = barcodeList[barcodeKey]
            
            revcomBarcode = Seq(barcodeText).reverse_complement()
            #print ('barcode:' + barcodeText)
            #print ('reverse complement barcode:' + revcomBarcode)
            # Here is where I check if the barcode exists in the reads.  
            if barcodeText in rec:
                #print('I found ' + barcodeKey + ' in the record with id:' + rec.id)
                allBarcodeCollections[barcodeKey].readCollection.append(rec)
                barcodeFound = True
                 
            elif revcomBarcode in rec:
                #beforeQualities = rec.letter_annotations["phred_quality"] 
                #beforeSequence = rec.seq
                #print('I found ' + barcodeKey + ' (REVCOM) in the record with id:' + rec.id)

                barcodeFound = True
                allBarcodeCollections[barcodeKey].readCollection.append(rec)
                # Ok I can revcom the sequence, but I won't right now.
                # Barcodes can be on the forward or reverse sequence
                # It is not worth doing the calculation.
                # The punkin-chunker utility will face reads forward relative to 5'->3' HLA
                #rec = rec.reverse_complement(id=rec.id+"_reverse_complement", name=True, description=True)
                #SeqIO.write([rec], allReadFileOutputs[barcodeKey], 'fastq')
                
            else:
                #Don't write the file here. You're too deep in the loop.
                #SeqIO.write([rec], unbarcodedOutput, 'fastq')
                
                pass
   
        if(barcodeFound):
            # Great.  Do nothing.
            pass
        else:
            unbarcodedReadCollection.readCollection.append(rec)
            #SeqIO.write([rec], unbarcodedOutput, 'fastq')
            
    #close output files. Delete empties and generate scatterplots
    for barcodeKey in barcodeList.keys():
        #print('How many reads for barcode:' + str(barcodeKey) + ':' + str(len(allBarcodeCollections[barcodeKey].readCollection)))
        
        # if barcode file has entries
        if (len(allBarcodeCollections[barcodeKey].readCollection) > 0 ):
            # calculate stats
            allBarcodeCollections[barcodeKey].calculateReadStats()
            
            # write fastq
            barcodeWriter = open(join(outputDirectory,  str(sampleID) + '_' + barcodeKey + '_Pass.fastq'), 'w')
            SeqIO.write(allBarcodeCollections[barcodeKey].readCollection, barcodeWriter, 'fastq')
            barcodeWriter.close()
                        
            # write scatterplot
            minion_read_collection.createScatterPlot("Barcode " + str(barcodeKey) + " Pass Reads"
                , allBarcodeCollections[barcodeKey].readLengths
                , allBarcodeCollections[barcodeKey].readAvgPhredQualities 
                , "Read Lengths"
                , "Avg Read Quality(Phred)"
                , join(outputDirectory,  str(sampleID) + '_' + barcodeKey + '_Pass'))
        else:
            #print('Barcode:' + str(barcodeKey) + ' has no reads. Nothing to do.')
            pass
        
    # write the unbarcoded reads
    unbarcodedReadCollection.calculateReadStats()
    
    unbarcodedWriter = open(join(outputDirectory,  str(sampleID) + '_Unbarcoded_Pass.fastq'), 'w')
    SeqIO.write(unbarcodedReadCollection.readCollection, unbarcodedWriter, 'fastq')
    unbarcodedWriter.close()
    
    minion_read_collection.createScatterPlot("Unbarcoded Pass Reads"
        , unbarcodedReadCollection.readLengths
        , unbarcodedReadCollection.readAvgPhredQualities 
        , "Read Lengths"
        , "Avg Read Quality(Phred)"
        , join(outputDirectory,  str(sampleID) + '_Unbarcoded_Pass'))
        
    #unbarcodedOutput.close()

    print('Done.')



