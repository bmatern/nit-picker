#import sys
import os
import numpy
from Bio import SeqIO
import minion_read_collection
from minion_read_collection import *#minionReadCollection
import barcode_demultiplexer
from barcode_demultiplexer import *
from os.path import split, join, isdir, isfile
from os import makedirs
#from Bio.Seq import Seq



#nit_picker.prepareReads(readInput, outputResultDirectory, sampleID, barcodeFileLocation, minimumReadLength, maximumReadLength, minimumQuality, maximumQuality )
#inputReads and outputDirectory are necessary.  The rest can be "None" if you would like.            
def prepareReads(inputReads, outputDirectory, sampleID, barcodeFileLocation, minimumReadLength, maximumReadLength, minimumReadQuality, maximumReadQuality ):
    print ('Preparing Reads')
    
    #TODO: Create an output file with read stats, like the read extractor did.
    
    allReads = minionReadCollection([])
    
    # Default sample id
    if (sampleID is None):
        sampleID = 'minion_reads'

    # Determine if input is a file, or a directory
    if (isfile(inputReads)):
        print ('Read input is a file that exists.')
        allReads = createCollectionFromReadFile(inputReads)
    elif (isdir(inputReads)):
        print ('Read input is a directory that exists.')
        for currentInputReadFile in os.listdir(inputReads):
            if (".fasta" == currentInputReadFile[-6:] or ".fa" == currentInputReadFile[-3:] or ".fastq"== currentInputReadFile[-6:] or ".fq" == currentInputReadFile[-3:]):
                print ('loading Reads from:' + str(join(inputReads,currentInputReadFile)))
                newReads = createCollectionFromReadFile(join(inputReads,currentInputReadFile))
                allReads.concatenate(newReads)
    else :
        print('I expect a .fasta or .fastq format for read input. Alternatively, specify a directory containing read inputs. Please check your input.')
        raise Exception('Bad Read Input Format')
    
    print ('Total # of allReads found = ' + str(len(allReads.readCollection)))

    passReads=[]
    lengthRejectReads=[]
    qualityRejectReads=[]

    # Iterate all reads
    print ('Rejecting reads for Length and Quality...')
    for currentRead in allReads.readCollection:
        
        phredQualities = currentRead.letter_annotations["phred_quality"]                    
        currentSeqLength = len(currentRead)
        currentAvgPhredQuality = numpy.mean(phredQualities)
        
        # Reject allReads that are wrong length
        if( (minimumReadLength is not None and currentSeqLength < minimumReadLength)
            or (maximumReadLength is not None and currentSeqLength > maximumReadLength)
            ):
            lengthRejectReads.append(currentRead)

        # Reject allReads that have the wrong quality
        elif( (minimumReadQuality is not None and currentAvgPhredQuality < minimumReadQuality)
            or (maximumReadQuality is not None and currentAvgPhredQuality > maximumReadQuality)
            ):
            qualityRejectReads.append(currentRead)
      
        # This read is okay.
        else:
            passReads.append(currentRead)
           
    passReadCollection=minionReadCollection(passReads)
    lengthRejectReadCollection=minionReadCollection(lengthRejectReads)
    qualityRejectReadCollection=minionReadCollection(qualityRejectReads)

    # Create Scatterplots

    if(len(allReads.readCollection)>0):
        
        allReads.calculateReadStats()        
        minion_read_collection.createScatterPlot("All Reads"
            , allReads.readLengths
            , allReads.readAvgPhredQualities 
            , "Read Lengths"
            , "Avg Read Quality(Phred)"
            , join(outputDirectory,  str(sampleID) + '_Unfiltered_Reads'))
    if(len(lengthRejectReadCollection.readCollection)>0):
        minion_read_collection.createScatterPlot("Length-Rejected Reads"
            , lengthRejectReadCollection.readLengths
            , lengthRejectReadCollection.readAvgPhredQualities 
            , "Read Lengths"
            , "Avg Read Quality(Phred)"
            , join(outputDirectory, str(sampleID) + '_Length_Rejected_Reads'))
    if(len(qualityRejectReadCollection.readCollection)>0):
        minion_read_collection.createScatterPlot("Quality-Rejected Reads"
            , qualityRejectReadCollection.readLengths
            , qualityRejectReadCollection.readAvgPhredQualities 
            , "Read Lengths"
            , "Avg Read Quality(Phred)"
            , join(outputDirectory, str(sampleID) + '_Quality_Rejected_Reads'))
            
    # Print fastq for Rejected Reads and all reads    
    readFormat = allReads.readInputFormat
    
    allReadOutputFile = createOutputFile(join(outputDirectory, str(sampleID) + '_Unfiltered_Reads.' + readFormat))
    SeqIO.write(allReads.readCollection, allReadOutputFile, readFormat)
    allReadOutputFile.close()
    
    lengthRejectReadOutputFile = createOutputFile(join(outputDirectory, str(sampleID) + '_Length_Rejected_Reads.' + readFormat))
    SeqIO.write(lengthRejectReadCollection.readCollection, lengthRejectReadOutputFile, readFormat)
    lengthRejectReadOutputFile.close()
    
    qualityRejectReadOutputFile = createOutputFile(join(outputDirectory, str(sampleID) + '_Quality_Rejected_Reads.' + readFormat))
    SeqIO.write(qualityRejectReadCollection.readCollection, qualityRejectReadOutputFile, readFormat)
    qualityRejectReadOutputFile.close()
    

    if(barcodeFileLocation is None):
        print('No barcode file was provided. I will not attempt to de-multiplex the allReads.')
        
        
        if(len(passReadCollection.readCollection)>0):
            minion_read_collection.createScatterPlot("Pass Reads"
                , passReadCollection.readLengths
                , passReadCollection.readAvgPhredQualities 
                , "Read Lengths"
                , "Avg Read Quality(Phred)"
                , join(outputDirectory,  str(sampleID) + '_PassReads'))
            
        passReadOutputFile = createOutputFile(join(outputDirectory, str(sampleID) + '_Pass.' + readFormat))
        SeqIO.write(passReadCollection.readCollection, passReadOutputFile, readFormat)
        passReadOutputFile.close()
        
        #qualityReadOutput
        # Print Pass

    else:
        print('Provided Barcode File:' + str(barcodeFileLocation) + '\nI will de-multiplex the allReads.')
        #print ('Just kidding, you still need to add this code.  Do it now.')
        # Debarcode pass reads
        
        #def splitByBarcode(readCollection, outputDirectory, barcodeFileNameWithPath):
        splitByBarcode(passReadCollection, outputDirectory, barcodeFileLocation, sampleID)
    
    



         
# TODO: Am i using this?
# Short method to print out summary stats for a group of allReads
def writeReadStatsSub(readType, barcode, readStats):
    
    statsSubText=''        
    
    if(len(readStats) > 0):
        #readLengths = readStats[:,0]
        #readQualities = readStats[:,1]
        
        statsSubText += ('\n' + barcode + ' ' + readType + ' Reads Summary:\n' 
            + 'Sequences Extracted,' + str(len(readStats)) + '\n'
            + 'Minimum Length,' + str(int(numpy.amin(readStats[:,0]))) + '\n'
            + 'Maximum Length,' + str(int(numpy.amax(readStats[:,0]))) + '\n'
            + 'Mean Length,' + str(numpy.mean(readStats[:,0])) + '\n'
            + 'Mean Quality,' + str(numpy.mean(readStats[:,1])) + '\n'
        )
        
    else:
        statsSubText += ('\n' + barcode + ' ' + readType + ' Reads Summary:\n' 
            + 'Sequences Extracted,0\n'

        )

    return statsSubText
