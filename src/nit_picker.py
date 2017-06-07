import sys
import os
import numpy
from Bio import SeqIO
import pylab
#from Bio.Seq import Seq




def prepareReads(inputReads, outputDirectory, barcoding):
    print ('Preparing Reads')
    
    reads = []

    # Determine if input is a file, or a directory
    if (os.path.isfile(inputReads)):
        print ('Read input is a file that exists.')
        reads = loadReadFile(inputReads)
    elif (os.path.isdir(inputReads)):
        print ('Read input is a directory that exists.')
        for currentInputReadFile in os.listdir(inputReads):
            if (".fasta" == currentInputReadFile[-6:] or ".fa" == currentInputReadFile[-3:] or ".fastq"== currentInputReadFile[-6:] or ".fq" == currentInputReadFile[-3:]):
                print ('loading reads from:' + str(os.path.join(inputReads,currentInputReadFile)))
                reads = reads + loadReadFile(os.path.join(inputReads,currentInputReadFile))
    else :
        print('I expect a .fasta or .fastq format for read input. Alternatively, specify a directory containing read inputs. Please check your input.')
        raise Exception('Bad Read Input Format')
    
    print ('Total # of reads found = ' + str(len(reads)))
    
    readLengths = []
    readAvgQualities = []
    
    for currentRead in reads:
        #currentID = currentRead.id
        phredQualities = currentRead.letter_annotations["phred_quality"]                    
        currentSeqLength = len(currentRead)
        currentAvgQuality = numpy.mean(phredQualities)
        
        readLengths.append(currentSeqLength)
        readAvgQualities.append(currentAvgQuality)
        
        #print('Analyzing read:' + str(currentID))
    
    createScatterPlot("All Reads", readLengths, readAvgQualities, "Read Lengths", "Avg Read Quality(Phred)", os.path.join(outputDirectory, 'ReadDistribution'))
    
    
    
    
def loadReadFile(readFile):
    global readInputFormat
    #print ('loading reads from:' + readFile)
    
    # Determine Input File Type
    if (".fasta" == readFile[-6:] or ".fa" == readFile[-3:]):
        readInputFormat = "fasta"
    elif (".fastq"== readFile[-6:] or ".fq" == readFile[-3:]):
        readInputFormat = "fastq"
    else:
        print('I expect a .fasta or .fastq format for read input. Alternatively, specify a directory containing read inputs. Please check your input.')
        raise Exception('Bad Read Input Format')
        
    return list(SeqIO.parse(readFile, readInputFormat))    
        
def createScatterPlot(graphTitleText, xValues, yValues, xAxisName, yAxisName, outputFileName):
    print('Creating a Scatter Plot: ' + outputFileName)
      
    #Clear the figure and start anew.
    pylab.clf()
    
    # K is black, we need to repeat N times.
    colors = ['K'] * len(xValues)

    pylab.scatter(xValues, yValues, s=1, c=colors, marker='.')
    
    pylab.xlabel(xAxisName)
    pylab.ylabel(yAxisName)
    pylab.title(graphTitleText)
    
    pylab.savefig(outputFileName)
         

# Short method to print out summary stats for a group of reads
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
