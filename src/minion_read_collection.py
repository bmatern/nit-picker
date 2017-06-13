from Bio import SeqIO
import numpy
import pylab
from os.path import split, isdir
from os import makedirs

# A class defining a collection of minION reads
class minionReadCollection:

    def __init__(self, readArray):
        self.readCollection=readArray
        self.readInputFormat='fastq'
        self.readLengths = []
        self.readAvgPhredQualities = []
        
        self.calculateReadStats()

    def calculateReadStats(self):
        #print('calculating read stats')
        
        self.readLengths = []
        self.readAvgPhredQualities = []
        
        for currentRead in self.readCollection:
            phredQualities = currentRead.letter_annotations["phred_quality"]                    
            currentSeqLength = len(currentRead)
            currentAvgPhredQuality = numpy.mean(phredQualities)
            
            self.readLengths.append(currentSeqLength)
            self.readAvgPhredQualities.append(currentAvgPhredQuality)   
        
    def concatenate(self, otherCollection): 
        self.readCollection = self.readCollection + otherCollection.readCollection
        self.readLengths = self.readLengths + otherCollection.readLengths
        self.readAvgPhredQualities = self.readAvgPhredQualities + otherCollection.readAvgPhredQualities
        #self.calculateReadStats()
        
        
def createCollectionFromReadFile(readFile):
    global readInputFormat
    #print ('loading reads from:' + readFile)
    
    # Determine Input File Type
    if (".fasta" == readFile[-6:] or ".fa" == readFile[-3:]):
        readInputFormat = "fasta"
        raise Exception('Fasta files are not supported.  You can try to add this functionality if you want.')
    elif (".fastq"== readFile[-6:] or ".fq" == readFile[-3:]):
        readInputFormat = "fastq"
    else:
        print('I expect a .fasta or .fastq format for read input. Alternatively, specify a directory containing read inputs. Please check your input.')
        raise Exception('Bad Read Input Format')
    
    newCollectionObject = minionReadCollection(list(SeqIO.parse(readFile, readInputFormat)))
    newCollectionObject.readInputFormat = readInputFormat
    
    return newCollectionObject    


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
    
# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput      
