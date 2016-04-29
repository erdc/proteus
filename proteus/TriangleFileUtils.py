#!/usr/bin/env python
#import standard Python modules
"""
A collection of functions for reading and manipulating files so that I
can try to read and write files in the format that Triangle expects

.. inheritance-diagram:: proteus.TriangleFileUtils
   :parts: 1
"""
import sys,os
import numpy


#setup basic script for doing simulations
def checkFileExists(filename):
    if not os.path.isfile(filename):
        print '%s is not a file!' % filename
        #sys.exit(1)
        return False
    #end
    else:
        return True
#end checkFileExists

def checkDataRead(dataArray,notReadVal):
    """
    just go through data variable and check to make sure
    no values are equal to a special value to which data
    should have been initialized before it was read

    returns boolean for ok or not
    """
    ok = True
    for d in dataArray.flat:
        ok = (float(math.abs(d-notReadVal)) > 1.0e-6)
    # end d
    return ok
#end checkDataRead

def generateDefaultFormatPattern(recInfo,verbose=0):
    """
    generate the initial line format for the Triangle data files
    based on the variable types specified in recInfo

    This is in most cases a line that looks like
       Nentries NumEntriesInRec0 NumEntriesInRec1 ....

    In some cases (e.g., segments) the number of entries in a segment
    is omitted so that the default format doesn't follow

    """
    #number of records looking for
    nrecs = recInfo['nRecords']
    if verbose > 2:
        print 'looking for pattern with up to ',nrecs,' records'
    #end if
    if verbose > 2:
        for i in range(nrecs):
            print 'record type[',i,']= ',recInfo['types'][i]
    #end verbose
    #format line is a sequence of unsigned integers
    unint = r'(\d+)'
    formatpattern  = r'^\s*'+unint
    for i in range(nrecs):
        formatpattern+=r'\s+'+unint
    #end i
    #try to force match for whole string, so that don't just match
    #the beginning?
    formatpattern += r'\s*$' #don't need this to get match at beg. of line
    if verbose > 2:
        print 'trying to match pattern ',formatpattern
    #end if

    return formatpattern
#end def

def findInitialFormatLine(lines,recInfo,formatpattern,
                          commchar='#',nbase=0,verbose=0):
    """
    find the line for the Triangle data files, specified in
    formatpattern.

    Generate dataInfo dictionary that should specify

      dataInfo['nEntries']    : number of entries in the file
      dataInfo['recordSizes'] : size of each record per entry (can be 0)
      dataInfo['dataStartLineNumber']    : line where actual data entries may
                                           start
      dataInfo['recordLocationsPerLine'] : array of indeces saying
                                           where each record
                                           should start on a valid entry line
                                           record I is in entries
                                           recStart[I]:recStart[I+1]

    """

    import re
    #mwf debug
    #verbose=6
    #number of records looking for
    nrecs = recInfo['nRecords']
    if verbose > 2:
        print 'looking for file with up to ',nrecs,' records'
    #end if
    if verbose > 2:
        for i in range(nrecs):
            print 'record type[',i,']= ',recInfo['types'][i]
    #end verbose
    if verbose > 2:
        print 'trying to match pattern ',formatpattern
    #end if
    #record information about what is read
    dataInfo = {}
    dataInfo['nEntries']=0        #number of entries read
    dataInfo['recordSizes'] = []  #fields for each record per entry
    dataInfo['recordSizes'] = [0 for i in range(nrecs)]
    if verbose > 9:
        print 'recordSizes= ',dataInfo['recordSizes']
    #end verbose
    done = False
    iline= 0
    while not done:
        #read next line
        line = lines[iline]
        if verbose > 5:
            print 'read line= ',line
        #throw away anything after # character
        dline= line.split(commchar)[0]
        if verbose > 5:
            print 'leading portion= ',dline
        match =  re.search(formatpattern,dline)
        if match:
            foundFirstLine = True
            done = True
            if verbose > 5:
                print 'first line match.groups()= ',match.groups()
            #end if
            dataInfo['nEntries']    = int(match.group(1))
            dataInfo['recordSizes'] = [int(match.group(i)) for i in range(2,nrecs+2)]
            iline = iline+1 #move on past this line
        else:
            iline = iline+1
            if iline == len(lines):
                done = True
                foundFirstLine = False
            #end reached the end
        #end else found match
    #end while
    #now proceed through the rest of the file and get the data
    nentries= dataInfo['nEntries']
    nValuesPerLine = sum(dataInfo['recordSizes'])+1 #include entry number
    #should have record I contained in entries recStart[I]:recStart[I+1]
    #also include first value as entry number
    recStart = numpy.zeros((nrecs+2,),'i')
    recStart[0] = 0
    recStart[1] = 1 #first value is entry number
    for i in range(nrecs):
        recStart[2+i] = recStart[i+1]+dataInfo['recordSizes'][i]
    #end i
    if verbose > 2:
        print 'recStart= ',recStart
    #end verbose
    dataInfo['recordLocationsPerLine'] = recStart
    dataInfo['dataStartLineNumber']    = int(iline)
    if verbose > 2:
        print 'end of findInitialFormatLine, foundFirstLine= ',foundFirstLine
        print 'next line is ',iline
        print """nEntries= %d nValuesPerLine= %d """ % (dataInfo['nEntries'],
                                                        nValuesPerLine)
    #end verbose
    return dataInfo

#end def
def readSimpleFormattedDataEntries(lines,recInfo,dataInfo,
                                   commchar='#',nbase=0,verbose=0):
    """
    Try to read some data consisting of records
    stored in lines. The records are described in recInfo.
    The basic thing I need is the possible number of records
    per line and the type for each one.

    I'll go ahead and assume the file looks like a series of blank
    lines or commented lines and then the first meaningful line which
    gives the file format that should be the same as what I'm
    expecting ...

    Nentries NumEntriesInRec0 NumEntriesInRec1 ....

    The regex format for this line is given in formatpattern.

    After that comes Nentry values with the basic format

    Entry Num. Rec0 Rec1 ...

    If NumEntriesInRecI is 0 then that record shouldn't appear

    returns
       data     : what was read
       ilast    : index of last line read

    """
    #
    if verbose > 2:
        print "rec info in read simple....",recInfo
    nrecs = recInfo['nRecords']
    #proceed through the rest of the file and get the data
    nentries= dataInfo['nEntries']
    nValuesPerLine = sum(dataInfo['recordSizes'])+1 #include entry number
    #should have record I contained in entries recStart[I]:recStart[I+1]
    #also include first value as entry number
    recStart       = dataInfo['recordLocationsPerLine']
    iline          = 0#assume start at beginning of lines
    if nentries == 0:
        ilast = iline
        return None,ilast
    #end no entries
    #describes where
    #setup default data read as something obvious if it's not
    #initialized
    data = []
    for i in range(nrecs):
        if dataInfo['recordSizes'][i] > 1:
            data.append(numpy.ones((nentries,
                                    dataInfo['recordSizes'][i]),
                                   dtype=recInfo['types'][i]))
            data[-1][:]=recInfo['defaultValues'][i]
        elif dataInfo['recordSizes'][i] == 1:
            data.append(numpy.ones((nentries,),
                                   dtype=recInfo['types'][i]))
            data[-1][:]=recInfo['defaultValues'][i]
        else:
            data.append(None)
        #end else
    #end i

    #now loop through the following lines and read in data
    done = (iline == len(lines))
    ientry = 0 #record number of entries read
    while not done:
        #read next line
        line = lines[iline]
        if verbose > 5:
            print 'read line= ',line
        #throw away anything after # character
        dline= line.split(commchar)[0]
        if verbose > 5:
            print 'leading portion= ',dline
        #end if
        #skip empty or all white space
        if (len(dline) > 0 and not dline.isspace()):
            entry = dline.split() #separate by white space
            if not (len(entry) == nValuesPerLine):
                print 'PROBLEM readFormatted info non-trivial line wrong length, quitting!'
                print 'len(entry)= ',len(entry),' problem nValuesPerLine= ',nValuesPerLine
                print 'entry = ',entry
                raise IOError, "file format in TriangleFileUtils.readSimpleFormattedDataEntries"
                return data,iline
            #end if
            entryNum = entry[recStart[0]:recStart[1]]
            if verbose > 5:
                print 'entryNum= ',entryNum
            index    = int(entryNum[0])-nbase #could be base 1
            if verbose > 5:
                print 'index= ',index

            for i in range(nrecs): #go through each record and insert if it is in file
                #get the substring for this record
                subentry = entry[recStart[i+1]:recStart[i+2]]
                if verbose > 5:
                    print """entry num %d rec %d subentry= %s """ % (ientry,i,subentry)
                    print """conversions[i]= %s """ % recInfo['conversions'][i]
                if dataInfo['recordSizes'][i] > 1:
                    for j,s in enumerate(subentry):
                        data[i][index,j] = recInfo['conversions'][i](s)
                    #end loop
                elif dataInfo['recordSizes'][i] == 1:
                    data[i][index] = recInfo['conversions'][i](subentry[0])
                #end 1d array
            #end i
            ientry = ientry+1 #added another entry
        #end nonempty line
        iline = iline+1
        done = (ientry == nentries) or (iline == len(lines))
    #end not done
    if verbose > 3:
        print """leaving readInfo iline= %d ientry= %d """ % (iline,ientry)
    #end
    ilast = iline
    return data,ilast
#end readFormatted data
def readSimpleFormattedDataEntriesLastOptional(lines,recInfo,dataInfo,
                                               commchar='#',nbase=0,verbose=0):
    """
    Try to read some data consisting of records
    stored in lines. The records are described in recInfo.
    The basic thing I need is the possible number of records
    per line and the type for each one.

    This version allows the last data entry to be omitted so that I
      can read the regions list. I need to generalize to allow
      a given entry to not be found I guess

    I'll go ahead and assume the file looks like a series of blank
    lines or commented lines and then the first meaningful line which
    gives the file format that should be the same as what I'm
    expecting ...

    Nentries NumEntriesInRec0 NumEntriesInRec1 ....

    The regex format for this line is given in formatpattern.

    After that comes Nentry values with the basic format

    Entry Num. Rec0 Rec1 ...

    If NumEntriesInRecI is 0 then that record shouldn't appear

    returns
       data     : what was read
       ilast    : index of last line read

    """
    #
    nrecs = recInfo['nRecords']
    #proceed through the rest of the file and get the data
    nentries= dataInfo['nEntries']
    nValuesPerLine = sum(dataInfo['recordSizes'])+1 #include entry number
    #should have record I contained in entries recStart[I]:recStart[I+1]
    #also include first value as entry number
    recStart       = dataInfo['recordLocationsPerLine']
    iline          = 0#assume start at beginning of lines
    if nentries == 0:
        ilast = iline
        return None,ilast
    #end no entries
    #describes where
    #setup default data read as something obvious if it's not
    #initialized
    data = []
    for i in range(nrecs):
        if dataInfo['recordSizes'][i] > 1:
            data.append(numpy.ones((nentries,
                                    dataInfo['recordSizes'][i]),
                                   dtype=recInfo['types'][i]))
            data[-1][:]=recInfo['defaultValues'][i]
        elif dataInfo['recordSizes'][i] == 1:
            data.append(numpy.ones((nentries,),
                                   dtype=recInfo['types'][i]))
            data[-1][:]=recInfo['defaultValues'][i]
        else:
            data.append(None)
        #end else
    #end i
    #now loop through the following lines and read in data
    done = (iline == len(lines))
    ientry = 0 #record number of entries read
    while not done:
        #read next line
        line = lines[iline]
        if verbose > 5:
            print 'read line= ',line
        #throw away anything after # character
        dline= line.split(commchar)[0]
        if verbose > 5:
            print 'leading portion= ',dline
        #end if
        #skip empty or all white space
        if (len(dline) > 0 and not dline.isspace()):
            entry = dline.split() #separate by white space
            if not (len(entry) == nValuesPerLine or
                    len(entry) == nValuesPerLine-1):
                print 'PROBLEM readFormatted info non-trivial line wrong length, quitting!'
                print 'len(entry)= ',len(entry),' problem nValuesPerLine= ',nValuesPerLine
                print 'entry = ',entry
                raise IOError, "file format in TriangleFileUtils.readSimpleFormattedDataEntries"
                return data,iline
            #end if
            entryNum = entry[recStart[0]:recStart[1]]
            if verbose > 5:
                print 'entryNum= ',entryNum
            index    = int(entryNum[0])-nbase #could be base 1
            if verbose > 5:
                print 'index= ',index

            for i in range(nrecs): #go through each record and insert if it is in file
                #get the substring for this record
                padEnd = False
                lastInd= recStart[i+2]
                if lastInd > len(entry) and dataInfo['recordSizes'][i] > 1:
                    padEnd = True
                    lastInd= len(entry)
                #end if
                subentry = entry[recStart[i+1]:lastInd]
                if verbose > 5:
                    print """entry num %d rec %d subentry= %s """ % (ientry,i,subentry)
                    #print """conversions[i]= %s """ % recInfo['conversions'][i]
                if dataInfo['recordSizes'][i] > 1:
                    for j,s in enumerate(subentry):
                        data[i][index,j] = recInfo['conversions'][i](s)
                    #end loop
                    if padEnd:
                        data[i][index,len(subentry):]=data[i][index,len(subentry)-1]
                    #end padEnd
                elif dataInfo['recordSizes'][i] == 1:
                    data[i][index] = recInfo['conversions'][i](subentry[0])
                #end 1d array
            #end i
            ientry = ientry+1 #added another entry
        #end nonempty line
        iline = iline+1
        done = (ientry == nentries) or (iline == len(lines))
    #end not done
    if verbose > 3:
        print """leaving readInfo iline= %d ientry= %d """ % (iline,ientry)
    #end
    ilast = iline
    return data,ilast
#end readFormatted data
