def getSeq(dir='',file='',seqName='',accurateMatching = True):
    fileIn = open(dir+file,'r')
    text = fileIn.readlines()
    fileIn.close()
    seq = ""
    seqNameInDB=""
    flag = 0
    isFound = 0
    seqName = seqName.strip("\n")
    if len(seqName)<1:
        return seqName,"SeqNameError"
    for eachLine in text:
        eachLine=eachLine.rstrip("\n")
        if accurateMatching==True:
            if len(eachLine)>0 and eachLine[0]==">" and eachLine==seqName:
                flag = 1
                seqNameInDB = eachLine
                # print(eachLine)
                isFound = 1
                continue
            if flag == 1:
                if eachLine[0]==">":
                    break
                else:
                    seq = seq + eachLine
        else:
            if len(eachLine)>0 and  eachLine[0]==">" and eachLine.find(seqName)!=-1:
                flag = 1
                seqNameInDB = eachLine
                # print(eachLine)
                isFound = 2
                continue
            if flag == 1:
                if eachLine[0]==">":
                    break
                else:
                    seq = seq + eachLine
    if isFound==0:
        return seqName,"\tcanNotFound"
    else:
        return seqNameInDB,seq

# a=getSeq("H2.fa","/home/think/18Mid/GisAidSeq/","111\n",True)
# print(a)