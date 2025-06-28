

def getAlignedSite(file = "", dir = ""):
    print("getAlignedSite ...")
    fileIn = open(dir + file, "r")
    text = fileIn.readlines()

    flag = 0
    name = []
    alignedSeqBeforeFilledStarKey = []  
    eachSeqStartAlignedPos = []  
    for each in text:
        if len(each) > 2:
            flag = 1
        if flag == 1 and len(each) > 2:
            name.append(each.split()[0])
            alignedSeqBeforeFilledStarKey.append(each.split("\t")[2])
            eachSeqStartAlignedPos.append(each.split()[1])
    alignedSeq = []
    flag1 = 0
    
    for each in alignedSeqBeforeFilledStarKey:
        if flag1 == 0:
            flag1 = 1
            querySeqLength = len(each)
        each = str(each).rstrip("\n") + "*" * (querySeqLength - len(each)) + "\n"
        alignedSeq.append(each)
    seqNum = 2
    seqLenth = len(alignedSeq[0]) - 1  
    seqArray = [([" "] * (seqLenth)) for i in range(seqNum)]  
    for i in range(0, seqNum):
        for j in range(0, seqLenth):
            Column = alignedSeq[i][j]
            seqArray[i][j] = (Column)
    
    seqArrayTransposed = zip(*seqArray)

    fileOut = open(dir + file + ".transposed", "w")
    fileOut.write("\t".join(name) + "\n")
    fileOut.write("\t".join(eachSeqStartAlignedPos) + "\n")  
    seqArrayTransposed = list(seqArrayTransposed)

    for l in seqArrayTransposed:
        fileOut.write("\t".join(l) + "\n")  
    fileOut.close()
    fileIn.close()

    fileOutName = file + ".transposed"
    return fileOutName, dir





