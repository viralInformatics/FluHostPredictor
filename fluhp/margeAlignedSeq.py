def MargeAlignedSeq(file="",dir=""):
    fileIn = open (dir+file,"r")
    text = fileIn.readlines()
    # print(text)
    flag = 0
    n = 0
    count = 100
    seqName = ["Query_1"]
    for each in text:
        n = n + 1
        if each[0:6]=="Query=":
            each = each.split("=")[-1].lstrip(" ")
            queryName = each.rstrip('\n')
            flag = 1
            count = 0
            print("\n"+queryName)
        if flag == 1:
            count = count+1
            if count >=7:
                if len(each)<2:
                    flag = 0
                    print(n,count)
                    break
                # temp = []
                temp = each.split()[0]
                seqName.append(temp)
    print(seqName,len(seqName))                          #得到序列名字
    seqNumber = len(seqName)                #包括Query_1的序列的个数


    LengthFirstLinequerySeq = 0
    for each in text:
        if len(each)>2:
            if each[0:6]=="Query_":
                querySeq = each.split()[-2]
                # print(each.split())
                pos1 = each.find(querySeq)
                pos2 = pos1 + len(querySeq)
                # print(querySeq,pos1,pos2)
                break

    seqNameAndPos=[]
    startPos = []
    flag = 0
    for name in seqName:
        print(1)
        for each in text:
            print(2)
            if each[0:len(name)] == name:
                startPos.append(each.split()[1])
                print(name)
                break
        seqNameAndPos = list(zip(seqName,startPos))
    filemarge = open(dir+file+".margeseq","w")
    for eachName in seqNameAndPos:
        n = 0
        # print(eachName)
        filemarge.write(eachName[0]+'\t'+eachName[1]+"\t")
        for eachLine in text:
            n = n + 1
            lineWrite = ""
            if eachLine[0:len(eachName[0])]==eachName[0] and n>seqNumber+20:
                    eachLine = eachLine[pos1:]
                    print(eachLine)
                    # eachLine = "".join(list(eachLine))
                    eachLine = "".join(list(filter(lambda x:x not in '0123456789',eachLine)))
                    eachLine = eachLine.rstrip("\n")[:-2]
                    filemarge.write(eachLine)
        filemarge.write("\n")
    filemarge.close()

    fileIn.close()

    return (str(file)+".margeseq")
# MargeAlignedSeq("query___H6_standard_KJ200805___In___H6___DB", "/home/think/18Mid/standard_seq/out/")