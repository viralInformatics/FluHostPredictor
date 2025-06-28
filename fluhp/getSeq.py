# -*- coding:UTF-8 -*-

def get(fileName="test",dir="",seqName=""):
    file = open(dir+"/"+fileName,"r")
    text = file.readlines()
    flag = 0
    seq = ""
    for eachline in text:
        if eachline.find(">")!=-1 and eachline.find(seqName)!=-1:
            flag = flag + 1
            fileOutName= str(eachline.replace(">","").replace("\n",""))
            fileOut = open(dir+fileOutName,"w")
            fileOut.write(eachline)
            seq = seq+eachline
        if flag == 1 and eachline.find(">")==-1:
            fileOut.write(eachline)
            seq = seq+eachline
            fileOut.close()
            flag = 0
            break
    file.close()
    return seq,fileOutName
# get(fileName = "test_query_seq2.fas",dir = "//home//think//18Mid//standard_seq/out//",seqName = "H2_standard_KC899733")
