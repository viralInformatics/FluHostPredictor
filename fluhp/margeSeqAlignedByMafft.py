# -*- coding:UTF-8 -*-

def margeSeqAlignedSeq(querySeq,queryFileName,standardDBDir,standardDB,mafftDir,outDir):
    print("AlignSeq ...")
    fileStandardSeq =  open(standardDBDir+standardDB,'r')
    text = fileStandardSeq.read()
    fileMargeQueryAndStandard = open(outDir+queryFileName+".MargeQueryAndStandard.fas",'w')
    margeFileName = queryFileName+".MargeQueryAndStandard.fas"
    fileMargeQueryAndStandard.write(querySeq)
    fileMargeQueryAndStandard.write(text)
    fileMargeQueryAndStandard.close()
    fileStandardSeq.close()


    import os
    c = os.system(mafftDir+"mafft --quiet "+outDir+margeFileName+" > "+outDir+margeFileName+".mafft")         #使用mafft进行序列比对
    # print(outDir+margeFileName+".mafft")
    margeFileName = margeFileName+".mafft"                                               #   比对后的输出文件名
    fileText = open(outDir+margeFileName,"r")
    text = fileText.readlines()
    fileTextInOneLine = open(outDir+margeFileName+".oneLine",'w')              #将序列合并成一行，便于操作
    fileOut = open(outDir+margeFileName+".out.fas","w")

    for each in text:
        if each.find(">")==-1:
            each = each.replace("\n","")
        else:
            each = each.replace("\n","\t__\t")
            each = each.replace(">","\n>")
        fileTextInOneLine.write(each)
    fileTextInOneLine.close()
    fileTextInOneLine = open(outDir+margeFileName+".oneLine",'r')
    textInOneLine = fileTextInOneLine.readlines()
    n=0
    for eachLine in textInOneLine:
        if len(eachLine)>2:
            hyphenNumber = 0
            for eachOne in eachLine.split("\t")[2]:
                if eachOne == "-":
                    hyphenNumber = hyphenNumber + 1
                else:
                    eachLine = eachLine.split("\t")[0]+"\t"+eachLine.split("\t")[1]+"\t"+"$"*hyphenNumber+eachLine.split("\t")[2][hyphenNumber:]
                    eachLine = eachLine.replace("__", str(-hyphenNumber))
                    if n==0:
                        eachLine=eachLine.replace("$","-")          #使得第一行查询序列的“-”符号数目不变，其他的替换为$符号
                    fileOut.write(eachLine)
                    # print(eachLine)
                    break

    fileOut.close()
    fileText.close()
    fileTextInOneLine.close()
    return margeFileName+".out.fas"

# margeSeqAlignedSeq("H6PlusQuerySeq.fas.fftns2", "/home/think/18Mid/standard_seq/out/")

