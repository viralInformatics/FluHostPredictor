def getSite(fileName,dir,seqName,site):
    file = open(dir+"/"+fileName,"r")
    text = file.readlines()
    # print(text)
    for each in text :
        ColumnNumber = each.split().index(seqName)
        print(each.split().index(seqName))
        break
    line = -2
    for each in text :
        line = line + 1
        if line==site:
            print(each.split()[ColumnNumber],end="")




# getSite("H6PlusQuerySeq.fas.fftns2.out.fas.transposed","/home/think/18Mid/standard_seq/out//",">standard_CY110730",1000)