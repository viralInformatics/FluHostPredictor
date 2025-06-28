def getDBSiteInfo(file,dir,proteinType):
    fileIn = open(dir+file,"r",encoding="utf-8")
    text = fileIn.readlines()
    fileIn.close()
    for each in text:
        if len(each)>1 and each.split(">")[0]==proteinType:
            # print(each.split(">")[1])
            return each.split(">")[1]


# getDbSiteInfo("DBSiteInfo","/home/think/18Mid/standard_seq/out/outout/","PB2")
