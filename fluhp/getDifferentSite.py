def getdifferentSite(file,fileDir,DBSite):
    print("getdifferentSite ...")
    dicDB = eval(DBSite)
    DBSiteNum = ""
    for eachDBNum in dicDB:
        DBSiteNum =DBSiteNum+"&"+eachDBNum
    fileIn = open(fileDir+file)
    fileOut = open(fileDir+file+".siteInfo","w")
    text = fileIn.readlines()
    line = 0
    siteInfo = ""
    hyphenNumber = {}           
    sifferentSiteDic = {}
    DBsiteDic = {}
    for each in text:
        line = line + 1
        if line == 1:
            
            name = each.split()
            
            continue
        if line == 2:
            startSite = each.split("\t")
            for i in range(0,len(name)):
                
                hyphenNumber[i]= 0
                sifferentSiteDic.setdefault(i,[])
                DBsiteDic.setdefault(i,[])
            
            continue
        
        colonm = 0
        for eachSite in each.strip("\n").replace("\t",''):          
            if colonm == 0:                                         
                firstSite = eachSite
            
            if eachSite == "-":                                     
                hyphenNumber[colonm] +=1
            if str(line-hyphenNumber[colonm]-2) in DBSiteNum.split("&"):
                DBsiteDic[colonm].append(firstSite + "__" + eachSite + "__" + str(line - hyphenNumber[0]-2) + "__" + str(line-hyphenNumber[colonm]-2))
                
                
            colonm = colonm + 1
    fileIn.close()
    fileOut.write(str(name)+str(DBsiteDic).replace("{","\n{"))
    fileOut.close()
    differentSiteFileName = file + ".siteInfo"
    return differentSiteFileName,sifferentSiteDic , name , DBsiteDic
