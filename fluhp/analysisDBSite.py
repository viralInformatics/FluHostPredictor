def analysisDBSite(file,dir,DBSite,type):
    if type=="Virulence":
        print("analysisSite ...")
        import re
        fileQuerySeqSite =  open(dir+file,'r')
        text = fileQuerySeqSite.readlines()
        n = 0
        querySeqSite = ''
        for each in text :
            n = n + 1
            if n>1:
                querySeqSite = querySeqSite + each.strip("\n")
        fileQuerySeqSite.close()

        dic = eval(querySeqSite)
        # print(dic)
        # print(DBSite)
        dicDB = eval(DBSite)
        # print(dicDB)
        dicSingle = {}
        dicLong = {}
        outInfo = ""
        for numInDBSite in dicDB:
            if len(str(numInDBSite).split("&"))<=1:                     
                for eachQuerySite in dic[1]:
                    numInStandardSeq = (str(eachQuerySite).split("__")[-1])     
                    AAInQuerySeq = str(eachQuerySite).split("__")[0]           
                    if numInStandardSeq == numInDBSite:
                        dicSingle[numInDBSite] = AAInQuerySeq+"<>"+dicDB[numInDBSite]
                        # print(dicSingle)
            else:
                seq = ""
                for each in (numInDBSite).split("&"):               
                    for eachQuerySite in dic[1]:
                        numInStandardSeq = (str(eachQuerySite).split("__")[-1])
                        AAInQuerySeq = str(eachQuerySite).split("__")[0]
                        if numInStandardSeq == each:
                            seq = seq + AAInQuerySeq
                            dicLong[numInDBSite] = seq+"<>"+dicDB[numInDBSite]
                            # print(dicLong)
        dicMerged2 = dict(dicSingle, **dicLong)  
        # print(dicMerged2)

        for eachmaker in sorted(dicMerged2.items(), key=lambda d:int(d[0].split('&')[0])):
            eachmaker = eachmaker[0]
            querySeq = dicMerged2[eachmaker].split("<>")[0]
            for string in dicMerged2[eachmaker].split('<>')[1].split('+'):
                # print(string)
                DBSeq = string.split("~")[0]
                reference = '~'.join(string.split("~")[1:-1])
                describe = string.split("~")[-1]
                zztmp = eachmaker.split("~")[0]
                #print(zztmp)
                #modified by zyq 09032018 begin
                #eachmakerSplited = eachmaker.split("&")
                eachmakerSplited = zztmp.split("&")
                #print(eachmaker)
                #modified by zyq 09032018 end
                if len(eachmakerSplited)>1:
                    # print(eachmakerSplited[-1])
                    if int(eachmakerSplited[-1])-int(eachmakerSplited[0]) == len(eachmakerSplited)-1:
                        eachmaker = eachmakerSplited[0]+'-'+eachmakerSplited[-1]


                deletionState = "123456789"
                if "AnyDeletion" in DBSeq or "NoDeletion" in DBSeq or("CompleteDeletion" in DBSeq):
                    if  querySeq.find("-")==-1:
                        # print("NoDeletion")
                        deletionState = "NoDeletion"
                    elif querySeq.replace("-","")=="":
                        # print("CompleteDeletion")
                        deletionState = "CompleteDeletion"
                    if querySeq.find("-") != -1 and querySeq.replace("-", "") != "":
                        # print("AnyDeletion")
                        deletionState = "AnyDeletion"

                if re.match(DBSeq.split("|")[1],querySeq) or re.match(DBSeq.split("|")[1],deletionState):
                    # print(eachmaker,querySeq,DBSeq,"Virulent",reference)
                    outInfo = outInfo + "\t" + eachmaker + querySeq + "\tVirulent\t" + reference + '\t' + describe + "\n"
                elif re.match(DBSeq.split("|")[0],querySeq) or DBSeq.split("|")[0]=="other" or re.match(DBSeq.split("|")[0],deletionState):
                    # print(eachmaker,querySeq,DBSeq,"Non Virulent",reference)
                    s=1
                    # outInfo = outInfo+"\t"+eachmaker+"\t"+querySeq+"\t"+DBSeq+"\tNon Virulent\t"+ reference+"\n"
                elif DBSeq.split("|")[1]=="other":
                    # print(eachmaker,querySeq,DBSeq,"Virulent",reference)
                    outInfo = outInfo + "\t" + eachmaker + querySeq + "\tVirulent\t" + reference + '\t' + describe + "\n"
                else:
                    # outInfo = outInfo+"\t"+eachmaker+"\t"+querySeq+"\t"+DBSeq+"\tMismatch\t"+ reference+"\n"
                    s=1

        # print(outInfo)
        return outInfo
    else:
        print("analysisDBSite ...")
        import re
        fileQuerySeqSite = open(dir + file, 'r')
        text = fileQuerySeqSite.readlines()
        n = 0
        querySeqSite = ''
        for each in text:
            n = n + 1
            if n > 1:
                querySeqSite = querySeqSite + each.strip("\n")
        fileQuerySeqSite.close()

        dic = eval(querySeqSite)
        # print(dic)
        # print(DBSite)
        dicDB = eval(DBSite)
        # print(dicDB)
        dicSingle = {}
        dicLong = {}
        outInfo = ""
        for numInDBSite in dicDB:
            if len(str(numInDBSite).split("&")) <= 1: 
                for eachQuerySite in dic[1]:
                    numInStandardSeq = (str(eachQuerySite).split("__")[-1])  
                    AAInQuerySeq = str(eachQuerySite).split("__")[0]  
                    if numInStandardSeq == numInDBSite:
                        dicSingle[numInDBSite] = AAInQuerySeq + "<>" + dicDB[
                            numInDBSite]  
                        # print(dicSingle)
            else:
                seq = ""
                for each in (numInDBSite).split("&"):  
                    for eachQuerySite in dic[1]:
                        numInStandardSeq = (str(eachQuerySite).split("__")[-1])
                        AAInQuerySeq = str(eachQuerySite).split("__")[0]
                        if numInStandardSeq == each:
                            seq = seq + AAInQuerySeq
                            dicLong[numInDBSite] = seq + "<>" + dicDB[numInDBSite]
                            # print(dicLong)
        dicMerged2 = dict(dicSingle, **dicLong)  
        # print(dicMerged2)

        for eachmaker in sorted(dicMerged2.items(), key=lambda d:int(d[0].split('&')[0])):
            eachmaker = eachmaker[0]  
            querySeq = dicMerged2[eachmaker].split("<>")[0]
            for string in dicMerged2[eachmaker].split('<>')[1].split('+'):
                DBSeq = string.split("~")[0]
                ICFold = DBSeq.split("|")[2]
                reference = string.split("~",1)[1]
                eachmakerSplited = eachmaker.split("&")
                if len(eachmakerSplited)>1:
                    if int(eachmakerSplited[-1])-int(eachmakerSplited[0]) == len(eachmakerSplited)-1:
                        eachmaker = eachmakerSplited[0]+'-'+eachmakerSplited[-1]


                deletionState = "123456789"
                if "AnyDeletion" in DBSeq or "NoDeletion" in DBSeq or ("CompleteDeletion" in DBSeq):
                    if querySeq.find("-") == -1:
                        # print("NoDeletion")
                        deletionState = "NoDeletion"
                    elif querySeq.replace("-", "") == "":
                        # print("CompleteDeletion")
                        deletionState = "CompleteDeletion"
                    if querySeq.find("-") != -1 and querySeq.replace("-", "") != "":
                        # print("AnyDeletion")
                        deletionState = "AnyDeletion"

                if re.match(DBSeq.split("|")[1], querySeq) or re.match(DBSeq.split("|")[1], deletionState):
                    # print(eachmaker, querySeq, DBSeq, "Resistance", reference)
                    outInfo = outInfo + "\t" + eachmaker + querySeq +"\t"+ICFold+"\tResistant\t" + reference + "\n"
                elif re.match(DBSeq.split("|")[0], querySeq) or DBSeq.split("|")[0] == "other" or re.match(
                        DBSeq.split("|")[0], deletionState):
                    # print(eachmaker, querySeq, DBSeq, "sensitive", reference)
                    # outInfo = outInfo + "\t" + eachmaker + "\t" + querySeq + "\t" + DBSeq + "\tsensitive\t" + reference + "\n"
                    s=1
                elif DBSeq.split("|")[1] == "other":
                    # print(eachmaker, querySeq, DBSeq, "Resistance", reference)
                    outInfo = outInfo + "\t" + eachmaker + querySeq +"\t"+ICFold+ "\tResistant\t" + reference + "\n"
                else:
                    # print(eachmaker, querySeq, DBSeq, "Mismatch", reference)
                    # outInfo = outInfo + "\t" + eachmaker + "\t" + querySeq + "\t" + DBSeq + "\tMismatch\t" + reference + "\n"
                    s=1

        # print(outInfo)
        return outInfo
