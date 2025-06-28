def markHost(fileHost,fileSeq,dir):
    seq = open(dir+fileSeq,'r')
    host = open(dir+fileHost,'r')
    fileOut = open(dir+fileSeq+".markedHost.fna.fas",'w')
    textSeq = seq.readlines()
    textHost = host.readlines()
    dic = {}
    for eachHost in textHost:
        eachHost = eachHost.split("\t")
        dic[eachHost[0]] = eachHost[1]
        # print(dic)


    for each in textSeq:
        if each[0]==">":
            gb = each.split("|")[3]
            # each = ">"+dic[gb]+"gb"+gb+"\n"
            # each = each.replace(" ","_")
            # each = dic[gb] + "\n"
            # name = str(dic[gb]).replace(" ","_")
            # print(name)
            each = each.replace(">",str(">"+(str(dic[gb]).replace(" ","-")+"_")))

            print(each)
        fileOut.write(each)

    seq.close()
    host.close()
    fileOut.close()
# markHost("influenza_na.dat","influenza.fna","/home/think/markedHost/Influenza_Virus_Sequence_Database/")