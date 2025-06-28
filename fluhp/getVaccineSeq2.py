from getSingleFastaSeq import getSeq
fileSeqDir = "/home/think/18Mid/GisAidSeq/"
fileSeqName= "margeNCBIGISAID.fa"
fileQueryDir = "/home/think/18Mid/GisAidSeq/"
fileQueryName = "Virus_seqName.txt"
fileName = open(fileQueryDir+fileQueryName,'r')
text = fileName.readlines()
fileName.close()
for each in text :
    if len(each.strip())>1 and len(each.strip())<=3:
        print(each.strip())
        fileOut = open("/home/think/18Mid/GisAidSeq/aa/"+ each.strip(), 'w')
        # fileOut.write("@@"+each.strip()+"\n")
        # continue

    if len(each.split('\t'))>1 and len(each.split("\t")[1])>1:
        print(getSeq(fileSeqName,fileSeqDir,each.split("\t")[1]))
        seqNameInDB, seq = getSeq(fileSeqName, fileSeqDir, each.split("\t")[1])
        fileOut.write(">"+each.split("\t")[0]+"\n"+seq+"\n")

# fileOut.close()