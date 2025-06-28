def translateSeqInThisFile():
    import os
    import  subprocess
    for root ,dirs,files in os.walk("/home/think/18Mid/translatePerl/standard_seq"):
        for file in files:
            # print(file,root)
            DNA2protein6Dir = "/home/think/18Mid/translatePerl/translate//"
            DNASeqDir = "/home/think/18Mid/translatePerl/standard_seq/"
            DNASeq = file
            blastBinDir = "/home/think/app/blast/bin//"
            forblastDir = "/home/think/18Mid/translatePerl/translate/forblast//"
            dataDir = "/home/think/18Mid/translatePerl/data//"
            tempDir = "/home/think/18Mid/temp//"
            outPutDir = "/home/think/18Mid/translatePerl/standard_seq_protein//"
            # outputName = DNASeq.replace(".fas",".trans2protein.fas")
            outputName = DNASeq
            # print("perl "+DNA2protein6Dir+"DNA2protein6.pl "+DNASeqDir+DNASeq+" "+blastBinDir+" "+forblastDir+" "+tempDir+" "+dataDir+" "+outPutDir+outputName)
            result = os.system("perl "+DNA2protein6Dir+"DNA2protein6.pl "+DNASeqDir+DNASeq+" "+blastBinDir+" "+forblastDir+" "+tempDir+" "+dataDir+" "+outPutDir+outputName)
            # print(result)


# translateSeqInThisFile()

