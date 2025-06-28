filein = open("/home/think/platform/18Mid/translatePerl/standard_seq_protein/PA.fas")
fileOut1 = open("/home/think/platform/18Mid/translatePerl/standard_seq_protein/PA1.fas",'w')
fileOut2 = open("/home/think/platform/18Mid/translatePerl/standard_seq_protein/PA-X.fas",'w')
import getSingleFastaSeq
text = filein.readlines()
for each in text:
    if each[0]=='>'and "_PA("in each:
        # print(each)
        seq = getSingleFastaSeq.getSeq("/home/think/platform/18Mid/translatePerl/standard_seq_protein/","PA.fas",each,False)
        fileOut1.write(seq[0]+"\n"+seq[1]+'\n')
    elif each[0]=='>'and "_PA-X("in each:
        print(each)
        seq = getSingleFastaSeq.getSeq("/home/think/platform/18Mid/translatePerl/standard_seq_protein/","PA.fas",each,False)
        fileOut2.write(seq[0]+"\n"+seq[1]+'\n')


fileOut2.close()
fileOut1.close()
filein.close()
