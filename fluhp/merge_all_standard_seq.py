import os
fileMarge = open("/home/think/18Mid/translatePerl/fileMargeProtein.fas","a")
for root ,dir ,file in os.walk("/home/think/18Mid/translatePerl/standard_seq_protein/",topdown=False):
    print(root)
    print(file)
    for fileName in file:
        fileTemp = open(root+fileName,"r")
        text = fileTemp.read()
        text = text.replace(">",str(">"+fileName.replace(".fas","_")))
        fileMarge.write(text)

        fileTemp.close()

fileMarge.close()

# makeblastdb -in a.fa -dbtype prot -parse_seqids -out protseq