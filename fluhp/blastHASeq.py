import os


def renameFastaSeq():
    for i in range(1, 17):
        fileIn = open('../18Mid/antigen/HA_blast/HA_protein/H' + str(i) + '.faa', 'r')
        textIn = fileIn.readlines()
        fileIn.close()

        fileOut = open('../18Mid/antigen/HA_blast/H' + str(i), 'w')
        for eachLine in textIn:
            eachLine = eachLine.strip('\n')
            if '>' in eachLine:
                print(eachLine.split(' | '))
                eachLine = '>' + eachLine.split(' | ')[1] + '____GISAID_' + eachLine.split(' | ')[0].strip(
                    '>') + '____PMID_0\n'
            else:
                eachLine = eachLine + '\n'
            fileOut.write(eachLine)
        fileOut.close()


def makeSeqDic():
    for i in range(1, 17):
        In = '../18Mid/antigen/HA_blast/HA_protein/H' + str(i)
        fileIn = open(In, 'r')
        textIn = fileIn.read()
        fileIn.close()
        dic = {}
        fileOut = open(In + '.dic', 'w')
        for eachSeq in textIn.split('>'):
            key = eachSeq.split('\n')[0].split('____')[0]
            value = '>' + eachSeq
            dic[key] = value
        fileOut.write(str(dic))
        fileOut.close()
        # break


def removeDup():
    for i in range(1, 17):
        Out = '../18Mid/antigen/HA_blast/HA_protein/H' + str(i)
        fileIn = open(Out + '.dic', 'r')
        dic = eval(fileIn.read())
        fileIn.close()

        fileOut = open(Out, 'w')
        for each in dic:
            fileOut.write(dic[each])
        fileOut.close()


def makeBLastDB():
    import os
    for i in range(1, 17):
        In = '../18Mid/antigen/HA_blast/HA_protein/H' + str(i)
        Out = '../18Mid/antigen/HA_blast/HA_DB/H' + str(i)
        os.system('../app/blast+/bin/makeblastdb -in ' + In + ' -dbtype prot -parse_seqids -out ' + Out + '/H' + str(i))

def blastHASeq(prefixDir, blastQuerySeqHA, HAType, tempDir, mafftDir):
    # prefixDir = "D:/Learning/edegfile/think/platform"
    print("blastHASeq")
    with open(f'{prefixDir}/18Mid/antigen/HA_blast/HA_protein/{HAType}.dic', 'r') as fileDic:
        dic = eval(fileDic.read())

    os.system(
        # D:/user/app/blast-BLAST_VERSION+/bin/blastp
        f'blastp -db {prefixDir}/18Mid/antigen/HA_blast/HA_DB/{HAType}/{HAType} -query {tempDir}{blastQuerySeqHA} -out {tempDir}{blastQuerySeqHA}.blastHAOut -evalue 1e-5 -num_threads 1 -outfmt 7 \n')

    # 处理BLAST结果
    with open(f'{tempDir}{blastQuerySeqHA}.blastHAOut', 'r') as fileIn:
        textIn = fileIn.read()
    blastResult = textIn.split(' hits found\n')[1].split('# BLAST processed')[0].strip('\n')

    with open(f'{tempDir}HAHitSeq', 'w') as fileHitSeq:
        lines = blastResult.split('\n')
        for eachLine in lines[:min(5, len(lines))]:
            key = eachLine.split('\t')[1].split('____')[0]
            fileHitSeq.write(dic[key])

    os.system(
        f"perl {prefixDir}/18Mid/antigen/predAV.pl {HAType} {tempDir}{blastQuerySeqHA}"
        f" {tempDir}HAHitSeq {tempDir}{blastQuerySeqHA}.antigenInfo.blast {prefixDir}/18Mid/antigen/modelData/ {tempDir} ")
    # print( f"perl {prefixDir}/18Mid/antigen/predAV.pl {HAType} {tempDir}{blastQuerySeqHA}"
    #     f" {tempDir}HAHitSeq {tempDir}{blastQuerySeqHA}.antigenInfo.blast {prefixDir}/18Mid/antigen/modelData/ {tempDir} ")


# c = subprocess.getoutput(dirBlast+blastType+" -db "+DBDir+dataBaseName+" -query "+querySeqFileDir+querySeqFile+" -out "+outFileDir+outName+" -evalue "+eValue+" -num_threads 2 -outfmt "+outfmt+" -num_descriptions 5 -num_alignments 5\n"+"cd "+outFileDir+"\n")

if __name__ == "__main__":
    # renameFastaSeq()
    # makeBLastDB()
    # removeDup()
    # makeSeqDic()
    # blastHASeq(blastQuerySeqHA = 'querySeq1_HA',HAType = 'H5',tempDir = '../temp/',mafftDir = '../app/mafft/mafft-7.158-without-extensions/scripts/')
    blastHASeq(blastQuerySeqHA = 'querySeq1', HAType = 'H1', tempDir = '/home/think/platform/temp/',
               mafftDir = '/home/think/platform/app/mafft/mafft-7.158-without-extensions/scripts/')
