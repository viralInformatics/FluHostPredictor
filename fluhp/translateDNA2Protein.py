# !/usr/bin/python
# -*- coding: utf-8 -*-

import os
# chmod +x /usr/local/lib/python3.6/site-packages/app/blast/bin/blastall
# chmod +x /usr/local/lib/python3.6/site-packages/18Mid/translatePerl/translate/DNA2protein6.pl

def translate(file, DNASeqDir, DIR, outFileDir):

    DNA2protein6Dir = os.path.join(DIR, "18Mid/translatePerl/translate/")
    blastBinDir = os.path.join(DIR, "app/blast/bin/")
    forblastDir = os.path.join(DIR, "18Mid/translatePerl/translate/forblast/")
    dataDir = os.path.join(DIR, "18Mid/translatePerl/data/")
    mafftDir = os.path.join(DIR, "app/mafft/mafft-7.158-without-extensions/scripts/")
    outputName = file.replace(".fasta", ".trans2protein.fas")
    commands = [
        f"chmod +x {os.path.join(blastBinDir, 'blastall')}",
        f"chmod +x {os.path.join(DNA2protein6Dir, 'DNA2protein6.pl')}"
    ]

    for command in commands:
        os.system(command)

    command = f"perl {os.path.join(DNA2protein6Dir, 'DNA2protein6.pl')} {os.path.join(DNASeqDir, file)} {blastBinDir} {forblastDir} {outFileDir} {dataDir} {os.path.join(DIR,outFileDir, outputName)} {mafftDir}"

    # print(command)
    os.system(command)

    dicCDS = {}
    with open(os.path.join(DIR,outFileDir, outputName), 'r') as fileIn, open(
            os.path.join(outFileDir, outputName.replace(".trans2protein.fas",".trans2protein.fasta")), 'w') as fileOut:
        text = fileIn.readlines()
        for each in text:
            if each.startswith('>'):
                key = each.split('_', 1)[0].lstrip('>')
                value = each.split('_', 1)[1]
                dicCDS[key] = dicCDS.get(key, "") + ',' + value.strip() if key in dicCDS else value.strip()
                each = each.split("(")[0] + "\n"
            fileOut.write(each)

    for eachKey in dicCDS:
        dicCDS[eachKey] = dicCDS[eachKey].strip() if dicCDS[eachKey].strip() else 'Unknown'
    dicCDS = sorted(dicCDS.items(), key = lambda d: d[0])

    return outputName.replace(".trans2protein.fas",".trans2protein.fasta"), outFileDir, dicCDS


def makeProteinFileForDownload(dirUserTemp, file, dirUserOut, dicOriginalName, listProteinType):

    dicProteinType = {f'>{k[0]}': k[1] for k in listProteinType}

    with open(os.path.join(dirUserTemp, file), 'r') as fileIn, open(os.path.join(dirUserOut, file + '.annotation.fa'),
                                                                    'w') as fileOut:
        for eachLine in fileIn:
            if '>' in eachLine:
                seqSTDname = eachLine.strip().split('_')[0]
                annotation = f'_{dicProteinType[seqSTDname]}|{dicOriginalName[seqSTDname].strip(">")}' if seqSTDname in dicOriginalName else ''
                fileOut.write(eachLine.strip() + annotation + '\n')
            else:
                fileOut.write(eachLine)
    return os.path.join(dirUserOut, file + '.annotation.fa')

if __name__ == '__main__':
    pass
