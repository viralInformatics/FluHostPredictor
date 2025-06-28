import os
import subprocess
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
    dirUser = '/data/FluPhenotype/temp/PROT/'
    file = 'Isolate1.stdName'
    dirUserTemp = dirUser + '/temp/'
    dirUserOut = dirUser + '/result/'
    dicOriginalName = eval(
        "{'>querySeq1': '>querySeq3_NA111111111111111111', '>querySeq5': '>querySeq1_NS2111111111111111111', '>querySeq6': '>querySeq2_M1111111111111111111', '>querySeq10': '>querySeq7_PB1111111111111111111', '>querySeq9': '>querySeq4_NP111111111111111111', '>querySeq4': '>querySeq1_NS1111111111111111111', '>querySeq8': '>querySeq5_HA111111111111111111', '>querySeq2': '>querySeq6_PA111111111111111111', '>querySeq11': '>querySeq7_PB1-F2111111111111111111', '>querySeq12': '>querySeq8_PB2111111111111111111', '>querySeq3': '>querySeq6_PA-X111111111111111111', '>querySeq7': '>querySeq2_M2111111111111111111'}")
    listProteinType = eval(
        "[('querySeq1', 'N6'), ('querySeq2', 'PA'), ('querySeq3', 'PA-X'), ('querySeq4', 'NS1'), ('querySeq5', 'NS2'), ('querySeq6', 'M1'), ('querySeq7', 'M2'), ('querySeq8', 'H1'), ('querySeq9', 'NP'), ('querySeq10', 'PB1'), ('querySeq11', 'PB1-F2'), ('querySeq12', 'PB2')]")

    makeProteinFileForDownload(dirUserTemp, file, dirUserOut, dicOriginalName, listProteinType)

    dirUser = '/data/FluPhenotype/temp/DNA/'
    file = 'Isolate1.stdName'
    dirUserTemp = dirUser + '/temp/'
    dirUserOut = dirUser + '/result/'
    dicOriginalName = eval(
        "{'>querySeq1': '>test_NS', '>querySeq5': '>test_HA', '>querySeq4': '>test_NP', '>querySeq6': '>test_PA', '>querySeq2': '>test_MP', '>querySeq8': '>test_PB2', '>querySeq3': '>test_NA', '>querySeq7': '>test_PB1'}")
    listProteinType = eval(
        "[('querySeq1', 'NS'), ('querySeq2', 'MP'), ('querySeq3', 'N6'), ('querySeq4', 'NP'), ('querySeq5', 'H1'), ('querySeq6', 'PA'), ('querySeq7', 'PB1'), ('querySeq8', 'PB2')]")

    makeProteinFileForDownload(dirUserTemp, file, dirUserOut, dicOriginalName, listProteinType)
