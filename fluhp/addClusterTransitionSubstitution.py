# -*- coding:UTF-8 -*-
import re
import os
import csv
import pandas as pd

def extract_year(protein_str):
    match = re.search(r'/(\d{2,4})$', protein_str)
    if match:
        year_str = match.group(1)
        if len(year_str) == 2:
            year = int(year_str)
            if year > 22:
                return 1900 + year
            else:
                return 2000 + year
        elif len(year_str) == 4:
            return int(year_str)
    return None

def sort_by_year(df, column_name):
    df['Year'] = df[column_name].apply(lambda x: extract_year(x.split(',')[-1].strip()) if ',' in x else extract_year(x.strip()))

    df_sorted = df.sort_values(by='Year', ascending=False).reset_index(drop=True)

    df_sorted.drop(columns=['Year'], inplace=True)

    return df_sorted

def addClusterTransitionSubstitution(resultFileDir, resultFileName, DBFileName, dir, proteinTypeList):
    dic = {}
    dicFile = open(dir + '/' + DBFileName, 'r', encoding = 'utf-8')
    for eachLine in dicFile:
        dic[eachLine.split('>')[0]] = eval(eachLine.split('>')[1].strip('\n'))
    dicFile.close()

    dicDBAllProteinType = {}
    fileInDic = open(dir + '/DBSiteInfoClusterTransition_deleteionSite', 'r')
    for eachLine in fileInDic:
        dicDBAllProteinType[eachLine.split('>')[0]] = eachLine.split('>', 1)[1].strip('\n')
    fileInDic.close()
    pattern = re.compile(r'^querySeq\d+')
    isolate_name = resultFileName.split(".")[0]
    if os.path.getsize(resultFileDir + '/'+ resultFileName) != 0:
        fileIn = open(resultFileDir + '/'+ resultFileName)
        textIn = fileIn.readlines()
        fileIn.close()
        # fileOut = open(resultFileDir + '/'+ resultFileName, 'w')
        with open(resultFileDir + '/' + resultFileName, 'w', encoding = 'utf-8', newline = '') as fileOut:
            writer = csv.writer(fileOut, quoting=csv.QUOTE_MINIMAL)
            if "blast" in resultFileName:
                writer.writerow([
                    "Query sequence",
                    "Top 5 most similar HA proteins",
                    "Antigenic relationship",
                    "Predicted odds ratio",
                    "Genetic difference",
                    "Antigenic cluster-transition substitution"
                ])
            else:
                writer.writerow([
                    "Query sequence",
                    "Vaccine or reference strain",
                    "Antigenic relationship",
                    "Predicted odds ratio",
                    "Genetic difference",
                    "Antigenic cluster-transition substitution"
                ])
            for eachLine in textIn:
                # print(eachLine)
                eachLine = eachLine.replace('"','')
                clusterTransitionInfo_deletionSite = ""
                if eachLine.startswith('querySeq'):
                    proteinType = eachLine.split(',')[0].split(" ")[-1].strip()
                    if proteinType not in dic.keys() and proteinType not in dicDBAllProteinType.keys():
                        eachLine = pattern.sub(isolate_name, eachLine)
                        writer.writerow([eachLine.strip('\n')])
                        continue
                    dicSite = dic[proteinType]

                    if proteinType in dicDBAllProteinType:
                        variantSiteList = eachLine.split(',')[4].split()
                        clusterTransitionInfo_deletionSite = getDeletionSiteClusterTransitionSubstitution(
                            eval(dicDBAllProteinType[proteinType]), variantSiteList)
                    clusterTransitionInfo = analysesSite(eachLine, dicSite, proteinType)
                    if clusterTransitionInfo and clusterTransitionInfo_deletionSite:
                        clusterTransitionInfo = clusterTransitionInfo + "&&" + clusterTransitionInfo_deletionSite
                    else:
                        clusterTransitionInfo = clusterTransitionInfo + clusterTransitionInfo_deletionSite
                    # clusterTransitionInfo = clusterTransitionInfo + "\t" + clusterTransitionInfo_deletionSite
                    eachLine = pattern.sub(isolate_name, eachLine)
                    row_data = eachLine.strip('\n').split(',') + [clusterTransitionInfo]
                    writer.writerow(row_data)
        with open(resultFileDir + '/' + resultFileName, 'r', encoding = 'utf-8') as file:
            lines = file.readlines()

        with open(resultFileDir + '/' + resultFileName, 'w', encoding = 'utf-8', newline = '') as file:
            for line in lines:
                if line.startswith('"') and line.endswith('"\n'):
                    line = line[1:-2] + '\n'
                file.write(line)

        df = pd.read_csv(resultFileDir + '/' + resultFileName)
        if 'Top 5 most similar HA proteins' in df.columns:
            sorted_df = sort_by_year(df, 'Top 5 most similar HA proteins')
        else:
            sorted_df = sort_by_year(df, 'Vaccine or reference strain')

        sorted_df.to_csv(resultFileDir + '/' + resultFileName, index = False, encoding = 'utf-8')


def process_dbsite(dbsite, site_info_dic, dic_site, out_info):
    AAChange = ['', '']
    DBseq, reference, describe = dic_site[dbsite].split('~')[0], '~'.join(dic_site[dbsite].split("~")[1:-1]), \
                                 dic_site[dbsite].strip('\n').split("~")[-1]

    if '&' in dbsite:
        for eachDBSite in dbsite.split('&'):
            if eachDBSite not in site_info_dic:
                return out_info
            AAChange[0] += site_info_dic[eachDBSite].split('__')[1]
            AAChange[1] += site_info_dic[eachDBSite].split('__')[2]
    else:
        if dbsite not in site_info_dic:
            return out_info
        AAChange[0] = site_info_dic[dbsite].split('__')[1]
        AAChange[1] = site_info_dic[dbsite].split('__')[2]

    text = AAChange[0] + dbsite + AAChange[1] + '+' + reference + '+' + describe
    if (re.match(DBseq.split("|")[0], AAChange[0]) and re.match(DBseq.split("|")[1], AAChange[1])) or \
            (re.match(DBseq.split("|")[0], AAChange[1]) and re.match(DBseq.split("|")[1], AAChange[0])):
        if out_info:
            out_info += '|' + text
        else:
            out_info = text

    return out_info


def analysesSite(eachLine, dicSite, proteinType):
    siteInfoDic = {}
    for eachSite in eachLine.strip('\n').split(',')[4].split():
        groups = re.match("([A-Z])(\d{1,3})([A-Z])", eachSite)
        if groups:
            site = groups.group(2)
            siteInfo = site + '__' + groups.group(1) + '__' + groups.group(3)
            siteInfoDic[site] = siteInfo

    outInfo = ''
    for DBSite in dicSite:
        outInfo = process_dbsite(DBSite, siteInfoDic, dicSite, outInfo)

    if proteinType == 'H3':
        DicSiteH3 = {
            '156': 'K|E~PMID:24264991~The single cluster-transition substitution for the HK68 to EN72 cluster transition was 155TY and 156KE alone was responsible for the antigenic difference between the TX77 and BK79 clusters.'
        }
        for DBSite in DicSiteH3:
            if DBSite in siteInfoDic:
                outInfo = process_dbsite(DBSite, siteInfoDic, DicSiteH3, outInfo)

    return outInfo


def getDeletionSiteClusterTransitionSubstitution(dicDB, variantSiteList):
    siteInfoDic = {}
    for eachSite in variantSiteList:
        site = ''
        AAChange = []
        for eachLetter in eachSite:
            if eachLetter in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                site += eachLetter
            else:
                AAChange.append(eachLetter)

        if site == '130' and (AAChange[1] not in ['-', 'X', 'x', 'K', 'k']):
            siteInfo = (AAChange[0], '#')
        else:
            siteInfo = (AAChange[0], AAChange[1])
        siteInfoDic[site] = siteInfo

    outInfo = ''
    for each in dicDB:
        flag = 0
        for eachLetter in each.split('&'):
            if (eachLetter in siteInfoDic) and (siteInfoDic[eachLetter][1] in ['#']):
                continue
            if eachLetter not in siteInfoDic:
                continue
            if siteInfoDic[eachLetter][0] in ['-', 'X', 'x']:
                if (eachLetter in siteInfoDic) and (siteInfoDic[eachLetter][1] not in ['-', 'X', 'x']):
                    flag = 1
            else:
                if (eachLetter in siteInfoDic) and (siteInfoDic[eachLetter][1] in ['-', 'X', 'x']):
                    flag = 1

        if flag == 1:
            DBseq = dicDB[each].split('~')[0]
            reference = '~'.join(dicDB[each].split("~")[1:-1])
            describe = dicDB[each].strip('\n').split("~")[-1]

            if each == '206&207&208&209&210&211&212&213&214&215': each = '220-loop '
            if outInfo:
                outInfo = outInfo + '|' + DBseq.split('|')[0] + each + DBseq.split('|')[
                    1] + '+' + reference + '+' + describe
            else:
                outInfo = DBseq.split('|')[0] + each + DBseq.split('|')[
                    1] + '+' + reference + '+' + describe
    return outInfo


if __name__ == '__main__':
    dirUser = '../'
    # addClusterTransitionSubstitution(resultFileDir=dirUser+'/result/',resultFileName = 'example'+'.result',DBFileName="DBSiteInfoClusterTransition",dir="/home/think/platform/18Mid/",proteinTypeList=eval("[('querySeq1', 'NS'), ('querySeq2', 'MP'), ('querySeq3', 'N6'), ('querySeq4', 'NP'), ('querySeq5', 'H1'), ('querySeq6', 'PA'), ('querySeq7', 'PB1'), ('querySeq8', 'PB2')]"))
    # addClusterTransitionSubstitution(resultFileDir = dirUser + '/test/result/',
    #                                  resultFileName = 'isolate1.fasta.antigenInfo.blast.csv',
    #                                  DBFileName = "DBSiteInfoClusterTransition", dir = "../18Mid/",
    #                                  proteinTypeList = eval(
    #                                      "[('querySeq1_NS1', 'NS1'), ('querySeq1_NS2', 'NS2'), ('querySeq2_M1', 'M1'), ('querySeq2_M2', 'M2'), ('querySeq3', 'N6'), ('querySeq4', 'NP'), ('querySeq5', 'H1'), ('querySeq6_PA', 'PA'), ('querySeq6_PA-X', 'PA-X'), ('querySeq7_PB1', 'PB1'), ('querySeq7_PB1-F2', 'PB1-F2'), ('querySeq8', 'PB2')]"))
    addClusterTransitionSubstitution(resultFileDir = dirUser + '/test/result/',
                                     resultFileName = 'isolate1.fasta.antigenInfo.csv',
                                     DBFileName = "DBSiteInfoClusterTransition", dir = "../18Mid/",
                                     proteinTypeList = eval(
                                         "[('querySeq1_NS1', 'NS1'), ('querySeq1_NS2', 'NS2'), ('querySeq2_M1', 'M1'), ('querySeq2_M2', 'M2'), ('querySeq3', 'N6'), ('querySeq4', 'NP'), ('querySeq5', 'H1'), ('querySeq6_PA', 'PA'), ('querySeq6_PA-X', 'PA-X'), ('querySeq7_PB1', 'PB1'), ('querySeq7_PB1-F2', 'PB1-F2'), ('querySeq8', 'PB2')]"))
    # 示例数据
    data = {
        'Protein': ['ABC/21', 'DEF/1998', 'GHI/01', 'JKL/2023', 'MNO/30', 'PQR/22'],
        'Value': [1, 2, 3, 4, 5, 6]
    }
    df = pd.DataFrame(data)

    # 对 'Protein' 列按年份排序
    sorted_df = sort_by_year(df, 'Protein')
    print(sorted_df)
