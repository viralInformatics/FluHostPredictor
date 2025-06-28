def sortResult(fileDir,fileName):
    fileIn = open(fileDir+fileName,'r', encoding='utf8')
    textIn = fileIn.readlines()
    fileIn.close()
    fileIn = open(fileDir+fileName,'r', encoding='utf8')
    textInRead = fileIn.read()
    fileIn.close()

    textOneLine  = textInRead
    fileOut = open(fileDir+fileName,'w',encoding='utf-8')
    if textInRead.count('ProteinType->')>1:
        textFixed = textOneLine.split(']\n\n')[0]+']\n\n'
        fileOut.write(textFixed)
        textOneLine = textOneLine.split(']\n\n')[-1]
    else:
        textFixed = textOneLine.split(')\n\n')[0] + ')\n\n'
        fileOut.write(textFixed)
        textOneLine = textOneLine.split(')\n\n')[-1]
    textFixedEnd = '#'+textOneLine.split('#',1)[-1].split('\n')[0]+'\n'
    textHostMarker = textOneLine.split('down\n')[1]
    textOneLine = textOneLine.split('#',1)[0]
    textOneLineSplit = textOneLine.split('\n')

    flag = 0
    textOneLine = ''
    for each in textOneLineSplit:
        if (each[0:len('querySeq')] == 'querySeq') and ('>' in each):
            flag = 1
        if len(each.strip('\n'))==0:
            flag = 0
            each = '\n'
        if flag == 1:
            each = each+'&&&&'
        textOneLine = textOneLine + each
    textOneLine = textOneLine.split('\n')
    for eachType in ['PB2','PB1','PB1-F2','PA','PA-X','NP','M1','M2','NS1','NS2','H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12','H13','H14','H15','H16','N1','N2','N3','N4','N5','N6','N7','N8','N9','Unknown']:
        for eachOneLine in textOneLine:
            eachOneLine = eachOneLine+'\n'
            if len(eachOneLine.strip())>1 and eachOneLine.split('\t')[2].split('&&&&')[0] == eachType:
                # print(eachOneLine)
                fileOut.write(eachOneLine.replace('&&&&','\n'))
            #     fileOut.write(eachOneLine)
    fileOut.write(textFixedEnd)

    flag = 0
    textOneLine = ''
    textOneLineSplit = textHostMarker.split('\n')
    for each in textOneLineSplit:
        if (each[0:len('querySeq')] == 'querySeq') and ('>' in each):
            flag = 1
        if len(each.strip('\n'))==0:
            flag = 0
            each = '\n'
        if flag == 1:
            each = each+'&&&&'
        textOneLine = textOneLine + each
    textOneLine = textOneLine.split('\n')
    for eachType in ['PB2','PB1','PB1-F2','PA','PA-X','NP','M1','M2','NS1','NS2','H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12','H13','H14','H15','H16','N1','N2','N3','N4','N5','N6','N7','N8','N9','Unknown']:
        for eachOneLine in textOneLine:
            eachOneLine = eachOneLine+'\n'
            if len(eachOneLine.strip())>1 and eachOneLine.split('\t')[2].split('&&&&')[0] == eachType:
                # print(eachOneLine)
                fileOut.write(eachOneLine.replace('&&&&','\n'))
            #     fileOut.write(eachOneLine)


    fileOut.close()

