# !/usr/bin/python
# -*- coding: utf-8 -*-

def refreshStandardName(predictedProteinType, dic):
    predictedProteinTypeNew = []
    for eachType in predictedProteinType:
        if eachType[1] == "NS":
            predictedProteinTypeNew.append((eachType[0] + "_NS1", "NS1"))
            predictedProteinTypeNew.append((eachType[0] + "_NS2", "NS2"))
            dic[">" + eachType[0] + "_NS1"] = dic[">" + eachType[0]]
            dic[">" + eachType[0] + "_NS2"] = dic[">" + eachType[0]]
            dic.pop(">" + eachType[0])
        elif eachType[1] == "MP":
            predictedProteinTypeNew.append((eachType[0] + "_M1", "M1"))
            predictedProteinTypeNew.append((eachType[0] + "_M2", "M2"))
            dic[">" + eachType[0] + "_M1"] = dic[">" + eachType[0]]
            dic[">" + eachType[0] + "_M2"] = dic[">" + eachType[0]]
            dic.pop(">" + eachType[0])
        elif eachType[1] == "PB1":
            predictedProteinTypeNew.append((eachType[0] + "_PB1", "PB1"))
            predictedProteinTypeNew.append((eachType[0] + "_PB1-F2", "PB1-F2"))
            dic[">" + eachType[0] + "_PB1"] = dic[">" + eachType[0]]
            dic[">" + eachType[0] + "_PB1-F2"] = dic[">" + eachType[0]]
            dic.pop(">" + eachType[0])
        elif eachType[1] == "PA":
            predictedProteinTypeNew.append((eachType[0] + "_PA", "PA"))
            predictedProteinTypeNew.append((eachType[0] + "_PA-X", "PA-X"))
            dic[">" + eachType[0] + "_PA"] = dic[">" + eachType[0]]
            dic[">" + eachType[0] + "_PA-X"] = dic[">" + eachType[0]]
            dic.pop(">" + eachType[0])
        else:
            predictedProteinTypeNew.append(eachType)
    return predictedProteinTypeNew
