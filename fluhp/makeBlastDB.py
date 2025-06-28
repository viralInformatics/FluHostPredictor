def makedb(dir):
    import  subprocess
    fileName = subprocess.getoutput("cd "+str(dir)+'\n'+'ls')
    print(fileName)
    fileName =fileName.split("\n")
    for each in fileName:
        print(each)
    a = subprocess.getoutput('ls')
    # print(a)
    return




makedb("/home/think/Downloads")