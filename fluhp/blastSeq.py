# !/usr/bin/python
# -*- coding: utf-8 -*-
# Author:lihuiru

import os
import subprocess
import py7zr

def blastSeq(querySeqFile, seqType = "protein", DBDir = "", dataBaseName = "", eValue = "1e-5", outfmt = "6",
             outName = "", querySeqFileDir = "", outFileDir = ""):
    print("blastSeq ...")
    if os.path.exists(DBDir + dataBaseName + ".7z") and not os.path.exists(DBDir + dataBaseName + ".dmnd"):
        with py7zr.SevenZipFile(DBDir + dataBaseName + ".7z", mode = 'r') as z:
            z.extractall(DBDir + "/")
    cmd = None
    if seqType == "protein":
        diamondType = "blastp"
        cmd = (
            f"diamond {diamondType} -d {DBDir}{dataBaseName} -q {querySeqFileDir}{querySeqFile} -o {outFileDir}{outName} "
            f"-e {eValue} -p 1 --outfmt {outfmt} -k 5"
        )
    elif seqType == "nucleo":
        diamondType = "blastn"
        cmd = (
            f"{diamondType} -db {DBDir}{dataBaseName} -query {querySeqFileDir}{querySeqFile} -out {outFileDir}{outName} "
            f"-evalue {eValue} -num_threads 1 -outfmt {outfmt} -max_target_seqs 5"
        )
    subprocess.getoutput(cmd)

