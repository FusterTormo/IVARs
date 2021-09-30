#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import shlex
import subprocess
import sys

def findVariantFiles(folder) :
    if os.path.isdir(folder) :
        filename = input("INPUT: Tell me the raw annotated variants name (possible solutions: filtro0.allInfo, raw.reanno.tsv, filtro0.hg19_mutianno.txt) ")
        cmd = "find {} -name {}".format(folder, filename)
        args = shlex.split(cmd)
        p = subprocess.Popen(args, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
        out, err = p.communicate()
        files = out.decode().strip().split("\n")
        if len(files) > 0 :
            for f in files :
                with open(f, "r") as fi :
                    header = []
                    tmp = {}
                    for l in fi :
                        aux = l.strip().split("\t")
                        if len(header) == 0 :
                            header = aux
                            # Buscar si todas las columnas para la base de datos tiene el nombre esperado
                            if "Chr" not in header :
                                print("ERROR: Chromosome column not found")
                        else :
                            chr = aux[header.index("Chr")]
        else :
            print("ERROR: No files found with the name {} in {}".format(filename, folder))
    else :
        print("ERROR: {} folder not found".format(folder))


if __name__ == "__main__" :
    if len(sys.argv) < 2 :
        folder = input("INPUT: Where can I find the NGS run? ")
        findVariantFiles(folder)
    else :
        findVariantFiles(sys.argv[1])
