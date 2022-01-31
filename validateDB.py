#!/usr/bin/python3
# -*- coding: utf-8 -*-

import mysql.connector
import os
import shlex
import subprocess
import sys

"""
Validar que las variantes introducidas en la base de datos son las mismas que hay en los archivos originales
"""
path = input("INPUT: Folder you want to validate: ")
if os.path.isdir(path) :
    # Buscar cuantas muestras hay en la carpeta que se ha introducido
    cmd = "find {} -name raw.hg19_multianno.txt".format(path)
    args = shlex.split(cmd)
    p = subprocess.Popen(args, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
    out, err = p.communicate()
    arx = out.decode().strip().split("\n")
    if arx[0] == "" :
        print("ERROR: No samples found")
        sys.exit()
    else :
        print("INFO: {} samples found".format(len(arx)))
        # Contar el numero de variantes que tiene cada gen en cada muestra
        # Comprobar si en la base de datos los numeros son iguales. Imprimir si hay algo que no cuadra
        for a in arx :
            # Contar el numero total de variantes que tiene la muestra
            cmd = "wc -l {}".format(a)
            args = shlex.split(cmd)
            p = subprocess.Popen(args, stdout = subprocess.PIPE)
            out, err = p.communicate()
            # Leer el nombre de la muestra
            muestra = arx.split("/")[-3]
            print("INFO: Getting data from {}".format(muestra))
            total = out.decode().strip().split(" ")[0]
            # Contar el numero de genes distintos con variante que tiene la muestra
            genes = {}
            cmd = "cut -f 7 {}".format(a)
            args = shlex.split(cmd)
            p = subprocess.Popen(args, stdout = subprocess.PIPE)
            out, err = p.communicate()
            for b in out.decode().strip().split("\n") :
                if b in genes :
                    genes[b] += 1
                else :
                    genes[b] = 1
            print(genes)

else :
    print("ERROR: Cannot find {} folder".format(path))
    sys.exit(1)
