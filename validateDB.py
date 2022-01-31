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
        for a in arx :
            # Contar el numero total de variantes que tiene la muestra
            cmd = "wc -l {}".format(a)
            args = shlex.split(cmd)
            p = subprocess.Popen(args, stdout = subprocess.PIPE)
            out, err = p.communicate()
            # Leer el nombre de la muestra
            muestra = a.split("/")[-3]
            print("INFO: Getting data from {}".format(muestra))
            total = int(out.decode().strip().split(" ")[0]) - 1 # El total de filas que tiene el archivo, menos la primera linea, que es la cabecera
            # Contar el numero de genes distintos con variante que tiene la muestra (en la clave del diccionario)
            # Contar el numero de variantes que tiene cada gen en cada muestra (el contenido del diccionario)
            genes = {}
            different = False
            cmd = "cut -f 7 {}".format(a)
            args = shlex.split(cmd)
            p = subprocess.Popen(args, stdout = subprocess.PIPE)
            out, err = p.communicate()
            for b in out.decode().strip().split("\n") :
                if b != "Gene.refGene" :
                    if b in genes :
                        genes[b] += 1
                    else :
                        genes[b] = 1
            # Comprobar si en la base de datos los numeros son iguales. Imprimir si hay algo que no cuadra
            print("INFO: Comparing the information obtained with the stored in the database")
            sql = "SELECT count(*) FROM run WHERE id_mostra='{}'".format(muestra)
            dbcon = mysql.connector.connect(host="localhost", user="ffuster", password="Aetaeb6e", database="MDSvar")
            with dbcon as con:
                with con.cursor() as cur :
                    cur.execute(sql)
                    res = cur.fetchone()[0]
                    if int(res) != total :
                        print("WARNING: Total variants stored are different: {} found in database, {} found in file".format(res, total))
                        different = True
                    for k,v in genes.items() :
                        sql = "SELECT count(*) FROM run WHERE id_mostra='{samp}' AND id_variant IN (SELECT id FROM variant WHERE gen='{gen}')".format(samp = muestra, gen = k)
                        cur.execute(sql)
                        res = cur.fetchone()[0]
                        if int(res) != v :
                            print("WARNING: Variants in gene {} are different. Found in database {}, in file {}".format(k, res, v))
                            different = True

                        if not different :
                            print("OK")


else :
    print("ERROR: Cannot find {} folder".format(path))
    sys.exit(1)
