#!/usr/bin/python3
# -*- coding: utf-8 -*-

import mysql.connector
import os
import shlex
import subprocess

"""
CONSTANTES
"""
pop_keys = ["1000g2015aug_all", "1000g2015aug_afr", "1000g2015aug_amr", "1000g2015aug_eas", "1000g2015aug_eur", "1000g2015aug_sas", "ExAC_ALL", "ExAC_AFR", "ExAC_AMR",
"ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "AF", "AF_popmax", "AF_male", "AF_female", "AF_raw", "AF_afr", "AF_sas", "AF_amr", "AF_eas",
"AF_nfe", "AF_fin", "AF_asj", "AF_oth", "non_topmed_AF_popmax", "non_neuro_AF_popmax", "non_cancer_AF_popmax", "controls_AF_popmax", "AF", "AF_popmax",
"AF_male", "AF_female", "AF_raw", "AF_afr", "AF_sas", "AF_amr", "AF_eas", "AF_nfe", "AF_fin", "AF_asj", "AF_oth", "non_topmed_AF_popmax", "non_neuro_AF_popmax",
"non_cancer_AF_popmax", "controls_AF_popmax", "esp6500siv2_all", "esp6500siv2_ea", "esp6500siv2_aa"]
pred_keys = ["SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "PROVEAN_pred",
"MetaSVM_pred", "MetaLR_pred", "DANN_score"]

def summaryPredictors(keys, values) :
    deleterious = 0
    tolerated = 0
    unknown = 0
    temp = {}
    it = 0
    idx = keys.index("SIFT_pred")
    if values[idx] == 'D' :
        deleterious += 1
    elif values[idx] == 'T' :
        tolerated += 1
    else :
        unknown += 1
    idx = keys.index("Polyphen2_HDIV_pred")
    if values[idx] == 'D' or values[idx] == 'P' :
        deleterious += 1
    elif values[idx] == 'B' :
        tolerated += 1
    else :
        unknown += 1
    idx = keys.index("Polyphen2_HVAR_pred")
    if values[idx] == 'D' or values[idx] == 'P' :
        deleterious += 1
    elif values[idx] == 'B' :
        tolerated += 1
    else :
        unknown += 1
    idx = keys.index("LRT_pred")
    if values[idx] == 'D' :
        deleterious += 1
    elif values[idx] == 'N' :
        tolerated += 1
    else :
        unknown += 1
    idx = keys.index("MutationTaster_pred")
    if values[idx] == 'A' or values[idx] == 'D' :
        deleterious += 1
    elif values[idx] == 'N' or values[idx] == 'P' :
        tolerated += 1
    else :
        unknown += 1
    idx = keys.index("MutationAssessor_pred")
    if values[idx] == 'H' or values[idx] == 'M' :
        deleterious += 1
    elif values[idx] == 'L' or values[idx] == 'N' :
        tolerated += 1
    else :
        unknown += 1
    idx = keys.index("FATHMM_pred")
    if values[idx] == 'D' :
        deleterious += 1
    elif values[idx] == 'T' :
        tolerated += 1
    else :
        unknown += 1
    idx = keys.index("PROVEAN_pred")
    if values[idx] == 'D' :
        deleterious += 1
    elif values[idx] == 'T' :
        tolerated += 1
    else :
        unknown += 1
    idx = keys.index("MetaSVM_pred")
    if values[idx] == 'D' :
        deleterious += 1
    elif values[idx] == 'T' :
        tolerated += 1
    else :
        unknown += 1
    idx = keys.index("MetaLR_pred")
    if values[idx] == 'D' :
        deleterious += 1
    elif values[idx] == 'T' :
        tolerated += 1
    else :
        unknown += 1
    idx = keys.index("DANN_score")
    if values[idx] == "NA" or values[idx] == "." :
        unknown += 1
    elif float(values[idx]) > 0.9 :
        deleterious += 1
    elif float(values[idx]) <= 0.9 :
        tolerated += 1

    return "{}D, {}T, {}U".format(deleterious, tolerated, unknown)

def getPopMax(keys, values) :
    """Buscar la poblacion que tiene la MAF mayor y devolver el valor, junto con la poblacion reportada"""
    pop = ""
    maf = -1
    temp = {}
    it = 0
    for k in keys :
        temp[k] = values[it]
        it += 1
    for k,v in temp.items() :
        if v != "NA" and v != ".":
            val = float(v)
            if val > maf :
                maf = val
                pop = k
    if maf == -1 :
        maf = "NA"

    return pop, maf

def filtrarVariante(type, exoType, maf, vaf) :
    filtro = "raw_candidate"
    # Comprobar si la variante es codificante
    if type == "exonic" or type == "splicing" :
        if exoType != "synonymous SNV" :
            filtro = "consequence"
    # Comprobar si tiene una MAF alta y si tiene una VAF baja
    if filtro == "consequence" :
        if maf == "NA" :
            auxMaf = -1
        else :
            auxMaf = float(maf)
        if auxMaf >= 0.01 :
            filtro = "population_SNP"
        else :
            try :
                if float(vaf) < 5 :
                    filtro = "low_VAF"
                else :
                    filtro = "candidate_variant"
            except ValueError :
                print("Invalid VAF: {}".format(vaf))
    return filtro

def anotarVariante(varFile) :
    """Anota el archivo vcf pasado por parametro usando ANNOVAR. Crea dos archivos: raw.av y raw.hg19_mutianno.txt. Los archivos se guardan en la misma carpeta donde esta el vcf
    Devuelve el nombre del archivo anotado"""
    vDir = os.path.dirname(varFile)
    av = "{}/raw.av".format(vDir)
    txt = "{}/raw.hg19_multianno.txt".format(vDir)
    conv = "/opt/annovar20200607/convert2annovar.pl -format vcf4 -outfile {input} -includeinfo {fic}".format(input = av, fic = varFile)
    anno = "/opt/annovar20200607/table_annovar.pl {input} /home/ffuster/share/biodata/Indexes/ANNOVAR/humandb -buildver hg19 -out {dir}/raw -remove --protocol refGene,avsnp150,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,gnomad211_exome,gnomad211_genome,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,clinvar_20190305,cosmic70,dbnsfp35a --operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --nastring NA --otherinfo".format(dir = vDir, input = av)
    if not os.path.isfile(av) :
        args = shlex.split(conv)
        p = subprocess.Popen(args, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
        out, err = p.communicate()
        if p.returncode != 0 :
            print("ERROR: While executing ANNOVAR. Check below possible errors. Command executed:\n{}".format(anno))
            sys.exit(1)
    if not os.path.isfile(txt) :
        args = shlex.split(anno)
        p = subprocess.Popen(args, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
        out, err = p.communicate()
        if p.returncode != 0 :
            print("ERROR: While executing ANNOVAR. Check below possible errors. Command executed:\n{}".format(anno))
            sys.exit(1)

    return txt

def getCommonData(aux, header) :
    var = {}
    auxChr = aux[header.index("Chr")]
    if auxChr.startswith("chr") :
        var["chr"] = auxChr
    else :
        var["chr"] =  "chr{}".format(auxChr)
    var["start"] = aux[header.index("Start")]
    var["end"] = aux[header.index("End")]
    var["ref"] = aux[header.index("Ref")]
    var["alt"] = aux[header.index("Alt")]
    var["refGen"] = "hg19"
    var["gene"] = aux[header.index("Gene.refGene")]
    var["type"] = aux[header.index("Func.refGene")]
    var["exoType"] = aux[header.index("ExonicFunc.refGene")]
    temp = aux[header.index("AAChange.refGene")]
    # Recoger la informacion de HGVS para cDNA y proteina
    if temp != "NA" :
        # ANNOVAR reporta la informacion sobre los distintos transcritos separados por comas. Recoger el primero, suponiendo que es el canonico
        transcritos = temp.split(",")[0]
        try :
            genAux, transAux, exon, cdna, prot = transcritos.split(":")
        except ValueError :
            genAux = "NA"
            transAux = "NA"
            exon = "NA"
            cdna = "NA"
            prot = "NA"
    else :
        genAux = transAux = exon = cdna = prot = "NA"
    var["transcript"] = transAux
    var["exon"] = exon
    var["cdna"] = cdna
    var["protein"] = prot
    var["dbsnp"] = aux[header.index("avsnp150")]
    var["clinvar"] = aux[header.index("CLNSIG")]
    var["cosmic"] = aux[header.index("cosmic70")]
    if var["cosmic"] != "NA" :
        temp = var["cosmic"].split(";")[0]
        var["cosmic"] = temp.split(",")[0].strip("ID=") # Recoger el primer identificador de COSMIC que reporta ANNOVAR

    # Guardar el resumen de predictores
    pred_vals = []
    for p in pred_keys :
        pred_vals.append(aux[header.index(p)])
    var["predictors"] = summaryPredictors(pred_keys, pred_vals)

    # Guardar las MAFs en una lista para sacar la mayor
    pop_vals = []
    for p in pop_keys :
        pop_vals.append(aux[header.index(p)])

    var["popMax"], var["maxMaf"] = getPopMax(pop_keys, pop_vals)

    # Como las columnas del vcf no estan anotadas el filtro0.hg19_mutianno.txt, se usa el valor directamente
    var["filter"] = aux[144]

    return var

def varscan2db(line, header) :
    aux = line.strip().split("\t")
    var = getCommonData(aux, header)

    format = aux[146].split(":")
    vals = aux[147].split(":")

    # Comprobar que los numeros de columna son correctos
    if "RD" in format and "AD" in format and "RDF" in format and "RDR" in format and "ADF" in format and "ADR" in format and "FREQ" in format :
        var["covRef"] = vals[format.index("RD")]
        var["covAlt"] = vals[format.index("AD")]
        var["refFw"] = vals[format.index("RDF")]
        var["refRv"] = vals[format.index("RDR")]
        var["altFw"] = vals[format.index("ADF")]
        var["altRv"] = vals[format.index("ADR")]
        var["vaf"] = vals[format.index("FREQ")].replace("%", "").replace(",", ".")
        if var["filter"] == "PASS" :
            var["filter"] = filtrarVariante(var["type"], var["exoType"], var["maxMaf"], var["vaf"])
    else :
        raise ValueError("ERROR: Mandatory columns not found in FORMAT VCF data. Error column: {}".format(line))

    return var

def mutect2db(line, header) :
    aux = line.strip().split("\t")
    var = getCommonData(aux, header)
    format = aux[146].split(":")
    vals = aux[147].split(":")

    # PARCHE! Los filtros de Mutect2 son acumulativos y no creo que haga falta ponerlos todos. Como el maximo de
    if len(var["filter"]) > 60 :
        aux = var["filter"].split(";")
        var["filter"] = ";".join(aux[0:3])

    if "SB" in format and "AD" in format :
        lst = vals[format.index("AD")].split(",")
        var["covRef"] = int(lst[0])
        var["covAlt"] = int(lst[1])
        lst = vals[format.index("SB")].split(",")
        var["refFw"] = int(lst[0])
        var["refRv"] = int(lst[1])
        var["altFw"] = int(lst[2])
        var["altRv"] = int(lst[3])
        var["vaf"] = round(100*(float(var["covAlt"])/float(var["covRef"] + var["covAlt"])), 3)
    else :
        raise ValueError("ERROR: Mandatory columns (SB, AD) not found in VCF. Data found {} - {}".format(format, vals))

    return var

def saveInDB(variant, dbName) :
    # Comprobar si la variante se ha guardado previamente en la base de datos. Recoger el identificador de la variante
    dbcon = mysql.connector.connect(host="localhost", user="ffuster", password="Aetaeb6e", database=dbName)
    with dbcon as con :
        query = "SELECT id FROM variant WHERE cromosoma='{chr}' AND inicio={sta} AND observado='{alt}' AND genoma_ref='{gref}'".format(
        chr = variant["chr"], sta = variant["start"], alt = variant["alt"], gref = variant["refGen"])
        qvar = "INSERT INTO variant(cromosoma,inicio,fin,referencia,observado,genoma_ref,gen,tipo_var,tipo_ex,transcrito,hgvs_cDNA,hgvs_prot,exon,dbsnp,clinvar,cosmic,sum_pred,maf,pop_maf) "
        with con.cursor() as cur :
            cur.execute(query)
            res = cur.fetchone()
            # La variante no esta guardada en la base de datos
            if res != None :
                variant["id_variant"] = res[0]
            else :
                # Guardar los datos de una nueva variante
                if variant["maxMaf"] == "NA" :
                    variant["maxMaf"] = "NULL"
                qvar += "VALUES('{chr}',{sta},{end},'{ref}','{alt}','{gRF}','{gen}','{typ}','{eTP}','{trans}','{cdna}','{prot}','{exon}','{dbsnp}','{clnvar}','{cosmic}','{preds}',{maf},'{popMaf}')".format(
                chr = variant["chr"], sta = variant["start"], end = variant["end"], ref = variant["ref"], alt = variant["alt"], gRF = variant["refGen"], gen = variant["gene"],
                typ = variant["type"], eTP = variant["exoType"], trans = variant["transcript"], cdna = variant["cdna"], prot = variant["protein"], exon = variant["exon"], dbsnp = variant["dbsnp"],
                clnvar = variant["clinvar"], cosmic = variant["cosmic"], preds = variant["predictors"], maf = variant["maxMaf"], popMaf = variant["popMax"])
                try :
                    cur.execute(qvar)
                    con.commit()
                    variant["id_variant"] = cur.lastrowid
                except mysql.connector.Error as err:
                    print("ERROR: Executing \n{}\nDescription".format(qvar))
                    print("{} - {}".format(err.errno, err.msg))
                    cur.execute("DELETE FROM mostra WHERE id='{}'".format(variant["id_mostra"]))
                    con.commit()
                    sys.exit(1)

            # Guardar los datos del run
            qrun = "INSERT INTO run(id_variant,id_mostra,coverage,cov_ref,cov_alt,reads_FW_ref,reads_RV_ref,reads_FW_alt,reads_RV_alt,vaf,filtro) VALUES("
            qrun += "{idv},'{idm}',{cov},{cref},{calt},{rFr},{rRr},{rFa},{rRa},{vaf},'{filt}')".format(idv = variant["id_variant"], idm = variant["id_mostra"],
                cov = int(variant["covRef"]) + int(variant["covAlt"]), cref = variant["covRef"], calt = variant["covAlt"], rFr = variant["refFw"], rRr = variant["refRv"],
                rFa = variant["altFw"], rRa = variant["altRv"], vaf = variant["vaf"], filt = variant["filter"])
            try :
                cur.execute(qrun)
                con.commit()
            except mysql.connector.Error as err:
                print("ERROR: Executing \n{}\nDescription".format(qrun))
                print("{} - {}".format(err.errno, err.msg))
                cur.execute("DELETE FROM mostra WHERE id='{}'".format(variant["id_mostra"]))
                con.commit()
                sys.exit(1)
