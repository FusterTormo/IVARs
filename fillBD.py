#!/usr/bin/python3
# -*- coding: utf-8 -*-

import mysql.connector
import os
import shlex
import shutil
import subprocess
import sys

import updateBD as bd

"""
Columnes a guardar en la base de dades

-------------
Taula variant
-------------
cromosoma
inici
final
referencia
observat
tipus
tipus_exonic
HGVS_cDNA
HGVS_proteina
genoma_referencia
dbsnp (actualitzable)
clinvar (actualitzable) (per ara nomes guardar la significancia, perque annovar no anota be l'identificador)
cosmic (actualitzable)
resum_predictors (actualitzable)
max_MAF (actualitzable)
pop_max_MAF (actualitzable)
anotacions (log amb els canvis que han rebut cadascuna de les columnes actualitzables i altres anotacions)

----------
Taula run
----------
id_variant (autonumeric de la variant)
idmostra
Coverage
reads_referencia
reads_observat
reads_FW_referencia (?)
reads_FW_observat (?)
reads_RV_referencia (?)
reads_RV_observat (?)
VAF
filtre_variant_caller (o manual)
info_del_run (log amb comentaris que puga escriure l'usuari)

----------
Taula mostra
----------
id (identificador de la mostra)
[...] (es podrien afegir mes dades segons vulguen els usuaris)
"""

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
    vDir = os.path.dirname(varFile)
    conv = "/opt/annovar20200607/convert2annovar.pl -format vcf4 -outfile {dir}/raw.av -includeinfo {fic}".format(dir = vDir, fic = varFile)
    anno = "/opt/annovar20200607/table_annovar.pl raw.av /home/ffuster/share/biodata/Indexes/ANNOVAR/humandb -buildver hg19 -out raw -remove --protocol refGene,avsnp150,1000g2015aug_all,\
    1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,gnomad211_exome,gnomad211_genome,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,clinvar_20190305,\
    cosmic70,dbnsfp35a --operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --nastring NA --otherinfo"
    args = shlex.split(conv)
    p = subprocess.Popen(args, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0 :
        print(err)
        print(conv)
    sys.exit()
    shutil.move("raw.hg19_mutianno.txt", "{}/raw.hg19_mutianno.txt".format(vDir))

def fillALLdb(filename) :
    """
    # Columnas si el archivo es tipo filtro0.hg19_mutianno.txt
    Chr	Start	End	Ref	Alt	Func.refGene	Gene.refGene	GeneDetail.refGene	ExonicFunc.refGene	AAChange.refGene	avsnp150	1000g2015aug_all	1000g2015aug_afr	1000g2015aug_amr	1000g2015aug_eas	1000g2015aug_eur	1000g2015aug_sas	ExAC_ALL	ExAC_AFR	ExAC_AMR	ExAC_EAS	ExAC_FIN	ExAC_NFE	ExAC_OTH	ExAC_SAS	AF	AF_popmax	AF_male	AF_female	AF_raw	AF_afr	AF_sas	AF_amr	AF_eas	AF_nfe	AF_fin	AF_asj	AF_oth	non_topmed_AF_popmax	non_neuro_AF_popmax	non_cancer_AF_popmax	controls_AF_popmax	AF	AF_popmax	AF_male	AF_female	AF_raw	AF_afr	AF_sas	AF_amr	AF_eas	AF_nfe	AF_fin	AF_asj	AF_oth	non_topmed_AF_popmax	non_neuro_AF_popmax	non_cancer_AF_popmaxcontrols_AF_popmax	esp6500siv2_all	esp6500siv2_ea	esp6500siv2_aa	CLNALLELEID	CLNDN	CLNDISDB	CLNREVSTAT	CLNSIG	cosmic70	SIFT_score	SIFT_converted_rankscore	SIFT_pred	Polyphen2_HDIV_score	Polyphen2_HDIV_rankscore	Polyphen2_HDIV_pred	Polyphen2_HVAR_score	Polyphen2_HVAR_rankscore	Polyphen2_HVAR_pred	LRT_score	LRT_converted_rankscore	LRT_pred	MutationTaster_score	MutationTaster_converted_rankscore	MutationTaster_pred	MutationAssessor_score	MutationAssessor_score_rankscore	MutationAssessor_pred	FATHMM_score	FATHMM_converted_rankscore	FATHMM_pred	PROVEAN_score	PROVEAN_converted_rankscore	PROVEAN_pred	VEST3_score	VEST3_rankscore	MetaSVM_score	MetaSVM_rankscore	MetaSVM_pred	MetaLR_score	MetaLR_rankscore	MetaLR_pred	M-CAP_score	M-CAP_rankscore	M-CAP_pred	REVEL_score	REVEL_rankscore	MutPred_score	MutPred_rankscore	CADD_raw	CADD_raw_rankscore	CADD_phred	DANN_score	DANN_rankscore	fathmm-MKL_coding_score	fathmm-MKL_coding_rankscore	fathmm-MKL_coding_pred	Eigen_coding_or_noncoding	Eigen-raw	Eigen-PC-raw	GenoCanyon_score	GenoCanyon_score_rankscore	integrated_fitCons_score	integrated_fitCons_score_rankscore	integrated_confidence_value	GERP++_RS	GERP++_RS_rankscore	phyloP100way_vertebrate	phyloP100way_vertebrate_rankscore	phyloP20way_mammalian	phyloP20way_mammalian_rankscore	phastCons100way_vertebrate	phastCons100way_vertebrate_rankscore	phastCons20way_mammalian	phastCons20way_mammalian_rankscore	SiPhy_29way_logOdds	SiPhy_29way_logOdds_rankscore	Interpro_domain	GTEx_V6p_gene	GTEx_V6p_tissue	Otherinfo
    """
    pop_keys = ["1000g2015aug_all", "1000g2015aug_afr", "1000g2015aug_amr", "1000g2015aug_eas", "1000g2015aug_eur", "1000g2015aug_sas", "ExAC_ALL", "ExAC_AFR", "ExAC_AMR",
    "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "AF", "AF_popmax", "AF_male", "AF_female", "AF_raw", "AF_afr", "AF_sas", "AF_amr", "AF_eas",
    "AF_nfe", "AF_fin", "AF_asj", "AF_oth", "non_topmed_AF_popmax", "non_neuro_AF_popmax", "non_cancer_AF_popmax", "controls_AF_popmax", "AF", "AF_popmax",
    "AF_male", "AF_female", "AF_raw", "AF_afr", "AF_sas", "AF_amr", "AF_eas", "AF_nfe", "AF_fin", "AF_asj", "AF_oth", "non_topmed_AF_popmax", "non_neuro_AF_popmax",
    "non_cancer_AF_popmax", "controls_AF_popmax", "esp6500siv2_all", "esp6500siv2_ea", "esp6500siv2_aa"]
    pred_keys = ["SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "PROVEAN_pred",
    "MetaSVM_pred", "MetaLR_pred", "DANN_score"]
    header = []
    mostra = ""
    files = os.listdir(os.path.dirname(filename))
    backup = ""
    for f in files :
        if f.endswith(".xlsx") :
            mostra = f.replace(".xlsx", "")
            break
    # Comprobar si la muestra ya se ha guardado en la base de datos previamente
    dbcon = mysql.connector.connect(host="localhost", user="ffuster", password="Aetaeb6e", database="ALLvar")
    with dbcon as con :
        query = "SELECT id FROM mostra WHERE id='{}'".format(mostra);
        with con.cursor() as cur :
            cur.execute(query)
            res = cur.fetchall()
    if len(res) > 0 :
        print("WARNING: Sample {} already stored in the database. Variants will not be stored".format(mostra))
    else :
        # Guardar la muestra en la base de datos
        dbcon = mysql.connector.connect(host="localhost", user="ffuster", password="Aetaeb6e", database="ALLvar")
        with dbcon as con :
            query = "INSERT INTO mostra(id) VALUES('{}')".format(mostra);
            with con.cursor() as cur :
                cur.execute(query)
                con.commit()
                backup = backup + query + ";\n"

        with open(filename, "r") as fi :
            for l in fi :
                aux = l.strip().split("\t")
                if len(header) == 0 :
                    header = aux
                    # Buscar si todas las columnas para la base de datos tiene el nombre esperado
                    if "Chr" not in header :
                        raise ValueError("ERROR: Chromosome column not found")
                else :
                    # Comprobar si la muestra ya se ha guardado en la base de datos previamente
                    chr = aux[header.index("Chr")]
                    start = aux[header.index("Start")]
                    end = aux[header.index("End")]
                    ref = aux[header.index("Ref")]
                    alt = aux[header.index("Alt")]
                    refGenome = "hg19"
                    gen = aux[header.index("Gene.refGene")]
                    typ = aux[header.index("Func.refGene")]
                    exo = aux[header.index("ExonicFunc.refGene")]
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
                    dbsnp = aux[header.index("avsnp150")]
                    sgClinvar = aux[header.index("CLNSIG")]
                    cosmic = aux[header.index("cosmic70")]
                    if cosmic != "NA" :
                        temp = cosmic.split(";")[0]
                        cosmic = temp.split(",")[0].strip("ID=") # Recoger el primer identificador de COSMIC que reporta ANNOVAR

                    # Guardar el resumen de predictores
                    pred_vals = []
                    for p in pred_keys :
                        pred_vals.append(aux[header.index(p)])
                    sumPreds = summaryPredictors(pred_keys, pred_vals)

                    # Guardar las MAFs en una lista para sacar la mayor
                    pop_vals = []
                    for p in pop_keys :
                        pop_vals.append(aux[header.index(p)])

                    popMax, maxMaf = getPopMax(pop_keys, pop_vals)

                    # Como las columnas del vcf no estan anotadas el filtro0.hg19_mutianno.txt, se usa el valor directamente
                    filter = aux[144]
                    format = aux[146].split(":")
                    vals = aux[147].split(":")
                    # Comprobar que los numeros de columna son correctos
                    if format[0] == "GT" :
                        covRef = vals[format.index("RD")]
                        covAlt = vals[format.index("AD")]
                        refFw = vals[format.index("RDF")]
                        refRv = vals[format.index("RDR")]
                        altFw = vals[format.index("ADF")]
                        altRv = vals[format.index("ADR")]
                        vaf = vals[format.index("FREQ")].replace("%", "").replace(",", ".")
                        if filter == "PASS" :
                            filter = filtrarVariante(typ, exo, maxMaf, vaf)
                    else :
                        print("ERROR: Invalid column found when getting the FORMAT VCF data. Please, change that value in 'format = aux[146].split(...' line")
                        sys.exit()

                    # Comprobar si la variante existe en la base de datos. Guardar la variante en caso de que no exista
                    dbcon = mysql.connector.connect(host="localhost", user="ffuster", password="Aetaeb6e", database="ALLvar")
                    qvar = "INSERT INTO variant(cromosoma,inicio,fin,referencia,observado,genoma_ref,gen,tipo_var,tipo_ex,hgvs_cDNA,hgvs_prot,exon,dbsnp,clinvar,cosmic,sum_pred,maf,pop_maf) "
                    with dbcon as con :
                        query = "SELECT id FROM variant WHERE cromosoma='{chr}' AND inicio={sta} AND observado='{alt}' AND genoma_ref='{gref}'".format(chr = chr, sta = start, alt = alt, gref = refGenome)
                        with con.cursor() as cur :
                            cur.execute(query)
                            res = cur.fetchall()
                            if len(res) == 0:
                                if maxMaf == "NA" :
                                    maxMaf = "NULL"
                                qvar += "VALUES('{chr}',{sta},{end},'{ref}','{alt}','{genRef}','{gen}','{typ}','{exType}','{cdna}','{prot}','{exon}','{snp}','{cli}','{cosm}','{pred}',{maf},'{pop}')".format(
                                chr=chr, sta=start, end=end, ref=ref, alt=alt, genRef=refGenome, gen=gen, typ=typ, exType=exo, cdna=cdna, prot=prot, exon=exon, snp=dbsnp, cli=sgClinvar, cosm=cosmic, pred=sumPreds,
                                maf=maxMaf, pop=popMax)
                                try :
                                    cur.execute(qvar)
                                    con.commit()
                                except mysql.connector.Error as err:
                                    print("ERROR: ejecutando \n{}\nDescripcion".format(qvar))
                                    print("{} - {}".format(err.errno, err.msg))
                                    cur.execute("DELETE FROM mostra WHERE id='{}'".format(mostra))
                                    con.commit()
                                    sys.exit(1)
                                idVariant = cur.lastrowid
                                backup = backup + qvar + ";\n"
                            # Guardar el identificador de la variante para ponerlo en los datos de la tabla run
                            else :
                                idVariant = res[0][0]
                    # Guardar los datos del run en la base de datos
                    qrun = "INSERT INTO run(id_variant, id_mostra, coverage, cov_ref, cov_alt, reads_FW_ref, reads_FW_alt, reads_RV_ref, reads_RV_alt, vaf, filtro) "
                    dbcon = mysql.connector.connect(host="localhost", user="ffuster", password="Aetaeb6e", database="ALLvar")
                    with dbcon as con :
                        qrun += "VALUES({var},'{samp}',{cov},{rcov},{acov},{rfcov},{afcov},{rrcov},{arcov},{vaf},'{filt}')".format(
                        var=idVariant, samp=mostra, cov=int(covRef)+int(covAlt), rcov=covRef, acov=covAlt, rfcov=refFw, afcov=altFw, rrcov=refRv, arcov=altRv, vaf=vaf, filt=filter)
                        with con.cursor() as cur :
                            cur.execute(qrun)
                            con.commit()
                            backup = backup + qrun + ";\n"


    return backup


def fillVHdb(filename) :
    """
    # Columnas si el archivo es tipo raw.anno.tsv

    Chr	Start	End	Ref	Alt	Func.refGene	Gene.refGene	GeneDetail.refGene	ExonicFunc.refGene	AAChange.refGene	avsnp150	1000g2015aug_all	1000g2015aug_afr	1000g2015aug_amr	1000g2015aug_eas	1000g2015aug_eur	1000g2015aug_sas	ExAC_ALL	ExAC_AFR	ExAC_AMR	ExAC_EAS	ExAC_FIN	ExAC_NFE	ExAC_OTH	ExAC_SAS	gnomad_exome_AF	gnomad_exome_AF_popmax	gnomad_exome_AF_male	gnomad_exome_AF_female	gnomad_exome_AF_raw	gnomad_exome_AF_afr	gnomad_exome_AF_sas	gnomad_exome_AF_amr	gnomad_exome_AF_eas	gnomad_exome_AF_nfe	gnomad_exome_AF_fin	gnomad_exome_AF_asj	gnomad_exome_AF_oth	gnomad_exome_non_topmed_AF_popmax	gnomad_exome_non_neuro_AF_popmax	gnomad_exome_non_cancer_AF_popmax	gnomad_exome_controls_AF_popmax	gnomad_genome_AF	gnomad_genome_AF_popmax	gnomad_genome_AF_male	gnomad_genome_AF_female	gnomad_genome_AF_raw	gnomad_genome_AF_afr	gnomad_genome_AF_sas	gnomad_genome_AF_amr	gnomad_genome_AF_eas	gnomad_genome_AF_nfe	gnomad_genome_AF_fin	gnomad_genome_AF_asj	gnomad_genome_AF_oth	gnomad_genome_non_topmed_AF_popmax	gnomad_genome_non_neuro_AF_popmax	gnomad_genome_non_cancer_AF_popmax	gnomad_genome_controls_AF_popmax	esp6500siv2_all	esp6500siv2_ea	esp6500siv2_aa	CLNALLELEID	CLNDN	CLNDISDB	CLNREVSTAT	CLNSIG	cosmic70	SIFT_score	SIFT_converted_rankscore	SIFT_pred	Polyphen2_HDIV_score	Polyphen2_HDIV_rankscore	Polyphen2_HDIV_pred	Polyphen2_HVAR_score	Polyphen2_HVAR_rankscore	Polyphen2_HVAR_pred	LRT_score	LRT_converted_rankscore	LRT_pred	MutationTaster_score	MutationTaster_converted_rankscore	MutationTaster_pred	MutationAssessor_score	MutationAssessor_score_rankscore	MutationAssessor_pred	FATHMM_score	FATHMM_converted_rankscore	FATHMM_pred	PROVEAN_score	PROVEAN_converted_rankscore	PROVEAN_predVEST3_score	VEST3_rankscore	MetaSVM_score	MetaSVM_rankscore	MetaSVM_pred	MetaLR_score	MetaLR_rankscore	MetaLR_pred	M-CAP_score	M-CAP_rankscore	M-CAP_pred	REVEL_scoreREVEL_rankscore	MutPred_score	MutPred_rankscore	CADD_raw	CADD_raw_rankscore	CADD_phred	DANN_score	DANN_rankscore	fathmm-MKL_coding_score	fathmm-MKL_coding_rankscore	fathmm-MKL_coding_pred	Eigen_coding_or_noncoding	Eigen-raw	Eigen-PC-raw	GenoCanyon_score	GenoCanyon_score_rankscore	integrated_fitCons_score	integrated_fitCons_score_rankscore	integrated_confidence_value	GERP++_RS	GERP++_RS_rankscore	phyloP100way_vertebrate	phyloP100way_vertebrate_rankscore	phyloP20way_mammalian	phyloP20way_mammalian_rankscore	phastCons100way_vertebrate	phastCons100way_vertebrate_rankscore	phastCons20way_mammalian	phastCons20way_mammalian_rankscore	SiPhy_29way_logOdds	SiPhy_29way_logOdds_rankscore	Interpro_domain	GTEx_V6p_gene	GTEx_V6p_tissue	CHROM	POS	ID	REFERENCE	ALTERATED	QUAL	FILTER	INFO	FORMAT	SAMPLE	END	BLOCKAVG	SNVHPOL	CIGAR	RU	REFPREP	IDREP	MQ	GT	GQ	GQX	DP	DPF	MIN_DP	AD	ADF	ADR	FT	DPI	PL	PS	SB	population_max	population_max_name	predictor_summary	Strand_bias_score	Ref_depth	Alt_depth	VAF	IGV_link	sample	CADD_1000g_all	CADD_1000g_afr	CADD_1000g_amr	CADD_1000g_eur	CADD_1000g_eas	CADD_1000g_sas	dbNSFP_1000g_all	dbNSFP_1000g_afr	dbNSFP_1000g_amr	dbNSFP_1000g_eur	dbNSFP_1000g_eas	dbNSFP_1000g_sas	CADD_ESP6500_all	CADD_ESP6500_ea	CADD_ESP6500_aa	dbNSFP_esp6500_all	dbNSFP_esp6500_ea	dbNSFP_esp6500_aa	ExAC_ExAC_all	ExAC_ExAC_afr	ExAC_ExAC_amr	ExAC_ExAC_eas	ExAC_ExAC_fin	ExAC_ExAC_nfe	ExAC_ExAC_oth	ExAC_ExAC_sas	dbNSFP_ExAC_all	dbNSFP_ExAC_afr	dbNSFP_ExAC_amr	dbNSFP_ExAC_eas	dbNSFP_ExAC_fin	dbNSFP_ExAC_nfe	dbNSFP_ExAC_oth	dbNSFP_ExAC_sas	gNOMAD_Exome_all	gNOMAD_Exome_afr	gNOMAD_Exome_amr	gNOMAD_Exome_asj	gNOMAD_Exome_eas	gNOMAD_Exome_fin	gNOMAD_Exome_nfe	gNOMAD_Exome_oth	gNOMAD_Exome_popmax	gNOMAD_Exome_raw	gNOMAD_Exome_sas	gNOMAD_Genome_all	gNOMAD_Genome_afr	gNOMAD_Genome_amr	gNOMAD_Genome_asj	gNOMAD_Genome_eas	gNOMAD_Genome_fin	gNOMAD_Genome_nfe	gNOMAD_Genome_oth	gNOMAD_Genome_popmax	gNOMAD_Genome_raw	dbSNP_MAF

    """
    pass

def findMDSfiles(dir) :
    cmd = "find {} -name *vcf".format(dir)
    args = shlex.split(cmd)
    p = subprocess.Popen(args, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    out, err = p.communicate()
    arx = out.decode().strip().split("\n")
    for a in arx :
        if a.endswith("vcf") :
            vc = os.path.basename(a)
            if vc == "varscan.vcf" :
                print("Trobat un varscan")
            elif vc == "mutect.filtered.vcf" :
                ## TODO: Convertir el mutect.filtered.vcf en mutect.revised.vcf per evitar variants multiples (executar AUP/revisaVcf.py)
                newf = a.replace("mutect.filtered.vcf", "mutect.revised.vcf")
                # if not os.path.isfile(newf)
            else :
                print("Trobat format estrany".format(a))


def findVariantFiles(folder, db = None) :
    filename = "" # Nombre del archivo que se va a buscar. Depende de la base de datos
    files = [] # Lista con  los archivos encontrados por el comando find. Con estos archivos se rellenara la base de datos
    sqlcode = "" # Datos para guardar en la base de datos. Se guardaran en un archivo sql aparte como copia de seguridad
    if os.path.isdir(folder) :
        # Preguntar por la base de datos para saber que archivo de variantes hay que buscar
        if db == None :
            db = input("INPUT: Which database do you want to fill (ALL, MDS, VH)? ")

        if db == "ALL" :
            filename = "filtro0.hg19_multianno.txt"
        elif db == "VH" :
            filename = "raw.reanno.tsv"
        elif db == "MDS" :
            filename = findMDSfiles(folder)
            sys.exit(0)
        else :
            raise ValueError("ERROR: Unknown database")

        # Buscar los archivos de variantes con los que se rellenara la base de datos
        cmd = "find {} -name {}".format(folder, filename)
        args = shlex.split(cmd)
        p = subprocess.Popen(args, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
        out, err = p.communicate()
        files = out.decode().strip().split("\n")
        if len(files) > 0 and files[0] != "":
            print("INFO: {} files found with name = {}".format(len(files), filename))
            for f in files :
                sampId = f.split("/")[-3]
                print("INFO: Getting {} variants".format(sampId))
                with open(f, "r") as fi :
                    header = []
                    tmp = {}
                    if db == "ALL" :
                        sqlcode += fillALLdb(f)

            with open("sql/{}.sql".format(db), "a") as fi :
                fi.write(sqlcode)
                print("INFO: SQL backup file stored in sql/{}.sql".format(db))
        else :
            # No se ha encontrado ningun archivo con el nombre que se buscaba
            raise TypeError("ERROR: No files found with the name {} in {}".format(filename, folder))
    else :
        # No se ha encontrado la carpeta
        raise TypeError("ERROR: {} folder not found".format(folder))


if __name__ == "__main__" :
    if len(sys.argv) < 2 :
        folder = input("INPUT: Where can I find the NGS run? ")
        findVariantFiles(folder)
    elif len(sys.argv) == 2 :
        findVariantFiles(sys.argv[1])
    else :
        findVariantFiles(sys.argv[1], sys.argv[2])
