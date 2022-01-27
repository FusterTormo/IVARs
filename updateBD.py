#!/usr/bin/python3
# -*- coding: utf-8 -*-

import myvariant

def actualizarVariante(chr, start, end, ref, alt, refGenome) :
    snv = ['A', 'C', 'G', 'T']
    campos = ["dbsnp.rsid", "cosmic.cosmic_id", "clinvar", "dbnsfp"]
    if len(ref) == 1 and len(alt) == 1 and ref in snv and alt in snv :
        if chr.startswith("chr") :
            hgvs = "{chr}:g.{pos}{ref}>{alt}".format(chr = chr, pos = start, ref = ref, alt = alt)
        else :
            hgvs = "chr{chr}:g.{pos}{ref}>{alt}".format(chr = chr, pos = start, ref = ref, alt = alt)

        cnt = myvariant.MyVariantInfo()
        json = cnt.getvariant(hgvs, fields=campos)
        print(json)
# dbsnp, sgClinvar, cosmic, predResum, maxMaf, popMax, anotacions = bd.actualizarVariante(chr, start, end, ref, alt, refGenome)
