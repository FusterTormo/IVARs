# Crear una base de datos MySQL de variantes

Lee, filtra y crea una base de datos con los datos de los paneles analizados en los grupos [MDS](https://www.carrerasresearch.org/en/Myelodysplastic_Syndromes) y [ALL](https://www.carrerasresearch.org/en/Acute_Lymphoblastic_Leukemia_(ALL)), del [Instituto de investigaci&oacute;n contra la leucemia Josep Carreras](https://www.carrerasresearch.org).

Esta base de datos estar&aacute; accesible solo para los miembros de dichos grupos.

## Estructura de la base de datos

### Tabla variante

* <u>id</u>
* cromosoma
* inicio
* final
* referencia
* observado
* tipos
* tipo_exonico
* HGVS_cDNA
* HGVS_proteina
* genoma_referencia
* dbsnp (actualizable)
* clinvar (actualizable) (por ahora solo guardar la significancia, porque ANNOVAR no anota bien el identificador de variante)
* cosmic (actualizable)
* resumen_predictores (actualizable)
* max_MAF (actualizable)
* pop_max_MAF (actualizable)
* anotaciones (log con los cambios que ha recibido cada una de las columnas actualizables, entre otras anotaciones)

### Tabla run

* <u>id_variant</u> (autonumerico de la variante)
* <u>idmostra</u>
* coverage
* reads_referencia
* reads_observado
* reads_FW_referencia
* reads_FW_observado
* reads_RV_referencia
* reads_RV_observado
* VAF
* filtro_variant_caller (o manual)
* info_del_run (log con comentarios que puede escribir el usuario)

### Tabla muestra

* <u>id</u> (identificador de la muestra)
* [...] (se puede poner mas datos segun quieran los usuarios)
