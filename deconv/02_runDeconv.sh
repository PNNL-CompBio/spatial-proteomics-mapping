#!/bin/bash

cwltool https://raw.githubusercontent.com/PNNL-CompBio/decomprolute/main/run-deconv-on-file.cwl --data spatialProtMat.tsv --protAlg cibersort --signature LM7c

cwltool https://raw.githubusercontent.com/PNNL-CompBio/decomprolute/main/run-deconv-on-file.cwl --data spatialProtMat.tsv --protAlg xcell --signature LM7c

cwltool https://raw.githubusercontent.com/PNNL-CompBio/decomprolute/main/run-deconv-on-file.cwl --data spatialProtMat.tsv --protAlg mcpcounter --signature LM7c

cwltool https://raw.githubusercontent.com/PNNL-CompBio/decomprolute/main/run-deconv-on-file.cwl --data globalProtMat.tsv --protAlg cibersort --signature LM7c

cwltool https://raw.githubusercontent.com/PNNL-CompBio/decomprolute/main/run-deconv-on-file.cwl --data globalProtMat.tsv --protAlg xcell --signature LM7c

cwltool https://raw.githubusercontent.com/PNNL-CompBio/decomprolute/main/run-deconv-on-file.cwl --data globalProtMat.tsv --protAlg mcpcounter --signature LM7c
