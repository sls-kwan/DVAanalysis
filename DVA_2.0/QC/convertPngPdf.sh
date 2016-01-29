#!/bin/bash

cd $1

pdfExt=".pdf"
pngExt=".png"
lsLook="*.png"
convert +append $1$lsLook $1ALL$pngExt
convert $1$lsLook $1ALL$pdfExt
