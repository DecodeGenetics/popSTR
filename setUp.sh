#!/bin/bash
echo "Unzipping and relocating markerInfo files"
mkdir markerInfo
for i in {1..22}; do gunzip -c chr${i}markerInfo.gz > markerInfo/chr${i}markerInfo; rm chr${i}markerInfo.gz; done
mv longRepeats markerInfo/
mv defaultModel markerInfo/
echo "Unpacking kernel"
gunzip kernel.tar.gz
tar -xvf kernel.tar
echo "Unpacking panel"
gunzip panelMarkerInfo.tar.gz
tar -xvf panelMarkerInfo.tar