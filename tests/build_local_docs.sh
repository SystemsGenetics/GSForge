#!/usr/bin/env bash

OUTPUT_DIRECTORY="~/gsforge_local_docs"

echo "generating module rst containers..."
nbsite_generate_modules.py GSForge -d ./doc/Reference_Manual -n GSForge
echo "generating notebook rst containers..."
nbsite generate-rst --org SystemsGenetics --project-name GSForge
echo "building docs..."
nbsite build --what=html --output=$OUTPUT_DIRECTORY

