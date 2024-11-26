#!/bin/bash

# Directorio que contiene los archivos
DIR="/home/alvaro/geometric_dis/pruebas/"

# Iterar sobre los archivos que coinciden con el patrón B6_mask_c4.###.fits
#~ for file in "$DIR"/B6_mask_c4.[0-9][0-9][0-9].fits; do
    #~ # Verificar si el archivo existe
    #~ if [[ -e "$file" ]]; then
        #~ # Extraer el número de tres dígitos
        #~ num=$(echo "$file" | grep -oP '(?<=B6_mask_c4\.)\d{3}(?=\.fits)')
        
        #~ # Formatear el número a cuatro dígitos
        #~ new_num=$(printf "%04d" "$num")
        
        #~ # Construir el nuevo nombre de archivo
        #~ new_file="${file/B6_mask_c4.${num}.fits/B6_mask_c4.${new_num}.fits}"
        
        #~ # Renombrar el archivo
        #~ mv "$file" "$new_file"
        #~ echo "Renombrado: $file -> $new_file"
    #~ fi
#~ done

#~ echo "Proceso de renombrado completado."


for file in "$DIR"/B6_mask_c4.[0-9][0-9][0-9].fits; do
    # Verificar si el archivo existe
    if [[ -e "$file" ]]; then
        # Extraer el número de tres dígitos
        num=$(echo "$file" | grep -oP '(?<=B6_image_c4\.)\d{3}(?=\.fits)')
        
        # Formatear el número a cuatro dígitos
        new_num=$(printf "%04d" "$num")
        
        # Construir el nuevo nombre de archivo
        new_file="${file/B6_mask_c4.${num}.fits/B6_mask_c4.${new_num}.fits}"
        
        # Renombrar el archivo
        mv "$file" "$new_file"
        echo "Renombrado: $file -> $new_file"
    fi
done

echo "Proceso de renombrado completado."
