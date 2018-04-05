#!/bin/bash

(
cd out

for FILE in *.txt; do
  gnuplot <<- EOF
    set terminal png
    set output "${FILE}.png"
    plot "${FILE}" notitle pointtype 6
EOF
done

ffmpeg -r 60 -f image2 -s 1080x960 -i 'fout-%05d.txt.png' -pix_fmt yuv420p test.mp4

rm *.png
)
