#!/bin/zsh




# Minimal example
python3 main.py -gfa xxx.gfa -gaf xxx.gaf -maternal HG00421:0 -paternal HG00422:1 -out ./result.tsv

echo "\n"

# Full example
python3 main.py -gfa xxx.gfa -gaf xxx.gaf -maternal HG00421:0 -paternal HG00422:1 -out ./result.tsv -ignore_sex True



