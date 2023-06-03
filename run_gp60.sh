#! /bin/bash -l

module load ncbi-blast+/LATEST
module load perl/5.30.1-mt

perl /scicomp/home-pure/qxu5/CE_Mock/data-example/ceroot/bin/Crypto_gp60Typing/gp60Typer_json_encrypt_mod4.pl 

module unload ncbi-blast+/LATEST

