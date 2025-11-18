#FUS intron/exon Neurolincs
./coverage_calculator.sh -g FUS_intron -d neurolincs -r FUS_exon_intron.bed -c 8 -m 64G -t 48:00:00

#FUS intron/exon NYGC 
./coverage_calculator.sh -g FUS_intron -d nygc -r FUS_exon_intron.bed -c 8 -m 96G -t 48:00:00

#FUS intron/exon ANSWER ALS 
./coverage_calculator.sh -g FUS_intron -d answerals -r FUS_exon_intron.bed -c 8 -m 96G -t 48:00:00