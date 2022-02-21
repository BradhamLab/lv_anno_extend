#!/bin/sh
# properties = {"type": "single", "rule": "format_rRNA", "local": false, "input": ["/projectnb/bradham/data/ReferenceSequences/wray-genome/rRNA.gff"], "output": ["/projectnb/bradham/data/ReferenceSequences/modified-wray-extended/rRNA.gff"], "wildcards": {}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 4, "cluster": {}}
cd /projectnb2/bradham/workflows/lvar-anno && \
/projectnb/bradham/dyh0110/.conda/envs/alignment/bin/python3.8 \
-m snakemake /projectnb/bradham/data/ReferenceSequences/modified-wray-extended/rRNA.gff --snakefile /projectnb2/bradham/workflows/lvar-anno/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /projectnb2/bradham/workflows/lvar-anno/.snakemake/tmp.53hpx4qr /projectnb/bradham/data/ReferenceSequences/wray-genome/rRNA.gff --latency-wait 30 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules format_rRNA --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/projectnb2/bradham/workflows/lvar-anno/.snakemake/tmp.53hpx4qr/4.jobfinished" || (touch "/projectnb2/bradham/workflows/lvar-anno/.snakemake/tmp.53hpx4qr/4.jobfailed"; exit 1)

