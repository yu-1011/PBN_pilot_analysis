#hail dataproc
hailctl dataproc start ychen-wes \
--project 'pbn-analysis'  \
--vep GRCh38 --requester-pays-allow-all --requester-pays-allow-annotation-db \
--packages 'ipython<8.22',gnomad \
--region "us-central1" \
--no-off-heap-memory \
--subnet default  \
--max-idle=20m
--master-machine-type e2-standard-8 --worker-machine-type e2-standard-8 --num-workers 2 \


hailctl dataproc connect ychen-wes notebook
hailctl dataproc list
hailctl dataproc stop ychen-wes
