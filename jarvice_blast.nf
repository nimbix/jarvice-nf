process stage1 {
    executor 'jarvice'

    container "us-docker.pkg.dev/jarvice/images/nf-blast:latest"
    debug true
    machineType 'n1'
    cpus 4

    input:
    path infile

    output:
    path "query.out"

    """
    #!/bin/bash
    echo Making DB
    makeblastdb -in $infile -out $infile-db -dbtype prot -title "Test" -parse_seqids
    echo

    echo Querying sequences
    time blastp -query $infile -db $infile-db -max_target_seqs 5 -max_hsps 1 -evalue 1e-6 -outfmt '7 qseqid sseqid length qlen slen qstart qend sstart send evalue' -out ./query.out -num_threads $cpus
    echo
    echo  ============================================================
    """
}


workflow {
  Channel.fromPath( '/data/blast/*.fasta' ) | stage1
}
