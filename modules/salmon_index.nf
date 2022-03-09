process INDEX {

    input:
    path fasta

    output:
    path 'index'

    """
    salmon index \\
        -t "${fasta}" \\
        -i index \\
    """
}

