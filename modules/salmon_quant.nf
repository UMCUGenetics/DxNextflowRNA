process QUANT {

    input:
    path index
    tuple val(pair_id),rg_id, path(reads) 

    output:
    path (pair_id)

    script:
    def barcode = rg_id.split('_')[1]

    """
    salmon quant \\
        --index $index\\
        --libType=U  \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o $pair_id
    """
}

