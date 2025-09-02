include { OUTRIDER      } from '../../../modules/local/outrider/main'


workflow GENE_EXON_OUTRIDER {
    take:
    ch_gene_counts // gene expression counts
    ch_exon_counts // exon expression counts
    ch_gtf

    main:

    reference_panels = Channel.from(params.outrider_refs
        .collectMany { level, treatments ->
            treatments
                .collectMany { treatment, panels ->
                    panels.collect { p ->
                        def meta = [level: level, treatment: treatment, panel:
                                        (p.id as String)]
                        tuple(meta, (p.path ?: p.glob))
                    }
                }
        }
    )
    .map { meta, path ->
        [meta, files("${path}/*", checkIfExists: true)]
    }


    // create a regex to match a treatment condition out of multiple possible options
    // example params.conditions and condition_patterns:
    //  conditions = ['chx', 'ctrl']
        // condition_patterns = [
            // chx  : ['pluschx','chx'],
            // ctrl : ['control','cntrl','ctrl']
        // ]
    // so any file containg pluschx or chx will return chx, etc.
    def CONDITION_REGEX = params.conditions.collectEntries { c ->
        def pats = (params.condition_patterns[c] ?: [])
        def quoted = pats.collect { java.util.regex.Pattern.quote(it) }
        // case-insensitive, match word borders so (for example) 'ctrl' wil not match 'control'
        def rx = java.util.regex.Pattern.compile("(?i)\\b(?:${quoted.join('|')})\\b")
        [(c): rx]
    }

    // check the file name and match that with the regex pattern
    def detect_condition = { String fname ->
        def hits = CONDITION_REGEX.findAll { cond, rx -> rx.matcher(fname).find() }
            .collect { it.key }
        if (hits.size() == 1) return hits[0]
        if (hits.isEmpty())   error "No condition match for file: ${fname} (patterns: ${condition_patterns})"
        error "Ambiguous condition for file: ${fname}; matches: ${hits.join(', ')}"
    }

    // Annotate meta field with treatment ----
    def add_treatment = { meta, file ->
        def name = (file instanceof java.nio.file.Path) ? file.getFileName().toString() : file.toString()
        def cond = detect_condition(name)
        return [[ *: meta, treatment: cond ], file]
    }


    ch_genes = ch_gene_counts.map(add_treatment)
        .map{meta, file -> [[*: meta, level: "gene"], file]}
    ch_exons = ch_exon_counts.map(add_treatment)
        .map{meta, file -> [[*: meta, level: "exon"], file]}


    ch_samples = ch_genes.concat(ch_exons)

    // channel: [meta, query, [ref_files]]
    // meta: [id: <sample_name> , level: exon|gene,  treatment: cntrl|..., panel: <name>]
    ch_outrider_inputs = ch_samples
        .combine(reference_panels)
        .filter { sample_meta, sample_file, ref_meta, ref_files ->
            sample_meta.level == ref_meta.level &&
            sample_meta.treatment == ref_meta.treatment
        }.map{ sample_meta, sample_file, ref_meta, ref_files ->
        [[*:sample_meta, panel: ref_meta.panel], sample_file, ref_files]
        }



    OUTRIDER(ch_outrider_inputs, ch_gtf)


    ch_versions = OUTRIDER.out.versions

    emit:
    tsv_full   = OUTRIDER.out.tsv_full ?: [:] // channel: [ val(meta), path(*_outrider_results.tsv) ]
    tsv_signif = OUTRIDER.out.tsv_signif ?: [:]// channel: [ val(meta), path(*_outrider_results.tsv) ]
    versions        = ch_versions ?: [:]



}
