rule all:
    input:
        auspice_json = "auspice/sars-like-cov.json",

rule files:
    params:
        input_sequences = "data/sequences.fasta",
        input_metadata = "data/metadata.tsv",
        include = "config/include.txt",
        exclude = "config/exclude.txt",
        reference = "config/sars-like-cov_reference.gb",
        clades = "config/clades.tsv",
        auspice_config = "config/auspice_config.json",
        colors = "config/colors.tsv",
        description = "config/description.md"

files = rules.files.params

rule filter:
    message:
        """
        Filtering to
          - excluding strains in {input.exclude}
          - minimum genome length of {params.min_length}
        """
    input:
        sequences = files.input_sequences,
        metadata = files.input_metadata,
        include = files.include,
        exclude = files.exclude
    output:
        sequences = "results/filtered.fasta"
    params:
        group_by = "year",
        sequences_per_group = 5,
        min_length = 5000
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --include {input.include} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-length {params.min_length}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference
        """

rule mask:
    message:
        """
        Mask bases in alignment
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        """
    input:
        alignment = rules.align.output.alignment
    output:
        alignment = "results/masked.fasta"
    params:
        mask_from_beginning = 15,
        mask_from_end = 15,
        mask_sites = 18460
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --output {output.alignment}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.mask.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.mask.output,
        metadata = files.input_metadata
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        root = "BtKY72"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --root {params.root} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """

rule prune_outgroup:
    message: "Pruning the outgroup from the tree"
    input:
        tree = rules.refine.output.tree
    output:
        tree = "results/tree_pruned.nwk"
    params:
        root = "BtKY72"
    run:
        from Bio import Phylo
        T = Phylo.read(input[0], "newick")
        outgroup = [c for c in T.find_clades() if str(c.name) == params[0]][0]
        T.prune(outgroup)
        Phylo.write(T, output[0], "newick")

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.prune_outgroup.output.tree,
        alignment = rules.mask.output
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.prune_outgroup.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule clades:
    message: " Labeling clades as specified in config/clades.tsv"
    input:
        tree = rules.prune_outgroup.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = files.clades
    output:
        node_data = "results/clades.json"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output {output.node_data}
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.prune_outgroup.output.tree,
        metadata = files.input_metadata,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        clades = rules.clades.output.node_data,
        auspice_config = files.auspice_config,
        colors = files.colors,
        description = files.description
    output:
        auspice_json = rules.all.input.auspice_json
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.clades} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --description {input.description} \
            --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
