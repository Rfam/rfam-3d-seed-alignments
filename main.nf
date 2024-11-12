nextflow.preview.output = true

process fetch_svn_files {
  tag { "$accession" }
  // Try not to overwhelm the SVN server, but retry when it gets overwhelmed
  maxForks 10
  errorStrategy 'retry'
  maxRetries 5

  input:
  tuple val(url), val(accession)

  output:
  tuple val(accession), path("${accession}.seed"), path("${accession}.cm"), emit: svn_files
  path "${accession}-info.jsonl", emit: alignment_info

  // TODO Reformat alignment so all sequences are on one line?
  """
  curl ${url}/${accession}/SEED > ${accession}.seed
  curl ${url}/${accession}/CM > ${accession}.cm
  cmbuild --hand -F ${accession}.cm ${accession}.seed
  rfam_3d parse-alignment ${accession} ${accession}.seed ${accession}-info.jsonl
  """
}

process merge_info {
  input:
  path("alignment-info*.jsonl")

  output:
  path("all-alignments.jsonl")

  """
  find . -name 'alignment-info*.jsonl' | xargs cat > all-alignments.jsonl
  """
}

process find_actions {
  // If starting without anything cached this can be very slow to fetch all
  // data.
  time '12h'

  input:
  tuple path(matches_file), path(disallow), path(info)

  output:
  path "missing-sequences/*.fa", emit: sequences
  path "pdb-info/*.json", emit: pdb_info
  path "report.txt", emit: report

  """
  rfam_3d compute-actions --disallow-file ${disallow} ${matches_file} ${info} report.txt missing-sequences/ pdb-info/
  """
}

process update_alignment {
  memory 5.GB
  tag { "$accession" }

  input:
  tuple val(accession), path(seed), path(cm), path(sequences)

  output:
  tuple val(accession), path("${accession}.extended.sto")

  """
  update-alignment "${seed}" "${cm}" "${sequences}" "${accession}.extended.sto"
  """
}

process add_structure_info {
  tag { "$accession" }
  errorStrategy 'ignore'

  input:
  tuple val(accession), path("alignment.sto"), path(info)

  output:
  path "${accession}.updated.sto", emit: updated

  """
  esl-reformat --informat stockholm pfam alignment.sto > "${accession}.pfam.sto"
  rfam_3d add-structure-information "${accession}.pfam.sto" "$info" "${accession}.info.sto"
  esl-reformat --informat pfam stockholm "${accession}.info.sto" > "${accession}.updated.sto"
  """
}

workflow add_3d {
  take:
    matches
  main:
    Channel.fromPath('config/limits.yaml') | set { disallow }

    matches \
    | map { it.text } \
    | splitCsv(sep: '\t')
    | map { it[0] } \
    | toList \
    | flatMap { it.unique() } \
    | map { [params.svn, it] } \
    | fetch_svn_files

    fetch_svn_files.out.alignment_info \
    | collect \
    | merge_info \
    | set { info }

    matches \
    | combine(disallow) \
    | combine(info) \
    | find_actions

    find_actions.out.sequences | flatten | map { [it.baseName, it] } | set { sequences }
    find_actions.out.pdb_info | flatten | map { [it.baseName, it] } | set { structures }

    fetch_svn_files.out.svn_files \
    | join(sequences) \
    | update_alignment \
    | join(structures) \
    | add_structure_info

  emit:
    report = find_actions.out.report
    alignments = add_structure_info.out.updated
}

workflow {
  main:
    add_3d(Channel.fromPath("mapping.tsv"))
  publish:
    add_3d.out.report >> 'report'
    add_3d.out.alignments >> 'alignments'
}

output {
  directory "$launchDir/results"
  mode 'copy'
}
