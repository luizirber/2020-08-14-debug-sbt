from pathlib import Path

SBT_ZIP = Path("/home/luizirber/work/sourmash-bio/sourmash_resources_debug/genbank-bacteria-x1e6-k51.sbt.zip")

rule all:
  input:
    "outputs/signames_new_catalog.txt",
    "outputs/signames_new.txt"

def extract(inp, out):
    import sourmash

    index = sourmash.load_file_as_index(inp)
    with open(out, 'w') as out:
      for leaf in index.leaves():
        sig = leaf.data
        out.write(sig.name() + "\n")
        leaf.unload()

rule extract_sig_names:
  output: "outputs/signames.txt"
  input: SBT_ZIP
  run:
    extract(input[0], output[0])

rule extract_sig_names_new_sbt:
  output: "outputs/signames_new.txt"
  input: "outputs/newtree.sbt.zip"
  run:
    extract(input[0], output[0])

rule extract_sig_names_new_sbt_catalog:
  output: "outputs/signames_new_catalog.txt"
  input: "outputs/newtree_catalog.sbt.zip"
  run:
    extract(input[0], output[0])

rule extract_sig_names_subset:
  output: "outputs/signames_subset.txt"
  input: "outputs/duplicated_subset.sbt.zip"
  run:
    extract(input[0], output[0])


rule build_sbt_from_set:
  output: "outputs/newtree.sbt.zip"
  input: "outputs/signames.txt"
  run:
    import sourmash

    index = sourmash.create_sbt_index(bloom_filter_size=64)
    signames = set()
    with open(input[0], 'r') as f:
      for line in f:
        signames.add(line.strip())

    for i, signame in enumerate(signames):
      mh = sourmash.MinHash(0, 31, scaled=2000, mins=[i])
      new_sig = sourmash.SourmashSignature(mh, name=signame)
      index.insert(new_sig)
      if i % 1000 == 0:
        print(f"processed {i} sigs")

    index.save(output[0])

rule build_sbt_from_catalog:
  output: "outputs/newtree_catalog.sbt.zip"
  input: "outputs/catalog.txt"
  run:
    import sourmash

    index = sourmash.create_sbt_index(bloom_filter_size=64)
    signames = []
    with open(input[0], 'r') as f:
      for line in f:
        signames.append(line.strip())

    for i, signame in enumerate(signames):
      mh = sourmash.MinHash(0, 31, scaled=2000, mins=[i])
      new_sig = sourmash.SourmashSignature(mh, name=signame)
      index.insert(new_sig)
      if i % 1000 == 0:
        print(f"processed {i} sigs")

    index.save(output[0])

rule generate_catalog:
  output: "outputs/catalog.txt"
  input: "outputs/assembly_summary_genbank.txt.gz"
  run:
    import csv
    import gzip

    with open(output[0], 'w') as out:
      with gzip.open(input[0], 'rt') as fp:
        fp.readline() # skip first line
        fp.read(2) # skip initial comment in header
        data = csv.DictReader(fp, delimiter='\t')
        for row in data:
          name_parts = [row["assembly_accession"], " ", row['organism_name']]
          if row['infraspecific_name']:
            name_parts += [" ", row['infraspecific_name']]
          name_parts += [', ', row['asm_name']]
          out.write("".join(name_parts) + "\n")

rule duplicated_sigs:
  output: "outputs/duplicated.txt"
  input: "outputs/signames.txt"
  shell: """
    sort {input} | uniq -c | sort -nr | \
    grep -v '  1 ' | tr -s ' ' | cut -d ' ' -f3 > {output}
  """

rule new_subset_catalog:
  output: "outputs/duplicated_subset.txt"
  input:
    dup_accs = "outputs/duplicated.txt",
    sig_links = "/home/irber/sourmash_databases/outputs/catalog/bacteria/genbank.txt"
  run:
    accessions = set()
    with open(input.dup_accs, 'r') as f:
      for line in f:
        accessions.add(line.strip())

    with open(output[0], 'w') as out:
      with open(input.sig_links, 'r') as f:
        for line in f:
            sig_acc = line.strip().split("/")[-1][:-4]
            if sig_acc in accessions:
              out.write(line)

rule duplicated_sbt:
  output: "outputs/duplicated_subset.sbt.zip"
  input: "outputs/duplicated_subset.txt"
  params:
    abs_output = lambda w, input, output: os.path.abspath(output[0]),
    abs_input = lambda w, input, output: os.path.abspath(input[0]),
  shell: """
    cd ~/sourmash_databases
    sourmash index -k 51 -x 64 --from-file <(tail +2 {params.abs_input}) {params.abs_output} $(head -1 {params.abs_input})
  """
