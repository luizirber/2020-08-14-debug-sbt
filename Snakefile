from pathlib import Path

SBT_ZIP = Path("/home/luizirber/work/sourmash-bio/sourmash_resources_debug/genbank-bacteria-x1e6-k51.sbt.zip")

rule all:
  input: "outputs/newtree.sbt.zip"

rule extract_sig_names:
  output: "outputs/signames.txt"
  input: SBT_ZIP
  run:
    import sourmash

    index = sourmash.load_file_as_index(input[0])
    with open(output[0], 'w') as out:
      for leaf in index.leaves():
        sig = leaf.data
        out.write(sig.name() + "\n")
        leaf.unload()

rule build_sbt:
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
