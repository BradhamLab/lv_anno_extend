import gffutils

def write(db, fn):
    """Write gene, rRNA, and exon features to a new gff file."""
    with open(fn, 'w') as fout:
        for record in db.all_features():
            fout.write(str(record) + '\n')

if __name__ == "__main__":
    try:
        snakemake
    except NameErorr:
        snakemake = None
    if snakemake is not None:
        db = gffutils.create_db(snakemake.input['gff'],
                                ":memory:",
                                merge_strategy="create-unique")
        with open(snakemake.input['overlaps'], 'r') as f:
            for line in f:
                gene_id = line.split('\t')[3]
                print(gene_id)
                children = db.children(gene_id)
                db = db.delete(children)
        write(db, snakemake.output['gff'])
        