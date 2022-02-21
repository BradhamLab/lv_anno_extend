"""
Update Lvar 3.0 rRNA annotation gff to follow required gene-exon feature formats
for proper quantification.

@author: Dakota Hawkins
"""
import gffutils 

def sort_attributes(attrs):
    """Sort feature attributes putting ID an Parent keys first."""
    out = {}
    out['ID'] = attrs['ID']
    if 'Parent' in attrs:
        out['Parent'] = attrs['Parent']
    for key in attrs:
        if key not in ['ID', 'Parent']:
            out[key] = attrs[key]
    return out

def get_gene_and_exon_features(rRNA):
    """
    Create exon + gene annotations from written rRNA annotations.

    Updates IDs and adds parent relations as gene -> rRNA -> exon.

    Parameters
    ----------
        rRNA : gffutils.Feature
            rRNA feature from gff input.
    Returns
    -------
    (gffutils.Feature, gffutils.Feature)
        Gene and exon features, respectively. 
    """
    attrs = dict(rRNA.attributes.items())
    gene_id = "{}-{}".format(rRNA.id, 'gene')
    attrs['ID'] = [gene_id]
    rRNA.attributes['Parent'] = ['gene_id']
    rRNA.attributes = sort_attributes(rRNA.attributes)
    gene = gffutils.Feature(seqid=rRNA.chrom,
                            source=rRNA.source,
                            featuretype='gene',
                            start=rRNA.start,
                            end=rRNA.end,
                            id=gene_id,
                            strand=rRNA.strand,
                            attributes=sort_attributes(attrs))
    exon_id = "{}-{}".format(rRNA.id, 'exon')
    attrs['ID'] = [exon_id]
    attrs['Parent'] = [gene_id]
    exon = gffutils.Feature(seqid=rRNA.chrom,
                            source=rRNA.source,
                            featuretype='exon',
                            start=rRNA.start,
                            end=rRNA.end,
                            id=exon_id,
                            strand=rRNA.strand,
                            attributes=sort_attributes(attrs))
    return gene, exon

def main(db):
    """Create gene + exon annotations for each rRNA feature in the gff file."""
    out = []
    for rRNA in db.features_of_type('rRNA'):
        gene, exon = get_gene_and_exon_features(rRNA)
        out += [gene, rRNA, exon]
    return out

def write(outlines, fn):
    """Write gene, rRNA, and exon features to a new gff file."""
    with open(fn, 'w') as fout:
        for record in outlines:
            fout.write(str(record) + '\n')

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        db = gffutils.create_db(snakemake.input['gff'], ':memory:', 
                                merge_strategy='create_unique')
        outlines = main(db)
        write(outlines, snakemake.output['gff'])
