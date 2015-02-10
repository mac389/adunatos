cd src
python differential-pathway-analysis.py
echo 'Extracted regions of interest by z-score, deduplicated'
python add-uniprot.py
echo 'Found Uniprot alias for genes'
python add-GO.py
echo 'Mapped genes to gene ontology using Uniprot alias'
python compare-expression.py
echo 'Made figure comparing expression, printed list of upregulated genes to'